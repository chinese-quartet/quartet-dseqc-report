(ns quartet-dseqc-report.dseqc
  "A wrapper for dseqc tool."
  (:require [tservice-core.plugins.env :refer [get-context-path add-env-to-path]]
            [local-fs.core :as fs-lib]
            [clojure.java.shell :as shell :refer [sh]]
            [clj-yaml.core :as yaml]
            [clojure.data.csv :as csv]
            [clojure.string :as clj-str]
            [clojure.java.io :as io]
            [remote-fs.core :as remote-fs]
            [tservice-core.plugins.util :refer [call-command!]]
            [quartet-dseqc-report.version :as v]
            [clojure.tools.logging :as log])
  (:import [org.apache.commons.io.input BOMInputStream]))

(defn sort-exp-data
  [coll]
  (sort-by :GENE_ID coll))

(defn csv-data->maps [csv-data]
  (map zipmap
       (->> (first csv-data) ;; First row is the header
            (map keyword)    ;; Drop if you want string keys instead
            repeat)
       (rest csv-data)))

(defn bom-reader
  "Remove `Byte Order Mark` and return reader"
  [filepath]
  (-> filepath
      io/input-stream
      BOMInputStream.
      io/reader))

(defn guess-separator
  [filepath]
  (with-open [reader (bom-reader filepath)]
    (let [header (first (line-seq reader))
          seps [\tab \, \; \space]
          sep-map (->> (map #(hash-map % (count (clj-str/split header (re-pattern (str %))))) seps)
                       (into {}))]
      (key (apply max-key val sep-map)))))

(defn read-csv
  [^String file]
  (when (.isFile (io/file file))
    (with-open
     [reader (io/reader file)]
      (doall
       (->> (csv/read-csv reader :separator (guess-separator file))
            csv-data->maps)))))

(defn vec-remove
  "Remove elem in coll"
  [pos coll]
  (vec (concat (subvec coll 0 pos) (subvec coll (inc pos)))))

(defn write-csv!
  "Write row-data to a csv file, row-data is a vector that each element is a map."
  [path row-data]
  (let [columns (keys (first row-data))
        headers (map name columns)
        rows (mapv #(mapv % columns) row-data)]
    (with-open [file (io/writer path)]
      (csv/write-csv file (cons headers rows) :separator \tab))))

(defn write-csv-by-cols! [path row-data columns]
  (let [headers (map name columns)
        rows (mapv #(mapv % columns) row-data)]
    (with-open [file (io/writer path)]
      (csv/write-csv file (cons headers rows)))))

(defn read-csvs
  [files]
  (map #(sort-exp-data (read-csv %)) files))

(defn reorder
  [data]
  (let [cols (vec (sort (keys (first data))))]
    (cons :GENE_ID (vec-remove (.indexOf cols :GENE_ID) cols))))

(defn merge-exp
  "[[{:GENE_ID 'XXX0' :YYY0 1.2} {:GENE_ID 'XXX1' :YYY1 1.3}]
    [{:GENE_ID 'XXX0' :YYY2 1.2} {:GENE_ID 'XXX1' :YYY3 1.3}]]"
  [all-exp-data]
  (apply map merge all-exp-data))

(defn write-csv-by-ordered-cols!
  [path row-data]
  (let [cols (reorder row-data)]
    (write-csv-by-cols! path row-data cols)))

(defn merge-exp-files!
  "Assumption: all files have the same GENE_ID list, no matter what order."
  [files path]
  (->> (read-csvs files)
       (merge-exp)
       (write-csv-by-ordered-cols! path)))

(defn gen-multiqc-config
  [config-file]
  (let [config (yaml/generate-string
                {:run_modules ["dnaseq_data_generation_information"
                               "dnaseq_performance_assessment"
                               "dnaseq_raw_qc"
                               "dnaseq_post_alignment_qc"
                               "dnaseq_quantification_qc"
                               "dnaseq_supplementary"]
                 :skip_generalstats true}
                :dumper-options {:flow-style :block})]
    (spit config-file config)))

(defn decompression-tar
  [filepath]
  (shell/with-sh-env {:PATH   (add-env-to-path "quartet-dnaseq-report")
                      :LC_ALL "en_US.utf-8"
                      :LANG   "en_US.utf-8"}
    (let [command ["bash" "-c"
                   (format "tar -xvf %s -C %s" filepath (fs-lib/parent-path filepath))]
          result  (apply sh command)
          status (if (= (:exit result) 0) "Success" "Error")
          msg (str (:out result) "\n" (:err result))]
      {:status status
       :msg msg})))

(defn fs-service?
  [filepath]
  (re-matches #"^[a-zA-Z0-9]+:\/\/.*" filepath))

(defn not-empty?
  [coll]
  ((complement empty?) coll))

(defn parse-path
  "Parse path and extract protocol, bucket, and object path."
  [path]
  (let [path-lst (rest (re-find (re-matcher #"^([a-zA-Z0-9]+):\/\/([^\/]+)\/(.*)" path)))]
    (when (not-empty? path-lst)
      {:protocol (first path-lst)
       :bucket (second path-lst)
       :prefix (nth path-lst 2)})))

(defn make-link
  [protocol bucket object]
  (str protocol "://" bucket "/" object))

(defn filter-remote-fn
  [mode]
  (cond
    (= mode "directory") #".*\/"
    (= mode "file") #".*[^\/]$"
    :else #".*"))

(defn filter-local-fn
  [mode]
  (cond
    (= mode "directory") (fn [file] (.isDirectory (io/file file)))
    (= mode "file") (fn [file] (.isFile (io/file file)))
    :else (fn [file] file)))

(defn objects->links
  [protocol bucket objects mode]
  (->> objects
       (map #(make-link protocol bucket (:key %)))
       (filter #(re-matches (filter-remote-fn mode) %))))

(defn list-files
  "Local path must not contain file:// prefix. options - directory | file"
  [path & options]
  (let [{:keys [mode]} (first options)
        {:keys [protocol bucket prefix]} (parse-path path)
        is-service? (fs-service? path)]
    (if is-service?
      (->> (concat
            (map (fn [item] (remote-fs/with-conn protocol (remote-fs/list-objects bucket (:key item))))
                 (remote-fs/with-conn protocol (remote-fs/list-objects bucket prefix true))))
           (flatten)
           (map #(make-link protocol bucket (:key %)))
           (filter #(re-matches (filter-remote-fn mode) %)))
      (->> (io/file path)
           file-seq
           (map #(.getAbsolutePath %))
           (filter #((filter-local-fn mode) %))))))

(defn make-pattern-fn
  [patterns]
  (map #(re-pattern %) patterns))

(defn filter-files
  [all-files pattern]
  (filter #(re-matches (re-pattern pattern) %) all-files))

(defn batch-filter-files
  [path patterns]
  (-> (map #(filter-files (list-files path) %)
           (make-pattern-fn patterns))
      flatten
      dedupe))

(defn copy-local-files!
  ":replace-existing, :copy-attributes, :nofollow-links"
  [files dest-dir options]
  (doseq [file-path files]
    (let [file (io/file file-path)
          dest (fs-lib/join-paths dest-dir (fs-lib/base-name file-path))]
      (if (.isFile file)
        (fs-lib/copy file-path dest options)
        (fs-lib/copy-recursively file-path dest options)))))

(defn copy-local-file!
  [file-path dest-dir options]
  (let [file (io/file file-path)
        dest (fs-lib/join-paths dest-dir (fs-lib/base-name file-path))]
    (if (.isFile file)
      (fs-lib/copy file-path dest options)
      (fs-lib/copy-recursively file-path dest options))))

(defn basename
  [path]
  (let [groups (re-find #".*\/(.*\/?)" path)]
    (if groups
      (first (rest groups))
      nil)))

(defn dirname
  [filepath]
  (let [groups (re-find #".*\/(.*)\/$" filepath)]
    (if groups
      (first (rest groups))
      nil)))

(defn copy-remote-file!
  [file-path dest-dir options]
  ;; TODO: Any options can help to improve the action.
  (let [is-dir? (re-matches #".*\/" file-path)
        filename (basename file-path)
        {:keys [protocol bucket prefix]} (parse-path file-path)]
    (if is-dir?
      (let [dest-dir (fs-lib/join-paths dest-dir (dirname file-path))]
        (fs-lib/create-directory dest-dir)
        (map #(copy-remote-file! % dest-dir options)
             (list-files file-path)))
      (remote-fs/with-conn protocol (remote-fs/download-object bucket prefix (fs-lib/join-paths dest-dir filename))))))

(defn copy-files!
  ":replace-existing, :copy-attributes, :nofollow-links"
  [files dest-dir options]
  (doseq [file-path files]
    (if (fs-service? file-path)
      (copy-remote-file! file-path dest-dir options)
      (copy-local-file! file-path dest-dir options))))

(defn multiqc
  "A multiqc wrapper for generating multiqc report:
   TODO: set the absolute path of multiqc binary instead of environment variable

  Required:
  analysis-dir: Analysis directory, e.g. data directory from project
  outdir: Create report in the specified output directory.

  Options:
  | key                | description |
  | -------------------|-------------|
  | :dry-run?          | Dry run mode |
  | :filename          | Report filename. Use 'stdout' to print to standard out. |
  | :comment           | Custom comment, will be printed at the top of the report. |
  | :title             | Report title. Printed as page header, used for filename if not otherwise specified. |
  | :force?            | Overwrite any existing reports |
  | :prepend-dirs?     | Prepend directory to sample names |
  | :template          | default, other custom template    |
  | :config            | Where is the config file          |
  | :env               | An environemnt map for running multiqc, such as {:PATH (get-path-variable)} |

  Example:
  (multiqc 'XXX' 'YYY' {:filename       'ZZZ'
                        :comment        ''
                        :title          ''
                        :force?         true
                        :prepend-dirs?  true})"
  [analysis-dir outdir {:keys [dry-run? filename comment title force? prepend-dirs? template config env]
                        :or   {dry-run?      false
                               force?        true
                               prepend-dirs? false
                               filename      "multiqc_report.html"
                               comment       ""
                               template      "default"
                               title         "iSEQ Analyzer Report"}}]
  (let [force-arg   (if force? "--force" "")
        dirs-arg    (if prepend-dirs? "--dirs" "")
        config-arg  (if config (str "-c " config) "")
        multiqc-command (filter #(> (count %) 0) ["multiqc"
                                                  force-arg dirs-arg config-arg
                                                  "--title" (format "'%s'" title)
                                                  "--comment" (format "'%s'" comment)
                                                  "--filename" filename
                                                  "--outdir" outdir
                                                  "-t" template
                                                  analysis-dir])
        command (clj-str/join " " multiqc-command)]
    (if dry-run?
      (log/info command)
      (if env
        (call-command! command env)
        (call-command! command)))))

(defn is-localpath?
  [filepath]
  (re-matches #"^file:\/\/.*" filepath))

(defn correct-filepath
  [filepath]
  (if (is-localpath? filepath)
    (clj-str/replace filepath #"^file:\/\/" "")
    filepath))
