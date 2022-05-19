(ns quartet-dseqc-report.task
  (:require [quartet-dseqc-report.dseqc :as dseqc]
            [local-fs.core :as fs-lib]
            [tservice-core.plugins.env :refer [get-workdir make-remote-link add-env-to-path create-task! update-task!]]
            [tservice-core.plugins.util :as util]
            [clojure.string :as clj-str]
            [clojure.data.json :as json]
            [clojure.tools.logging :as log]
            [tservice-core.tasks.async :refer [publish-event! make-events-init]]
            [quartet-dseqc-report.version :as v]))

(defn date
  []
  (.format (java.text.SimpleDateFormat. "yyyy-MM-dd")
           (new java.util.Date)))

(defn update-process!
  [^String task-id ^Integer percentage]
  (let [record (cond
                 (= percentage 100) {:status "Finished"
                                     :percentage 100
                                     :finished_time (util/time->int (util/now))}
                 (= percentage -1) {:status "Failed"
                                    :finished_time (util/time->int (util/now))}
                 :else {:percentage percentage})
        record (merge {:id task-id} record)]
    (update-task! record)))

(defn update-log-process!
  "Update message into log file and process into database.
   log-coll: {:msg \"xxx\" :status \"Failed\"}"
  [log-path log-coll task-id process]
  (log/info (:status log-coll)
            (format "%s... (More details on %s)"
                    (apply str (take 100 (:msg log-coll)))
                    log-path))
  (let [exists? (fs-lib/exists? log-path)
        results (if exists?
                  [(json/read-str (slurp log-path) :key-fn keyword) log-coll]
                  [log-coll])
        status (:status (last results))
        msg (apply str (map :msg results))]
    ;; (println results status msg)
    (spit log-path (json/write-str {:status status :msg msg}))
    (update-process! task-id process)))

(defn post-handler
  [{:keys [body owner plugin-context uuid workdir]
    :as payload}]
  (log/info (format "Create a report with %s" payload))
  (let [{:keys [name filepath description]
         :or {description (format "Quality control report for %s" name)}} body
        payload (merge {:description description} (:body payload))
        ;; slash in the end of filepath is crucial when you want to use it to filter oss files
        data-dir (str (clj-str/replace (dseqc/correct-filepath filepath) #"/$" "") "/")
        log-path (fs-lib/join-paths workdir "log")
        response {:report (make-remote-link (format "%s/multiqc_report.html" workdir))
                  :log (make-remote-link log-path)}
        task-id (create-task! {:id             uuid
                               :name           name
                               :description    description
                               :payload        payload
                               :owner          owner
                               :plugin-name    v/plugin-name
                               :plugin-type    "ReportPlugin"
                               :plugin-version (:plugin-version plugin-context)
                               :response       response})
        result-dir (fs-lib/join-paths workdir "results")]
    (fs-lib/create-directories! result-dir)
    (spit log-path (json/write-str {:status "Running"
                                    :msg ""}))
    (update-process! task-id 0)
    (publish-event! "quartet_dseqc_report"
                    {:data-dir data-dir
                     :dest-dir workdir
                     :task-id task-id
                     :parameters {:name name
                                  :description description
                                  :plugin-name v/plugin-name
                                  :plutin-type "ReportPlugin"
                                  :plugin-version (:plugin-version plugin-context)}})
    response))

(defn- filter-mkdir-copy
  [fmc-datadir fmc-patterns fmc-destdir fmc-newdir]
  (let [files-keep (dseqc/batch-filter-files fmc-datadir fmc-patterns)
        files-keep-dir (fs-lib/join-paths fmc-destdir fmc-newdir)]
    (fs-lib/create-directories! files-keep-dir)
    (if (empty? files-keep)
      (log/warn (format "Cannot find any files with pattern %s, please check your data." fmc-patterns))
      (dseqc/copy-files! files-keep files-keep-dir {:replace-existing true}))))

(defn copy-files-to-dir
  [data-dir dest-dir]
  (let [basename (fs-lib/base-name data-dir)
        dest-dir (fs-lib/join-paths dest-dir basename)]
    (filter-mkdir-copy (format "%s%s" data-dir "call-extract_tables") [".*.txt"] dest-dir "call-extract_tables")
    (filter-mkdir-copy (format "%s%s" data-dir "call-qualimap_D5") [".*zip"] dest-dir "call-qualimap_D5")
    (filter-mkdir-copy (format "%s%s" data-dir "call-qualimap_D6") [".*zip"] dest-dir "call-qualimap_D6")
    (filter-mkdir-copy (format "%s%s" data-dir "call-qualimap_F7") [".*zip"] dest-dir "call-qualimap_F7")
    (filter-mkdir-copy (format "%s%s" data-dir "call-qualimap_M8") [".*zip"] dest-dir "call-qualimap_M8")
    (filter-mkdir-copy (format "%s%s" data-dir "call-fastqc_D5") [".*.(zip|html)"] dest-dir "call-fastqc_D5")
    (filter-mkdir-copy (format "%s%s" data-dir "call-fastqc_D6") [".*.(zip|html)"] dest-dir "call-fastqc_D6")
    (filter-mkdir-copy (format "%s%s" data-dir "call-fastqc_F7") [".*.(zip|html)"] dest-dir "call-fastqc_F7")
    (filter-mkdir-copy (format "%s%s" data-dir "call-fastqc_M8") [".*.(zip|html)"] dest-dir "call-fastqc_M8")
    (filter-mkdir-copy (format "%s%s" data-dir "call-merge_mendelian") [".*.summary.txt"] dest-dir "call-merge_mendelian")))

(defn make-report!
  "Chaining Pipeline: filter-files -> copy-files -> multiqc."
  [{:keys [data-dir parameters dest-dir task-id]}]
  (log/info "Generate quartet dnaseq report: " data-dir parameters dest-dir)
  (let [parameters-file (fs-lib/join-paths dest-dir "general_information.json")
        log-path (fs-lib/join-paths dest-dir "log")
        subdirs (dseqc/list-dirs data-dir)]
    (log/info "List subdirs: " subdirs)
    (try
      (doseq [subdir subdirs]
        (copy-files-to-dir subdir dest-dir))
      (update-log-process! log-path {:status "Running" :msg "Download all files sucessfully.\n"} task-id 10)
      (spit parameters-file (json/write-str parameters))
      (doseq [files-qualimap-tar (dseqc/batch-filter-files dest-dir [".*qualimap.zip"])]
        (dseqc/decompression-tar files-qualimap-tar))
      (update-log-process! log-path {:status "Running" :msg "Prepare results successfully.\n"} task-id 50)
      (update-process! task-id 50)
      (spit parameters-file (json/write-str {"Report Name" (or (:name parameters) "Quartet QC Report for DNA-Seq")
                                             "Description" (or (:description parameters) "Visualizes Quality Control(QC) Results for Quartet DNA-Seq Data.")
                                             "Report Tool" (format "%s-%s"
                                                                   (:plugin-name parameters)
                                                                   (:plugin-version parameters))
                                             "Team" "Quartet Team"
                                             "Date" (date)}))
      (let [result (dseqc/multiqc dest-dir dest-dir {:template "report_templates"
                                                     :title "Quartet DNA report"
                                                     :env {:PATH (add-env-to-path "quartet-dseqc-report")}})]
        (if (= (:status result) "Error")
          (throw (Exception. (:msg result)))
          (update-log-process! log-path result task-id 100)))
      (catch Exception e
        (update-log-process! log-path {:status "Error" :msg (.toString e)} task-id -1)))))

(def events-init
  "Automatically called during startup; start event listener for quartet_dseqc_report events."
  (make-events-init "quartet_dseqc_report" make-report!))
