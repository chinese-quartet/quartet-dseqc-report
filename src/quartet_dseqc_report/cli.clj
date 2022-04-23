(ns quartet-dseqc-report.cli
  (:gen-class)
  (:require [quartet-dseqc-report.task :refer [make-report!]]
            [local-fs.core :refer [file? directory?]]
            [clojure.string :as clj-str]
            [clojure.tools.cli :refer [parse-opts]]
            [quartet-dseqc-report.version :refer [version]]))

(def cli-options
  [["-d" "--data PATH" "Data Directory"
    :validate [#(directory? %) "Must be a valid directory"]]
   ["-o" "--output PATH" "Result Directory"
    :validate [#(directory? %) "Must be a valid directory."]]
   ["-n" "--name NAME" "Report name"
    :default "report"]
   ["-D" "--description DESC" "Report Description"
    :default "Quality control report"]
   ["-v" "--version" "Show version" :default false]
   ["-h" "--help"]])

(defn usage [options-summary]
  (->> ["DSeQC - Visualizes Quality Control(QC) results for Quartet Project."
        ""
        "Usage: dseqc [options]"
        ""
        "Options:"
        options-summary
        ""
        "Please refer to the manual page for more information."]
       (clj-str/join \newline)))

(defn error-msg [errors]
  (str "The following errors occurred while parsing your command:\n\n"
       (clj-str/join \newline errors)))

(defn validate-args
  "Validate command line arguments. Either return a map indicating the program
  should exit (with an error message, and optional ok status), or a map
  indicating the action the program should take and the options provided."
  [args]
  (let [{:keys [options arguments errors summary]} (parse-opts args cli-options)]
    (cond
      (:help options) ; help => exit OK with usage summary
      {:exit-message (usage summary) :ok? true}

      errors ; errors => exit with description of errors
      {:exit-message (error-msg errors)}

      ;; custom validation on arguments
      (:version options)
      {:exit-message (format "v%s" version)}

      (nil? (:data options))
      {:exit-message "You need to specified -d/--data argument."}

      (nil? (:output options))
      {:exit-message "You need to specified -o/--output argument."}

      (and (:data options) (:name options) (:output options))
      {:options options}

      :else ; failed custom validation => exit with usage summary
      {:exit-message (usage summary)})))

(defn exit [status msg]
  (println msg)
  (System/exit status))

(defn -main
  "Generate DSeQC report for quartet project."
  [& args]
  (let [{:keys [options exit-message ok?]} (validate-args args)]
    (if exit-message
      (exit (if ok? 0 1) exit-message)
      (make-report! {:data-dir (:data options)
                     :dest-dir (:output options)
                     :parameters {:name (:name options)
                                  :description (:description options)
                                  :plugin-name "quartet-dseqc-report"
                                  :plutin-type "ReportPlugin"
                                  :plugin-version version}
                     :task-id nil}))
    (shutdown-agents)))
