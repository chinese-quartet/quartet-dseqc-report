(ns quartet-dseqc-report.core
  (:require [tservice-core.tasks.http :as http-task]
            [quartet-dseqc-report.spec :as spec]
            [quartet-dseqc-report.task :as task]))

(def metadata
  (http-task/make-routes "quartet-dseqc-report" :ReportPlugin
                         {:method-type :post
                          :endpoint "quartet-dseqc-report"
                          :summary "Generate the QC Report for Quartet Genomics Data."
                          :body-schema spec/quartet-dseqc-report-params-body
                          :response-schema any?
                          :handler task/post-handler}))

(def events-init task/events-init)
