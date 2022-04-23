(ns quartet-dseqc-report.spec
  (:require [clojure.spec.alpha :as s]
            [spec-tools.core :as st]))

(s/def ::filepath
  (st/spec
   {:spec                (s/and string? #(re-matches #"^[a-zA-Z0-9]+:\/\/(\/|\.\/)[a-zA-Z0-9_]+.*" %))
    :type                :string
    :description         "File path for genomics profiled data, such as file:///xxx/xxx/data.csv"
    :swagger/default     nil
    :reason              "The filepath must be string."}))

(s/def ::name
  (st/spec
   {:spec                string?
    :type                :string
    :description         "The name of the report"
    :swagger/default     ""
    :reason              "Not a valid name"}))

(s/def ::description
  (st/spec
   {:spec                string?
    :type                :string
    :description         "Description of the report"
    :swagger/default     ""
    :reason              "Not a valid description."}))

(def quartet-dseqc-report-params-body
  "A spec for the body parameters."
  (s/keys :req-un [::name ::filepath]
          :opt-un [::description]))
