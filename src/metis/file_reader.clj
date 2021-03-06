(ns metis.file-reader
  (:use [clojure.string]))


(defn read-file
  [path-to-file]
  (slurp path-to-file))


(defn- extract-topics
  "Extracts the topics (e.g. the column names) of the csv file"
  [raw-data]
  (rest (split (first (split raw-data #"\n")) #",")))


; only for test purposes
(defn format-data
  "Formats data into hashmap from testfile"
  [data]
  (apply hash-map
         (split (reduce #(str %1 " " %2)
                        (split data #"\n"))#"\s")))


(defn create-country-entry
  "Format big csv data"
  [raw-entry topics]
  {:country (first raw-entry) :data ()})


(defn parse-number
  "Reads a number from a string. Returns nil if not a number."
  [s]
  (if (re-find #"^-?\d+\.?\d*$" s)
    (read-string s)))
