(ns metis.numeric-test
  (:require [clojure.test :refer :all]
            [midje.sweet :refer [facts fact]]
            [metis.numeric :refer :all]))


(def f (partial polynomial [2 1 3])) ; f(x) = 2x^2 + x + 3
(def f-prime (partial derivative [2 1 3])) ; f'(x) = 4x + 1
(def F (partial antiderivative [2 1 3] 5)) ; F(x) = 2/3*x^3 + 1/2*x^2 + 3x + 5


(fact "'defactor' succeeds"
      (defactor [2 3]) => [1 -5 6])

(fact "'lagrange-weights' succeeds"
      (lagrange-weights 3) => [1/3 4/3 1/3])

(fact "'polynomial-integral' succeeds"
      (polynomial-integral [2 1 3] 0 2) => (F 2))

(fact "'newton-cotes' succeeds"
      (newton-cotes f 0 2 3) => (F 2))

(fact "'newton-cotes-summed' succeds"
      (newton-cotes-summed f 0 2 3 5) => (F 2))
