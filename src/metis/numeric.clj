(ns metis.numeric)


;; Following functions are taken from:
;; http://java.ociweb.com/mark/clojure/article.html
(defn polynomial
  "computes the value of a polynomial
   with the given coefficients for a given value x"
  [coefs x]
  (reduce #(+ (* x %1) %2) coefs))

(defn derivative
  "computes the value of the derivative of a polynomial
   with the given coefficients for a given value x"
  [coefs x]
  (let [exponents (reverse (range (count coefs)))
        derivative-coefs (map #(* %1 %2) (butlast coefs) exponents)]
    (polynomial derivative-coefs x)))

;; so the antiderivative should be:
(defn antiderivative
  "computes the value of the antiderivative of a polynomial
   with the given coefficients for a given value x"
  [coefs c x]
  (let [exponents (reverse (range 1 (inc (count coefs))))
        ad-coefs (conj (vec (map #(/ %1 %2) coefs exponents)) 0)]
    (polynomial ad-coefs x)))

;; integration
(defn polynomial-integral
  "computes the value of the integral for a polynomial"
  [coefs a b]
  (- (antiderivative coefs 0 b) (antiderivative coefs 0 a)))

(defn defactor
  "defactors a polynomial and returns the coefficients:
   (x-a)(x-b) -> (defactor [a b]) -> [1 -a-b a*b]"
  [v]
  (let [n (count v)
        x (map unchecked-negate v)]
    (loop [i 1
           coefs [1 (nth x 0)]]
      (if (< i n)
        (recur (inc i) (defac-rec coefs (nth x i)))
        coefs)
      )))

(defn defac-rec
  "returns the new coefficients of a polynomial with
   the given coefficients multiplied by (x - c):
   [a b] c -> (ax + b)(x - c) -> [a b-ac -bc]"
  [coefs c]
  (let [x (cons 0 (conj coefs 0))]
    (vec (map + (next x) (map #(* % c) x)))))


(defn lagrange-weights
  "computes the Lagrange weights for n points (not n-th degree);
   (Lagrange 3)*((b-a)/2) -> weights for Simpson's rule"
  [n]
  (let [coefs (map #(into (range %) (range (inc %) n))
                   (range n))
        div (map #(reduce * (map - (repeat %) (nth coefs %)))
                 (range n))
        integral (map #(polynomial-integral (defactor %) 0 (dec n)) coefs)]
      (map / integral div)))


(defn newton-cotes
  "computes the integral between a and b for a function f
   with the closed Newton-Cotes formula of degree d"
  ;; d < 7 recommended (else: rounding errors and overflow)
  [f a b d]
  (let [h (/ (- b a) d)
        w (lagrange-weights (inc d))
        y (map f (range a (inc b) d))
        res (reduce + (map * w y))]
    (* h res)))

;; likely to prove inefficient
(defn newton-cotes-summed
  "computes the integral between a and b for a function f
   with the closed Newton-Cotes formula of degree d applied
   to n partitions independently"
  [f a b d n]
  (let [h (/ (- b a) n)
        parts (range a (inc b) h)]
    (reduce + (map #(newton-cotes f %1 %2 d) parts (next parts)))))
