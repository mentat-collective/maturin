(require 'maturin.better)
(in-ns 'maturin.better)

;; create xs
(defn gen-xs [n]
  (take n (repeatedly #(rand 10))))

;; create linear regression mock data
(defn gen-ys [xs slope intercept]
  (map #(+ (d/sample (d/normal 0 1)) intercept (* % slope)) xs))

;; return linear function
(defn linear-func [slope intercept]
  (fn [x] (+ intercept (* slope x))))

(defn line->offsets [xs ys]
  (fn [[slope intercept]]
    (let [preds (map (linear-func slope intercept) xs)]
      (mapv - preds ys))))

(defn L-springs [k]
  (fn [diffs]
    (- (* 1/2 k (square diffs)))))

(defn L-regression-composed [xs ys k]
  (compose (L-springs k) (F->C (line->offsets xs ys))))

;; lagrangian for the springs-rod system, regression case:
;; springs constrained to have the same x-coordiate on both axes = minimize MSE
(defn L-regression
  [xs ys k]
  (fn [[_ [slope intercept] _]]
    (let [preds (map (linear-func slope intercept) xs)
          diffs (mapv - preds ys)]
      (- (* 1/2 k (square diffs))))))

(defn distance-point-line [slope intercept]
  (fn [x y]
    (/ (Math/abs (+ (* slope x) (- y) intercept)) (Math/sqrt (square [slope 1])))))

;; lagrangian for the springs-rod system, PCA case:
;; springs can "slide" on the bar's end = minimize point-line distance
(defn L-PCA
  [xs ys k]
  (fn [[_ [slope intercept] _]]
    (let [diffs (mapv (distance-point-line slope intercept) xs ys)]
      (- (* 1/2 k (square diffs))))))


;; trying to draw something,
;; prettify later


(defn draw-line [convert]
  (fn [{:keys [state color]}]
    (q/background 100)
    (q/fill color 255 255)
    (let [[slope intercept] (convert state)]
      (attach [0 intercept] [10 ((linear-func slope intercept) 10)])
      (bob [0 intercept])
      (bob [10 ((linear-func slope intercept) 10)]))))


;; fake data        


(def x (gen-xs 10))
(def y (gen-ys x 2 3))

;; animate
(let [L (L-regression-composed x y 1)
      initial-state (up 1
                        (up 1 1)
                        1)
      built (build L initial-state)]
  (q/defsketch regression
    :title "Regression line"
    :size [500 500]
    :setup (:setup built)
    :update (:update built)
    :draw (draw-line (:xform built))
    :features [:keep-on-top]
    :middleware [m/fun-mode m/navigation-2d]))