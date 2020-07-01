(ns maturin.better
  "Next, I think we want to transform this into a Middleware for Quil: https://github.com/quil/quil/wiki/Middleware"
  (:refer-clojure :exclude [partial zero? + - * / ref])
  (:require [maturin.ode :as o]
            [sicmutils.env :refer :all]
            [sicmutils.structure :as struct]
            [quil.core :as q]
            [quil.middleware :as m]))

(defn timestamp
  "TODO update this to a cljc compatible catch."
  []
  (try (/ (q/millis) 1000.0)
       (catch Exception _ 0)))

(defn setup-fn
  "I THINK this can stay the same for any of the Lagrangians."
  [initial-state]
  (fn []
    (q/frame-rate 60)
    ;; Set color mode to HSB (HSV) instead of default RGB.
    (q/color-mode :hsb)
    {:state initial-state
     :time (timestamp)
     :color 0
     :tick 0
     :navigation-2d {:zoom 4}}))

(defn Lagrangian-updater
  "Returns an update function that uses the supplied Lagrangian to tick forward."
  [L initial-state]
  (let [state-derivative (Lagrangian->state-derivative L)
        {:keys [integrator equations dimension] :as m}
        ((o/little-integrator (constantly state-derivative) [])
         initial-state
         :epsilon 1e-6
         :compile true)
        buffer (double-array dimension)]
    (fn [{:keys [state time color tick] :as m}]
      (let [s (double-array (flatten state))
            t2 (timestamp)]
        (.integrate integrator equations time s t2 buffer)
        (merge m {:color (mod (inc color) 255)
                  :state (struct/unflatten buffer state)
                  :time t2
                  :tick (inc tick)})))))

;; # API Attempt

(defn init
  "Generates an initial state dictionary."
  [L]
  {:lagrangian L
   :transforms []})

(defn transform
  "Takes a state dictionary and stores the coordinate transformation, AND composes
  it with the Lagrangian."
  [m f]
  (let [xform (F->C f)]
    (update-in m [:transforms] conj xform)))

(defn build
  "Returns the final keys for the sketch."
  [{:keys [lagrangian transforms]} initial-state]
  (let [xform (apply compose (reverse transforms))
        L (compose lagrangian xform)]
    {:setup (setup-fn initial-state)
     :update (Lagrangian-updater L initial-state)
     :xform (compose coordinate xform)}))

;; Next, let's do the double pendulum. This is going to require some manual work
;; to get the coordinate changes working.

(defn L-rectangular
  "Accepts:

  - a down tuple of the masses of each particle
  - a potential function of all coordinates

  Returns a function of local tuple -> the Lagrangian for N cartesian
  coordinates."
  [masses U]
  (fn [[_ q qdot]]
    (- (* 1/2 masses (mapv square qdot))
       (U q))))

(defn U-uniform-gravity
  "Accepts:

  - a down tuple of each particle's mass
  - g, the gravitational constant
  - optionally, the cartesian coordinate index of the vertical component. You'd
    need to update this if you went into 3d.

  Returns a function of the generalized coordinates to a uniform vertical
  gravitational potential."
  ([masses g]
   (U-uniform-gravity masses g 1))
  ([masses g vertical-coordinate]
   (let [vertical #(nth % vertical-coordinate)]
     (fn [q]
       (* g masses (mapv vertical q))))))

(defn attach-pendulum
  "Replaces an angle at the ith index with a pendulum. The supplied attachment
  can be either a:

  - map of the form {:coordinate idx}, in which case the pendulum will attach
    there
  - A 2-tuple, representing a 2d attachment point

  If the attachment index doesn't exist, attaches to the origin.
  "
  [l idx attachment]
  (fn [[_ q]]
    (let [v (structure->vector q)
          [x y] (if (vector? attachment)
                  attachment
                  (get v (:coordinate attachment) [0 0]))
          update (fn [angle]
                   (up (+ x (* l (sin angle)))
                       (- y (* l (cos angle)))))]
      (vector->up
       (update-in v [idx] update)))))

(defn double-pendulum->rect
  "Convert to rectangular coordinates from the angles. This is a simpler, explicit
  version of a transformation that attaches a double pendulum."
  [l1 l2]
  (fn [[_ [theta phi]]]
    (let [x1 (* l1 (sin theta))
          y1 (* -1 l1 (cos theta))]
      (up (up x1 y1)
          (up (+ x1 (* l2 (sin phi)))
              (-  y1 (* l2 (cos phi))))))))


(defn L-chain
  "Lagrangian for a chain of pendulums under the influence of a uniform
  gravitational pull in the -y direction."
  [m lengths g]
  (let [U (U-uniform-gravity m g)]
    (reduce transform
            (init (L-rectangular m U))
            (map-indexed
             (fn [i l]
               (attach-pendulum l i (if (zero? i)
                                      [0 0]
                                      {:coordinate (dec i)})))
             lengths))))

;; Drawing API

(defn attach
  "Attach two points with a line."
  [[x1 y1] [x2 y2]]
  (q/line x1 (- y1) x2 (- y2)))

(defn bob
  "Generate a bob at that point."
  [[x y]]
  (q/ellipse x (- y) 2 2))

(defn draw-chain [convert]
  (fn [{:keys [state color]}]
    ;; Clear the sketch by filling it with light-grey color.
    (q/background 100)

    ;; Set a fill color to use for subsequent shapes.
    (q/fill color 255 255)

    ;; Calculate x and y coordinates of the circle.
    (let [[b1 b2 b3] (convert state)]
      ;; Move origin point to the center of the sketch.
      (q/with-translation [(/ (q/width) 2)
                           (/ (q/height) 2)]
        (attach [0 0] b1)
        (attach b1 b2)
        (attach b2 b3)
        (bob b1)
        (bob b2)
        (bob b3)))))

;; Notes on where I'm leaving this. I really want to get this so I can actually
;; stack and reuse the transformations, as that's what's so interesting here.
;;
;; Right now, there is still some duplication with how the transformations get
;; stacked internally. We also have duplication on the initial state, and how
;; many components we have there.
;;
;; That needs to get cleaned up, AND we can switch to the interface above where
;; we

(let [g 9.8
      m (down 1 1 1)
      lengths [4 8 12]
      L (L-chain m lengths g)
      initial-state (up 0
                        (up (/ pi 2) (/ pi 2) (/ pi 2))
                        (up 0 0 0))
      built (build L initial-state)]
  (q/defsketch triple-pendulum
    :title "Triple pendulum"
    :size [500 500]
    :setup (:setup built)
    :update (:update built)
    :draw (draw-chain (:xform built))
    :features [:keep-on-top]
    :middleware [m/fun-mode m/navigation-2d]))

;; Remaining intro stuff

(defn draw-particle [convert]
  (fn [{:keys [state color]}]
    ;; Clear the sketch by filling it with light-grey color.
    (q/background 100)

    ;; Set a fill color to use for subsequent shapes.
    (q/fill color 255 255)

    ;; Calculate x and y coordinates of the circle.
    (let [[[x y]] (convert state)]
      ;; Move origin point to the center of the sketch.
      (q/with-translation [(/ (q/width) 2)
                           (/ (q/height) 2)]
        ;; Draw the circle.
        (q/ellipse x (- y) 5 5)))))

(let [m 1
      g 9.8
      L (init (L-rectangular (down m) (fn [q] 0)))
      initial-state (up 0 (up (up 5 5)) (up (up 4 10)))
      built (build L initial-state)]
  #_(q/defsketch uniform-particle
      :title "Particle in uniform gravity"
      :size [500 500]
      ;; setup function called only once, during sketch initialization.
      :setup (:setup built)
      :update (:update built)
      :draw (draw-particle (:xform built))
      :features [:keep-on-top]
      :middleware [m/fun-mode m/navigation-2d]))


;; (let [m 1
;;       k 9.8
;;       initial-state (up 0 (up 5 5) (up 4 10))]
;;   (q/defsketch harmonic-oscillator-sketch
;;     :title "Harmonic oscillator"
;;     :size [500 500]
;;     ;; setup function called only once, during sketch initialization.
;;     :setup (setup-fn initial-state)
;;     :update (Lagrangian-updater (L-harmonic m k) initial-state)
;;     :draw draw-state
;;     :features [:keep-on-top]
;;     :middleware [m/fun-mode m/navigation-2d]))

;; ;; # Driven Pendulum

;; (defn L-pend-rect
;;   "Single particle to start, under some potential. This is finicky because we have
;;   to ignore the y coordinate here.

;;   It WOULD be cool to color the rod based on its stress."
;;   [m g]
;;   (fn [[_ [_ y] [xdot ydot]]]
;;     (- (* 1/2 m (+ (square xdot)
;;                    (square ydot)))
;;        (* m g y))))

;; (defn driven-pend->rect
;;   "Convert to rectangular coordinates from a single angle."
;;   [l yfn]
;;   (fn [[t [theta]]]
;;     (up (* l (sin theta))
;;         (- (yfn t)
;;            (* l (cos theta))))))

;; (defn L-driven
;;   "Lagrangian for the double pendulum."
;;   [m l g yfn]
;;   (compose (L-pend-rect m g)
;;            (F->C (driven-pend->rect l yfn))))

;; (defn draw-driven [l yfn]
;;   (let [convert (driven-pend->rect l yfn)]
;;     (fn [{:keys [state color time tick]}]
;;       ;; Clear the sketch by filling it with light-grey color.
;;       (q/background 100)

;;       ;; Set a fill color to use for subsequent shapes.
;;       (q/fill color 255 255)

;;       ;; Calculate x and y coordinates of the circle.
;;       (let [[x y] (convert state)]
;;         ;; Move origin point to the center of the sketch.
;;         (q/with-translation [(/ (q/width) 2)
;;                              (/ (q/height) 2)]
;;           (q/line 0 (- (yfn (state->t state))) x (- y))
;;           (q/ellipse x (- y) 2 2))))))

;; (let [m 1
;;       l 6
;;       g 9.8
;;       yfn (fn [t]
;;             (* 2 (cos (* t (sqrt (/ g l))))))
;;       L (L-driven m l g yfn)
;;       initial-state (up 0
;;                         (up (/ pi 4))
;;                         (up 0))]
;;   (q/defsketch driven-pendulum
;;     :title "Driven pendulum"
;;     :size [500 500]
;;     ;; setup function called only once, during sketch initialization.
;;     :setup (setup-fn initial-state)
;;     :update (Lagrangian-updater L initial-state)
;;     :draw (draw-driven l yfn)
;;     :features [:keep-on-top]
;;     :middleware [m/fun-mode m/navigation-2d]))
