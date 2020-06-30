(ns maturin.better
  (:refer-clojure :exclude [partial zero? + - * / ref])
  (:require [maturin.ode :as o]
            [sicmutils.env :refer :all]
            [sicmutils.structure :as struct]
            [quil.core :as q]
            [quil.middleware :as m]))

;; # Double Pendulum

(defn timestamp
  "TODO update this to a cljc compatible catch."
  []
  (try (/ (q/millis) 1000.0)
       (catch Exception _ 0)))

(defn setup-fn
  "I THINK this can stay the same for any of the Lagrangians."
  [initial-state]
  (fn []
    ;; Set frame rate to 30 frames per second.
    (q/frame-rate 60)
    ;; Set color mode to HSB (HSV) instead of default RGB.
    (q/color-mode :hsb)

    ;; setup function returns initial state. It contains
    ;; circle color and position.
    {:state initial-state
     :time (timestamp)
     :color 0
     :tick 0
     :navigation-2d {:zoom 4}}))

(defn Lagrangian-updater
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

(defn draw-state [{:keys [state color time tick]}]
  ;; Clear the sketch by filling it with light-grey color.
  (q/background 100)

  ;; Set a fill color to use for subsequent shapes.
  (q/fill color 255 255)

  ;; Calculate x and y coordinates of the circle.
  (let [[x y] (coordinate state)]
    ;; Move origin point to the center of the sketch.
    (q/with-translation [(/ (q/width) 2)
                         (/ (q/height) 2)]
      ;; Draw the circle.
      (q/ellipse x (- y) 5 5))))

;; # Attempt at an API

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
    (-> m
        (update-in :transforms assoc xform)
        (update-in :lagrangian compose xform))))

(defn build
  "Returns the final keys for the sketch."
  [{:keys [lagrangian transforms]} initial-state]
  {:setup (setup-fn initial-state)
   :update (Lagrangian-updater lagrangian initial-state)})

(let [m 1
      g 9.8
      initial-state (up 0 (up 5 5) (up 4 10))]
  (q/defsketch uniform-particle
    :title "Particle in uniform gravity"
    :size [500 500]
    ;; setup function called only once, during sketch initialization.
    :setup (setup-fn initial-state)
    :update (Lagrangian-updater (L-particle m g) initial-state)
    :draw draw-state
    :features [:keep-on-top]
    :middleware [m/fun-mode m/navigation-2d]))


(let [m 1
      k 9.8
      initial-state (up 0 (up 5 5) (up 4 10))]
  (q/defsketch harmonic-oscillator-sketch
    :title "Harmonic oscillator"
    :size [500 500]
    ;; setup function called only once, during sketch initialization.
    :setup (setup-fn initial-state)
    :update (Lagrangian-updater (L-harmonic m k) initial-state)
    :draw draw-state
    :features [:keep-on-top]
    :middleware [m/fun-mode m/navigation-2d]))

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

(defn attach-pendulum [l idx]
  (fn [[_ q]]
    (let [v (structure->vector q)
          [x y] (get v (dec idx) [0 0])
          update (fn [angle]
                   (up (+ x (* l (sin angle)))
                       (- y (* l (cos angle)))))]
      (vector->up
       (update-in v [idx] update)))))

(defn double-pendulum->rect
  "Convert to rectangular coordinates from the angles."
  [l1 l2]
  (fn [[_ [theta phi]]]
    (let [x1 (* l1 (sin theta))
          y1 (* l1 (cos theta))]
      (up (up x1 (- y1))
          (up (+ x1 (* l2 (sin phi)))
              (-  y1 (* l2 (cos phi))))))))

(defn L-double-pendulum
  "Lagrangian for the double pendulum.

  TODO: Note that you can SEPARATELY compose the potential with various
  transformations."
  [m1 m2 l1 l2 g]
  (let [m (down m1 m2 m1 m2)
        U (U-uniform-gravity m g)]
    (compose (L-rectangular m U)
             (F->C (attach-pendulum l1 3))
             (F->C (attach-pendulum l1 2))
             (F->C (attach-pendulum l2 1))
             (F->C (attach-pendulum l1 0)))))

(defn attach [[x1 y1] [x2 y2]]
  (q/line x1 (- y1) x2 (- y2)))

(defn bob [[x y]]
  (q/ellipse x (- y) 2 2))

(defn draw-double [convert]
  (fn [{:keys [state color]}]
    ;; Clear the sketch by filling it with light-grey color.
    (q/background 100)

    ;; Set a fill color to use for subsequent shapes.
    (q/fill color 255 255)

    ;; Calculate x and y coordinates of the circle.
    (let [[b1 b2 b3 b4] (convert state)]
      ;; Move origin point to the center of the sketch.
      (q/with-translation [(/ (q/width) 2)
                           (/ (q/height) 2)]
        (attach [0 0] b1)
        (attach b1 b2)
        (attach b2 b3)
        (attach b3 b4)
        (bob b1)
        (bob b2)
        (bob b3)
        (bob b4)
        ))))

(let [m1 1 m2 1
      l1 6 l2 6
      g 9.8
      L (L-double-pendulum m1 m2 l1 l2 g)
      initial-state (up 0
                        (up (/ pi 2) (/ pi 2) (/ pi 2) (/ pi 2))
                        (up 0 0 0 0))]
  (q/defsketch double-pendulum
    :title "Double pendulum"
    :size [500 500]
    ;; setup function called only once, during sketch initialization.
    :setup (setup-fn initial-state)
    :update (Lagrangian-updater L initial-state)
    :draw (draw-double (compose coordinate
                                (F->C (attach-pendulum l1 3))
                                (F->C (attach-pendulum l1 2))
                                (F->C (attach-pendulum l2 1))
                                (F->C (attach-pendulum l1 0))))
    :features [:keep-on-top]
    :middleware [m/fun-mode m/navigation-2d]))

;; Driven Pendulum

(defn L-pend-rect
  "Single particle to start, under some potential. This is finicky because we have
  to ignore the y coordinate here.

  It WOULD be cool to color the rod based on its stress."
  [m g]
  (fn [[_ [_ y] [xdot ydot]]]
    (- (* 1/2 m (+ (square xdot)
                   (square ydot)))
       (* m g y))))

(defn driven-pend->rect
  "Convert to rectangular coordinates from a single angle."
  [l yfn]
  (fn [[t [theta]]]
    (up (* l (sin theta))
        (- (yfn t)
           (* l (cos theta))))))

(defn L-driven
  "Lagrangian for the double pendulum."
  [m l g yfn]
  (compose (L-pend-rect m g)
           (F->C (driven-pend->rect l yfn))))

(defn draw-driven [l yfn]
  (let [convert (driven-pend->rect l yfn)]
    (fn [{:keys [state color time tick]}]
      ;; Clear the sketch by filling it with light-grey color.
      (q/background 100)

      ;; Set a fill color to use for subsequent shapes.
      (q/fill color 255 255)

      ;; Calculate x and y coordinates of the circle.
      (let [[x y] (convert state)]
        ;; Move origin point to the center of the sketch.
        (q/with-translation [(/ (q/width) 2)
                             (/ (q/height) 2)]
          (q/line 0 (- (yfn (state->t state))) x (- y))
          (q/ellipse x (- y) 2 2))))))

(let [m 1
      l 6
      g 9.8
      yfn (fn [t]
            (* 2 (cos (* t (sqrt (/ g l))))))
      L (L-driven m l g yfn)
      initial-state (up 0
                        (up (/ pi 4))
                        (up 0))]
  (q/defsketch driven-pendulum
    :title "Driven pendulum"
    :size [500 500]
    ;; setup function called only once, during sketch initialization.
    :setup (setup-fn initial-state)
    :update (Lagrangian-updater L initial-state)

    :draw (draw-driven l yfn)
    :features [:keep-on-top]
    :middleware [m/fun-mode m/navigation-2d]))
