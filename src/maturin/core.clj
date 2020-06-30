(ns maturin.core
  (:refer-clojure :exclude [partial zero? + - * / ref])
  (:require [maturin.ode :as o]
            [sicmutils.env :refer :all]
            [sicmutils.structure :as struct]
            [quil.core :as q]
            [quil.middleware :as m])
  (:import (com.google.common.base Stopwatch)))

;; # Double Pendulum

(def tex tex$$)

(def s->infix (compose ->infix simplify))

(defn L-particle
  "Single particle to start, under some potential."
  [m g]
  (fn [[_ [_ y] qdot]]
    (- (* 1/2 m (square qdot))
       (* m g y))))

;; Check the equations of motion.
(s->infix
 ((L-particle 'm 'g)
  (up 't
      (up 'x 'y)
      (up 'xdot 'ydot))))

(defn particle-state-derivative
  "Returns a function that can generate the derivative of the state for us."
  [m g]
  (Lagrangian->state-derivative
   (L-particle m g)))

(def particle-equations-of-motion
  "Equations of motion for a particle under the influence of a uniform
  gravitational field."
  ((particle-state-derivative 'm 'g)
   (up 't
       (up 'x 'y)
       (up 'xdot 'ydot))))

;; # Harmonic Oscillator

(defn L-harmonic
  "Lagrangian for a harmonic oscillator"
  [m k]
  (fn [[_ q qdot]]
    (- (* 1/2 m (square qdot))
       (* 1/2 k (square q)))))

(def harmonic-state-derivative
  (compose Lagrangian->state-derivative
           L-harmonic))

(def harmonic-equations-of-motion
  "Equations of motion for a particle under the influence of a uniform
  gravitational field."
  ((harmonic-state-derivative 'm 'k)
   (up 't
       (up 'x 'y)
       (up 'xdot 'ydot))))

;; # Animation Code

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

(defn L-two-particles
  "I know we can do better here, but start with two particles, explicitly
  flattened out."
  [m1 m2 g]
  (fn [[_ [_ y1 _ y2] [xdot1 ydot1 xdot2 ydot2]]]
    (+
     (- (* 1/2 m1 (+ (square xdot1)
                     (square ydot1)))
        (* m1 g y1))

     (- (* 1/2 m2 (+ (square xdot2)
                     (square ydot2)))
        (* m2 g y2)))))

(defn double-pend->rect
  "Convert to rectangular coordinates from the angles."
  [l1 l2]
  (fn [[_ [theta phi]]]
    (let [x1 (* l1 (sin theta))
          y1 (* l1 (cos theta))]
      (up x1
          y1
          (+ x1 (* l2 (sin phi)))
          (+ y1 (* l2 (cos phi)))))))

(defn L-double
  "Lagrangian for the double pendulum."
  [m1 m2 l1 l2 g]
  (compose (L-two-particles m1 m2 g)
           (F->C (double-pend->rect l1 l2))))

(defn draw-double [l1 l2]
  (let [convert (double-pend->rect l1 l2)]
    (fn [{:keys [state color time tick]}]
      ;; Clear the sketch by filling it with light-grey color.
      (q/background 100)

      ;; Set a fill color to use for subsequent shapes.
      (q/fill color 255 255)

      ;; Calculate x and y coordinates of the circle.
      (let [[x1 y1 x2 y2] (convert state)]
        ;; Move origin point to the center of the sketch.
        (q/with-translation [(/ (q/width) 2)
                             (/ (q/height) 2)]
          (q/line 0 0 x1 (- y1))
          (q/line x1 (- y1) x2 (- y2))
          (q/ellipse x2 (- y2) 2 2)
          (q/ellipse x1 (- y1) 2 2)
          )))))

(let [m1 4 m2 1
      l1 6 l2 6
      g 9.8
      L (L-double m1 m2 l1 l2 g)
      initial-state (up 0
                        (up (/ pi 4) (/ pi 8))
                        (up 0 0))]
  (q/defsketch double-pendulum
    :title "Double pendulum"
    :size [500 500]
    ;; setup function called only once, during sketch initialization.
    :setup (setup-fn initial-state)
    :update (Lagrangian-updater L initial-state)
    :draw (draw-double l1 l2)
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
