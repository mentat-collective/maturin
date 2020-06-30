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

(defn V-double
  "Potential for the double pendulum"
  [m_1 m_2 l_1 l_2 g]
  (fn [[_ [theta phi] _]]
    (let [y_1 (- (* l_1 (cos theta)))
          y_2 (- y_1 (* l_2 (cos phi)))]
      (+ (* m_1 g y_1)
         (* m_2 g y_2)))))

(defn T-double
  "Kinetic energy for the double pendulum"
  [m_1 m_2 l_1 l_2 _]
  (fn [[_ [theta phi] [thetadot phidot]]]
    (let [v1sq (* (square l_1) (square thetadot))
          v2sq (* (square l_2) (square phidot))]
      (+ (* 1/2 m_1 v1sq)
         (* 1/2 m_2
            (+ v1sq
               v2sq
               (* 2 l_1 l_2 thetadot phidot (cos (- theta phi)))))))))

(def L-double
  "The Lagrangian!"
  (- T-double V-double))

;; #

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

(def equations-of-motion
  ((particle-state-derivative 'm 'g)
   (up 't
       (up 'x 'y)
       (up 'xdot 'ydot))))

(def out (atom []))

(defn observe [t [_ [x y] :as v]]
  (swap! out conj [t x y]))

(defn big-test []
  (let [m 1
        g 9.8]
    (s->infix
     ((o/make-integrator particle-state-derivative [m g])
      (up 0 (up 0 0) (up 1 1)) ;initial state
      observe ;callback
      1/60 ; step size
      5 ; initial time
      1e-6 ; allowed absolute and relative error
      :compile true))))

(defn little-test []
  (let [m 1
        g 9.8
        {:keys [integrator equations]}
        ((o/little-integrator particle-state-derivative [m g])
         (up 0 (up 0 0) (up 1 1));initial state
         0
         1e-6
         :compile true)
        ]
    ;; this needs to be equal to the FULL dimension, time included.
    (let [dimension 5
          out (double-array dimension)]
      (fn [state tick]
        (let [s (double-array (flatten state))]
          evaluation-time (Stopwatch/createUnstarted)
          (.start evaluation-time)
          (.integrate integrator equations 0 s tick out)
          (.stop evaluation-time)

          (.reset evaluation-time))
        (struct/unflatten out state)))))

;; # Notes: I think we need to have a framerate, so we render every n steps, AND
;; # some internal rate where we make the simulation proceed.
;;
;; Step one is to get anything working and displaying.
;;
;;

(defn timestamp
  "TODO update this to a cljc compatible catch."
  []
  (try (/ (q/millis) 1000.0)
       (catch Exception _ 0)))

;; # Animation Code
(defn setup-fn [initial-state]
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
     :tick 0}))

(defn grav-particle-updater
  "Returns an updater for the particle affected by particle."
  [m g initial-state]
  (let [{:keys [integrator equations dimension] :as m}
        ((o/little-integrator particle-state-derivative [m g])
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
  (let [[x y] (coordinate state)
        zoom 2]
    (q/scale zoom)
    ;; Move origin point to the center of the sketch.
    (q/with-translation [(/ (q/width) 2 zoom)
                         (/ (q/height) 2 zoom)]
      ;; Draw the circle.
      (q/ellipse x (- y) 5 5))))

(let [m 1
      g 9.8
      initial-state (up 0 (up 0 0) (up 10 10))]
  (q/defsketch lagrange
    :title "You spin my circle right round"
    :size [500 500]
    ;; setup function called only once, during sketch initialization.
    :setup (setup-fn initial-state)
    :update (grav-particle-updater m g initial-state)
    :draw draw-state
    :features [:keep-on-top]
    :middleware [m/fun-mode m/navigation-2d]))
