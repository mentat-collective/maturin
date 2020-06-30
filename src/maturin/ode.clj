(ns maturin.ode
  (:require [clojure.tools.logging :as log]
            [sicmutils
             [structure :as struct]
             [simplify]]
            [sicmutils.numerical.compile :refer :all])
  (:import (org.apache.commons.math3.ode.nonstiff GraggBulirschStoerIntegrator.)
           (org.apache.commons.math3.ode FirstOrderDifferentialEquations)
           (org.apache.commons.math3.ode.sampling StepHandler)
           (com.google.common.base Stopwatch)))

(defn attach-handler
  " We implement the observation callback by adding a StepHandler to the
  integration. The StepHandler is not invoked at every grid point; rather, it is
  invoked once in a while over a range of time within which the integrated
  function may be accurately evaluated. The handler we install does this,
  invoking the callback for each requested grid point within the valid range,
  ensuring that we also invoke the callback for the final point."
  [integrator observe step-size initial-state]
  (let [array->state #(struct/unflatten % initial-state)]
    (.addStepHandler
     integrator
     (reify StepHandler
       (handleStep
           [_ interpolator is-last]
         (let [it0 (.getPreviousTime interpolator)
               it1 (.getCurrentTime interpolator)
               adjust (mod it0 step-size)
               t0 (if (> adjust 0) (+ (- it0 adjust) step-size) it0)
               last-state (when is-last (double-array (.getInterpolatedState interpolator)))]
           (doseq [t (range t0 it1 step-size)]
             (.setInterpolatedTime interpolator t)
             (observe t (-> interpolator .getInterpolatedState array->state)))
           (when is-last
             (observe it1 (array->state last-state)))))
       (init [_ _ _ _])))))

(defn little-integrator
  "Trimming down the big guy. Let's see where I get.

  TODO get the equations out of here. We can use this with different
  integrations. The equations are what matters."
  [state-derivative derivative-args]
  (fn [initial-state & {:keys [compile epsilon]}]
    (let [state->array #(-> % flatten double-array)
          array->state #(struct/unflatten % initial-state)
          initial-state-array (doubles (state->array initial-state))
          dimension (alength initial-state-array)
          eps (double epsilon)

          derivative-fn (if compile
                          (compile-state-function state-derivative derivative-args initial-state)
                          (do (log/warn "Not compiling function for ODE analysis")
                              (let [d:dt (apply state-derivative derivative-args)]
                                #(-> % array->state d:dt))))
          ;; This uses the derivative fn to calculate the next step. I don't get
          ;; why we're concatenating the derivative args back on.
          equations (reify FirstOrderDifferentialEquations
                      (computeDerivatives [_ _ y out]
                        (let [y' (doubles (-> (concat y derivative-args)
                                              derivative-fn
                                              state->array))]
                          (System/arraycopy y' 0 out 0 (alength y'))))
                      (getDimension [_] dimension))
          integrator (GraggBulirschStoerIntegrator. 0. 1. eps eps)]
      {:integrator integrator
       :equations equations
       :dimension dimension})))

(defn make-integrator
  "make-integrator takes a state derivative function (which in this
  system is assumed to be a map from a structure to a structure of the
  same shape, as differentiating a function does not change its
  shape), and returns an integrator, which is a function of several
  arguments: the initial state, an intermediate-state observation
  function, the step size desired, the final time to seek, and an
  error tolerance. If the observe function is not nil, it will be
  invoked with the time as first argument and integrated state as the
  second, at each intermediate step."
  [state-derivative derivative-args]
  (fn [initial-state observe step-size t ε & {:keys [compile]}]
    (let [total-time (Stopwatch/createStarted)
          evaluation-count (atom 0)
          evaluation-time (Stopwatch/createUnstarted)
          state->array #(-> % flatten double-array)
          array->state #(struct/unflatten % initial-state)
          initial-state-array (doubles (state->array initial-state))

          derivative-fn (if compile
                          (compile-state-function state-derivative derivative-args initial-state)
                          (let [d:dt (apply state-derivative derivative-args)]
                            #(-> % array->state d:dt)))

          dimension (alength initial-state-array)

          integrator (GraggBulirschStoerIntegrator. 0. 1. (double ε) (double ε))

          ;; This uses the derivative fn to calculate the next step. I don't get
          ;; why we're concatenating the derivative args back on.
          equations (reify FirstOrderDifferentialEquations
                      (computeDerivatives [_ _ y out]
                        (.start evaluation-time)
                        (swap! evaluation-count inc)
                        (let [y' (doubles (-> (concat y derivative-args)
                                              derivative-fn
                                              state->array))]
                          (System/arraycopy y' 0 out 0 (alength y')))
                        (.stop evaluation-time))
                      (getDimension [_] dimension))
          out (double-array dimension)]

      (when-not compile
        (log/warn "Not compiling function for ODE analysis"))
      (when observe
        ;; We implement the observation callback by adding a StepHandler
        ;; to the integration. The StepHandler is not invoked at every grid
        ;; point; rather, it is invoked once in a while over a range of time
        ;; within which the integrated function may be accurately evaluated.
        ;; The handler we install does this, invoking the callback for
        ;; each requested grid point within the valid range, ensuring that we
        ;; also invoke the callback for the final point.
        (.addStepHandler
         integrator
         (reify StepHandler
           (handleStep
               [_ interpolator is-last]
             (let [it0 (.getPreviousTime interpolator)
                   it1 (.getCurrentTime interpolator)
                   adjust (mod it0 step-size)
                   t0 (if (> adjust 0) (+ (- it0 adjust) step-size) it0)
                   last-state (when is-last (double-array (.getInterpolatedState interpolator)))]
               (doseq [t (range t0 it1 step-size)]
                 (.setInterpolatedTime interpolator t)
                 (observe t (-> interpolator .getInterpolatedState array->state)))
               (when is-last
                 (observe it1 (array->state last-state)))))
           (init [_ _ _ _]))))



      (.integrate integrator equations 0 initial-state-array t out)
      (log/info "#" @evaluation-count "total" (str total-time) "f" (str evaluation-time))
      (array->state out))))

(defn state-advancer
  "state-advancer takes a state derivative function constructor
  followed by the arguments to construct it with. The state derivative
  function is constructed and an integrator is produced which takes
  the initial state, target time, and error tolerance as
  arguments. The final state is returned. The state derivative is
  expected to map a structure to a structure of the same shape,
  and is required to have the time parameter as the first element."
  [state-derivative & state-derivative-args]
  (let [I (make-integrator state-derivative state-derivative-args)]
    (fn [initial-state t ε & options]
      (apply I initial-state nil 0 t ε options))))

(defn evolve
  "evolve takes a state derivative function constructor and its
  arguments, and returns an integrator via make-integrator. In
  particular, the returned function accepts a callback function which
  will be invoked at intermediate grid points of the integration."
  [state-derivative & state-derivative-args]
  (make-integrator state-derivative state-derivative-args))

(defn integrate-state-derivative
  "A wrapper for evolve, which is more convenient when you just
  want a vector of (time, state) pairs over the integration interval
  instead of having to deal with a callback. Integrates the supplied
  state derivative (and its argument package) from [0 to t1] in steps
  of size dt"
  [state-derivative state-derivative-args initial-state t1 dt]
  (let [I (make-integrator state-derivative state-derivative-args)
        out (atom [])
        collector (fn [t state]
                    (swap! out conj state))]
    (I initial-state collector dt t1 1e-6 :compile true)
    @out))
