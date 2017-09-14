# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------
# Copyright 2017 Jürgen Probst
#
# This file is part of pyMPB.
#
# pyMPB is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMPB is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyMPB.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------


# the main .ctl template:

template = """%(initcode)s

(set! geometry-lattice %(lattice)s)

(set! resolution %(resolution)s)

(set! mesh-size %(meshsize)s)

(set! num-bands %(numbands)s)

(set! k-points %(kspace)s)

(set! geometry (list %(geometry)s))

%(runcode)s

%(postcode)s

(display-eigensolver-stats)"""


# template used in simulations that use a function of one variable for
# epsilon:
template_epsilon_function = r"""
; ######################################################
; ########## epsilon function of one variable ##########
; ######################################################

; the knots of the cubic spline describing the epsilon function:
(define epsknots (list %(epsknots)s))

; the coefficients of the cubic spline describing the epsilon function:
(define epscoeffs #2(%(epscoeffs)s))

; a node for a segment tree realized by a simple list:
(define (segtree-node boundary left right) (list boundary left right))

(define (build-segment-tree xlist . ilist)
    ;Return the root segtree-node to a segment tree of a partitioned interval.
    ;Built using xlist as list of segment boundaries.
    ;The branches are segtree-nodes, the leaves are indexes (optionally
    ;given by ilist, must have same length than xlist).
    (let ( (count (length xlist)) )
        (if (null? ilist)
            ; initialize list of indexes, usually on first call:
            (set! ilist (arith-sequence 0 1 count))
        )
        (case count
            ((0) 0) ; should never get here
            ((1) ; arrived at a single interval, i.e. a leave;
                 ; Don't return a segtree-node, just the index of the interval:
                (car ilist)
            )
            (else ; split the list of interval boundaries in half:
                (let ( (halfway (quotient count 2)) )
                    (segtree-node
                        (list-ref xlist halfway)
                        (apply build-segment-tree
                            (list-head xlist halfway)
                            (list-head ilist halfway))
                        (apply build-segment-tree
                            (list-tail xlist halfway)
                            (list-tail ilist halfway))
                    )
                ) ; let
            ) ; else
        ) ; case
    ) ; let
) ; define (build-segment-tree )

(define (find-index x segtree)
    ;Return the index of the segment where x is located.
    ;segtree is the root of a segment tree of segtree-nodes.
    ;The segments are interpreted this way:
    ;(note that x0 is ignored)
    ;(-inf, x1) [x1 x2) [x2 x3) [x3 x4) [x4 +inf)
    (if (list? segtree)
        (if (< x (car segtree))
            (find-index x (cadr segtree))
            (find-index x (caddr segtree))
        )
        segtree
    )
)

; Build the segment tree from epsknots:
; It is important to drop the last limit, since the last segment goes to +inf.
; Don't drop the first limit, though.
(define segtree-root
    (build-segment-tree (list-head epsknots (- (length epsknots) 1))) )

; the epsilon function of one variable:
(define (epsfunc v)
    (let (
            ; find the segment index where v is located:
            (ind (find-index v segtree-root))
            ; third degree polygon:
            (k 3) )

        (if (or (< v (car epsknots))
                (> v (list-ref epsknots (- (length epsknots) 1))) )
            (print
                "sim-info: WARNING: extrapolated epsilon for " v "\n")
        )
        ; evaluate and return the polygon:
        (apply + (map
            (lambda (m)
                (*
                    (array-ref epscoeffs m ind)
                    (expt ; power
                        (- v (list-ref epsknots ind))
                        (- k m)
                    )
                )
            ) ; lambda
            (arith-sequence 0 1 (+ k 1))
        )); map ; apply
    )
)

; #######################################################
"""

# template for runcode used for simulations that use a
# frequency dependent epsilon:
template_runcode_freq_dependent_epsilon = r"""

; #######################################################
; ################### run functions #####################
; #######################################################

; Run an MPB simulation at kvec, with default-material eps. Return the
; frequency of band bandnum:
(define (simulate-at-eps eps kvec bandnum p reset-fields?)
    (print "\nsim-info: Setting default-material epsilon: " eps "\n")
    (set! default-material (make dielectric (epsilon eps)) )
    (init-params p reset-fields?)
    (set! current-k kvec)
    (solve-kpoint kvec)
    (print "sim-info: frequency: " (list-ref freqs (- bandnum 1)) "\n")
    (list-ref freqs (- bandnum 1))
)


(define (run-sim-%(mode_lower)s
            eps-func kvec bandnum init-freq tolerance p . band-functions)
    (let (  (left-f init-freq)
            (right-f init-freq)
            ; frequency steps to determine initial bracketing:
            (init-freq-steps (* 0.01 init-freq))
            (kpoint-index (+ (list-index k-points kvec) 1))
            (temp "")
            (result 0)
            ; At a single given bandnum n and kvec k, the MPB simulation
            ; returns a frequency at a given effective eps, f_kn(eps).
            ; But, eps is a function of frequency: eps = eps(f^)
            ; So, to get the frequency at the proper eps, we must find the
            ; root of the function y = (x - f_kn(eps(x))).
            ; This is this function:
            (ffunc (lambda (f)
                (let ( (eps (eps-func f)) )
                    (- f (simulate-at-eps eps kvec bandnum p false))
                ) ; let
            ))
         )

        (print "sim-info: running simulation at k = " kvec
               " for band " bandnum "\n")

        ; don't be interactive if we call (run-sim-%(mode_lower)s)
        (set! interactive? false)

        ; ##### bracket freq of interest #####
        (cond
            ((< (ffunc init-freq) 0)
                (set! right-f
                    (do ((f (+ init-freq init-freq-steps)
                            (+ f init-freq-steps)))
                        ((> (ffunc f) 0) f)
                        (set! left-f f)
                    ) ;do
                ); set right-f
            )
            (else
                (set! left-f
                    (do ((f (- init-freq init-freq-steps)
                            (- f init-freq-steps)))
                        ((< (ffunc f) 0) f)
                        (set! right-f f)
                    ) ;do
                ); set left-f
            )
        ) ; cond

        ;(print "sim-debug: init: " init-freq " -> " (ffunc init-freq) "\n")
        ;(print "sim-debug: left: " left-f " -> " (ffunc left-f) "\n")
        ;(print "sim-debug: right: " right-f " -> " (ffunc right-f) "\n")
        (print "sim-info: bracketing end\n")


        ; ########### find root ###########
        (set! result (find-root ffunc tolerance left-f right-f))

        ; run simulation one last time with result, which was not run
        ; in find-root before:
        (set! result (simulate-at-eps (eps-func result) kvec bandnum p false))

        ; set output variables:
        (set! all-freqs (cons freqs all-freqs))
        (set! band-range-data
            (update-band-range-data band-range-data freqs kvec))
        (set! eigensolver-iters
            (append eigensolver-iters
                (list (/ iterations num-bands))))

        ; print results:

        ; if this is the first k point, print out a header line for
        ; the frequency grep data:
        (if (eqv? (car k-points) kvec)
            (begin
                (print "sim-" parity "freqs:, k index, k1, k2, k3, kmag/2pi")
                (do ( (i 0 (+ i 1)) )
                    (( = i num-bands) (print "\n"))
                    (print ", band " (+ i 1)) ))
        ) ; if
        (print "sim-" parity "freqs:, " kpoint-index ", "
            (vector3-x kvec) ", " (vector3-y kvec) ", " (vector3-z kvec) ", "
            (vector3-norm kvec))
        (do ( (i 0 (+ i 1)) )
            (( = i num-bands) (print "\n"))
            (print ", " (list-ref freqs i)) )

        ; if this is the first k point, print out a header line for
        ; the sim-result grep data:
        (if (eqv? (car k-points) kvec)
            (begin
                (print "sim-" parity "result:, k index, k1, k2, k3, kmag/2pi, "
                       "band number, frequency, y-parity, z-parity, "
                       "velocity-x, velocity-y, velocity-z\n"))
        ) ; if
        (print "sim-" parity "result:, " kpoint-index ", "
            (vector3-x kvec) ", " (vector3-y kvec) ", " (vector3-z kvec) ", "
            (vector3-norm kvec) ", " bandnum ", " result
            ", " (list-ref (compute-yparities) (- bandnum 1))
            ", " (list-ref (compute-zparities) (- bandnum 1))
            ", " (compute-1-group-velocity-component
                  (cartesian->reciprocal (vector3 1 0 0)) bandnum)
            ", " (compute-1-group-velocity-component
                  (cartesian->reciprocal (vector3 0 1 0)) bandnum)
            ", " (compute-1-group-velocity-component
                  (cartesian->reciprocal (vector3 0 0 1)) bandnum) "\n")

        ; As we only calculate one kpoint after we call init-params,
        ; the mpb-internal kpoint-index is always 1. This affects the
        ; output functions, which output to files with the kpoint-index
        ; in their names. So we temporarily correct this index:
        (set-kpoint-index kpoint-index)

        ; We also correct the mode name for the display-functions,
        ; so the output can be grepped like the frequencies with the "sim-"
        ; prefix:
        (set! temp parity)
        (set! parity (string-append "sim-" parity))

        ; run band functions (only for bandnum):
        (map (lambda (fun)
                 (if (zero? (procedure-num-args fun))
                     (fun)
                     (fun bandnum)))
            band-functions)

        ; Change it back so it does not break things:
        (set-kpoint-index 1)
        (set! parity temp)

    ) ; let
) ; define


; #######################################################
; ################## run simulation #####################
; #######################################################

; output epsilon (with approx. background):
(set! default-material (make dielectric (epsilon (epsfunc init-freq))))
(init-params %(mode_upper)s false)
(output-epsilon)

(set! total-run-time (+ total-run-time
    (begin-time "sim-info: total elapsed time for run: "
        (map (lambda (kvec)
            (begin-time "sim-info: elapsed time for k point: "
                %(preparation)s
                (run-sim-%(mode_lower)s
                    epsfunc kvec bandnum init-freq 1e-4 %(mode_upper)s
                    %(bandfuncs)s
                )
                ; update init-freq for the next kvec
                (set! init-freq (list-ref freqs (- bandnum 1)) )
            )
        ) k-points)
    )))

; put all-freqs in the right order:
(set! all-freqs (reverse all-freqs))

(print-dos 0 1.2 121 "sim-")

(print "done.\n")
"""


# template for runcode used for simulations that use a
# wavenumber k dependent epsilon:
template_runcode_k_dependent_epsilon = r"""

; #######################################################
; ################### run functions #####################
; #######################################################

(define (run-sim-%(mode_lower)s
            eps-func kvec init-bandnum bandnumfunc p . band-functions)
    (let (
            (kpoint-index (+ (list-index k-points kvec) 1))
            (temp "")
            (bandnum 0)
            (eps (eps-func (vector3-norm kvec)) )
         )

        (print "sim-info: running simulation at k = " kvec "\n")

        ; don't be interactive if we call (run-sim-%(mode_lower)s)
        (set! interactive? false)

        ; prepare simulation at kvec:
        (print "\nsim-info: Setting default-material epsilon: " eps "\n")
        (set! default-material (make dielectric (epsilon eps)) )
        (init-params p false)
        (set! current-k kvec)

        ; run simulation:
        (solve-kpoint kvec)

        ; set output variables:
        (set! all-freqs (cons freqs all-freqs))
        (set! band-range-data
            (update-band-range-data band-range-data freqs kvec))
        (set! eigensolver-iters
            (append eigensolver-iters
                (list (/ iterations num-bands))))

        ; now we have access to the simulation results, including frequencies,
        ; y-parities and group velocities, to determine the correct band number:
        (set! bandnum (bandnumfunc init-bandnum))

        ; print results:

        ; if this is the first k point, print out a header line for
        ; the frequency grep data:
        (if (eqv? (car k-points) kvec)
            (begin
                (print "sim-" parity "freqs:, k index, k1, k2, k3, kmag/2pi")
                (do ( (i 0 (+ i 1)) )
                    (( = i num-bands) (print "\n"))
                    (print ", band " (+ i 1)) ))
        ) ; if
        (print "sim-" parity "freqs:, " kpoint-index ", "
            (vector3-x kvec) ", " (vector3-y kvec) ", " (vector3-z kvec) ", "
            (vector3-norm kvec))
        (do ( (i 0 (+ i 1)) )
            (( = i num-bands) (print "\n"))
            (print ", " (list-ref freqs i)) )

        ; if this is the first k point, print out a header line for
        ; the sim-result grep data:
        (if (eqv? (car k-points) kvec)
            (begin
                (print "sim-" parity "result:, k index, k1, k2, k3, kmag/2pi, "
                       "band number, frequency, y-parity, z-parity, "
                       "velocity-x, velocity-y, velocity-z\n"))
        ) ; if
        (print "sim-" parity "result:, " kpoint-index ", "
            (vector3-x kvec) ", " (vector3-y kvec) ", " (vector3-z kvec) ", "
            (vector3-norm kvec) ", " bandnum ", "
            (list-ref freqs (- bandnum 1))
            ", " (list-ref (compute-yparities) (- bandnum 1))
            ", " (list-ref (compute-zparities) (- bandnum 1))
            ", " (compute-1-group-velocity-component
                  (cartesian->reciprocal (vector3 1 0 0)) bandnum)
            ", " (compute-1-group-velocity-component
                  (cartesian->reciprocal (vector3 0 1 0)) bandnum)
            ", " (compute-1-group-velocity-component
                  (cartesian->reciprocal (vector3 0 0 1)) bandnum) "\n")

        ; As we only calculate one kpoint after we call init-params,
        ; the mpb-internal kpoint-index is always 1. This affects the
        ; output functions, which output to files with the kpoint-index
        ; in their names. So we temporarily correct this index:
        (set-kpoint-index kpoint-index)

        ; We also correct the mode name for the display-functions,
        ; so the output can be grepped like the frequencies with the "sim-"
        ; prefix:
        (set! temp parity)
        (set! parity (string-append "sim-" parity))

        ; run band functions (only for bandnum):
        (map (lambda (fun)
                 (if (zero? (procedure-num-args fun))
                     (fun)
                     (fun bandnum)))
            band-functions)

        ; Change it back so it does not break things:
        (set-kpoint-index 1)
        (set! parity temp)

    ) ; let
) ; define


; #######################################################
; ################## run simulation #####################
; #######################################################

; output epsilon (with background at k=0):
(set! default-material (make dielectric (epsilon (epsfunc 0))))
(init-params %(mode_upper)s false)
(output-epsilon)

(set! total-run-time (+ total-run-time
    (begin-time "sim-info: total elapsed time for run: "
        (map (lambda (kvec)
            (begin-time "sim-info: elapsed time for k point: "
                %(preparation)s
                (run-sim-%(mode_lower)s
                    epsfunc kvec bandnum %(bandnumfunc)s %(mode_upper)s
                    %(bandfuncs)s
                )
            )
        ) k-points)
    )))

; put all-freqs in the right order:
(set! all-freqs (reverse all-freqs))

(print-dos 0 1.2 121 "sim-")

(print "done.\n")
"""

template_y_odd_bandnum = r"""
(define (get-bandnum-for-y-odd-parity init-bandnum)
    (let ( (pars (compute-yparities))
           (velos (compute-group-velocity-component
                      (cartesian->reciprocal (vector3 1 0 0))) )
           (last-yodd 0)
           (bi (- init-bandnum 1))
         )
        (if (< (list-ref pars bi) -0.5)
            ; init-bandnum is already y-odd, return it:
            init-bandnum
            ; look for first y-odd band with negative velocity:
            (begin
                (while (or (> (list-ref pars bi) 0) (> (list-ref velos bi) 0) )
                    (if (< (list-ref pars bi) -0.5) (set! last-yodd (+ bi 1)))
                    (set! bi (+ bi 1))
                    (if (>= bi num-bands) (break))
                ) ; while
                (if (= bi num-bands)
                   ; reached last band, return last y-odd band encountered:
                   last-yodd
                   ; otherwise:
                   ; directly following bi, more y-odd bands might have
                   ; negative velocities; find the lowest thereof:
                   (let ( (bi2 bi) )
                      (while (or (> (list-ref pars bi2) 0)
                                 (<= (list-ref velos bi2) (list-ref velos bi)))
                          (if (< (list-ref pars bi2) -0.5) (set! bi bi2))
                          (set! bi2 (+ bi2 1))
                          (if (>= bi2 num-bands) (break))
                      ) ; while
                      (+ bi 1)
                   ) ; let
                ); if
            ); begin
        ) ; if
    ) ; let
) ; define


"""

