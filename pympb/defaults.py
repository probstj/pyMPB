    #Copyright 2009-2016 Seyed Hessam Moosavi Mehr, Juergen Probst
    #This program is free software; you can redistribute it and/or modify
    #it under the terms of the GNU General Public License as published by
    #the Free Software Foundation; either version 3 of the License, or
    #(at your option) any later version.

    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.

    #You should have received a copy of the GNU General Public License
    #along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division, print_function
from subprocess import check_output
import re
import numpy as np


#mpb_call = 'mpb'
mpb_call = 'mpirun -np %(num_procs)s mpbi-mpi'

# use -T if we run the simulation with mpb-mpi:
mpbdata_call = ('mpb-data -T -rn%(resolution)s '
                '-x%(number_of_tiles_to_output)s '
                '-y%(number_of_tiles_to_output)s '
                '-o%(output_file)s '
                '%(h5_file)s')
epsh5topng_call_2D = 'h5topng -S3 -Zrcbluered -oepsilon.png %(h5_file)s'
epsh5topng_call_3D = 'h5topng -0z0 -S3 -Zrcbluered -oepsilon.png %(h5_file)s'
epsh5topng_call_3D_cross_sect = ('h5topng -0x0 -S3 -Zrcbluered '
                              '-oepsilonslab.png %(h5_file)s')

fieldh5topng_call_2D = ('h5topng -S3 -Zcdkbluered -C%(eps_file)s '
                        '-o%(output_file)s %(h5_file)s')
fieldh5topng_call_2D_no_ovl = ('h5topng -S3 -Zcdkbluered '
                        '-o%(output_file_no_ovl)s %(h5_file)s')

fieldh5topng_call_3D = ('h5topng -0z0 -S3 -Zcdkbluered -C%(eps_file)s '
                        '-o%(output_file)s %(h5_file)s')
fieldh5topng_call_3D_no_ovl = ('h5topng -0z0 -S3 -Zcdkbluered '
                        '-o%(output_file_no_ovl)s %(h5_file)s')
display_png_call = 'display  %(files)s'

# get mpb version:
mpbversion = 'n/a'
for mpb in ['mbp', 'mpbi', 'mpb-mpi', 'mpbi-mpi']:
    try:
        mpbversionline = check_output(
            [mpb, '--version'], universal_newlines=True)
        # MPB made it hard to check the version. The line even changed
        # in version 1.5. Look for first non-alpha part, this might be
        # what we are looking for:
        try:
            mpbversion = re.search(
                '\s([0-9.]*)[,\s]',
                mpbversionline).groups()[0]
        except AttributeError:
            # did not find anything:
            mpbversion = 'n/a'
        break
    except OSError:
        pass

newmpb = mpbversion >= '1.5'

default_resolution = 32
default_mesh_size = 3
default_numbands = 8
# the number of bands to calculate if calculation is only supposed to be used
# for projection of bands:
num_projected_bands = 4

default_k_interpolation = 3
k_interpolation_function = 'interpolate'
#if newmpb:
#    k_uniform_interpolation_function = 'kinterpolate-uniform'
#else:
#    k_uniform_interpolation_function = 'interpolate'
k_uniform_interpolation_function = 'interpolate'

default_initcode = (
    ';load module for calculating density of states:\n'
    '(define dosmodule (%search-load-path "dosv2.scm"))\n'
    '(if dosmodule\n'
    '    (include dosmodule)\n'
    '    (throw \'error "dos.scm not found"))\n\n'
    ';remove the default filename-prefix:\n'
    ';before MPB 1.5:\n' +
    ('{0[0]}(set! filename-prefix "")\n'
     ';MPB 1.5 and newer:\n'
     '{0[1]}(set! filename-prefix #f)\n\n').format(
        [';', ''] if newmpb else ['', ';'])
)


default_postcode = ''
default_runcode = '(run-te)'

# data descriptors to try and grep in MPB-output after simulation:
# (Simulation.postprocess() will grep for `<mode><dataname>:`, e.g. `tefreqs`
# and save grepped data to the file `<jobname>_<mode><dataname>.csv`)
grep_datanames = ['freqs', 'velocity', 'dos', 'yparity', 'zparity']

number_of_tiles_to_output = 3
# Field patterns transformed to PNG will be placed in subfolders named
# (field_output_folder_prefix + '_' + mode):
field_output_folder_prefix = 'pngs'
# specify wheter the rather big hdf5 files should be kept after they
# were converted to png files:
delete_h5_after_postprocessing = True


def default_band_func(poi, outputfunc):
    """Return a string which will be supplied to (run %s) as a bandfunction.

    poi: k-points of interest, list of 3-tuples.
    outputfunc: mpb outputfunction, e.g. 'output-efield-z'

    """
    return (
        '\n    display-group-velocities'
        '\n    display-zparities display-yparities' +
        ''.join(
            [
                '\n    (output-at-kpoint (vector3 {0}) {1})'.format(
                    ' '.join(str(c) for c in vec),
                    outputfunc
                )
                for vec in poi
            ]
        )
    )

output_funcs_te = ['fix-hfield-phase', 'output-hfield-z']
output_funcs_tm = ['fix-efield-phase', 'output-efield-z']
# these are used for (run) function without specific modes:
output_funcs_other = output_funcs_te + output_funcs_tm

temporary_epsh5 = './temporary_eps.h5'
temporary_h5 = './temporary.h5'
temporary_h5_folder = './patterns~/'

isQuiet = False

log_format = "%(asctime)s %(levelname)s: %(message)s"
log_datefmt = "%d.%m.%Y %H:%M:%S"

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

# template for initcode used for simulations that use a
# frequency dependent epsilon:
template_initcode_epsilon_function = r"""
; #######################################################
; ############ frequency dependend epsilon ##############
; #######################################################

; the knots of the cubic spline describing the
; frequency dependent epsilon function:
(define epsknots (list %(epsknots)s))

; the coefficients of the cubic spline describing the
; frequency dependent epsilon function:
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

; the frequency dependent epsilon function:
(define (epsfunc f)
    (let (
            ; find the segment index where f is located:
            (ind (find-index f segtree-root))
            ; third degree polygon:
            (k 3) )

        (if (or (< f (car epsknots))
                (> f (list-ref epsknots (- (length epsknots) 1))) )
            (print
                "sim-info: WARNING: extrapolated epsilon for frequency "
                f "\n")
        )
        ; evaluate and return the polygon:
        (apply + (map
            (lambda (m)
                (*
                    (array-ref epscoeffs m ind)
                    (expt ; power
                        (- f (list-ref epsknots ind))
                        (- k m)
                    )
                )
            ) ; lambda
            (arith-sequence 0 1 (+ k 1))
        )); map ; apply
    )
)


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

"""

# template for runcode used for simulations that use a
# frequency dependent epsilon:
template_runcode_epsilon_function = r"""
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





#####################################################
###          bandplotter defaults                 ###
#####################################################

fig_size = (12, 9)

draw_bands_formatstr = 'o-'
# keyword arguments for band diagram (matplotlib.lines.Line2D properties,
# see :class:`~matplotlib.lines.Line2D` for details.):
draw_bands_kwargs = {'linewidth' : 2}
hide_band_gap = False;

# default kwargs for the tick labels for the k-vec-axis of band diagrams:
# (will be forwarded to underlying matplotlib.text.Text objects)
xticklabels_kwargs={'rotation':0, 'horizontalalignment':'center'}
# xticklabels_kwargs used when one of the labels strings is longer than
# long_xticklabels_when_longer_than:
long_xticklabels_kwargs={'rotation':45, 'horizontalalignment':'right'}
long_xticklabels_when_longer_than = 12

# Text added to gaps drawn in the band diagrams,
# formatted with default_gaptext.format(gapsize_in_percent):
default_gaptext='gap size: {0:.2f}%'
# for locale-aware formatting e.g.:
#default_gaptext='gap size: {0:.4n}%'

# Minimum gap size (normalized, i.e. \Delta\omega/\omega_{center});
# gaps smaller than this will not be drawn in the band diagram:
min_gapsize = 0.0

default_x_axis_hint = 5 # 5 equally spaced ticks, labeled with k-vector
default_y_axis_label = r'frequency $\omega a/2\pi c$'
default_x_axis_label = 'wave vector {0}'
# the x_axis_label used when showing high symmetry point labels on the k
# axis: Note: I am not entirely satisfied with this title. How do you
# really call it? 'Brillouin zone symmetry points'? 'Wave vector
# direction'? (this last one is good, but we also see the magnitude,
# when e.g. going from Gamma to M etc.) 'Wave vector point in Brillouin
# zone'? (too long)
default_kspace_axis_label = 'wave vector'

default_kvecformatter_format_str = '({0:.2f}  {1:.2f}  {2:.2f})'
# make it locale-aware:
#default_kvecformatter_format_str = '({0:.2n}  {1:.2n}  {2:.2n})'
# other possibilities:
#default_kvecformatter_format_str = r'$\binom{{ {0} }}{{ {1} }}$'
# unfortunately, \stackrel[3]{}{}{} does not work, so it looks bad:
#default_kvecformatter_format_str = \
#    r'$\left(\stackrel{{ {0} }}{{ \stackrel{{ {1} }}{{ {2} }} }}\right)$'
#default_kvecformatter_format_str ='{0}\n{1}\n{2}'

# Show fractions in tick labels of k-axis instead of floating point values:
ticks_fractions = True
# Always show a floating point value if the resulting fraction's denominator
# is greater than:
tick_max_denominator = 1000

# If correct_x_axis is set to True, the bands are plotted versus
# x-values which are equidistant according to the Euclidean distance
# between the k-vectors. That way distortions are avoided which occur
# when plotting versus the k-index.
correct_x_axis = True

color_by_parity_marker_size = 60

add_epsilon_as_inset = False
# The valid location codes are:
#     'upper right'  : 1,
#     'upper left'   : 2,
#     'lower left'   : 3,
#     'lower right'  : 4,
#     'right'        : 5,
#     'center left'  : 6,
#     'center right' : 7,
#     'lower center' : 8,
#     'upper center' : 9,
#     'center'       : 10,
epsilon_inset_location = 4
epsilon_inset_zoom = 0.5
epsilon_inset_transpose = False

# add density of states to bandplot?
add_dos_to_bandplot=False


def default_onclick(event, bandplotter):
    """This is the default function called if the bands are plotted with a
    picker supplied and the user clicks on a vertex in the plot. It then just
    prints some information about the vertex(ices) clicked on to stdout,
    including the mode, the k-vector and -index and the frequency(ies).

    """
    try:
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        xaxisformatter = event.mouseevent.inaxes.xaxis.major.formatter
    except AttributeError:
        return

    print(thisline.get_label() + ' mode(s): ', end='')
    for i in ind:
        kindex = thisline.data[2, i]
        xaxispos = xdata[i]
        freq = ydata[i]
        bandindex = thisline.data[3, i]
        parity = thisline.data[4, i]
        s = 'bandnum={0}, k_index={1:.0f}, k_vec={2}, freq={3}'.format(
            bandindex + 1, kindex, xaxisformatter(xaxispos), freq)
        if np.isfinite(parity):
            s += ', parity={0}'.format(parity)
        print(s + '; ', end='')
    print()

    # Other idea (not implemented): display mode pattern if it was exported
    # or even calculate it if not exported yet, e.g.:
    ## Start interactive mode:
    #plt.ion()
    ## create a new popup figure:
    #fig = plt.figure('mode pattern', figsize=(6, 2))
    #ax = fig.add_subplot(111) #put mode pattern image here

# In the field pattern distribution plot, should the real and imaginary
# parts be on top of each other? Otherwise, they go next to each other:
field_dist_vertical_cmplx_comps=True
field_dist_filetype = 'pdf'


contour_lines = {'colors':'k',
                 'linestyles':['dashed','solid'],
                 'linewidths':1.0}
contour_plain = {'linewidths':1.0}
contour_filled = {}
colorbar_style = {'extend':'both','shrink':0.8}



# uncomment to use locale-aware formatting on the numbers along the y-axis:
#import matplotlib.pyplot as plt
#plt.rcParams['axes.formatter.use_locale'] = True

