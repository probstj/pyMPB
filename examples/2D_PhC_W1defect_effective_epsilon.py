# -*- coding:utf-8 -*-
# ----------------------------------------------------------------------
# Copyright 2017 Juergen Probst
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

from __future__ import division

import sys
from os import path

import numpy as np
import matplotlib.pyplot as plt

from pympb.phc_simulations import TriHoles2D_Waveguide_effective_epsilon
from pympb import log, defaults

def main():
    if len(sys.argv) > 1:
        mode=sys.argv[1]
    else:
        mode='sim'

    # monkey patching:
    # (don't output multiple tiles along y)
    # (x and y are switched because of -T)
    defaults.mpbdata_call = (
        'mpb-data -T -rn%(resolution)s '
        '-y%(number_of_tiles_to_output)s '
        '-o%(output_file)s '
        '%(h5_file)s')
    defaults.number_of_tiles_to_output = 5
    # defaults.add_epsilon_as_inset = True
    defaults.output_funcs_te = [
        'fix-hfield-phase', 'output-hfield-z',
        'output-dpwr', 'output-hpwr', 'output-tot-pwr']

    # the knots of the cubic spline describing the
    # frequency dependent epsilon function:
    fknots = np.array(
        [ 0.43606492,  0.43632591,  0.43699937,  0.43826107,  0.43987166,
          0.4405708 ,  0.44359142,  0.4540687 ,  0.46606201,  0.47865032,
          0.49139404,  0.50431241,  0.51692864,  0.53029386,  0.54260207,
          0.55531914,  0.56803985])
    # the coefficients of the cubic spline describing the
    # frequency dependent epsilon function:
    fcoeffs = np.array(
        [[  1.05325162e+07,  -7.53192904e+06,  -5.34491361e+05,
            7.35887805e+06,  -1.78975790e+07,   1.11124611e+06,
            3.50296450e+04,  -7.75123125e+03,   2.52814337e+03,
           -7.88904847e+02,  -1.88305373e+02,  -3.86288083e+02,
            8.09187014e+03,  -2.33155197e+04,   3.20396213e+04,
           -1.69215816e+04],
         [  2.72848411e-11,   8.24663379e+03,  -6.97063083e+03,
           -8.99373723e+03,   2.65628394e+04,  -1.09756189e+04,
           -9.05673628e+02,   1.95372830e+02,  -8.35160080e+01,
            1.19591551e+01,  -1.82015878e+01,  -2.54993856e+01,
           -4.01198768e+01,   2.84329111e+02,  -5.76587375e+02,
            6.45763294e+02],
         [  2.43974026e+01,   2.65496888e+01,   2.74090204e+01,
            7.26674591e+00,   3.55635018e+01,   4.64610694e+01,
            1.05722403e+01,   3.13021818e+00,   4.47175210e+00,
            3.57097220e+00,   3.49142039e+00,   2.92687487e+00,
            2.09900754e+00,   5.36291884e+00,   1.76574464e+00,
            2.64545990e+00],
         [  2.95227326e+00,   2.95882797e+00,   2.97814765e+00,
            3.00055966e+00,   3.01967836e+00,   3.05140955e+00,
            3.12223419e+00,   3.17387217e+00,   3.22614445e+00,
            3.27424502e+00,   3.32006195e+00,   3.36172191e+00,
            3.39381359e+00,   3.43401941e+00,   3.49962690e+00,
            3.49472848e+00]])

    ksteps = np.linspace(0., 0.5, num=17)
    supcellsize = 13

    # The bigger the unit cell, the more bands fold in below the
    # waveguide band:
    first_wg_band = int(3 + 2 * (supcellsize - 1))

    sim = TriHoles2D_Waveguide_effective_epsilon(
        epsilon_cubspline_knots=fknots,
        epsilon_cubspline_coeffs=fcoeffs,
        band_number=first_wg_band,
        extra_bands=4,
        init_frequency=0.5,
        radius=0.38,
        mode='te',
        k_steps=ksteps,
        supercell_size=supcellsize,
        resolution=16,
        mesh_size=7,
        ensure_y_parity='odd',
        runmode=mode,
        num_processors=2,
        save_field_patterns_kvecs=[
            (x, 0, 0) for x in ksteps],
        plot_crop_y=False,
        convert_field_patterns=True)

    if not sim:
        log.error('an error occurred during simulation. See the .out file')
        return

    log.info(' ##### success! #####\n\n')


if __name__ == '__main__':
    main()
