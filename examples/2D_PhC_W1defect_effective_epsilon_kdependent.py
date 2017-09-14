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

from pympb.phc_simulations import \
    TriHoles2D_Waveguide_effective_epsilon_k_dependent
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
    # k dependent epsilon function:
    fknots = np.array(
        [0.     ,  0.03125,  0.0625 ,  0.09375,  0.125  ,  0.15625,
         0.1875 ,  0.21875,  0.25   ,  0.28125,  0.3125 ,  0.34375,
         0.375  ,  0.40625,  0.4375 ,  0.46875,  0.5    ])

    # the coefficients of the cubic spline describing the
    # k dependent epsilon function:
    fcoeffs = np.array(
        [[-1.26479352e-02,   3.79438056e-02,  -1.39127287e-01,
           5.18565342e-01,  -1.93513408e+00,   7.22197099e+00,
          -2.69527499e+01,   6.47721598e+01,  -7.50651090e+01,
          -1.26125949e+01,   2.63548544e+02,  -2.26207933e+02,
          -5.13747147e+01,   1.41246898e+02,  -3.11230826e+01,
           4.84425567e+01],
        [  3.95247975e-04,  -7.90495950e-04,   2.76673582e-03,
          -1.02764473e-02,   3.83390535e-02,  -1.43079767e-01,
           5.33980013e-01,  -1.99284029e+00,   4.07954969e+00,
          -2.95780428e+00,  -4.14023504e+00,   2.05674410e+01,
          -6.39552710e-01,  -5.45593221e+00,   7.78596450e+00,
           4.86817550e+00],
        [ -1.70254597e+00,  -1.70255832e+00,  -1.70249656e+00,
          -1.70273124e+00,  -1.70185428e+00,  -1.70512743e+00,
          -1.69291180e+00,  -1.73850118e+00,  -1.67329151e+00,
          -1.63823697e+00,  -1.86005070e+00,  -1.34670051e+00,
          -7.23954002e-01,  -9.14437905e-01,  -8.41624396e-01,
          -4.46182521e-01],
        [  3.64596772e+00,   3.59276316e+00,   3.53955860e+00,
           3.48635403e+00,   3.43314947e+00,   3.37994491e+00,
           3.32674035e+00,   3.27353579e+00,   3.21923818e+00,
           3.16864095e+00,   3.11417266e+00,   3.06004574e+00,
           3.03114342e+00,   3.00632747e+00,   2.97673374e+00,
           2.95708665e+00]])

    ksteps = np.linspace(0., 0.5, num=17)
    supcellsize = 13

    # The bigger the unit cell, the more bands fold in below the
    # waveguide band:
    first_wg_band = int(3 + 2 * (supcellsize - 1))

    sim = TriHoles2D_Waveguide_effective_epsilon_k_dependent(
        epsilon_cubspline_knots=fknots,
        epsilon_cubspline_coeffs=fcoeffs,
        band_number=first_wg_band,
        extra_bands=10,
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
        convert_field_patterns=True,
        gap=(0.44090, 0.52120))

    if not sim:
        log.error('an error occurred during simulation. See the .out file')
        return

    log.info(' ##### success! #####\n\n')


if __name__ == '__main__':
    main()
