# energy_spectra

Takes a two-dimensional cross correlation function (CCF) and writes multiple 
energy spectra from it (per specified phase bin). Fits multiple energy spectra
simultaneously with XSPEC, and fits a function to those spectral
parameter variations. Please see [Stevens & Uttley 2016](https://ui.adsabs.harvard.edu/#abs/2016MNRAS.460.2796S/abstract)
for reference.

## Contents

### energyspec.py
Takes CCF amplitudes of a specific time bin and writes them to a file. Can make
only a mean spectrum (of just the mean count rate per channel of interest),
mean+ccf for a phase-resolved spectrum, or only the ccf for a deviation of
each phase-resolved spectrum from the mean spectrum. Units should be photon
count rate.

### multifit_plots.py
Reads the xspec log file to plot how the untied spectral parameters change with
QPO phase. Fits the changing spectral parameters with a sinusoid+harmonic
function and gets the 'phase' of each parameter variation.

### out_es/ccfonly.pco
A QDP style file for plotting the ccf deviations with stepped lines.

### out_es/ccfonly_points.pco
A QDP style file for plotting the ccf deviations with unconnected points.

### out_es/ccfwmean.pco
A QDP style file for plotting the mean+ccf (i.e. phase-resolved spectra) with
stepped lines.

### out_es/ccfwmean_points.pco
A QDP style file for plotting the mean+ccf (i.e. phase-resolved spectra) with
unconnected points.

### phasespec_fitting.sh
Simultaneously fits multiple energy spectra from different points in the
QPO phase. Makes the phase-resolved energy spectra, writes the xspec script, 
runs the xspec script, outputs to a log file.

### run_energyspec.sh
Bash script to run energyspec.py and make and run XSPEC scripts to make plots of
the spectra.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Authors
* Abigail Stevens (UvA)

## Copyright
 
All content © 2014-2017 the Authors, and is distributed under the MIT
license. See LICENSE.md for details.

If you use this code, please cite [Stevens & Uttley 2016](https://ui.adsabs.harvard.edu/#abs/2016MNRAS.460.2796S/abstract).

The functionality of this software will be folded into [Stingray](http://stingraysoftware.github.io/),
so please get involved over there!


