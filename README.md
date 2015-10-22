# energy_spectra

Takes a two-dimensional cross correlation function (CCF) and writes multiple 
spectral energy distributions (SEDs) from it (per specified phase bin). Fits 
multiple SEDs simultaneously with XSPEC, and fits a function to those SED 
parameter variations.

## Authors
Abigail Stevens (UvA)

## Contents

### energyspec.py
Takes CCF amplitudes of a specific time bin and writes them to a file. Can make
only a mean SED (of just the mean count rate per channel of interest), mean \+ 
ccf for a phase-resolved SED, or only the ccf for a deviation of each phase-
resolved SEDs from the mean SED. Units should be photon count rate.

### run_energyspec.sh
Bash script to run energyspec.py and make and run XSPEC scripts to make plots of
the SEDs.

### spectra/ccfonly.pco
A QDP style file for plotting the ccf deviations with stepped lines.

### spectra/ccfonly_points.pco
A QDP style file for plotting the ccf deviations with unconnected points.

### spectra/ccfwmean.pco
A QDP style file for plotting the mean+ccf (i.e. phase-resolved SED) with 
stepped lines.

### spectra/ccfwmean_points.pco
A QDP style file for plotting the mean+ccf (i.e. phase-resolved SED) with 
unconnected points.

### sed_fitting.sh
Simultaneously fits multiple energy spectra (SEDs) from different points in the 
QPO phase. Makes the phase-resolved energy spectra, writes the xspec script, 
runs the xspec script, outputs to a log file.

### multifit_plots.py
Reads the xspec log file to plot how the untied SED parameters change with QPO 
phase. Fits the changing parameters with a function and gets the 'phase' of 
each parameter variation.


## Copyright
 
All content © 2014-2015 the authors. The code is distributed under the MIT 
license. See LICENSE.md for details.

If you use this code, please cite the paper as well as the code's Zenodo DOI! 
Stevens et al. 2015

Pull requests are welcome! If you are interested in the further development of 
energy_spectra, please [get in touch via the issues](https://github.com/abigailStev/energy_spectra/issues)!


 [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/) 