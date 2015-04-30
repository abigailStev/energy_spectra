# energy_spectra

Takes a two-dimensional cross correlation function and writes multiple energy 
spectra from it (per specified phase bin). Look at the 1-D or 2-D CCF to see
which phase bins might be interesting to look at. Now also able to fit multiple
energy spectra simultaneously!

## Contents

### energyspec.py
Takes CCF amplitudes of a specific time bin and writes them to a file. Can make
only a mean energy spectrum (of just the mean count rate per channel of 
interest), mean \+ ccf for a phase-resolved energy spectrum, or only the ccf
for a deviation of each phase-resolved energy spectrum from the mean energy 
spectrum. Units should be photon count rate.

### run_energyspec.sh
Bash script to run energyspec.py and make and run XSPEC scripts to make plots of
the energy spectra.

### spectra/ccfonly.pco
A QDP style file for plotting the ccf deviations with stepped lines.

### spectra/ccfonly_points.pco
A QDP style file for plotting the ccf deviations with unconnected points.

### spectra/ccfwmean.pco
A QDP style file for plotting the mean+ccf (i.e. phase-resolved energy 
spectrum) with stepped lines.

### spectra/ccfwmean_points.pco
A QDP style file for plotting the mean+ccf (i.e. phase-resolved energy 
spectrum) with unconnected points.

### sed_fitting.sh
Simultaneously fits multiple energy spectra (SEDs) from different points in the 
QPO phase. Makes the phase-resolved energy spectra, writes the xspec script, 
runs the xspec script, outputs to a log file.

### multifit_plots.py
Reads the xspec log file to plot how the free parameters change with QPO phase.
Fits the changing parameters with a sine wave to get the phase of each.

##### Disclaimer: This code comes with no legal guarantees.


## Old
#### ratio_spectrum.py
#### mean_spectrum.py
#### mean_spectrum.sh

