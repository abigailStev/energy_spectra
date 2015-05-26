import argparse
import numpy as np
from astropy.io import fits
from tools import type_positive_int

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"

"""
        energyspec.py

Makes an energy spectrum from a cross-correlation function at a specific time
(or phase) bin. User can indicate whether to produce spectra that are mean+ccf,
ccf only (i.e., deviations from the mean), or mean only.

Written in Python 2.7, 2014-2015

"""

################################################################################
def output(out_file, phase_bin, detchans, amps, err):
    """
    Writes the energy spectrum to a file. This output file is then used as input
    for the FTOOLS script ascii2pha.

    """
    print "in output"
    with open(out_file, 'w') as out:
        for i in xrange(0, detchans):
            out.write("%d\t%.6e\t%.6e\n" % (i, amps[i], err[i]))
# 			out.write("%d\t%.6e\t%.6e\n" % (i, amps[i], amps[i]*.1))


################################################################################
def get_mean_count_rate(string):
    """
    Gets the mean count rate as an array from a string listing it with values
    separated by commas. For example, this is what you get from the 'RATE_CI'
    header value in the CCF fits file.

    """

    start_index = string.index('[')

    return np.asarray(string[start_index+1:-1].split(', '), dtype=np.float64)


################################################################################
def fits_in(in_file, phase_bin):
    """
    Gets CCF at a specific time (or phase) bin from the FITS file of CCF output.
    """
    file_hdu = fits.open(in_file)
    table = file_hdu[1].data
    obs_time = file_hdu[0].header['EXPOSURE']
    detchans = file_hdu[0].header['DETCHANS']
    mean_count_rate = get_mean_count_rate(file_hdu[0].header['RATE_CI'])
    mean_count_rate[np.where(mean_count_rate < 0.0)] = 0
    file_hdu.close()

    time_bin_mask = table.field('TIME_BIN') == phase_bin
    table_i = table[time_bin_mask]
    ccf_amps = table_i.field('CCF')
    ccf_err = table_i.field('ERROR')

    return ccf_amps, ccf_err, obs_time, mean_count_rate, detchans


################################################################################
def main(in_file, out_file, phase_bin, spec_type):
    """
    Finds the time bin desired, gets the CCF amplitudes, computes energy
    spectrum, sends to output.

    """

    ###################
    ## Initializations
    ###################

    amps = []
    err = []

    #####################################
    ## Reading in CCF based on file type
    #####################################

    assert in_file[-4:].lower() == 'fits', "ERROR: Input file must have "\
            "extension .fits."

    ccf_amps, ccf_err, obs_time, mean_count_rate, detchans = fits_in(in_file, \
            phase_bin)

    ##############################################################
    ## Computes the type of energy spectrum indicated by the user
    ##############################################################

    mean_err = np.sqrt(mean_count_rate * obs_time) / obs_time

    if spec_type == 0:
        amps = np.add(ccf_amps, mean_count_rate)
        err = np.sqrt(np.add(np.square(ccf_err), np.square(mean_err)))
    elif spec_type == 1:
        amps = ccf_amps
        err = ccf_err
    elif spec_type == 2:
        amps = mean_count_rate
        err = mean_err
    else:
        raise Exception("ERROR: Spectrum type is not a valid option.")

    ##########
    ## Output
    ##########

    output(out_file, phase_bin, detchans, amps, err)


################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python energyspec.py infile "\
            "outfile -b phase_bin [-s spectrum_type]", description="Makes "\
            "energy spectrum from a cross-correlation function at a specific "\
            "time bin.", epilog="For optional arguments, default values are "\
            "given in brackets at end of description.")

    parser.add_argument('infile', help="The full path of the (.fits or .dat) "\
            "input file listing the CCF amplitudes per energy bin.")

    parser.add_argument('outfile', help="The full path of the FITS output file"\
            " to write the energy spectrum to.")

    parser.add_argument('-b', '--bin', type=type_positive_int, dest='phase_bin', \
            default=0, help="The phase bin number of the CCF to make an energy"\
            " spectrum for. [0]")

    parser.add_argument('-s', '--spec', type=int, dest='spec_type',
            choices=range(0, 3), default=0, help="Indicating the type of "\
            "spectrum to produce. 0 for mean+ccf, 1 for ccf, 2 for mean. [0]")

    args = parser.parse_args()

    main(args.infile, args.outfile, args.phase_bin, args.spec_type)

################################################################################
