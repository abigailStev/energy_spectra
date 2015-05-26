import argparse
import numpy as np
from astropy.io import fits
from tools import read_obs_time, type_positive_int

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
    separated by commas.

    """

    start_index = string.index('[')

    return np.asarray(string[start_index+1:-1].split(', '), dtype=np.float64)


################################################################################
def dat_in(in_file, ccf_amps_and_err):
    """
    Gets CCF from a .dat file.

    """

    ccf_amps_and_err = np.zeros(128, dtype=np.float64)

    with open(in_file, 'r') as f:
        for line in f:
            if line[0].strip() != "#":
                line = line.strip().split()
                if int(line[0]) == phase_bin:
# 					print "Bin found!"
                    for i in xrange(0, 128):
                        ccf_amps_and_err[i] = float(line[i+1])
                    break
            else:
                if "Mean" in line.strip() and \
                    "count rate" in line.strip() and \
                    "ci" in line.strip():
                    mean_count_rate = get_mean_count_rate(line.strip())
        else:
            raise Exception("ERROR: Phase bin not found. Check that it is \
within the range of the file.")
    ## End of with-block

    ccf_amps = ccf_amps_and_err[0:64]
    ccf_err = ccf_amps_and_err[64:128]
    obs_time = read_obs_time(in_file)

    return ccf_amps, ccf_err, obs_time, mean_count_rate


################################################################################
def fits_in(in_file, phase_bin):
    """
    Gets CCF at a specific time (or phase) bin from the FITS file of CCF output.
    """
    print in_file
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

    ccf_amps, ccf_err, obs_time, mean_count_rate, detchans = fits_in(in_file)
    print "\t", phase_bin

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
        print amps
        print err
# 		amps = np.where(mean_count_rate != 0, ccf_amps / mean_count_rate, 0)
# 		err = np.where(mean_count_rate != 0, ccf_err / mean_count_rate, 0)
    elif spec_type == 2:
        amps = mean_count_rate
        err = mean_err
    else:
        raise Exception("ERROR: Spectrum type not a valid option.")
    ## End of if/else spectrum type

    ##########
    ## Output
    ##########
    print "still here!"
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
