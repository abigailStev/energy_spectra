#!/usr/bin/env python
"""
Makes a spectral energy distribution from a cross-correlation function at a
specific time bin. User can indicate whether to produce spectra that are
mean+ccf, ccf only (i.e., deviations from the mean), or mean only.

Enter    python energyspec.py -h    at the command line for help.
"""
import argparse
import numpy as np
from astropy.table import Table
__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2014-2015"


###############################################################################
def get_mean_count_rate(keyword):
    """
    Get the mean count rate as an array from a string listing it with values
    separated by commas. For example, this is what you get from the 'RATE_CI'
    header value in the CCF fits file.

    Parameters
    ----------
    keyword : str
        The FITS header keyword containing the mean count rate of the channels
        of interest as a string surrounded by [ ] and separated by commas.

    Returns
    -------
    np.array of floats
        The mean count rate of the channels of interest. Size = (detchans).
    """

    start_index = keyword.index('[')
    return np.asarray(keyword[start_index+1:-1].split(', '), dtype=np.float64)


################################################################################
def ccf_in(ccf_file, time_bin):
    """
    Get the ccf at a specific time bin from the cross_correlation/ccf.py output
    file in FITS format.

    Parameters
    ----------
    ccf_file : str
        The FITS file containing the ccf, saved as an astropy table.

    time_bin : int
        The time bin of the ccf to get the energy spectrum for. Must be in range
        0 to n_bins-1, inclusive.

    Returns
    -------
    ccf, error : np.arrays of floats
        2-D arrays of the cross-correlation function (ccf) and error on the ccf
        as computed in cross_correlation/ccf.py. Sizes = (n_bins, detchans).

    obs_time : float
        The exposure time of the observation(s) used to compute the ccf.

    ci_mean_rate : np.array of floats
        1-D array of the mean count rate in the channels of interest.
        Size = (detchans)

    detchans : int
        The number of detector energy channels for the data mode.

    Raises
    ------
    Exception if ccf file cannot be opened by astropy.table.Table.read.

    Exception if time_bin is out of bounds.

    """

    try:
        in_table = Table.read(ccf_file)
    except IOError:
        print("ERROR: File does not exist: %s" % ccf_file)
        exit()

    obs_time = in_table.meta['EXPOSURE']
    n_bins = in_table.meta['N_BINS']
    ci_mean_rate = get_mean_count_rate(in_table.meta['RATE_CI'])
    detchans = in_table.meta['DETCHANS']

    try:
        ccf = in_table['CCF'][time_bin, :]
        error = in_table['ERROR'][time_bin, :]
    except IndexError:
        print("ERROR: Time bin %d is out of bounds. Must be between 0 and %d." \
                % (time_bin, n_bins))
        exit()

    return ccf, error, obs_time, ci_mean_rate, detchans


################################################################################
def main(in_file, out_file, time_bin, spec_type):
    """
    Get the ccf and error at the desired time bin, compute desired type of
    spectral energy distribution (SED), write to an ascii table.

    Parameters
    ----------
    in_file : str
        The FITS file containing the ccf, saved as an astropy table.

    out_file : str
        Name of the ASCII output file to write the spectral energy distribution
        to.

    time_bin : int
        The time bin of the ccf to make a spectral energy distribution for.

    spec_type : { 0 | 1 | 2 }
        The type of spectrum to produce. 0 for mean+ccf, 1 for ccf, 2 for mean.

    """

    ###########################
    ## Read in ccf output file
    ###########################
    assert in_file[-4:].lower() == 'fits', "ERROR: Input file must have "\
            "extension .fits."

    ccf_amps, ccf_err, obs_time, ci_mean_rate, detchans = ccf_in(in_file, \
            time_bin)

    ##############################################################
    ## Compute the type of energy spectrum indicated by the user
    ##############################################################

    ci_mean_rate[ci_mean_rate < 0.0] = np.nan
    mean_err = np.sqrt(ci_mean_rate * obs_time) / obs_time

    if spec_type == 0:
        amps = np.add(ccf_amps, ci_mean_rate)
        err = np.sqrt(np.add(np.square(ccf_err), np.square(mean_err)))
    elif spec_type == 1:
        amps = ccf_amps
        err = ccf_err
    elif spec_type == 2:
        amps = ci_mean_rate
        err = mean_err
    else:
        raise Exception("Spectrum type is not a valid option. Must indicate "\
                "0, 1, or 2.")

    ##########
    ## Output
    ##########

    out_table = Table([np.arange(0,detchans), amps, err],
            names=("CHAN", "RATE", "ERROR"))
    # out_table.pprint()
    out_table.write(out_file, format='ascii.no_header')


################################################################################
if __name__ == "__main__":

    #########################################
    ## Parse input arguments and call 'main'
    #########################################

    parser = argparse.ArgumentParser(usage="python energyspec.py infile "\
            "outfile [OPTONAL ARGUMENTS]", description=__doc__,
            epilog="For optional arguments, default values are given in "\
            "brackets at end of description.")

    parser.add_argument('infile', help="Name of the .fits format "\
            "ccf file as an astropy table.")

    parser.add_argument('outfile', help="Name of the .dat output file to write"\
            " the spectral energy distribution (SED) to.")

    parser.add_argument('-b', '--bin', type=int, dest='time_bin', default=0,
            help="The time bin of the ccf to make an SED for. [0]")

    parser.add_argument('-s', '--spec', type=int, dest='spec_type',
            choices=range(0, 3), default=0, help="Indicating the type of "\
            "spectrum to produce. 0 for mean+ccf, 1 for ccf, 2 for mean. [0]")

    args = parser.parse_args()

    main(args.infile, args.outfile, args.time_bin, args.spec_type)

################################################################################
