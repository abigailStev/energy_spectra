#!/usr/bin/env

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import argparse
import tools
import subprocess
from scipy.optimize import leastsq

__author__ = "Abigail Stevens"

"""
multifit_plots.py

Reads an XSPEC log file and makes plots of varying fit parameters as a function
of QPO or pulse phase.

Abigail Stevens, A.L.Stevens at uva.nl, 2015

"""

class Parameter(object):
    def __init__(self, label):
        self.label = label
        self.value = np.asarray([])
        self.error = np.asarray([])
        self.lo_v = np.asarray([])
        self.hi_v = np.asarray([])
        self.pos_err = np.asarray([])
        self.neg_err = np.asarray([])
        self.par_num = np.asarray([])
        self.varying = False
        self.sinefit = None
        self.phase = None
        self.phase_err = None


class Phabs(object):
    def __init__(self):#, nH):
        self.nH = Parameter(r"phabs: nH ($\times 10^{22}$)")

class Simpl(object):
    def __init__(self):
        self.Gamma = Parameter(r"simpl: $\Gamma$")
        self.FracSctr = Parameter("simpl: FracSctr")
        self.UpScOnly = Parameter("simpl: UpScOnly")

class Nthcomp(object):
    def __init__(self):
        self.Gamma = Parameter(r"nthComp: $\Gamma$")
        self.kT_e = Parameter(r"nthComp: kT$_{e}$ (keV)")
        self.kT_bb = Parameter(r"nthComp: kT$_{bb}$ (keV)")
        self.inp_type = Parameter("nthComp: inp type")
        self.Redshift = Parameter("nthComp: Redshift")
        self.norm = Parameter("nthComp: norm.")

class Diskbb(object):
    def __init__(self):
        self.Tin = Parameter(r"diskBB: T$_{in}$ (kev)")
        self.norm = Parameter("diskBB: norm.")

class Bbodyrad(object):
    def __init__(self):
        self.kT = Parameter("bbodyrad: kT (keV)")
        self.norm = Parameter("bbodyrad: norm.")

class Gaussian(object):
    def __init__(self):
        self.LineE = Parameter("gaussian: LineE (keV)")
        self.Sigma = Parameter(r"gaussian: $\sigma$ (keV)")
        self.norm = Parameter("gaussian: norm.")


################################################################################
def read_log_file(log_file):
    """
    Reads the XSPEC log file and assigns parameters to Parameter objects, with
    the expectation that chatter is set to 4.

    Params
    ------
    log_file : string
        Full path of the XSPEC log file, with chatter set to 4. Assuming the
        spectral models listed below are the only ones used (or the only
        interesting ones).

    Returns
    -------
    Parameter object
        phabs
    Parameter object
        simpl
    Parameter object
        diskbb
    Parameter object
        gaussian
    int
        Number of spectra for one QPO phase.

    """

    phabs = Phabs()
    simpl = Simpl()
    diskbb = Diskbb()
    nth = Nthcomp()
    bbrad = Bbodyrad()
    gauss = Gaussian()

    ###########################
    ## Reading in the log file
    ###########################

    print ""
    with open(log_file, 'r') as f:
        for line in f:

            ## Reading in phabs parameter
            if "phabs" in line and "nH" in line:
                phabs.nH.value = float(line.split()[6])
                phabs.nH.par_num = np.append(phabs.nH.par_num, \
                        int(line.split()[1]))

            ## Reading simpl parameters
            elif "simpl" in line:
                if "Gamma" in line:
                    simpl.Gamma.value = np.append(simpl.Gamma.value, \
                            float(line.split()[5]))
                    simpl.Gamma.par_num = np.append(phabs.nH.par_num, \
                            int(line.split()[1]))

                elif "FracSctr" in line:
                    simpl.FracSctr.value = np.append(simpl.FracSctr.value, \
                            float(line.split()[5]))
                    simpl.FracSctr.par_num = np.append(simpl.FracSctr.par_num, \
                            int(line.split()[1]))

            ## Reading nthComp parameters
            elif "nthComp" in line:
                if "Gamma" in line:
                    nth.Gamma.value = np.append(nth.Gamma.value, \
                            float(line.split()[5]))
                    nth.Gamma.par_num = np.append(nth.Gamma.par_num, \
                            int(line.split()[1]))

                elif "norm" in line:
                    nth.norm.value = np.append(nth.norm.value, \
                            float(line.split()[5]))
                    nth.norm.par_num = np.append(nth.norm.par_num, \
                            int(line.split()[1]))

            ## Reading diskbb parameters
            elif "diskbb" in line:
                if "Tin" in line:
                    diskbb.Tin.value = np.append(diskbb.Tin.value, \
                            float(line.split()[6]))
                    diskbb.Tin.par_num = np.append(diskbb.Tin.par_num, \
                        int(line.split()[1]))

                elif "norm" in line:
                    diskbb.norm.value = np.append(diskbb.norm.value, \
                            float(line.split()[5]))
                    diskbb.norm.par_num = np.append(diskbb.norm.par_num, \
                            int(line.split()[1]))

            ## Reading bbodyrad parameters
            elif "bbodyrad" in line and " 4 " in line:
                if "kT" in line:
                    bbrad.kT.value = np.append(bbrad.kT.value, \
                            float(line.split()[6]))
                    bbrad.kT.par_num = np.append(bbrad.kT.par_num, \
                            int(line.split()[1]))
                elif "norm" in line:
                    bbrad.norm.value = np.append(bbrad.norm.value, \
                            float(line.split()[5]))
                    bbrad.norm.par_num = np.append(bbrad.norm.par_num, \
                            int(line.split()[1]))

            ## Reading gaussian parameters
            elif "gaussian" in line:
                if "LineE" in line:
                    gauss.LineE.value = np.append(gauss.LineE.value, \
                            float(line.split()[6]))
                    gauss.LineE.par_num = np.append(gauss.LineE.par_num, \
                            int(line.split()[1]))

                elif "Sigma" in line:
                    gauss.Sigma.value = np.append(gauss.Sigma.value, \
                            float(line.split()[6]))
                    gauss.Sigma.par_num = np.append(gauss.Sigma.par_num, \
                            int(line.split()[1]))

                elif "norm" in line:
                    gauss.norm.value = np.append(gauss.norm.value, \
                            float(line.split()[5]))
                    gauss.norm.par_num = np.append(gauss.norm.par_num, \
                            int(line.split()[1]))


    num_spectra = len(nth.Gamma.value)

    ###################################################################
    ## If parameter is tied across phases, assign it the correct error
    ###################################################################

    # if simpl.Gamma.value[0] != simpl.Gamma.value[1]:
    #     simpl.Gamma.varying = True
    #
    # if simpl.FracSctr.value[0] != simpl.FracSctr.value[1]:
    #     simpl.FracSctr.varying = True

    # if diskbb.Tin.value[0] != diskbb.Tin.value[1]:
    #     diskbb.Tin.varying = True
    #
    # if diskbb.norm.value[0] != diskbb.norm.value[1]:
    #     diskbb.norm.varying = True

    if nth.Gamma.value[0] != nth.Gamma.value[1]:
        nth.Gamma.varying = True

    if nth.norm.value[0] != nth.norm.value[1]:
        nth.norm.varying = True

    if bbrad.kT.value[0] != bbrad.kT.value[1]:
        bbrad.kT.varying = True

    if bbrad.norm.value[0] != bbrad.norm.value[1]:
        bbrad.norm.varying = True

    if gauss.LineE.value[0] != gauss.LineE.value[1]:
        gauss.LineE.varying = True

    if gauss.Sigma.value[0] != gauss.Sigma.value[1]:
        gauss.Sigma.varying = True

    if gauss.norm.value[0] == gauss.norm.value[1]:
        gauss.norm.varying = True

    ##################################
    ## Reading in errors from 'chain'
    ##################################

    # error_file = "/Users/abigailstevens/Dropbox/Research/energy_spectra/out_es/GX339-BQPO_NTH-2BB_errors.txt"
    error_file = "/Users/abigailstevens/Dropbox/Research/energy_spectra/out_es/GX339-BQPO-NTH-2BB-longwalker.txt"
    par_nums = np.asarray([])
    lo_v = np.asarray([])
    hi_v = np.asarray([])
    pos_err = np.asarray([])
    neg_err = np.asarray([])

    with open(error_file, 'r') as f:
        for line in f:
            par_nums = np.append(par_nums, int(line.split()[0]))
            lo_v = np.append(lo_v, float(line.split()[1]))
            hi_v = np.append(hi_v, float(line.split()[2]))
            tup = line.split()[-1].replace('(', '').replace(')', '').split(',')
            neg_err = np.append(neg_err, float(tup[0]))
            pos_err = np.append(pos_err, float(tup[1]))

    temp_mask = np.array([], dtype=bool)
    if nth.Gamma.varying:
        for elt in par_nums:
            if elt in nth.Gamma.par_num:
                temp_mask = np.append(temp_mask, True)
            else:
                temp_mask = np.append(temp_mask, False)
        nth.Gamma.pos_err = pos_err[temp_mask]
        nth.Gamma.neg_err = np.abs(neg_err[temp_mask])
        nth.Gamma.lo_v = lo_v[temp_mask]
        nth.Gamma.hi_v = hi_v[temp_mask]

    temp_mask = np.array([], dtype=bool)
    if nth.norm.varying:
        for elt in par_nums:
            if elt in nth.norm.par_num:
                temp_mask = np.append(temp_mask, True)
            else:
                temp_mask = np.append(temp_mask, False)
        nth.norm.pos_err = pos_err[temp_mask]
        nth.norm.neg_err = np.abs(neg_err[temp_mask])
        nth.norm.lo_v = lo_v[temp_mask]
        nth.norm.hi_v = hi_v[temp_mask]

    temp_mask = np.array([], dtype=bool)
    if bbrad.kT.varying:
        for elt in par_nums:
            if elt in bbrad.kT.par_num:
                temp_mask = np.append(temp_mask, True)
            else:
                temp_mask = np.append(temp_mask, False)
        bbrad.kT.pos_err = pos_err[temp_mask]
        bbrad.kT.neg_err = np.abs(neg_err[temp_mask])
        bbrad.kT.lo_v = lo_v[temp_mask]
        bbrad.kT.hi_v = hi_v[temp_mask]

    # return phabs, simpl, diskbb, gauss, num_spectra
    return phabs, nth, bbrad, diskbb, gauss, num_spectra


################################################################################
def make_plots_var3(plot_name, num_spectra, param1, param2, param3):
    """
    Making plots of fit parameters vs phase for three co-varying parameters.

    """

    font_prop = font_manager.FontProperties(size=20)
    xLocator = MultipleLocator(0.05)  ## loc of minor ticks on y-axis

    phase = np.arange(num_spectra) / 23.646776
    tinybins = np.arange(-0.01, 1.03, 0.01)
    param1_max = -1
    param2_max = -1
    param3_max = -1

    if param1.sinefit is not None:
        param1_max = tinybins[np.argmax(param1.sinefit)]
    if param2.sinefit is not None:
        param2_max = tinybins[np.argmax(param2.sinefit)]
    if param3.sinefit is not None:
        param3_max = tinybins[np.argmax(param3.sinefit)]


    print "Plot file: %s" % plot_name

    fig = plt.figure(figsize=(12, 10))

    ########################
    ## Plotting parameter 1
    ########################

    ax1 = fig.add_subplot(311)
    ax1.errorbar(phase, param1.value, yerr=[param1.neg_err, param1.pos_err], lw=0, ecolor='red', \
                 marker='.', ms=10, mec='red', mfc='red', elinewidth=2,
                 capsize=2)
    if param1.sinefit is not None:
        ax1.plot(tinybins, param1.sinefit, c='black', lw=2)
        ax1.vlines(param1_max, 0.06, 0.2, lw=1, linestyles='dashed')

    ax1.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=False, labeltop=False)
    ax1.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax1.set_ylabel(param1.label, fontproperties=font_prop)

    ax1.set_ylim(0.06, 0.2)

    ax1.set_xlim(-0.01 , 1.02)
    y_maj_loc = ax1.get_yticks()
    ax1.set_yticks(y_maj_loc[1:])
    yLocator1 = MultipleLocator(.01)  ## loc of minor ticks on y-axis
    ax1.yaxis.set_minor_locator(yLocator1)
    ax1.xaxis.set_minor_locator(xLocator)

    ########################
    ## Plotting parameter 2
    ########################

    ax2 = fig.add_subplot(312, sharex=ax1)
    ax2.errorbar(phase, param2.value, yerr=[param2.neg_err, param2.pos_err], lw=0, ecolor='green', \
                 marker='.', ms=10, mec='green', mfc='green', elinewidth=2,
                 capsize=2)
    if param2.sinefit is not None:
        ax2.plot(tinybins, param2.sinefit, c='black', lw=2)
        ax2.vlines(param2_max, 2.1, 2.7, lw=1, linestyles='dashed')

    ax2.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=False, labeltop=False)
    ax2.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax2.set_ylabel(param2.label, fontproperties=font_prop)

    ax1.set_ylim(0.06, 0.2)
    ax2.set_ylim(2.1, 2.7)

    ax2.set_xlim(-0.01, 1.02)
    y_maj_loc = ax2.get_yticks()
    ax2.set_yticks(y_maj_loc[1:-1])
    y_min_mult = 0.5 * (y_maj_loc[1] - y_maj_loc[0])
    yLocator2 = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax2.yaxis.set_minor_locator(yLocator2)
    ax2.set_xticks(np.arange(0, 1.05, 0.25))
    ax2.xaxis.set_minor_locator(xLocator)

    ##################
    ## Plotting parameter 3
    ##################

    ax3 = fig.add_subplot(313, sharex=ax1)

    ax3.errorbar(phase, param3.value, yerr=[param3.neg_err, param3.pos_err], lw=0, ecolor='blue', \
                 marker='.', ms=10, mec='blue', mfc='blue', elinewidth=2,
                 capsize=2)
    if param3.sinefit is not None:
        ax3.plot(tinybins, param3.sinefit, c='black', lw=2)
        ax3.vlines(param3_max, 0.565, 0.58, lw=1, linestyles='dashed')
    ax3.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=True, labeltop=False)
    ax3.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax3.set_ylabel(param3.label, fontproperties=font_prop)
    ax3.set_xlabel('Normalized QPO phase', fontproperties=font_prop)

    ax1.set_ylim(0.06, 0.2)
    ax2.set_ylim(2.1, 2.7)
    ax3.set_ylim(0.565, 0.58)

    # y_maj_loc = ax3.get_yticks()
    # ax3.set_yticks(y_maj_loc[0:-1])
    y_maj_loc = [0.565, 0.570, 0.575, 0.58]
    ax3.set_yticks(y_maj_loc)
    ax3.set_xlim(-0.01, 1.02)
    y_min_mult = 0.1 * (y_maj_loc[1] - y_maj_loc[0])
    yLocator3 = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax3.yaxis.set_minor_locator(yLocator3)
    ax3.set_xticks(np.arange(0, 1.05, 0.25))
    ax3.xaxis.set_minor_locator(xLocator)

    fig.subplots_adjust(hspace=0.00)
    # 	plt.show()
    plt.savefig(plot_name)
    plt.close()
    subprocess.call(['open', plot_name])


################################################################################
def make_plots_var2(plot_name, num_spectra, param1, param2):
    """
    Making plots of fit parameters vs phase for two co-varying parameters,
    fracsctr and another.

    """

    font_prop = font_manager.FontProperties(size=20)
    xLocator = MultipleLocator(0.05)  ## loc of minor ticks on y-axis

    phase = np.arange(num_spectra) / 23.646776
    tinybins = np.arange(-0.01, 1.03, 0.01)
    param1_max = -1
    param2_max = -1

    if param1.sinefit is not None:
        param1_max = tinybins[np.argmax(param1.sinefit)]
    if param2.sinefit is not None:
        param2_max = tinybins[np.argmax(param2.sinefit)]

    print "Plot file: %s" % plot_name

    fig = plt.figure(figsize=(12, 10))

    #####################
    ## Plotting FracSctr
    #####################

    ax1 = fig.add_subplot(211)
    ax1.errorbar(phase, param1.value, yerr=[param1.neg_err, param1.pos_err], \
            lw=0, ecolor='green', marker='.', ms=10, mec='green', mfc='green', \
            elinewidth=2, capsize=2)
    if param1.sinefit is not None:
        ax1.plot(tinybins, param1.sinefit, c='black', lw=2)
        ax1.vlines(param1_max, 0.06, 0.2, lw=1, linestyles='dashed')
    #   ax1.vlines(param1_max - 3.2181e-03, 0.12, 0.24, lw=1, linestyles='dotted')
    #   ax1.vlines(param1_max + 3.2181e-03, 0.12, 0.24, lw=1, linestyles='dotted')

    ax1.set_ylim(0.06, 0.2)
    ax1.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=False, labeltop=False)
    ax1.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax1.set_ylabel(param1.label, fontproperties=font_prop)

    y_maj_loc = ax1.get_yticks()
    ax1.set_yticks(y_maj_loc[1:])
    yLocator1 = MultipleLocator(.01)  ## loc of minor ticks on y-axis
    ax1.yaxis.set_minor_locator(yLocator1)
    ax1.set_xticks(np.arange(0, 1.05, 0.25))
    ax1.xaxis.set_minor_locator(xLocator)

    #######################################################
    ## Plotting other param (diskBB Tin, diskBB norm, ...)
    #######################################################

    ax2 = fig.add_subplot(212, sharex=ax1)
    ax2.errorbar(phase, param2.value, yerr=[param2.neg_err, param2.pos_err], \
            lw=0, ecolor='blue', marker='.', ms=10, mec='blue', mfc='blue', \
            elinewidth=2, capsize=2)
    # 	ax2.set_xlim(bins[0]-0.5,bins[-1]+0.5)
    if param2.sinefit is not None:
        ax2.plot(tinybins, param2.sinefit, c='black', lw=2)
        ax2.vlines(param2_max, 0.55, 0.75, lw=1, linestyles='dashed')
    #   ax2.vlines(param1_max - 3.2181e-03, 0.12, 0.24, lw=1, linestyles='dotted')
    #   ax2.vlines(param1_max + 3.2181e-03, 0.12, 0.24, lw=1, linestyles='dotted')
    ax2.set_ylim(0.55, 0.75)
    ax2.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=True, labeltop=False)
    ax2.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax2.set_ylabel(param2.label, fontproperties=font_prop)

    ax2.set_xlabel('Normalized QPO phase', fontproperties=font_prop)

    y_maj_loc = ax2.get_yticks()
    ax2.set_yticks(y_maj_loc[0:-1])
    yLocator2 = MultipleLocator(.001)  ## loc of minor ticks on y-axis
    ax2.yaxis.set_minor_locator(yLocator2)
    ax2.set_xticks(np.arange(0, 1.05, 0.25))
    ax2.xaxis.set_minor_locator(xLocator)

    fig.subplots_adjust(hspace=0.00)
    # 	plt.show()
    plt.savefig(plot_name)
    plt.close()
    subprocess.call(['open', plot_name])


################################################################################
def sinewave(t, p):
    """
    Computing a sine wave for values t with amplitude p[0], phase shift p[1],
    and mean p[2].

    Params
    ------
    t : np.array of floats
        Time steps for the sine wave.
    p : np.array of floats
        The sine wave parameters.

    Returns
    -------
    np.array of floats
        A sine wave at steps t with parameters p.
    """
    return p[0] * np.sin(2.0 * np.pi * t + p[1]) + p[2]


################################################################################
def sine_residuals(p, data, data_err, t):
    """
    Getting the residual of the data with the current fit sine wave. Dividing by
    error bar to weight it appropriately like in weighted least squares, e.g.
    S. Vaughan 2013 eqn 6.12 (modified -- not squaring because according to
    scipy.optimize.leastsq documentation, it will square for me to compute the
    real residual and covariance matrix (which will also make it look exactly
    like eqn 6.12))

    Params
    ------
    p : np.array of floats
        The sine wave parameters.
    data : np.array of floats
        The data we want to fit to; in this case, the list of fit parameters per
        energy spectra over QPO phase.
    data_err : np.array of floats
        The error on the data.
    t : np.array of floats
        Time steps for the fitting sine wave.

    Return
    ------
    np.array of floats
        A modified weighted least squared residual of the current sinewave fit
        with the data. From S. Vaughan 2013 eqn 6.12 and scipy.optimize.leastsq
        documentation.
    """
    residual = np.abs(data - sinewave(t, p)) / data_err
    return residual


################################################################################
def get_phase(sed_parameter, num_spectra):
    """
    Fitting a sine wave to an energy spectra fit parameter to determine the
    phase of the parameter changes.

    Params
    ------
    sed_parameter : Parameter object
        Description
    num_spectra : int
        The number of energy spectra in use, i.e. the length of
        sed_fit_parameters.

    Return
    ------
    Parameter
        Description.

    """
    t = np.arange(num_spectra) / 23.646776
    p = [1.0, 0.0, np.mean(sed_parameter.value)]  ## Amplitude, phase shift, mean

    sed_parameter.error = np.mean((sed_parameter.pos_err, sed_parameter.neg_err), axis=0)
    print np.shape(sed_parameter.error)

    p_best = leastsq(sine_residuals, p, args=(sed_parameter.value, sed_parameter.error, t), \
            full_output=1)
    # print "P best:", p_best
    best_fit = p_best[0]
    # print "Best fit:", best_fit

    plt.errorbar(t, sed_parameter.value, xerr=None, yerr=sed_parameter.error)
    plt.plot(t, sinewave(t, best_fit))
    plt.xlim(0,1)
    # plt.show()

    ## Error on phase from S. Vaughan 2013 p 168
    bonus_matrix = p_best[1]  ## A Jacobian approximation to the Hessian of the
            ## least squares objective function.
    resid_var = np.var(sine_residuals(best_fit, sed_parameter.value, \
            sed_parameter.error, t), ddof=1)
    ## As outlined in the scipy.optimize.leastsq documentation, multiply the
    ## bonus matrix by the variance of the residuals to get the covariance
    ## matrix.
    # print "Bonus matrix:", bonus_matrix
    # print "Resid var:", resid_var
    cov_matrix = bonus_matrix * resid_var

    sed_parameter.sinefit = sinewave(np.arange(-0.01, 1.03, 0.01), best_fit)
    sed_parameter.phase = best_fit[1] / (2.0 * np.pi)
    sed_parameter.phase_err = np.sqrt(cov_matrix[1][1]) / (2.0 * np.pi)

    return sed_parameter


################################################################################
def main(log_file):
    """
    Params
    ------
    log_file : string
        The XSPEC log file, with chatter set to 4, with extension '.log'.

    Returns
    -------
    nothing
    """
    ##########################################
    ## Reading in the log file to data arrays
    ##########################################

    # phabs, simpl, diskbb, gauss, num_spectra = read_log_file(log_file)
    phabs, nth, bbrad, diskbb, gauss, num_spectra = read_log_file(log_file)
    print "Number of spectra:", num_spectra

    # print "norm(BB) = ", bb_norm_1[1]
    # print np.mean(bb_norm_1)
    # print np.min(bb_norm_1)
    # print np.max(bb_norm_1)

    ######################################################################
    ## Computing the phase of the best-fit sine wave and phase difference
    ######################################################################

    # simpl.FracSctr = get_phase(simpl.FracSctr, num_spectra)
    #
    # print "simpl fracsctr phase: %.4e +- %.4e" % (simpl.FracSctr.phase, \
    #         simpl.FracSctr.phase_err)
    #
    # if diskbb.Tin.varying:
    #     diskbb.Tin = get_phase(diskbb.Tin, num_spectra)
    #     phase_2 = diskbb.Tin.phase
    #     print "diskBB Tin phase: %.4e +- %.4e" % (diskbb.Tin.phase, \
    #             diskbb.Tin.phase_err)
    #
    # if diskbb.norm.varying:
    #     print "\tIN HERE"
    #     diskbb.norm = get_phase(diskbb.norm, num_spectra)
    #     phase_2 = diskbb.norm.phase
    #     print "diskBB norm. phase: %.4e +- %.4e" % (diskbb.norm.phase, \
    #             diskbb.norm.phase_err)
    #
    # phase_diff = simpl.FracSctr.phase - phase_2

    nth.norm = get_phase(nth.norm, num_spectra)
    nth.Gamma = get_phase(nth.Gamma, num_spectra)
    bbrad.kT = get_phase(bbrad.kT, num_spectra)
    print nth.norm.phase
    print nth.Gamma.phase
    print bbrad.kT.phase


    # if bbrad.kT.varying:
    #     bbrad.kT = get_phase(bbrad.kT, num_spectra)
    #     phase_2 = bbrad.kT.phase
    #     print "BB T phase: %.4e +- %.4e" % (bbrad.kT.phase, \
    #             bbrad.kT.phase_err)
    # print bbrad.norm.varying
    #
    # if bbrad.norm.varying:
    #     bbrad.norm = get_phase(bbrad.norm, num_spectra)
    #     phase_2 = bbrad.norm.phase
    #     print "BB norm. phase: %.4e +- %.4e" % (bbrad.norm.phase, \
    #             bbrad.norm.phase_err)
    #
    # phase_diff = nth.norm.phase - phase_2
    #
    # if phase_diff < 0:
    #     phase_diff += 1.0
    # elif phase_diff > 1.0:
    #     phase_diff -= 1.0
    #
    # print "Phase diff: %.4f" % phase_diff

    ###################
    ## Make the plots!
    ###################

    plot_name = log_file.replace('.log', '.eps')

    # if simpl.FracSctr.varying and simpl.Gamma.varying and diskbb.Tin.varying:
    #     make_plots_var3(plot_name, num_spectra, simpl.FracSctr, diskbb.Tin, \
    #             simpl.Gamma)
    # elif simpl.FracSctr.varying and simpl.Gamma.varying and diskbb.norm.varying:
    #     make_plots_var3(plot_name, num_spectra, simpl.FracSctr, diskbb.norm, \
    #             simpl.Gamma)
    # elif simpl.FracSctr.varying and simpl.Gamma.varying:
    #     make_plots_var2(plot_name, num_spectra, simpl.FracSctr, simpl.Gamma)
    # elif simpl.FracSctr.varying and diskbb.Tin.varying:
    #     make_plots_var2(plot_name, num_spectra, simpl.FracSctr,diskbb.Tin)
    # elif simpl.FracSctr.varying and diskbb.norm.varying:
    #     make_plots_var2(plot_name, num_spectra, simpl.FracSctr, diskbb.norm)

    if nth.norm.varying and nth.Gamma.varying and bbrad.kT.varying:
        print "in here!!"
        make_plots_var3(plot_name, num_spectra, nth.norm, \
                nth.Gamma, bbrad.kT)
    elif nth.norm.varying and nth.Gamma.varying and bbrad.norm.varying:
        make_plots_var3(plot_name, num_spectra, nth.norm, bbrad.norm, \
                nth.Gamma)
    elif nth.norm.varying and nth.Gamma.varying:
        make_plots_var2(plot_name, num_spectra, nth.norm, nth.Gamma)
    elif nth.norm.varying and bbrad.kT.varying:
        make_plots_var2(plot_name, num_spectra, nth.norm, bbrad.kT)
    elif nth.norm.varying and bbrad.norm.varying:
        make_plots_var2(plot_name, num_spectra, nth.norm, bbrad.norm)


################################################################################
if __name__ == '__main__':

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python multifit_plots.py log_file",\
            description="Reads an XSPEC log file and makes plots of varying "\
            "fit parameters as a function of QPO or pulse phase.")

    parser.add_argument('log_file', help="The XSPEC log file, with chatter set"\
            "to 4.")

    args = parser.parse_args()

    main(args.log_file)

################################################################################
