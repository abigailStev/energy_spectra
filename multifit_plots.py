#!/usr/bin/env

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import argparse
import subprocess
from scipy.optimize import leastsq
from scipy.signal import sawtooth
import os.path

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2015"

"""
Reads an XSPEC log file and makes plots of varying fit parameters as a function
of QPO or pulse phase.

"""

class Parameter(object):
    def __init__(self, mod_name, label, par_name):
        """
        mod_name : str
            The exact string of the model of this parameter to search for and
            match for in an XSPEC log file.
        label : str
            What to label a graph (done up nice and pretty).
        par_name : str
            The exact string of this parameter to search for and match for in an
            XSPEC log file.
        """
        self.mod_name = mod_name
        self.label = label
        self.par_name = par_name
        self.value = np.array([])
        self.error = 0
        self.lo_v = 0
        self.hi_v = 0
        self.pos_err = 0
        self.neg_err = 0
        self.par_num = np.array([])
        self.varying = False
        self.funcfit = None
        self.best_fit = None
        self.phase = None
        self.phase_err = None

    def __str__(self):
        return "%s" % (self.label)

class Phabs(object):
    def __init__(self):
        self.mod_name = "phabs"
        self.nH = Parameter(self.mod_name, r"phabs: nH ($\times 10^{22}$)", "nH")

    def __str__(self):
        return "%s" % (self.mod_name)

class Simpl(object):
    def __init__(self):
        self.mod_name = "simpl "
        self.Gamma = Parameter(self.mod_name, "simpl: Gamma", "Gamma")
        self.FracSctr = Parameter(self.mod_name, "simpl: FracSctr", "FracSctr")
        self.UpScOnly = Parameter(self.mod_name, "simpl: UpScOnly", "UpScOnly")

    def __str__(self):
        return "%s" % (self.mod_name)

class Simpler(object):
    def __init__(self):
        self.mod_name = "simpler"
        self.Gamma = Parameter(self.mod_name, "simpler: Gamma", "Gamma")
        self.FracSctr = Parameter(self.mod_name, "simpler: FracSctr", "FracSctr")
        self.UpScOnly = Parameter(self.mod_name, "simpler: UpScOnly", "UpScOnly")

    def __str__(self):
        return "%s" % (self.mod_name)

class Nthcomp(object):
    def __init__(self):
        self.mod_name = "nthComp"
        self.Gamma = Parameter(self.mod_name, "nthComp: Gamma", "Gamma")
        self.kT_e = Parameter(self.mod_name, r"nthComp: kT$_{e}$ (keV)", "kT_e")
        self.kT_bb = Parameter(self.mod_name, r"nthComp: kT$_{bb}$ (keV)", "kT_bb")
        self.inp_type = Parameter(self.mod_name, "nthComp: inp type", "inp_type")
        self.Redshift = Parameter(self.mod_name, "nthComp: Redshift", "Redshift")
        self.norm = Parameter(self.mod_name, "nthComp: norm", "norm")

    def __str__(self):
        return "%s" % (self.mod_name)

class Diskbb(object):
    def __init__(self):
        self.mod_name = "diskbb"
        self.Tin = Parameter(self.mod_name, r"diskbb: T$_{in}$ (keV)", "Tin")
        self.norm = Parameter(self.mod_name, "diskbb: norm", "norm")

    def __str__(self):
        return "%s" % (self.mod_name)

class Diskpbb(object):
    def __init__(self):
        self.mod_name = "diskpbb"
        self.Tin = Parameter(self.mod_name, r"diskpbb: T$_{in}$ (keV)", "Tin")
        self.p = Parameter(self.mod_name, "diskpbb: p", " p ")
        self.norm = Parameter(self.mod_name, "diskpbb: norm", "norm")

    def __str__(self):
        return "%s" % (self.mod_name)

class Bbodyrad(object):
    def __init__(self):
        self.mod_name = "bbodyrad"
        self.kT = Parameter(self.mod_name, "bbodyrad: kT (keV)", "kT")
        self.norm = Parameter(self.mod_name, "bbodyrad: norm", "norm")

    def __str__(self):
        return "%s" % (self.mod_name)

class Gaussian(object):
    def __init__(self):
        self.mod_name = "gaussian"
        self.LineE = Parameter(self.mod_name, "gaussian: LineE (keV)", "LineE")
        self.Sigma = Parameter(self.mod_name, "gaussian: Sigma (keV)", "Sigma")
        self.norm = Parameter(self.mod_name, "gaussian: norm", "norm")

    def __str__(self):
        return "%s" % (self.mod_name)

class Diskline(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.LineE = Parameter(mod_name, r"diskline: LineE (keV)", "LineE")
        self.norm = Parameter(mod_name, r"diskline: norm", "norm")

    def __str__(self):
        return "%s" % (self.mod_name)

class Cutoffpl(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.PhoIndex = Parameter(mod_name, r"cutoffpl: PhoIndex", "PhoIndex")
        self.norm = Parameter(mod_name, r"cutoffpl: norm", "norm")

    def __str__(self):
        return "%s" % (self.mod_name)

################################################################################
def get_logfile_errors(log_file, num_spectra, mod_components):
    """
    Reads (fake-ish) errors on the varying parameters from the xspec log file.
    Only uses the first value and repeats along an array, since just need this
    for approximating when making first plots with a new fit.

    You need to run MCMC error estimation (and use get_chain_error) to get real
    errors on the varying parameters!!

    Parameters
    ----------
    log_file : string
        Full path of the XSPEC log file, with chatter set to 4.

    num_spectra : int
        Number of energy spectra being simultaneously fit in one QPO phase.

    mod_components : list of Parameter objects
        The model parameters used.

    Returns
    -------
    mod_components
        List of model parameters used, with the errors assigned.

    """

    with open(log_file, 'r') as f:
        for line in f:

            for component in mod_components:
                if component.mod_name in line and component.par_name in line \
                        and component.varying:
                    temp = line.split()[-1].strip()
                    if temp == "frozen":
                        component.pos_err = np.zeros(num_spectra)
                        component.neg_err = np.zeros(num_spectra)
                    else:
                        component.pos_err = np.repeat(float(temp), num_spectra)
                        component.neg_err = np.repeat(float(temp), num_spectra)

    return mod_components


################################################################################
def read_chain(error_file_o):
    """
    Gets the error on the varying parameters from an MCMC error estimation.

    Paramters
    ---------
    error_file_o : string
        The open file object of the log file with chain errors at the end.

    Returns
    -------
    par_nums
    lo_v
    hi_v
    neg_err
    pos_err
    """
    par_nums = np.array([])
    lo_v = np.array([])
    hi_v = np.array([])
    pos_err = np.array([])
    neg_err = np.array([])

    for line in error_file_o:
        if ("XSPEC: quit" not in line) and ("Spectrum" not in line) and \
                (len(line) > 2):
            # print line.split()
            # print len(line)
            par_nums = np.append(par_nums, int(line.split()[1]))
            lo_v = np.append(lo_v, float(line.split()[2]))
            hi_v = np.append(hi_v, float(line.split()[3]))
            tup = line.split()[-1].replace('(', '').replace(')', '').split(',')
            neg_err = np.append(neg_err, float(tup[0]))
            pos_err = np.append(pos_err, float(tup[1]))
        else:
            return par_nums, lo_v, hi_v, neg_err, pos_err


################################################################################
def get_chain_errors(mod_components, par_nums, lo_v, hi_v, neg_err, pos_err):
    """
    Gets the error on the varying parameters from an MCMC error estimation.

    Paramters
    ---------
    var_pars

    par_nums

    lo_v

    hi_v

    neg_err

    pos_err

    Returns
    -------
    var_pars


    """

    for component in mod_components:
        temp_mask = np.array([], dtype=bool)

        if component.varying:
            for elt in par_nums:
                if elt in component.par_num:
                    temp_mask = np.append(temp_mask, True)
                else:
                    temp_mask = np.append(temp_mask, False)
            component.pos_err = pos_err[temp_mask]
            component.neg_err = np.abs(neg_err[temp_mask])
            component.lo_v = lo_v[temp_mask]
            component.hi_v = hi_v[temp_mask]

    return mod_components


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
    var_pars : list of Parameter objects
        A list of the parameters that vary with QPO phase.

    int
        Number of spectra being simultaneously fit for one QPO phase.

    """

    if not os.path.isfile(log_file) or os.path.getsize(log_file) == 0:
        raise Exception("Log file does not exist or is empty.")

    chains = False
    mod_components = [Phabs().nH, Simpler().Gamma, Simpler().FracSctr, \
            Simpler().UpScOnly, Simpl().Gamma, Simpl().FracSctr, \
            Simpl().UpScOnly, Diskbb().Tin, Diskbb().norm, Diskpbb().Tin, \
            Diskpbb().p, Diskpbb().norm, Bbodyrad().kT, Bbodyrad().norm, \
            Gaussian().LineE, Gaussian().Sigma, Gaussian().norm]

    #################################################
    ## Reading in parameter values from the log file
    #################################################

    print("")
    with open(log_file, 'r') as f:
        for line in f:
            for component in mod_components:
                if component.mod_name in line and component.par_name in line:
                    if "frozen" in line:
                        component.value = np.append(component.value, \
                            float(line.split()[-2]))
                        component.par_num = np.append(component.par_num, \
                            int(line.split()[1]))
                    else:
                        component.value = np.append(component.value, \
                                float(line.split()[-3]))
                        component.par_num = np.append(component.par_num, \
                                int(line.split()[1]))

            if "Parameter" in line and "Confidence Range" in line:
                chains = True
                par_nums, lo_v, hi_v, neg_err, pos_err = read_chain(f)


    #############################################################
    ## Delete components if they're not used/present
    ## Determine if parameter varies across QPO phase or is tied
    ## Assign zero error to components
    #############################################################

    # num_spectra = np.amax([len(nth.Gamma.value), len(simpler.Gamma.value)])
    unused_components = []
    num_spectra = 1

    for component in mod_components:
        # print component.mod_name, component.par_name
        if len(component.value) > 1:
            # print component.value[0], component.value[1], component.value[3]
            if component.value[0] != component.value[1] or \
                    component.value[7] != component.value[0]:
                component.varying = True
                num_spectra = len(component.value)
        elif len(component.value) == 0:
            unused_components.append(component)

        # component.pos_err = np.zeros(len(component.value))
        # component.neg_err = np.zeros(len(component.value))
        # component.lo_v = component.value
        # component.hi_v = component.value

    for elt in unused_components:
        mod_components.remove(elt)

    ############################################################################
    ## Reading in errors from 'chain' part of log, or bad errors from normal log
    ############################################################################

    var_pars = 0
    # print "VarPars"
    for component in mod_components:
        if component.varying:
            var_pars += 1
            # print component.mod_name, component.par_name

    if var_pars == 0:
        raise Exception("No parameters vary with QPO phase in this log file.")
        exit()

    if chains:
        mod_components = get_chain_errors(mod_components, par_nums, lo_v, hi_v,\
                neg_err, pos_err)

    if not chains:
        print("Using fake errors from log file. Need to run error analysis!!")
        mod_components = get_logfile_errors(log_file, num_spectra, \
                mod_components)

    return mod_components, num_spectra


################################################################################
def make_plots_var3(plot_name, num_spectra, param1, param2, param3):
    """
    Making plots of fit parameters vs phase for three co-varying parameters.

    """

    font_prop = font_manager.FontProperties(size=18)
    xLocator = MultipleLocator(0.05)  ## loc of minor ticks on y-axis

    phase = np.arange(num_spectra) / 23.646776
    tinybins = np.arange(-0.02, 1.02, 0.01)

    ## So that the plotted x-value is the MIDDLE of the 'bin', with and error of
    ## the width of the bin.
    plusphase = (phase[1]-phase[0])/2.0
    phase_err = np.repeat(plusphase, num_spectra)
    phase += plusphase
    tinybins += plusphase

    param1_max = -1
    param2_max = -1
    param3_max = -1
    if param1.funcfit is not None:
        param1_max = tinybins[np.argmax(param1.funcfit)]
    if param2.funcfit is not None:
        param2_max = tinybins[np.argmax(param2.funcfit)]
    if param3.funcfit is not None:
        param3_max = tinybins[np.argmax(param3.funcfit)]

    print("Plot file: %s" % plot_name)

    fig = plt.figure(figsize=(10, 12), tight_layout=True, dpi=300)

    ########################
    ## Plotting parameter 1
    ########################

    ymax1 = np.max(param1.value+param1.pos_err)*1.1
    ymin1 = np.min(param1.value-param1.pos_err)*0.9

    ax1 = fig.add_subplot(311)
    ax1.errorbar(phase, param1.value, xerr=phase_err, yerr=[param1.neg_err, param1.pos_err],
            lw=2, color='red', drawstyle='steps-mid', marker='.', ms=10,
            mec='red', mfc='red', ecolor='red', elinewidth=2, capsize=0)
    if param1.funcfit is not None:
        # rect1 = patches.Rectangle((param1_max-param1.phase_err, ymin1),
        #         2*param1.phase_err, ymax1-ymin1, color='pink', ec="none",
        #         alpha=0.3)
        # ax1.add_patch(rect1)
        ax1.plot(tinybins, param1.funcfit, c='black', lw=2)
        ax1.vlines(param1_max, ymin1, ymax1, lw=2, color='gray',
                linestyles='dashed')

    ax1.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=False, labeltop=False)
    ax1.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax1.set_ylabel(param1.label, fontproperties=font_prop)

    ax1.set_ylim(ymin1, ymax1)
    ax1.set_xlim(0.0, 1.01)
    y_maj_loc = ax1.get_yticks()
    # ax1.set_yticks(y_maj_loc[1:])
    yLocator1 = MultipleLocator(.01)  ## loc of minor ticks on y-axis
    ax1.yaxis.set_minor_locator(yLocator1)
    ax1.xaxis.set_minor_locator(xLocator)

    ########################
    ## Plotting parameter 2
    ########################

    ymax2 = np.max(param2.value+param2.pos_err)*1.1
    ymin2 = np.min(param2.value-param2.pos_err)*0.9

    ax2 = fig.add_subplot(312, sharex=ax1)
    ax2.errorbar(phase, param2.value, xerr=phase_err, yerr=[param2.neg_err, param2.pos_err],
            lw=2, color='green', drawstyle='steps-mid', marker='.', ms=10,
            mec='green', mfc='green', ecolor='green', elinewidth=2, capsize=0)
    if param2.funcfit is not None:
        # rect2 = patches.Rectangle((param2_max-param2.phase_err, ymin2),
        #         2*param2.phase_err, ymax2-ymin2, color='pink', ec="none",
        #         alpha=0.3)
        # ax2.add_patch(rect2)
        ax2.plot(tinybins, param2.funcfit, c='black', lw=2)
        ax2.vlines(param2_max, ymin2, ymax2, lw=2, color='gray',
                linestyles='dashed')

    ax2.tick_params(axis='x', labelsize=18, bottom=True, top=True,
                    labelbottom=False, labeltop=False)
    ax2.tick_params(axis='y', labelsize=18, left=True, right=True,
                    labelleft=True, labelright=False)
    ax2.set_ylabel(param2.label, fontproperties=font_prop)

    ax2.set_ylim(ymin2, ymax2)
    ax2.set_xlim(0.0, 1.01)
    y_maj_loc = ax2.get_yticks()
    # ax2.set_yticks(y_maj_loc[1:-1])
    y_min_mult = 0.5 * (y_maj_loc[1] - y_maj_loc[0])
    yLocator2 = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax2.yaxis.set_minor_locator(yLocator2)
    ax2.set_xticks(np.arange(0, 1.05, 0.25))
    ax2.xaxis.set_minor_locator(xLocator)

    ########################
    ## Plotting parameter 3
    ########################

    ymax3 = np.max(param3.value+param3.pos_err)*1.1
    ymin3 = np.min(param3.value-param3.neg_err)*0.9

    ax3 = fig.add_subplot(313, sharex=ax1)
    ax3.errorbar(phase, param3.value, xerr=phase_err, yerr=[param3.neg_err, param3.pos_err],
                 lw=2, color='blue', drawstyle='steps-mid', marker='.', ms=10,
                 mec='blue', mfc='blue', ecolor='blue', elinewidth=2, capsize=0)
    if param3.funcfit is not None:
        # rect3 = patches.Rectangle((param3_max-param3.phase_err, ymin3),
        #         2*param3.phase_err, ymax3-ymin3, color='pink', ec="none",
        #         alpha=0.3)
        # ax3.add_patch(rect3)
        ax3.plot(tinybins, param3.funcfit, c='black', lw=2)
        ax3.vlines(param3_max, ymin3, ymax3, lw=4, color='gray',
                linestyles='dashed')

    ax3.tick_params(axis='x', labelsize=18, bottom=True, top=True,
                    labelbottom=True, labeltop=False)
    ax3.tick_params(axis='y', labelsize=18, left=True, right=True,
                    labelleft=True, labelright=False)
    ax3.set_ylabel(param3.label, fontproperties=font_prop)
    ax3.set_xlabel('Normalized QPO phase', fontproperties=font_prop)

    ax3.set_ylim(ymin3, ymax3)
    ax3.set_xlim(0.0, 1.01)
    y_maj_loc = ax3.get_yticks()
    # ax3.set_yticks(y_maj_loc[0:-1])
    # y_maj_loc = [0.460, 0.464, 0.468, 0.472]
    # ax3.set_yticks(y_maj_loc)
    y_min_mult = 0.25 * (y_maj_loc[1] - y_maj_loc[0])
    yLocator3 = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax3.yaxis.set_minor_locator(yLocator3)
    ax3.set_xticks(np.arange(0, 1.05, 0.25))
    ax3.xaxis.set_minor_locator(xLocator)

    fig.subplots_adjust(hspace=0.00)
    # 	plt.show()
    plt.savefig(plot_name)
    plt.close()
    subprocess.call(['open', plot_name])
    subprocess.call(['cp', plot_name, "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])


################################################################################
def make_plots_var2(plot_name, num_spectra, param1, param2):
    """
    Making plots of fit parameters vs phase for two co-varying parameters,
    fracsctr and another.

    """

    font_prop = font_manager.FontProperties(size=18)
    xLocator = MultipleLocator(0.05)  ## loc of minor ticks on y-axis

    phase = np.arange(num_spectra) / 23.646776
    tinybins = np.arange(-0.01, 1.03, 0.01)
    param1_max = -1
    param2_max = -1

    if param1.funcfit is not None:
        param1_max = tinybins[np.argmax(param1.funcfit)]
    if param2.funcfit is not None:
        param2_max = tinybins[np.argmax(param2.funcfit)]

    print "Plot file: %s" % plot_name

    fig = plt.figure(figsize=(12, 10))

    #####################
    ## Plotting FracSctr
    #####################

    ax1 = fig.add_subplot(211)
    ax1.errorbar(phase, param1.value, yerr=[param1.neg_err, param1.pos_err], \
            lw=0, ecolor='green', marker='.', ms=10, mec='green', mfc='green', \
            elinewidth=2, capsize=2)
    if param1.funcfit is not None:
        ax1.plot(tinybins, param1.funcfit, c='black', lw=2)
        ax1.vlines(param1_max, 0.1, 0.35, lw=1, linestyles='dashed')
    #   ax1.vlines(param1_max - 3.2181e-03, 0.12, 0.24, lw=1, linestyles='dotted')
    #   ax1.vlines(param1_max + 3.2181e-03, 0.12, 0.24, lw=1, linestyles='dotted')

    ax1.set_ylim(0.1, 0.35)
    ax1.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=False, labeltop=False)
    ax1.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax1.set_ylabel(param1.label, fontproperties=font_prop)

    y_maj_loc = ax1.get_yticks()
    ax1.set_yticks(y_maj_loc[1:])
    # yLocator1 = MultipleLocator(.01)  ## loc of minor ticks on y-axis
    # ax1.yaxis.set_minor_locator(yLocator1)
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
    if param2.funcfit is not None:
        ax2.plot(tinybins, param2.funcfit, c='black', lw=2)
        # ax2.vlines(param2_max, 2.2, 2.8, lw=1, linestyles='dashed')
    #   ax2.vlines(param1_max - 3.2181e-03, 0.12, 0.24, lw=1, linestyles='dotted')
    #   ax2.vlines(param1_max + 3.2181e-03, 0.12, 0.24, lw=1, linestyles='dotted')
    # ax2.set_ylim(2.2, 2.8)
    ax2.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=True, labeltop=False)
    ax2.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax2.set_ylabel(param2.label, fontproperties=font_prop)

    ax2.set_xlabel('Normalized QPO phase', fontproperties=font_prop)

    y_maj_loc = ax2.get_yticks()
    ax2.set_yticks(y_maj_loc[0:-1])
    # yLocator2 = MultipleLocator(.001)  ## loc of minor ticks on y-axis
    # ax2.yaxis.set_minor_locator(yLocator2)
    ax2.set_xticks(np.arange(0, 1.05, 0.25))
    ax2.xaxis.set_minor_locator(xLocator)

    fig.subplots_adjust(hspace=0.00)
    # 	plt.show()
    plt.savefig(plot_name)
    plt.close()
    subprocess.call(['open', plot_name])


################################################################################
def make_var_plots(plot_file, num_spectra, var_pars):
    """
    Making plots of fit parameters vs phase for mutiple co-varying parameters.

    """

    font_prop = font_manager.FontProperties(size=18)
    xLocator = MultipleLocator(0.05)  ## loc of minor ticks on y-axis

    phase = np.arange(num_spectra) / 23.646776
    tinybins = np.arange(-0.02, 1.02, 0.01)

    ## So that the plotted x-value is the MIDDLE of the 'bin', with and error of
    ## the width of the bin.
    plusphase = (phase[1]-phase[0])/2.0
    phase_err = np.repeat(plusphase, num_spectra)
    phase += plusphase
    tinybins += plusphase

    colours=['red', 'green', 'blue', 'orange', 'gray']

    # print("Plot file: %s" % plot_name)

    ax_list = []

    fig = plt.figure(figsize=(10, 12), tight_layout=True, dpi=300)
    i=1
    for param in var_pars:

        param_max = -1

        if param.funcfit is not None:
            param_max = tinybins[np.argmax(param.funcfit)]
            # print param.par_name
            # print "\tPhase of max parameter value:", param_max

        temp = (np.max(param.value)-np.min(param.value)) * 0.25
        ymax = np.max(param.value+param.pos_err) + temp
        ymin = np.min(param.value-param.pos_err) - temp

        if i == 1:
            ax = fig.add_subplot(len(var_pars), 1, i)
        else:
            ax = fig.add_subplot(len(var_pars), 1, i, sharex=ax_list[0])

        ax.errorbar(phase, param.value, xerr=phase_err, yerr=[param.neg_err,
                param.pos_err], lw=2, color=colours[i-1], drawstyle='steps-mid',
                marker='.', ms=10, mec=colours[i-1], mfc=colours[i-1],
                ecolor=colours[i-1], elinewidth=2, capsize=0)

        if param.funcfit is not None:
            # rect = patches.Rectangle((param_max-param.phase_err, ymin),
            #         2*param.phase_err, ymax-ymin, color='pink', ec="none",
            #         alpha=0.5)
            # ax.add_patch(rect)
            phase_width = np.round(2.0*param.phase_err/0.002, decimals=1)
            # print phase_width
            # print 2*param.phase_err
            ax.plot(tinybins, param.funcfit, c='black', lw=2)
            ax.vlines(param_max, ymin, ymax, lw=phase_width, color='gray',
                    linestyles='dashed')

        ax.tick_params(axis='x', labelsize=18, bottom=True, top=True,
                labelbottom=False, labeltop=False)
        ax.tick_params(axis='y', labelsize=18, left=True, right=True,
                labelleft=True, labelright=False)
        ax.set_ylabel(param.label, fontproperties=font_prop)

        ax.set_ylim(ymin, ymax)
        ax.set_xlim(0.0, 1.01)
        y_maj_loc = ax.get_yticks()
        y_min_mult = 0.25 * (y_maj_loc[1] - y_maj_loc[0])
        yLocator = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
        ax.yaxis.set_minor_locator(yLocator)

        ax_list.append(ax)
        i += 1

    ax_list[-1].set_xlabel('Normalized QPO phase', fontproperties=font_prop)
    ax_list[-1].set_xticks(np.arange(0, 1.05, 0.25))
    ax_list[-1].xaxis.set_minor_locator(xLocator)
    ax_list[-1].tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                labelbottom=True, labeltop=False)

    fig.subplots_adjust(hspace=0.00)
    # 	plt.show()
    plt.savefig(plot_file)
    plt.close()
    subprocess.call(['open', plot_file])
    # subprocess.call(['cp', plot_file, "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])


################################################################################
def fit_function(t, p):
    """
    Computing a function to fit to the SED parameter variations.

    Parameters
    ----------
    t : np.array of floats
        Time steps for the fit function.
    p : np.array of floats
        The function parameters.

    Returns
    -------
    np.array of floats
        A function fit to the data at steps t with parameters p.
    """
    return p[0] * np.sin(2.0 * np.pi * t + p[1]) + p[2] * \
           sawtooth(2.0 * np.pi * t + p[3]) + p[4]
           # scipy.signal.sawtooth(2.0 * np.pi * t + p[3]) + p[4]


################################################################################
def function_residuals(p, data, data_err, t):
    """
    Getting the residual of the data with the current fit function. Dividing by
    error bar to weight it appropriately like in weighted least squares, e.g.
    S. Vaughan 2013 eqn 6.12 (modified -- not squaring because according to
    scipy.optimize.leastsq documentation, it will square for me to compute the
    real residual and covariance matrix (which will also make it look exactly
    like eqn 6.12))

    Params
    ------
    p : np.array of floats
        The function parameters.
    data : np.array of floats
        The data we want to fit to; in this case, the list of fit parameters per
        energy spectra over QPO phase.
    data_err : np.array of floats
        The error on the data.
    t : np.array of floats
        Time steps for the fitting function.

    Return
    ------
    np.array of floats
        A modified weighted least squared residual of the current function fit
        with the data. From S. Vaughan 2013 eqn 6.12 and scipy.optimize.leastsq
        documentation.
    """
    residual = np.abs(data - fit_function(t, p)) / data_err
    return residual


################################################################################
def get_phase(parameter, num_spectra):
    """
    Fitting a function to an energy spectrum fit parameter to determine the
    phase of the parameter changes.

    Params
    ------
    parameter : Parameter object
        The spectral energy distribution parameter.
    num_spectra : int
        The number of energy spectra in use (the number of energy spectra per
        QPO phase).

    Return
    ------
    Parameter object
        The energy spectrum parameter, with funcfit, phase, and phase_err
        assigned.

    """
    t = np.arange(num_spectra) / 23.646776
    p = [1.0, 0.0, 1.0, 0.0, np.mean((np.min(parameter.value), \
            np.max(parameter.value)))]  ## Amplitude, phase shift, mean

    parameter.error = np.mean((parameter.pos_err, parameter.neg_err), axis=0)

    p_best = leastsq(function_residuals, p, args=(parameter.value, \
            parameter.error, t), full_output=1)
    # print "P best:", p_best
    best_fit = p_best[0]
    print("\tBest fit: %s" % str(best_fit))

    # plt.errorbar(t, parameter.value, xerr=None, yerr=parameter.error)
    # plt.plot(t, fit_function(t, best_fit))
    # plt.xlim(0,1)
    # plt.show()

    ## Error on phase from S. Vaughan 2013 p 168
    bonus_matrix = p_best[1]  ## A Jacobian approximation to the Hessian of the
            ## least squares objective function.
    resid_var = np.var(function_residuals(best_fit, parameter.value, \
            parameter.error, t), ddof=1)
    ## As outlined in the scipy.optimize.leastsq documentation, multiply the
    ## bonus matrix by the variance of the residuals to get the covariance
    ## matrix.

    # print("\t %s" % str(resid_var))
    # print("\t %s" % str(bonus_matrix))
    cov_matrix = bonus_matrix * resid_var

    parameter.best_fit = best_fit
    parameter.funcfit = fit_function(np.arange(-0.02, 1.02, 0.01), best_fit)
    parameter.phase = best_fit[1] / (2.0 * np.pi)
    parameter.phase_err = np.sqrt(cov_matrix[1][1]) / (2.0 * np.pi)

    return parameter


################################################################################
def main(log_file, write_func, mod_string):
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

    mod_components, num_spectra = read_log_file(log_file)

    # print "Number of spectra:", num_spectra

    ######################################################################
    ## Computing the phase of the best-fit function and phase difference
    ######################################################################

    var_pars = []

    for component in mod_components:
        print component.par_name, "mean:", np.mean(component.value)
        if component.varying:
            component = get_phase(component, num_spectra)
            var_pars.append(component)
        # print "%s %s phase: %.4f +- %.4f" % (parameter.mod_name, \
        #         parameter.par_name, parameter.phase, parameter.phase_err)

    if write_func != "":
        print("Writing function parameters to: %s" % write_func)
        with open(write_func, 'w') as out:
            out.write("%s    " % (mod_string))
            for component in mod_components:
                if component.varying:
                    # print component.best_fit
                    out.write("[%.4e,%.4e,%.4e,%.4e,%.4e]    " % \
                            (component.best_fit[0], component.best_fit[1], \
                             component.best_fit[2], component.best_fit[3], \
                             component.best_fit[4]))
                else:
                    out.write("%.4e    " % component.value[0])
            out.write("\n")


    #################################################################
    ## Make plot showing the varying parameters and print phase diff
    #################################################################

    if len(var_pars) >= 1:
        plot_name = log_file.replace('.log', '.eps')
        make_var_plots(plot_name, num_spectra, var_pars)
    else:
        print("\tNo parameters are varying. Nothing to plot.")


################################################################################
if __name__ == '__main__':

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python multifit_plots.py log_file",\
            description="Reads an XSPEC log file and makes plots of varying "\
            "spectral energy distribution parameters as a function of QPO or "\
            "pulse phase, and fits a function to those parameter variations.")

    parser.add_argument('log_file', help="The XSPEC log file, with chatter set"\
            " to 4.")

    parser.add_argument('--mod_string', dest='mod_string', default="", \
            help="The energy spectral model as a string with no spaces. []")

    parser.add_argument('-w', '-W', dest='write_func', default="",
            help="Specifies a text file to write the best-fitting function "\
            "parameters to. []")

    args = parser.parse_args()

    main(args.log_file, args.write_func, args.mod_string)

################################################################################
