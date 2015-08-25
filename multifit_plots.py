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
import os.path

__author__ = "Abigail Stevens, A.L.Stevens@uva.nl"

"""
Reads an XSPEC log file and makes plots of varying fit parameters as a function
of QPO or pulse phase.

2015

"""

class Parameter(object):
    def __init__(self, mod_name, label, par_name):
        """
        label : str
            What to label a graph (done up nice and pretty).
        par_name : str
            The exact string to search for and match in an xspec log file.
        """
        self.mod_name = mod_name
        self.label = label
        self.par_name = par_name
        self.value = np.asarray([])
        self.error = 0
        self.lo_v = 0
        self.hi_v = 0
        self.pos_err = 0
        self.neg_err = 0
        self.par_num = np.asarray([])
        self.varying = False
        self.sinefit = None
        self.phase = None
        self.phase_err = None


class Phabs(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.nH = Parameter(mod_name, r"phabs: nH ($\times 10^{22}$)", "nH")

class Simpler(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.Gamma = Parameter(mod_name, r"simpler: $\Gamma$", "Gamma")
        self.FracSctr = Parameter(mod_name, r"simpler: FracSctr", "FracSctr")
        self.UpScOnly = Parameter(mod_name, r"simpler: UpScOnly", "UpScOnly")

class Nthcomp(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.Gamma = Parameter(mod_name, r"nthComp: $\Gamma$", "Gamma")
        self.kT_e = Parameter(mod_name, r"nthComp: kT$_{e}$ (keV)", "kT_e")
        self.kT_bb = Parameter(mod_name, r"nthComp: kT$_{bb}$ (keV)", "kT_bb")
        self.inp_type = Parameter(mod_name, "nthComp: inp type", "inp_type")
        self.Redshift = Parameter(mod_name, "nthComp: Redshift", "Redshift")
        self.norm = Parameter(mod_name, "nthComp: norm", "norm")

class Diskbb(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.Tin = Parameter(mod_name, r"diskbb: T$_{in}$ (kev)", "Tin")
        self.norm = Parameter(mod_name, r"diskbb: norm", "norm")

class Diskpbb(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.Tin = Parameter(mod_name, r"diskpbb: T$_{in}$ (kev)", "Tin")
        self.p = Parameter(mod_name, r"diskpbb: p", " p ")
        self.norm = Parameter(mod_name, r"diskpbb: norm", "norm")

class Bbodyrad(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.kT = Parameter(mod_name, r"bbodyrad: kT (keV)", "4   bbodyrad   kT")
        self.norm = Parameter(mod_name, r"bbodyrad: norm", "4   bbodyrad   norm")

class Gaussian(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.LineE = Parameter(mod_name, r"gaussian: LineE (keV)", "LineE")
        self.Sigma = Parameter(mod_name, r"gaussian: $\sigma$ (keV)", "Sigma")
        self.norm = Parameter(mod_name, r"gaussian: norm", "norm")

class Diskline(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.LineE = Parameter(mod_name, r"diskline: LineE (keV)", "LineE")
        self.norm = Parameter(mod_name, r"diskline: norm", "norm")

class Cutoffpl(object):
    def __init__(self, mod_name):
        self.mod_name = mod_name
        self.PhoIndex = Parameter(mod_name, r"cutoffpl: PhoIndex", "PhoIndex")
        self.norm = Parameter(mod_name, r"cutoffpl: norm", "norm")

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
        List of model parameters used.

    """

    with open(log_file, 'r') as f:
        for line in f:

            for component in mod_components:
                if component.mod_name in line:
                    if component.par_name in line and component.varying:
                        temp = float(line.split()[-1])
                        component.pos_err = np.repeat(temp, num_spectra)
                        component.neg_err = np.repeat(temp, num_spectra)

            # ## Reading simpler parameters
            # if "simpler" in line:
            #     if "Gamma" in line and simpler.Gamma.varying:
            #         temp = float(line.split()[-1])
            #         simpler.Gamma.pos_err = np.repeat(temp, num_spectra)
            #         simpler.Gamma.neg_err = np.repeat(temp, num_spectra)
            #
            #     elif "FracSctr" in line and simpler.FracSctr.varying:
            #         temp = float(line.split()[-1])
            #         simpler.FracSctr.pos_err = np.repeat(temp, num_spectra)
            #         simpler.FracSctr.neg_err = np.repeat(temp, num_spectra)
            #
            # ## Reading in bbodyrad parameters
            # elif "bbodyrad" in line and " 4 " in line:
            #     if "kT" in line and bbrad.kT.varying:
            #         temp = float(line.split()[-1])
            #         bbrad.kT.pos_err = np.repeat(temp, num_spectra)
            #         bbrad.kT.neg_err = np.repeat(temp, num_spectra)
            #         if bbrad.norm.varying is False:
            #             break
            #
            #     elif "norm" in line and bbrad.norm.varying:
            #         temp = float(line.split()[-1])
            #         bbrad.norm.pos_err = np.repeat(temp, num_spectra)
            #         bbrad.norm.neg_err = np.repeat(temp, num_spectra)
            #         return simpler, bbrad
            #
            # elif "diskbb" in line:
            #     if "Tin" in line and diskbb.Tin.varying:
            #         temp = float(line.split()[-1])
            #         diskbb.Tin.pos_err = np.repeat(temp, num_spectra)
            #         diskbb.Tin.neg_err = np.repeat(temp, num_spectra)
            #     elif "norm" in line and diskbb.norm.varying:
            #         temp = float(line.split()[-1])
            #         diskbb.norm.pos_err = np.repeat(temp, num_spectra)
            #         diskbb.norm.neg_err = np.repeat(temp, num_spectra)
            #
            # elif "diskpbb" in line:
            #     if "Tin" in line and diskpbb.Tin.varying:
            #         temp = float(line.split()[-1])
            #         diskpbb.Tin.pos_err = np.repeat(temp, num_spectra)
            #         diskpbb.Tin.neg_err = np.repeat(temp, num_spectra)
            #     elif " p " in line and diskpbb.p.varying:
            #         temp = float(line.split()[-1])
            #         diskpbb.p.pos_err = np.repeat(temp, num_spectra)
            #         diskpbb.p.neg_err = np.repeat(temp, num_spectra)
            #     elif "norm" in line and diskpbb.norm.varying:
            #         temp = float(line.split()[-1])
            #         diskpbb.norm.pos_err = np.repeat(temp, num_spectra)
            #         diskpbb.norm.neg_err = np.repeat(temp, num_spectra)

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
    par_nums = np.asarray([])
    lo_v = np.asarray([])
    hi_v = np.asarray([])
    pos_err = np.asarray([])
    neg_err = np.asarray([])

    for line in error_file_o:
        if "XSPEC: quit" not in line:
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
    mod_components


    """

    # temp_mask = np.array([], dtype=bool)
    # if nth.Gamma.varying:
    #     for elt in par_nums:
    #         if elt in nth.Gamma.par_num:
    #             temp_mask = np.append(temp_mask, True)
    #         else:
    #             temp_mask = np.append(temp_mask, False)
    #     nth.Gamma.pos_err = pos_err[temp_mask]
    #     nth.Gamma.neg_err = np.abs(neg_err[temp_mask])
    #     nth.Gamma.lo_v = lo_v[temp_mask]
    #     nth.Gamma.hi_v = hi_v[temp_mask]
    #
    # temp_mask = np.array([], dtype=bool)
    # if nth.norm.varying:
    #     for elt in par_nums:
    #         if elt in nth.norm.par_num:
    #             temp_mask = np.append(temp_mask, True)
    #         else:
    #             temp_mask = np.append(temp_mask, False)
    #     nth.norm.pos_err = pos_err[temp_mask]
    #     nth.norm.neg_err = np.abs(neg_err[temp_mask])
    #     nth.norm.lo_v = lo_v[temp_mask]
    #     nth.norm.hi_v = hi_v[temp_mask]

    for component in mod_components:
        temp_mask = np.array([], dtype=bool)

        if component.varying:
            for elt in par_nums:
                if elt in component.par_num:
                    temp_mask = np.append(temp_mask, True)
                else:
                    temp_mask = np.append(temp_mask, False)
            component.pos_err = pos_err[temp_mask]
            component.Gamma.neg_err = np.abs(neg_err[temp_mask])
            component.Gamma.lo_v = lo_v[temp_mask]
            component.Gamma.hi_v = hi_v[temp_mask]

    # if simpler.Gamma.varying:
    #     for elt in par_nums:
    #         if elt in simpler.Gamma.par_num:
    #             temp_mask = np.append(temp_mask, True)
    #         else:
    #             temp_mask = np.append(temp_mask, False)
    #     simpler.Gamma.pos_err = pos_err[temp_mask]
    #     simpler.Gamma.neg_err = np.abs(neg_err[temp_mask])
    #     simpler.Gamma.lo_v = lo_v[temp_mask]
    #     simpler.Gamma.hi_v = hi_v[temp_mask]
    #
    # temp_mask = np.array([], dtype=bool)
    # if simpler.FracSctr.varying:
    #     for elt in par_nums:
    #         if elt in simpler.FracSctr.par_num:
    #             temp_mask = np.append(temp_mask, True)
    #         else:
    #             temp_mask = np.append(temp_mask, False)
    #     simpler.FracSctr.pos_err = pos_err[temp_mask]
    #     simpler.FracSctr.neg_err = np.abs(neg_err[temp_mask])
    #     simpler.FracSctr.lo_v = lo_v[temp_mask]
    #     simpler.FracSctr.hi_v = hi_v[temp_mask]
    #
    # temp_mask = np.array([], dtype=bool)
    # if bbrad.kT.varying:
    #     for elt in par_nums:
    #         if elt in bbrad.kT.par_num:
    #             temp_mask = np.append(temp_mask, True)
    #         else:
    #             temp_mask = np.append(temp_mask, False)
    #     bbrad.kT.pos_err = pos_err[temp_mask]
    #     bbrad.kT.neg_err = np.abs(neg_err[temp_mask])
    #     bbrad.kT.lo_v = lo_v[temp_mask]
    #     bbrad.kT.hi_v = hi_v[temp_mask]
    #
    # temp_mask = np.array([], dtype=bool)
    # if diskbb.Tin.varying:
    #     for elt in par_nums:
    #         if elt in diskbb.Tin.par_num:
    #             temp_mask = np.append(temp_mask, True)
    #         else:
    #             temp_mask = np.append(temp_mask, False)
    #     diskbb.Tin.pos_err = pos_err[temp_mask]
    #     diskbb.Tin.neg_err = np.abs(neg_err[temp_mask])
    #     diskbb.Tin.lo_v = lo_v[temp_mask]
    #     diskbb.Tin.hi_v = hi_v[temp_mask]
    #
    # temp_mask = np.array([], dtype=bool)
    # if diskbb.norm.varying:
    #     for elt in par_nums:
    #         if elt in diskbb.norm.par_num:
    #             temp_mask = np.append(temp_mask, True)
    #         else:
    #             temp_mask = np.append(temp_mask, False)
    #     diskbb.norm.pos_err = pos_err[temp_mask]
    #     diskbb.norm.neg_err = np.abs(neg_err[temp_mask])
    #     diskbb.norm.lo_v = lo_v[temp_mask]
    #     diskbb.norm.hi_v = hi_v[temp_mask]
    #
    # temp_mask = np.array([], dtype=bool)
    # if diskpbb.Tin.varying:
    #     for elt in par_nums:
    #         if elt in diskpbb.Tin.par_num:
    #             temp_mask = np.append(temp_mask, True)
    #         else:
    #             temp_mask = np.append(temp_mask, False)
    #     diskpbb.Tin.pos_err = pos_err[temp_mask]
    #     diskpbb.Tin.neg_err = np.abs(neg_err[temp_mask])
    #     diskpbb.Tin.lo_v = lo_v[temp_mask]
    #     diskpbb.Tin.hi_v = hi_v[temp_mask]
    #
    # temp_mask = np.array([], dtype=bool)
    # if diskpbb.p.varying:
    #     for elt in par_nums:
    #         if elt in diskpbb.p.par_num:
    #             temp_mask = np.append(temp_mask, True)
    #         else:
    #             temp_mask = np.append(temp_mask, False)
    #     diskpbb.p.pos_err = pos_err[temp_mask]
    #     diskpbb.p.neg_err = np.abs(neg_err[temp_mask])
    #     diskpbb.p.lo_v = lo_v[temp_mask]
    #     diskpbb.p.hi_v = hi_v[temp_mask]
    #
    # temp_mask = np.array([], dtype=bool)
    # if diskpbb.norm.varying:
    #     for elt in par_nums:
    #         if elt in diskpbb.norm.par_num:
    #             temp_mask = np.append(temp_mask, True)
    #         else:
    #             temp_mask = np.append(temp_mask, False)
    #     diskpbb.norm.pos_err = pos_err[temp_mask]
    #     diskpbb.norm.neg_err = np.abs(neg_err[temp_mask])
    #     diskpbb.norm.lo_v = lo_v[temp_mask]
    #     diskpbb.norm.hi_v = hi_v[temp_mask]

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
    Phabs object
        Description.

    Simpler object
        Description.

    Nthcomp object
        Description.

    Diskbb object
        Description.

    Bbodyrad object
        Description.

    Gaussian object
        Description.

    int
        Number of spectra being simultaneously fit for one QPO phase.

    """

    if not os.path.isfile(log_file) or os.path.getsize(log_file) == 0:
        raise Exception("Log file does not exist or is empty.")

    simpler = Simpler("simpler")
    nth = Nthcomp("nthComp")
    diskbb = Diskbb("diskbb")
    diskpbb = Diskpbb("diskpbb")
    bbrad = Bbodyrad("bbodyrad")
    gauss = Gaussian("gaussian")
    chains = False
    mod_components = [simpler.FracSctr, simpler.Gamma, diskbb.Tin, \
                      diskbb.norm, diskpbb.Tin, diskpbb.p, diskpbb.norm, \
                      bbrad.kT, bbrad.norm, nth.Gamma, nth.norm, gauss.LineE, \
                      gauss.Sigma, gauss.norm]

    #################################################
    ## Reading in parameter values from the log file
    #################################################

    print ""
    with open(log_file, 'r') as f:
        for line in f:

            for component in mod_components:
                if component.mod_name in line and component.par_name in line:
                    print line.split()
                    component.value = np.append(component.value, \
                            float(line.split()[-3]))
                    component.par_num = np.append(component.par_num, \
                            int(line.split()[1]))

            # ## Reading simpler parameters
            # if simpler.mod_name in line:
            #     if simpler.Gamma.par_name in line:
            #         print line.split()
            #         simpler.Gamma.value = np.append(simpler.Gamma.value, \
            #                 float(line.split()[-3]))
            #         simpler.Gamma.par_num = np.append(simpler.Gamma.par_num, \
            #                 int(line.split()[1]))
            #
            #     elif simpler.FracSctr.par_name in line:
            #         simpler.FracSctr.value = np.append(simpler.FracSctr.value, \
            #                 float(line.split()[-3]))
            #         simpler.FracSctr.par_num = np.append(simpler.FracSctr.par_num, \
            #                 int(line.split()[1]))
            #
            # ## Reading nthComp parameters
            # elif nth.mod_name in line:
            #     if nth.Gamma.par_name in line:
            #         nth.Gamma.value = np.append(nth.Gamma.value, \
            #                 float(line.split()[-3]))
            #         nth.Gamma.par_num = np.append(nth.Gamma.par_num, \
            #                 int(line.split()[1]))
            #
            #     elif nth.norm.par_name in line:
            #         nth.norm.value = np.append(nth.norm.value, \
            #                 float(line.split()[5]))
            #         nth.norm.par_num = np.append(nth.norm.par_num, \
            #                 int(line.split()[1]))
            #
            # ## Reading diskbb parameters
            # elif diskbb.mod_name in line:
            #     if diskbb.Tin.par_name in line:
            #         print line.split()
            #         diskbb.Tin.value = np.append(diskbb.Tin.value, \
            #                 float(line.split()[-3]))
            #         diskbb.Tin.par_num = np.append(diskbb.Tin.par_num, \
            #             int(line.split()[1]))
            #
            #     elif diskbb.norm.par_name in line:
            #         diskbb.norm.value = np.append(diskbb.norm.value, \
            #                 float(line.split()[-3]))
            #         diskbb.norm.par_num = np.append(diskbb.norm.par_num, \
            #                 int(line.split()[1]))
            #
            # ## Reading diskpbb parameters
            # elif diskpbb.mod_name in line:
            #     if diskpbb.Tin.par_name in line:
            #         diskpbb.Tin.value = np.append(diskpbb.Tin.value, \
            #                 float(line.split()[-3]))
            #         diskpbb.Tin.par_num = np.append(diskpbb.Tin.par_num, \
            #             int(line.split()[1]))
            #
            #     elif diskpbb.p.par_name in line:
            #         diskpbb.p.value = np.append(diskpbb.p.value, \
            #                 float(line.split()[-3]))
            #         diskpbb.p.par_num = np.append(diskpbb.p.par_num, \
            #                 int(line.split()[1]))
            #
            #     elif diskpbb.norm.par_name in line:
            #         diskpbb.norm.value = np.append(diskpbb.norm.value, \
            #                 float(line.split()[-3]))
            #         diskpbb.norm.par_num = np.append(diskpbb.norm.par_num, \
            #                 int(line.split()[1]))
            #
            # ## Reading bbodyrad parameters
            # elif bbrad.mod_name in line:
            #     if bbrad.kT.par_name in line:
            #         bbrad.kT.value = np.append(bbrad.kT.value, \
            #                 float(line.split()[-3]))
            #         bbrad.kT.par_num = np.append(bbrad.kT.par_num, \
            #                 int(line.split()[1]))
            #
            #     elif bbrad.norm.par_name in line:
            #         bbrad.norm.value = np.append(bbrad.norm.value, \
            #                 float(line.split()[-3]))
            #         bbrad.norm.par_num = np.append(bbrad.norm.par_num, \
            #                 int(line.split()[1]))
            #
            # ## Reading gaussian parameters
            # elif gauss.mod_name in line:
            #     if gauss.LineE.par_name in line:
            #         gauss.LineE.value = np.append(gauss.LineE.value, \
            #                 float(line.split()[-3]))
            #         gauss.LineE.par_num = np.append(gauss.LineE.par_num, \
            #                 int(line.split()[1]))
            #
            #     elif gauss.Sigma.par_name in line:
            #         gauss.Sigma.value = np.append(gauss.Sigma.value, \
            #                 float(line.split()[-3]))
            #         gauss.Sigma.par_num = np.append(gauss.Sigma.par_num, \
            #                 int(line.split()[1]))
            #
            #     elif gauss.norm.par_name in line:
            #         gauss.norm.value = np.append(gauss.norm.value, \
            #                 float(line.split()[-3]))
            #         gauss.norm.par_num = np.append(gauss.norm.par_num, \
            #                 int(line.split()[1]))

            if "Parameter" in line and "Confidence Range" in line:
                chains = True
                par_nums, lo_v, hi_v, neg_err, pos_err = read_chain(f)

    num_spectra = np.amax([len(nth.Gamma.value), len(simpler.Gamma.value)])

    #############################################################
    ## Determine if parameter varies across QPO phase or is tied
    #############################################################

    for component in mod_components:
        print component.mod_name, component.par_name
        if len(component.value) > 1:
            if component.value[0] != component.value[1]:
                component.varying = True

    # if len(simpler.Gamma.value) > 1:
    #     if simpler.Gamma.value[0] != simpler.Gamma.value[1]:
    #         simpler.Gamma.varying = True
    #
    # if len(simpler.FracSctr.value) > 1:
    #     if simpler.FracSctr.value[0] != simpler.FracSctr.value[1]:
    #         simpler.FracSctr.varying = True
    #
    # if len(nth.Gamma.value) > 1:
    #     if nth.Gamma.value[0] != nth.Gamma.value[1]:
    #         nth.Gamma.varying = True
    #
    # if len(nth.norm.value) > 1:
    #     if nth.norm.value[0] != nth.norm.value[1]:
    #         nth.norm.varying = True
    #
    # if len(diskbb.Tin.value) > 1:
    #     if diskbb.Tin.value[0] != diskbb.Tin.value[1]:
    #         diskbb.Tin.varying = True
    #
    # if len(diskbb.norm.value) > 1:
    #     if diskbb.norm.value[0] != diskbb.norm.value[1]:
    #         diskbb.norm.varying = True
    #
    # if len(bbrad.kT.value) > 1:
    #     if bbrad.kT.value[0] != bbrad.kT.value[1]:
    #         bbrad.kT.varying = True
    #
    # if len(bbrad.norm.value) > 1:
    #     if bbrad.norm.value[0] != bbrad.norm.value[1]:
    #         bbrad.norm.varying = True
    #
    # if len(gauss.LineE.value) > 1:
    #     if gauss.LineE.value[0] != gauss.LineE.value[1]:
    #         gauss.LineE.varying = True
    #
    # if len(gauss.Sigma.value) > 1:
    #     if gauss.Sigma.value[0] != gauss.Sigma.value[1]:
    #         gauss.Sigma.varying = True
    #
    # if len(gauss.norm.value) > 1:
    #     if gauss.norm.value[0] != gauss.norm.value[1]:
    #         gauss.norm.varying = True

    ############################################################################
    ## Reading in errors from 'chain' part of log, or bad errors from normal log
    ############################################################################

    var_pars = []
    print "VarPars"
    for component in mod_components:
        if component.varying:
            var_pars.append(component)
            print component.mod_name, component.par_name

    if len(var_pars) == 0:
        raise Exception("No parameters vary with QPO phase in this log file.")
        exit()

    if chains:
        var_pars = get_chain_errors(var_pars, par_nums, lo_v, hi_v, neg_err, \
                pos_err)
    if not chains:
        print "Using fake errors from log file. Need to run error analysis!!"
        var_pars = get_logfile_errors(log_file, num_spectra, \
                var_pars)

    return var_pars, num_spectra


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
    if param1.sinefit is not None:
        param1_max = tinybins[np.argmax(param1.sinefit)]
        print param1_max
    if param2.sinefit is not None:
        param2_max = tinybins[np.argmax(param2.sinefit)]
        print param2_max
    if param3.sinefit is not None:
        param3_max = tinybins[np.argmax(param3.sinefit)]
        print param3_max


    print "Plot file: %s" % plot_name

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
    if param1.sinefit is not None:
        # rect1 = patches.Rectangle((param1_max-param1.phase_err, ymin1),
        #         2*param1.phase_err, ymax1-ymin1, color='pink', ec="none",
        #         alpha=0.3)
        # ax1.add_patch(rect1)
        ax1.plot(tinybins, param1.sinefit, c='black', lw=2)
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
    if param2.sinefit is not None:
        # rect2 = patches.Rectangle((param2_max-param2.phase_err, ymin2),
        #         2*param2.phase_err, ymax2-ymin2, color='pink', ec="none",
        #         alpha=0.3)
        # ax2.add_patch(rect2)
        ax2.plot(tinybins, param2.sinefit, c='black', lw=2)
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
    if param3.sinefit is not None:
        # rect3 = patches.Rectangle((param3_max-param3.phase_err, ymin3),
        #         2*param3.phase_err, ymax3-ymin3, color='pink', ec="none",
        #         alpha=0.3)
        # ax3.add_patch(rect3)
        ax3.plot(tinybins, param3.sinefit, c='black', lw=2)
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
def make_var_plots(plot_name, num_spectra, var_pars):
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

    print "Plot file: %s" % plot_name

    ax_list = []

    fig = plt.figure(figsize=(10, 12), tight_layout=True, dpi=300)
    i=1
    for param in var_pars:
        param_max = -1

        if param.sinefit is not None:
            param_max = tinybins[np.argmax(param.sinefit)]
            print param.par_name
            print "\tPhase of max parameter value:", param_max

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

        if param.sinefit is not None:
            # rect = patches.Rectangle((param_max-param.phase_err, ymin),
            #         2*param.phase_err, ymax-ymin, color='pink', ec="none",
            #         alpha=0.5)
            # ax.add_patch(rect)
            phase_width = np.round(2.0*param.phase_err/0.002, decimals=1)
            print phase_width
            # print 2*param.phase_err
            ax.plot(tinybins, param.sinefit, c='black', lw=2)
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
    plt.savefig(plot_name)
    plt.close()
    subprocess.call(['open', plot_name])
    # subprocess.call(['cp', plot_name, "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])


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
    if param2.sinefit is not None:
        ax2.plot(tinybins, param2.sinefit, c='black', lw=2)
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
    Parameter object
        Description.

    """
    t = np.arange(num_spectra) / 23.646776
    p = [1.0, 0.0, np.mean((np.min(sed_parameter.value), \
            np.max(sed_parameter.value)))]  ## Amplitude, phase shift, mean

    sed_parameter.error = np.mean((sed_parameter.pos_err, \
            sed_parameter.neg_err), axis=0)

    p_best = leastsq(sine_residuals, p, args=(sed_parameter.value, \
            sed_parameter.error, t), full_output=1)
    # print "P best:", p_best
    best_fit = p_best[0]
    print "\tBest fit:", best_fit

    # plt.errorbar(t, sed_parameter.value, xerr=None, yerr=sed_parameter.error)
    # plt.plot(t, sinewave(t, best_fit))
    # plt.xlim(0,1)
    # plt.show()

    ## Error on phase from S. Vaughan 2013 p 168
    bonus_matrix = p_best[1]  ## A Jacobian approximation to the Hessian of the
            ## least squares objective function.
    resid_var = np.var(sine_residuals(best_fit, sed_parameter.value, \
            sed_parameter.error, t), ddof=1)
    ## As outlined in the scipy.optimize.leastsq documentation, multiply the
    ## bonus matrix by the variance of the residuals to get the covariance
    ## matrix.

    # print "\t", resid_var
    # print "\t", bonus_matrix
    cov_matrix = bonus_matrix * resid_var

    sed_parameter.sinefit = sinewave(np.arange(-0.02, 1.02, 0.01), best_fit)
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

    var_pars, num_spectra = read_log_file(log_file)

    print "Number of spectra:", num_spectra

    for param in var_pars:
        if param.varying:
            print param.par_name, "mean:", np.mean(param.value)
            param = get_phase(param, num_spectra)
            print "%s %s phase: %.4f +- %.4f" % (param.mod_name, \
                    param.par_name, param.phase, param.phase_err)

    ######################################################################
    ## Computing the phase of the best-fit sine wave and phase difference
    ######################################################################

    # if simpler.FracSctr.varying:
    #     simpler.FracSctr = get_phase(simpler.FracSctr, num_spectra)
    #     print "simpler FracSctr phase: %.4f +- %.4f" % (simpler.FracSctr.phase,\
    #             simpler.FracSctr.phase_err)
    #
    # if simpler.Gamma.varying:
    #     simpler.Gamma = get_phase(simpler.Gamma, num_spectra)
    #     print "simpler Gamma phase: %.4f +- %.4f" % (simpler.Gamma.phase, \
    #             simpler.Gamma.phase_err)
    #
    # if nth.norm.varying:
    #     nth.norm = get_phase(nth.norm, num_spectra)
    #     print "nth norm phase: %.4f +- %.4f" % (nth.norm.phase, \
    #             nth.norm.phase_err)
    #
    # if nth.Gamma.varying:
    #     nth.Gamma = get_phase(nth.Gamma, num_spectra)
    #     print "nth Gamma phase: %.4f +- %.4f" % (nth.Gamma.phase, \
    #             nth.Gamma.phase_err)
    #
    # if diskbb.Tin.varying:
    #     diskbb.Tin = get_phase(diskbb.Tin, num_spectra)
    #     print "diskbb Tin phase: %.4f +- %.4f" % (diskbb.Tin.phase, \
    #             diskbb.Tin.phase_err)
    #
    # if diskbb.norm.varying:
    #     diskbb.norm = get_phase(diskbb.norm, num_spectra)
    #     print "diskbb norm phase: %.4f +- %.4f" % (diskbb.norm.phase, \
    #             diskbb.norm.phase_err)
    #
    # if bbrad.kT.varying:
    #     bbrad.kT = get_phase(bbrad.kT, num_spectra)
    #     print "bbrad kT phase: %.4f +- %.4f" % (bbrad.kT.phase, \
    #             bbrad.kT.phase_err)
    #
    # if bbrad.norm.varying:
    #     bbrad.norm = get_phase(bbrad.norm, num_spectra)
    #     print "bbrad norm phase: %.4f +- %.4f" % (bbrad.norm.phase, \
    #             bbrad.norm.phase_err)

    #################################################################
    ## Make plot showing the varying parameters and print phase diff
    #################################################################

    plot_name = log_file.replace('.log', '.eps')

    make_var_plots(plot_name, num_spectra, var_pars)

    # if len(var_pars) == 1:
    #     pass
    # elif len(var_pars) == 2:
    #     make_plots_var2(plot_name, num_spectra, var_pars[0], var_pars[1])
    # elif len(var_pars) == 3:
    #     make_plots_var3(plot_name, num_spectra, var_pars[0], var_pars[1], \
    #             var_pars[2])

    # if simpler.FracSctr.varying and simpler.Gamma.varying and bbrad.kT.varying:
    #     print "Phase diff: %.4f" % \
    #             np.abs(simpler.FracSctr.phase - bbrad.kT.phase)
    #     make_plots_var3(plot_name, num_spectra, simpler.FracSctr, \
    #             simpler.Gamma, bbrad.kT)
    #
    # elif simpler.FracSctr.varying and simpler.Gamma.varying and \
    #         diskbb.Tin.varying:
    #     print "Phase diff: %.4f" % \
    #             np.abs(simpler.FracSctr.phase - diskbb.Tin.phase)
    #     make_plots_var3(plot_name, num_spectra, simpler.FracSctr, \
    #             simpler.Gamma, diskbb.Tin)
    #
    # elif simpler.FracSctr.varying and simpler.Gamma.varying and \
    #         bbrad.norm.varying:
    #     print "Phase diff: %.4f" % \
    #             np.abs(simpler.FracSctr.phase - bbrad.norm.phase)
    #     make_plots_var3(plot_name, num_spectra, simpler.FracSctr, simpler.Gamma,
    #             bbrad.norm)
    #
    # elif simpler.FracSctr.varying and bbrad.kT.varying and \
    #         bbrad.norm.varying:
    #     make_plots_var3(plot_name, num_spectra, simpler.FracSctr, bbrad.kT, \
    #             bbrad.norm)
    #
    # elif simpler.FracSctr.varying and simpler.Gamma.varying:
    #     print "Phase diff: %.4f" % \
    #             np.abs(simpler.FracSctr.phase - simpler.Gamma.phase)
    #     make_plots_var2(plot_name, num_spectra, simpler.FracSctr, simpler.Gamma)
    #
    # elif simpler.FracSctr.varying and bbrad.kT.varying:
    #     print "Phase diff: %.4f" % \
    #             np.abs(simpler.FracSctr.phase - bbrad.kT.phase)
    #     make_plots_var2(plot_name, num_spectra, simpler.FracSctr, bbrad.kT)
    #
    # elif simpler.FracSctr.varying and bbrad.norm.varying:
    #     make_plots_var2(plot_name, num_spectra, simpler.FracSctr, bbrad.norm)
    #
    # elif nth.norm.varying and nth.Gamma.varying and bbrad.kT.varying:
    #     make_plots_var3(plot_name, num_spectra, nth.norm, \
    #             nth.Gamma, bbrad.kT)
    #
    # elif nth.norm.varying and nth.Gamma.varying and bbrad.norm.varying:
    #     make_plots_var3(plot_name, num_spectra, nth.norm, bbrad.norm, \
    #             nth.Gamma)
    #
    # elif nth.norm.varying and nth.Gamma.varying:
    #     make_plots_var2(plot_name, num_spectra, nth.norm, nth.Gamma)
    #
    # elif nth.norm.varying and bbrad.kT.varying:
    #     make_plots_var2(plot_name, num_spectra, nth.norm, bbrad.kT)
    #
    # elif nth.norm.varying and bbrad.norm.varying:
    #     make_plots_var2(plot_name, num_spectra, nth.norm, bbrad.norm)


################################################################################
if __name__ == '__main__':

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python multifit_plots.py log_file",\
            description="Reads an XSPEC log file and makes plots of varying "\
            "fit parameters as a function of QPO or pulse phase. All arguments"\
            " are required.")

    parser.add_argument('log_file', help="The XSPEC log file, with chatter set"\
            "to 4.")

    args = parser.parse_args()

    main(args.log_file)

################################################################################
