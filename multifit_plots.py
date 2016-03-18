#!/usr/bin/env python
"""
Reads an XSPEC log file and makes plots of varying spectral parameters as a
function of QPO or pulse phase. Fits the changing spectral parameters with a
function and gets the 'phase' of each parameter variation.

"""
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
from datetime import datetime
import sed_pars
import matplotlib.colors as colors

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2015-2016"


################################################################################
def get_logfile_errors(log_file, num_spectra, mod_parameters):
    """
    Reads (fake) errors on the varying parameters from the XSPEC log file.
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

    mod_parameters : list of Parameter objects
        The model parameters used.

    Returns
    -------
    mod_parameters
        List of model parameters used, with the errors assigned.
    """

    with open(log_file, 'r') as f:
        for line in f:

            for parameter in mod_parameters:
                if parameter.mod_name in line and parameter.par_name in line \
                        and parameter.varying:
                    temp = line.split()[-1].strip()
                    if temp == "frozen":
                        parameter.pos_err = np.zeros(num_spectra)
                        parameter.neg_err = np.zeros(num_spectra)
                    else:
                        parameter.pos_err = np.repeat(float(temp), num_spectra)
                        parameter.neg_err = np.repeat(float(temp), num_spectra)

    return mod_parameters


################################################################################
def read_chain(error_file_o):
    """
    Reads in the MCMC chain from XSPEC's built-in MCMC error estimation.

    Paramters
    ---------
    error_file_o : str
        The open file object of the log file with chain errors at the end.

    Returns
    -------
    par_nums : np.array of ints
        1-D array of the XSPEC parameter number of each varying parameter.

    lo_v : np.array of floats
        1-D array of the parameter value representing the lower bound of the
        error bar.

    hi_v : np.array of floats
        1-D array of the parameter value representing the upper bound of the
        error bar.

    neg_err : np.array of floats
        1-D array of the negative error bars as deviations from the parameter
        values.

    pos_err : np.array of floats
        1-D array of the positive error bars as deviations from the parameter
        values.
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
def get_chain_errors(mod_parameters, par_nums, lo_v, hi_v, neg_err, pos_err):
    """
    Gets the error on the varying spectral parameters from the MCMC chain.

    Paramters
    ---------
    mod_parameters : list of Parameter objects
        The parameters of the spectral model.

    par_nums : np.array of int
        1-D array of the XSPEC parameter number of each varying parameter.

    lo_v : np.array of floats
        1-D array of the parameter value representing the lower bound of the
        error bar.

    hi_v : np.array of floats
        1-D array of the parameter value representing the upper bound of the
        error bar.

    neg_err : np.array of floats
        1-D array of the negative error bars as deviations from the parameter
        values.

    pos_err : np.array of floats
        1-D array of the positive error bars as deviations from the parameter
        values.

    Returns
    -------
    var_pars : list of Parameter objects
        A 1-D list of the untied spectral parameters that vary with QPO phase.
    """

    # print type(par_nums)
    # print type(mod_parameters)
    # print type(lo_v)
    # print type(hi_v)
    # print type(neg_err)
    # print type(pos_err)

    for parameter in mod_parameters:
        temp_mask = np.array([], dtype=bool)

        if parameter.varying:
            for elt in par_nums:
                if elt in parameter.par_num:
                    temp_mask = np.append(temp_mask, True)
                else:
                    temp_mask = np.append(temp_mask, False)
            parameter.pos_err = pos_err[temp_mask]
            parameter.neg_err = np.abs(neg_err[temp_mask])
            parameter.lo_v = lo_v[temp_mask]
            parameter.hi_v = hi_v[temp_mask]

    return mod_parameters


################################################################################
def read_log_file(log_file, quiet=True):
    """
    Reads the XSPEC log file and assigns parameters to Parameter objects, with
    the expectation that chatter is set to 4.

    Parameters
    ----------
    log_file : str
        Full path of the XSPEC log file, with chatter set to 4. Assuming the
        spectral models listed below are the only ones used (or the only
        interesting ones).

    quiet : bool
        If True, suppresses printing to the screen.

    Returns
    -------
    mod_parameters : list of Parameter objects
        A 1-D list of all spectral parameters for the model.

    num_spectra : int
        Number of spectra being simultaneously fit for one QPO phase.
    """

    if not os.path.isfile(log_file) or os.path.getsize(log_file) == 0:
        print log_file
        raise Exception("Log file does not exist or is empty.")

    chains = False
    mod_parameters = [sed_pars.Phabs().nH, sed_pars.Simpler().Gamma,
            sed_pars.Simpler().FracSctr, sed_pars.Simpler().UpScOnly,
            sed_pars.Simpl().Gamma, sed_pars.Simpl().FracSctr,
            sed_pars.Simpl().UpScOnly, sed_pars.Diskbb().Tin,
            sed_pars.Diskbb().norm, sed_pars.Diskpbb().Tin,
            sed_pars.Diskpbb().p, sed_pars.Diskpbb().norm,
            sed_pars.Bbodyrad().kT, sed_pars.Bbodyrad().norm,
            sed_pars.Gaussian().LineE, sed_pars.Gaussian().Sigma,
            sed_pars.Gaussian().norm]

    #################################################
    ## Reading in parameter values from the log file
    #################################################

    with open(log_file, 'r') as f:
        for line in f:
            for parameter in mod_parameters:
                if parameter.mod_name in line and parameter.par_name in line:
                    if "frozen" in line:
                        parameter.value = np.append(parameter.value,
                            float(line.split()[-2]))
                        parameter.par_num = np.append(parameter.par_num,
                            int(line.split()[1]))
                    else:
                        parameter.value = np.append(parameter.value,
                                float(line.split()[-3]))
                        parameter.par_num = np.append(parameter.par_num,
                                int(line.split()[1]))

            if "Parameter" in line and "Confidence Range" in line:
                chains = True
                par_nums, lo_v, hi_v, neg_err, pos_err = read_chain(f)


    #############################################################
    ## Delete parameters if they're not used/present
    ## Determine if parameter varies across QPO phase or is tied
    ## Assign zero error to parameters
    #############################################################

    # num_spectra = np.amax([len(nth.Gamma.value), len(simpler.Gamma.value)])
    unused_parameters = []
    num_spectra = 1

    for parameter in mod_parameters:
        # print parameter.mod_name, parameter.par_name
        if len(parameter.value) > 1:
            # print parameter.value[0], parameter.value[1], parameter.value[3]
            if parameter.value[0] != parameter.value[1] or \
                    parameter.value[7] != parameter.value[0]:
                parameter.varying = True
                num_spectra = len(parameter.value)
        elif len(parameter.value) == 0:
            unused_parameters.append(parameter)

        # parameter.pos_err = np.zeros(len(parameter.value))
        # parameter.neg_err = np.zeros(len(parameter.value))
        # parameter.lo_v = parameter.value
        # parameter.hi_v = parameter.value

    for elt in unused_parameters:
        mod_parameters.remove(elt)

    ############################################################################
    ## Reading in errors from 'chain' part of log, or bad errors from normal log
    ############################################################################

    var_pars = 0
    # print "VarPars"
    for parameter in mod_parameters:
        if parameter.varying:
            var_pars += 1
            # print parameter.mod_name, parameter.par_name

    if var_pars == 0:
        raise Exception("No parameters vary with QPO phase in this log file.")
        exit()

    if chains:
        mod_parameters = get_chain_errors(mod_parameters, par_nums, lo_v, hi_v,
                neg_err, pos_err)

    if not chains:
        if not quiet:
            print("Using fake errors from log file. Need to run error "\
                    "analysis!!")
        mod_parameters = get_logfile_errors(log_file, num_spectra, \
                mod_parameters)

    return mod_parameters, num_spectra


################################################################################
def make_var_plots(plot_file, num_spectra, var_pars, quiet=False, title=" "):
    """
    Making plots of spectral parameters vs phase for multiple co-varying
    parameters.

    Parameters
    ----------
    plot_file : str
        The full path of the file name to save the plot to.

    num_spectra : int
        The number of spectra that were co-fit.

    var_pars : list of Parameter objects
        A 1-D list of the spectral parameters that vary with QPO phase.

    quiet : bool
        If True, will not open the plot made. [False]
    """

    # title="phabs$\\times\\,$(simpler*diskbb+gauss)"

    font_prop = font_manager.FontProperties(size=18)
    xLocator = MultipleLocator(0.05)  ## loc of minor ticks on x-axis

    phase = np.arange(num_spectra) / 23.5
    if num_spectra == 24:
        tinybins = np.arange(-0.02, 1.02, 0.001)
    elif num_spectra == 47:
        tinybins = np.arange(-0.02, 2.02, 0.001)
    else:
        phase -= 1.0
        tinybins = np.arange(-1.02, 2.02, 0.001)

    ## So that the plotted x-value is the MIDDLE of the 'bin', with and error of
    ## the width of the bin.
    plusphase = (phase[1]-phase[0])/2.0
    phase_err = np.repeat(plusphase, num_spectra)
    phase += plusphase
    tinybins += plusphase

    colours = ['red', 'green', 'blue', colors.cnames['darkorange'],
            colors.cnames['darkorchid'], colors.cnames['darkslategray']]
    ax_list = []

    # fig = plt.figure(figsize=(14, 12), tight_layout=True, dpi=300)  ## for 2 pds
    fig = plt.figure(figsize=(10, 12), tight_layout=True, dpi=300)

    i = 1

    for param in var_pars:

        param_max = -1

        if param.phi_max is not None:
            param_max = param.phi_max
        elif param.funcfit is not None:
            param_max = tinybins[np.argmax(param.funcfit)]
            print param.par_name
            print "\tPhase of max parameter value:", param_max

        temp = (np.max(param.value)-np.min(param.value)) * 0.25
        print np.max(param.value+param.pos_err) + temp
        print np.min(param.value-param.pos_err) - temp

        if i == 1:
            ymax = 2.77
            ymin = 2.16
        elif i == 2:
            ymax = 0.28
            ymin = 0.07
        else:
            temp = (np.max(param.value)-np.min(param.value)) * 0.25
            ymax = np.max(param.value+param.pos_err) + temp
            ymin = np.min(param.value-param.pos_err) - temp
            # print ymax
            # print ymin

        if i == 1:
            ax = fig.add_subplot(len(var_pars), 1, i)
        else:
            ax = fig.add_subplot(len(var_pars), 1, i, sharex=ax_list[0])

        ax.errorbar(phase, param.value, yerr=[param.neg_err,
                param.pos_err], lw=2, color=colours[i-1], drawstyle='steps-mid',
                marker='.', ms=10, mec=colours[i-1], mfc=colours[i-1],
                ecolor=colours[i-1], elinewidth=2, capsize=0)

        if param.funcfit is not None and param.funcfit is not np.nan:
            if param.phi_max_err is not None:
                phase_width = np.round(2.0 * param.phi_max_err / 0.002,
                        decimals=1)
            else:
                phase_width = np.round(2.0 * param.phase_err / 0.002,
                        decimals=1)
            # print phase_width
            # print 2*param.phase_err
            ax.plot(tinybins, param.funcfit, c='black', lw=2)
            ax.vlines(param_max, ymin, ymax, lw=phase_width, color='gray',
                    linestyles='dashed')
            # ax.vlines(param_max+1, ymin, ymax, lw=phase_width, color='gray',
            #         linestyles='dashed')
            # ax.vlines(param_max+2, ymin, ymax, lw=phase_width, color='gray',
            #         linestyles='dashed')

        ax.tick_params(axis='x', labelsize=18, bottom=True, top=True,
                labelbottom=False, labeltop=False)
        ax.tick_params(axis='y', labelsize=18, left=True, right=True,
                labelleft=True, labelright=False)
        ax.set_ylabel(param.label, fontproperties=font_prop)

        ax.set_ylim(ymin, ymax)

        if num_spectra == 24:
            ax.set_xlim(0.0, 1.01)
        elif num_spectra == 47:
            ax.set_xlim(0.0, 2.01)
        else:
            ax.set_xlim(-1.0, 2.01)
            # ax.set_xlim(0, 3.01)

        y_maj_loc = ax.get_yticks()
        y_min_mult = 0.25 * (y_maj_loc[1] - y_maj_loc[0])
        yLocator = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
        ax.yaxis.set_minor_locator(yLocator)

        ax_list.append(ax)
        i += 1

    ax_list[-1].set_xlabel('Normalized QPO phase', fontproperties=font_prop)
    if num_spectra == 24:
        ax_list[-1].set_xticks(np.arange(0, 1.05, 0.25))
    elif num_spectra == 47:
        ax_list[-1].set_xticks(np.arange(0, 2.05, 0.25))
    else:
        ax_list[-1].set_xticks(np.arange(-1.00, 2.05, 0.25))
        # ax_list[-1].set_xticks(np.arange(0, 3.05, 0.25))


    ax_list[-1].xaxis.set_minor_locator(xLocator)
    ax_list[-1].tick_params(axis='x', labelsize=18, bottom=True, top=True, 
                labelbottom=True, labeltop=False)

    ax_list[0].set_title(r'%s' % title, fontproperties=font_prop)
    fig.subplots_adjust(hspace=0.00)
    # 	plt.show()
    plt.savefig(plot_file)
    plt.close()

    if not quiet:
        subprocess.call(['open', plot_file])
    # subprocess.call(['cp', plot_file,
    #         "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])


################################################################################
def fit_function(t, p):
    """
    Computing a function to fit to the spectral parameter variations.

    Parameters
    ----------
    t : np.array of floats
        1-D array of time steps for the fit function.
    p : np.array of floats
        1-D array of the function parameters.

    Returns
    -------
    np.array of floats
        1-D array of the function fit to the data at steps t with parameters p.
    """
    return p[0] * np.sin(2. * np.pi * t + p[1]) + \
           p[2] * np.sin(2. * 2. * np.pi * t + p[3]) + p[4]

################################################################################
def function_residuals(p, data, data_err, t):
    """
    Getting the residual of the data with the current fit function. Dividing by
    error bar to weight it appropriately like in weighted least squares, e.g.
    S. Vaughan 2013 eqn 6.12 (modified -- not squaring because according to
    scipy.optimize.leastsq documentation, it will square for me to compute the
    real residual and covariance matrix (which will also make it look exactly
    like eqn 6.12))

    Parameters
    ------
    p : np.array of floats
        1-D array of the function parameters.

    data : np.array of floats
        1-D array of the data we want to fit to; in this case, the list of
        spectral fit parameters over QPO phase.

    data_err : np.array of floats
        1-D array of the error on the data.

    t : np.array of floats
        1-D array of the time steps for the fitting function.

    Returns
    -------
    np.array of floats
        1-D array of a modified weighted least squared residual of the current
        function fit with the data. From S. Vaughan 2013 eqn 6.12 and
        scipy.optimize.leastsq documentation.
    """
    residual = np.abs(data - fit_function(t, p)) / data_err
    return residual


################################################################################
def get_phase(parameter, num_spectra, quiet):
    """
    Fitting a function to an energy spectrum fit parameter to determine the
    phase of the parameter changes.

    Parameters
    ----------
    parameter : Parameter object
        The spectral energy distribution parameter.

    num_spectra : int
        The number of energy spectra in use (the number of energy spectra per
        QPO phase).

    quiet : bool
        If True, suppresses printing to the screen.

    Returns
    -------
    Parameter object
        The energy spectrum parameter, with funcfit, phase, and phase_err
        assigned.

    """
    t = np.arange(num_spectra) / 23.5
    p = [1., 0., 1., 0., np.mean((np.min(parameter.value), \
            np.max(parameter.value)))]  ## Amplitude, phase shift, mean

    parameter.error = np.mean((parameter.pos_err, parameter.neg_err), axis=0)

    p_best = leastsq(function_residuals, p, args=(parameter.value, \
            parameter.error, t), full_output=1)
    # print "P best:", p_best
    best_fit = p_best[0]

    # if not quiet:
    #     print("\tBest fit: %s" % str(best_fit))

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

    try:
        cov_matrix = bonus_matrix * resid_var
    except TypeError:
        # print("\t %s" % str(resid_var))
        # print("\t %s" % str(bonus_matrix))
        parameter.best_fit = np.nan
        parameter.funcfit = np.nan
        parameter.phase = np.nan
        parameter.phase_err = np.nan
        return parameter

    parameter.best_fit = best_fit
    if num_spectra == 24:
        parameter.funcfit = fit_function(np.arange(-0.02, 1.02, 0.001),
                best_fit)
    elif num_spectra == 47:
        parameter.funcfit = fit_function(np.arange(-0.02, 2.02, 0.001),
                best_fit)
    else:
        parameter.funcfit = fit_function(np.arange(-1.02, 2.02, 0.001),
                best_fit)
    parameter.phase = best_fit[1] / (2.0 * np.pi)
    parameter.phase_err = np.sqrt(cov_matrix[1][1]) / (2.0 * np.pi)

    return parameter


################################################################################
def write_varpars(varying_params, fitfunc_file, num_spec=24):
    """
    NOT DOING THIS ANYMORE. Need to keep all segments so I can compute the
    variance.
    """
    varpar_fits_file = fitfunc_file.replace("_funcfit.txt", "_varpars.txt")

    this_boot = np.zeros(num_spec)
    for single_par in varying_params:
        this_boot = np.vstack((this_boot, single_par.value))
    this_boot = this_boot[1:,]
    # print this_boot
    # print "Boot shape:", np.shape(this_boot)

    if not os.path.isfile(varpar_fits_file):
        to_save = this_boot
    else:
        data = np.loadtxt(varpar_fits_file, delimiter='    ')
        # print "Data shape:", np.shape(data)
        to_save = data + this_boot
        # print to_save
    # print np.shape(to_save)

    np.savetxt(varpar_fits_file, to_save, fmt='%.6e', delimiter='    ')


################################################################################
def determine_varying_parameters(mod_parameters, n_spectra=24, quiet=False):
    """
    Determines which spectral parameters are varying with QPO phase, based on
    the log file from XSPEC.

    Parameters
    ----------
    mod_parameters : np.array of sed_pars.Parameter objects
        1-D array of the parameters of the whole energy spectral model.

    n_spectra : int
        The number of spectra in one QPO cycle being fit simultaneously. [24]

    quiet : bool
        Flag to print output or not; if True, will print each parameter name
        and its mean value. [False]

    Returns
    -------
    var_pars : np.array of sed_pars.Parameter objects
        1-D array of the spectral parameters that vary with QPO phase.

    """
    var_pars = np.array([])
    for parameter in mod_parameters:
        if not quiet:
            print parameter.par_name, "mean:", np.mean(parameter.value)
            print "\t", np.min(parameter.value), np.max(parameter.value)
        if parameter.varying:
            parameter = get_phase(parameter, n_spectra, quiet)
            var_pars = np.append(var_pars, parameter)
        # print "%s %s phase: %.4f +- %.4f" % (parameter.mod_name, \
        #         parameter.par_name, parameter.phase, parameter.phase_err)
    return var_pars


################################################################################
def main(log_file, mod_string="", write_func="", quiet=False):
    """
    Reads the XSPEC log file to get the parameter values, fits a function to the
    varying spectral parameters, writes those fit function parameters to a file
    (if specified), and plots the varying parameters with their best-fit
    function.

    Parameters
    ------
    log_file : str
        The full path of the XSPEC log file, with chatter set to 4, with
        extension '.log'.

    mod_string : str
        The energy spectral model, as a string with no spaces. Used when saving
        information to a table.

    write_func : str
        The full path of the text file to write the best-fitting function
        parameters from the spectral parameter variations to.

    quiet : bool
        If True, suppresses printing to the screen and will not open the plot
        made.

    """
    ##########################################
    ## Reading in the log file to data arrays
    ##########################################

    mod_parameters, num_spectra = read_log_file(log_file, quiet=quiet)

    ######################################################################
    ## Computing the phase of the best-fit function and phase difference
    ######################################################################

    var_pars = determine_varying_parameters(mod_parameters, \
            n_spectra=num_spectra, quiet=quiet)

    if write_func != "":
        # if not quiet:
        #     print("Writing function parameters to: %s" % write_func)
        with open(write_func, 'a') as out:
            out.write("%s    " % (mod_string))
            for parameter in mod_parameters:
                if parameter.varying:
                    # print parameter.best_fit
                    try:
                        out.write("[%.4e,%.4e,%.4e,%.4e,%.4e]    " % \
                                (parameter.best_fit[0], parameter.best_fit[1],
                                 parameter.best_fit[2], parameter.best_fit[3],
                                 parameter.best_fit[4]))
                    except TypeError:
                        out.write("[nan,nan,nan,nan,nan]    ")
                else:
                    out.write("%.4e    " % parameter.value[0])
            out.write("\n")

    #################################################################
    ## Make plot showing the varying parameters and print phase diff
    #################################################################

    if len(var_pars) >= 1:
        plot_name = log_file.replace('.log', '.eps')
        model_name = mod_string.upper()
        model_name = model_name.replace("PHABS*", "PHABS$\\times\\,$")
        model_name = model_name.replace("*", " * ")
        model_name = model_name.replace("+", " + ")
        make_var_plots(plot_name, num_spectra, var_pars, quiet=quiet,
                title=model_name)
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

    parser.add_argument('-q', '--quiet', dest='quiet', action='store_true',
            default=False, help="If present, quiets output and does not open "\
            "the plot.")

    args = parser.parse_args()

    main(args.log_file, args.mod_string, args.write_func, args.quiet)

################################################################################
