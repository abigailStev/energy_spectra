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

################################################################################
def read_log_file(log_file):
    """
    Reads the XSPEC log file, with the expectation that chatter is set to 4.

    Params
    ------
    log_file : string
        Description.

    Returns
    -------
    np.array of floats
        simpl Gamma
    np.array of floats
        Error on simpl Gamma
    np.array of floats
        simpl FracSctr
    np.array of floats
        Error on simpl FracSctr
    np.array of floats
        diskBB Tin
    np.array of floats
        Error on diskBB Tin
    np.array of floats
        diskBB normalization
    np.array of floats
        Error on diskBB normalization
    np.array of floats
        gaussian line E
    np.array of floats
        Error on gaussian line E
    np.array of floats
        gaussian Sigma
    np.array of floats
        Error on gaussian Sigma
    np.array of floats
        gaussian normalization
    np.array of floats
        Error on gaussian normalization
    int
        Number of spectra for one QPO phase.

    """
    gamma = np.asarray([])
    gamma_err = np.asarray([])
    fracsctr = np.asarray([])
    fracsctr_err = np.asarray([])
    bb_t = np.asarray([])
    bb_t_err = np.asarray([])
    bb_norm = np.asarray([])
    bb_norm_err = np.asarray([])
    line_e = np.asarray([])
    line_e_err = np.asarray([])
    sigma = np.asarray([])
    sigma_err = np.asarray([])
    e_norm = np.asarray([])
    e_norm_err = np.asarray([])

    ###########################
    ## Reading in the log file
    ###########################

    print ""
    with open(log_file, 'r') as f:
        for line in f:
            ## Reading simPL parameters
            if "simpl" in line:
                if "Gamma" in line:
                    gamma = np.append(gamma, float(line.split()[5]))
                    gamma_err = np.append(gamma_err, float(line.split()[7]))

                elif "FracSctr" in line:
                    fracsctr = np.append(fracsctr, float(line.split()[5]))
                    fracsctr_err = np.append(fracsctr_err,
                                             float(line.split()[7]))

            ## Reading diskBB parameters
            elif "diskbb" in line:
                if "Tin" in line:
                    bb_t = np.append(bb_t, float(line.split()[6]))
                    bb_t_err = np.append(bb_t_err, float(line.split()[8]))

                elif "norm" in line:
                    bb_norm = np.append(bb_norm, float(line.split()[5]))
                    # bb_norm_err = np.append(bb_norm_err, float(line.split()[7]))


            ## Reading gaussian parameters
            elif "gaussian" in line:
                if "LineE" in line:
                    line_e = np.append(line_e, float(line.split()[6]))
                    line_e_err = np.append(line_e_err, float(line.split()[8]))

                elif "Sigma" in line:
                    sigma = np.append(sigma, float(line.split()[6]))
                    sigma_err = np.append(sigma_err, float(line.split()[8]))

                elif "norm" in line:
                    e_norm = np.append(e_norm, float(line.split()[5]))
                    e_norm_err = np.append(e_norm_err, float(line.split()[7]))

    ###########################################################
    ## Assert statements, to ensure the log was read correctly
    ###########################################################

    assert len(gamma) == len(gamma_err), "ERROR: Issue with reading log file. "\
            "Length of Gamma != length of error on Gamma."
    assert len(fracsctr) == len(fracsctr_err), "ERROR: Issue with reading log "\
            "file. Length of FracSctr != length of error on FracSctr."
    assert len(bb_t) == len(bb_t_err), "ERROR: Issue with reading log file. "\
            "Length of Tin != length of error on Tin."
    # assert len(bb_norm) == len(bb_norm_err), "ERROR: Issue with reading log "\
    #         "file. Length of diskbb norm != length of error on diskbb norm."
    assert len(line_e) == len(line_e_err), "ERROR: Issue with reading log "\
            "file. Length of LineE != length of error on LineE."
    assert len(sigma) == len(sigma_err), "ERROR: Issue with reading log file. "\
            "Length of Sigma != length of error on Sigma."
    assert len(e_norm) == len(e_norm_err), "ERROR: Issue with reading log "\
            "file. Length of LineE norm != length of error on LineE norm."

    assert len(gamma) == len(fracsctr), "ERROR: Issue with reading log file. "\
            "Length of gamma != length of FracSctr."
    assert len(gamma) == len(bb_t), "ERROR: Issue with reading log file. "\
            "Length of gamma != length of Tin."
    assert len(gamma) == len(bb_norm), "ERROR: Issue with reading log file. "\
            "Length of gamma != length of diskbb norm."
    assert len(gamma) == len(line_e), "ERROR: Issue with reading log file. "\
            "Length of gamma != length of LineE."
    assert len(gamma) == len(sigma), "ERROR: Issue with reading log file. "\
            "Length of gamma != length of Sigma."
    assert len(gamma) == len(e_norm), "ERROR: Issue with reading log file. "\
            "Length of gamma != length of LineE norm."

    num_spectra = len(gamma)

    return gamma, gamma_err, fracsctr, fracsctr_err, bb_t, bb_t_err, bb_norm, \
        bb_norm_err, line_e, line_e_err, sigma, sigma_err, e_norm, e_norm_err, \
        num_spectra


################################################################################
def make_plots_var3(plot_name, num_spectra, fracsctr, fracsctr_err, gamma, \
                    gamma_err, param, param_err, fs_sinefit, param_sinefit,
                    param_title):
    """
    Making plots of fit parameters vs phase for three co-varying parameters.

    """

    font_prop = font_manager.FontProperties(size=20)
    xLocator = MultipleLocator(0.05)  ## loc of minor ticks on y-axis

    phase = np.arange(num_spectra) / 23.646776
    tinybins = np.arange(-0.01, 1.03, 0.01)
    fs_max = tinybins[np.argmax(fs_sinefit)]
    param_max = tinybins[np.argmax(param_sinefit)]


    print "Plot file: %s" % plot_name

    fig = plt.figure(figsize=(12, 10))

    #####################
    ## Plotting FracSctr
    #####################

    ax1 = fig.add_subplot(311)

    ax1.plot(tinybins, fs_sinefit, c='black', lw=2)
    ax1.errorbar(phase, fracsctr, yerr=fracsctr_err, lw=0, ecolor='red', \
                 marker='.', ms=10, mec='red', mfc='red', elinewidth=2,
                 capsize=2)
    ax1.vlines(fs_max, 0.12, 0.24, lw=1, linestyles='dashed')
    # ax1.vlines(fs_max - 3.2181e-03, 0.12, 0.24, lw=1, linestyles='dotted')
    # ax1.vlines(fs_max + 3.2181e-03, 0.12, 0.24, lw=1, linestyles='dotted')
    ax1.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=False, labeltop=False)
    ax1.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    # ax1.set_ylabel('Scattering \nfraction', fontproperties=font_prop)
    ax1.set_ylabel(r'simpl: FracSctr', fontproperties=font_prop)

    ax1.set_xlim(-0.01 , 1.02)
    y_maj_loc = ax1.get_yticks()
    ax1.set_yticks(y_maj_loc[1:])
    yLocator1 = MultipleLocator(.01)  ## loc of minor ticks on y-axis
    ax1.yaxis.set_minor_locator(yLocator1)
    ax1.xaxis.set_minor_locator(xLocator)

    #####################
    ## Plotting diskBB T
    #####################

    ax2 = fig.add_subplot(312, sharex=ax1)

    ax2.plot(tinybins, param_sinefit, c='black', lw=2)
    ax2.errorbar(phase, param, yerr=param_err, lw=0, ecolor='green', \
                 marker='.', ms=10, mec='green', mfc='green', elinewidth=2,
                 capsize=2)
    ax2.vlines(param_max, 0.798, 0.810, lw=1, linestyles='dashed')
    # ax2.vlines(param_max - 4.3787e-03, 0.798, 0.810, lw=1, linestyles='dotted')
    # ax2.vlines(param_max + 4.3787e-03, 0.798, 0.810, lw=1, linestyles='dotted')
    ax2.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=False, labeltop=False)
    ax2.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax2.set_ylabel(param_title, fontproperties=font_prop)

    ax2.set_xlim(-0.01, 1.02)
    y_maj_loc = ax2.get_yticks()
    ax2.set_yticks(y_maj_loc[1:-1])
    y_min_mult = 0.5 * (y_maj_loc[1] - y_maj_loc[0])
    yLocator2 = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax2.yaxis.set_minor_locator(yLocator2)
    ax2.set_xticks(np.arange(0, 1.05, 0.25))
    ax2.xaxis.set_minor_locator(xLocator)

    ##################
    ## Plotting Gamma
    ##################

    ax3 = fig.add_subplot(313, sharex=ax1)

    ax3.errorbar(phase, gamma, yerr=gamma_err, lw=0, ecolor='blue', \
                 marker='.', ms=10, mec='blue', mfc='blue', elinewidth=2,
                 capsize=2)
    ax3.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=True, labeltop=False)
    ax3.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax3.set_xlabel('Normalized QPO phase', fontproperties=font_prop)
    ax3.set_ylabel(r'simpl: $\Gamma$', fontproperties=font_prop)

    # 	ax3.set_xlim(bins[0]-0.5,bins[-1]+0.5)
    y_maj_loc = ax3.get_yticks()
    ax3.set_yticks(y_maj_loc[0:-1])
    ax3.set_xlim(-0.01, 1.02)
    y_min_mult = 0.2 * (y_maj_loc[1] - y_maj_loc[0])
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
def make_plots_var2(plot_name, num_spectra, fracsctr, fracsctr_err, param, \
                    param_err, param_title):
    """
    Making plots of fit parameters vs phase for two co-varying parameters,
    fracsctr and another.

    """

    font_prop = font_manager.FontProperties(size=20)
    xLocator = MultipleLocator(0.05)  ## loc of minor ticks on y-axis

    phase = np.arange(num_spectra) / 23.646776

    print "Plot file: %s" % plot_name

    fig = plt.figure(figsize=(12, 10))

    #####################
    ## Plotting FracSctr
    #####################

    ax1 = fig.add_subplot(211)
    ax1.errorbar(phase, fracsctr, yerr=fracsctr_err, lw=0, ecolor='red', \
                 marker='.', ms=10, mec='red', mfc='red', elinewidth=2,
                 capsize=2)
    ax1.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=False, labeltop=False)
    ax1.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    # ax1.set_ylabel('Scattering \nfraction', fontproperties=font_prop)
    ax1.set_ylabel(r'simpl: FracSctr', fontproperties=font_prop)

    # 	ax1.set_xlim(phase[0]-0.01, phase[-1]+0.01)
    y_maj_loc = ax1.get_yticks()
    ax1.set_yticks(y_maj_loc[1:])
    yLocator1 = MultipleLocator(.01)  ## loc of minor ticks on y-axis
    ax1.yaxis.set_minor_locator(yLocator1)
    ax1.set_xticks(np.arange(0, 1.05, 0.25))
    ax1.xaxis.set_minor_locator(xLocator)

    #####################
    ## Plotting diskBB T
    #####################

    ax2 = fig.add_subplot(212, sharex=ax1)
    ax2.errorbar(phase, param, yerr=param_err, lw=0, ecolor='green', \
                 marker='.', ms=10, mec='green', mfc='green', elinewidth=2,
                 capsize=2)
    # 	ax2.set_xlim(bins[0]-0.5,bins[-1]+0.5)
    ax2.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
                    labelbottom=True, labeltop=False)
    ax2.tick_params(axis='y', labelsize=18, left=True, right=True, \
                    labelleft=True, labelright=False)
    ax2.set_ylabel(param_title, fontproperties=font_prop)
    ax2.set_xlabel('Normalized phase', fontproperties=font_prop)

    # 	ax1.set_xlim(phase[0]-0.01, phase[-1]+0.01)
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
    # subprocess.call(['open', plot_name])


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
def get_phase(sed_parameters, sed_err, num_spectra):
    """
    Fitting a sine wave to an energy spectra fit parameter to determine the
    phase of the parameter changes.

    Params
    ------
    sed_parameters : np.array of floats
        The spectral energy distribution (energy spectra) parameters fitted by
        XSPEC over the QPO phase.
    sed_err : np.array of floats
        The error on each data point of sed_parameters.
    num_spectra : int
        The number of energy spectra in use, i.e. the length of
        sed_fit_parameters.

    Return
    ------
    float
        The phase of the best-fit sine wave.
    float
        The error on the phase of the best-fit sine wave.
    """
    t = np.arange(num_spectra) / 23.646776
    p = [1.0, 0.0, np.mean(sed_parameters)]  ## Amplitude, phase shift, mean

    p_best = leastsq(sine_residuals, p, args=(sed_parameters, sed_err, t), \
            full_output=1)
    print "P best:", p_best
    best_fit = p_best[0]
    print "Best fit:", best_fit

    plt.errorbar(t, sed_parameters, xerr=None, yerr=sed_err)
    plt.plot(t, sinewave(t, best_fit))
    plt.xlim(0,1)
    # plt.show()

    ## Error on phase from S. Vaughan 2013 p 168
    bonus_matrix = p_best[1]  ## A Jacobian approximation to the Hessian of the
            ## least squares objective function.
    resid_var = np.var(sine_residuals(best_fit, sed_parameters, sed_err, t), \
            ddof=1)
    ## As outlined in the scipy.optimize.leastsq documentation, multiply the
    ## bonus matrix by the variance of the residuals to get the covariance
    ## matrix.
    print "Bonus matrix:", bonus_matrix
    print "Resid var:", resid_var
    cov_matrix = bonus_matrix * resid_var

    tiny_bins = np.arange(-0.01, 1.03, 0.01)
    smooth_sine_fit = sinewave(tiny_bins, best_fit)

    return best_fit[1] / (2.0 * np.pi), np.sqrt(cov_matrix[1][1]) / \
            (2.0 * np.pi), smooth_sine_fit


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

    gamma, gamma_err, fracsctr, fracsctr_err, bb_t, bb_t_err, bb_norm, \
            bb_norm_err, line_e, line_e_err, sigma, sigma_err, e_norm, \
            e_norm_err, num_spectra = read_log_file(log_file)

    print "Number of spectra:", num_spectra

    # print "norm(BB) = ", bb_norm[1]
    # print np.mean(bb_norm)
    # print np.min(bb_norm)
    # print np.max(bb_norm)

    #####################################################
    ## If parameter is tied, assign it the correct error
    #####################################################

    var_gamma = True
    var_fracsctr = True
    var_bb_t = True
    var_bb_norm = True
    var_line_e = True
    var_sigma = True
    var_e_norm = True

    if gamma[0] == gamma[1] and gamma_err[0].is_integer():
        gamma_err = np.repeat(gamma_err[0], num_spectra)
        var_gamma = False

    if fracsctr[0] == fracsctr[1] and fracsctr_err[0].is_integer():
        fracsctr_err = np.repeat(fracsctr_err[0], num_spectra)
        var_fracsctr = False

    if (len(bb_t) == 1 or bb_t[0] == bb_t[1]) and bb_t[0].is_integer():
        bb_t_err = np.repeat(bb_t_err[0], num_spectra)
        var_bb_t = False

    if (bb_norm[0] == bb_norm[1] or len(bb_norm) == 1) and \
            (np.size(bb_norm_err) == 0 or bb_norm_err[0].is_integer()):
        if np.size(bb_norm_err) == 0:
            bb_norm_err = [0.0]
        bb_norm_err = np.repeat(bb_norm_err[0], num_spectra)
        var_bb_norm = False

    if line_e[0] == line_e[1] and line_e_err[0].is_integer():
        line_e_err = np.repeat(line_e_err[0], num_spectra)
        var_line_e = False

    if sigma[0] == sigma[1] and sigma_err[0].is_integer():
        sigma_err = np.repeat(sigma_err[0], num_spectra)
        var_sigma = False

    if e_norm[0] == e_norm[1] and e_norm_err[0].is_integer():
        e_norm_err = np.repeat(e_norm_err[0], num_spectra)
        var_e_norm = False

    ######################################################################
    ## Computing the phase of the best-fit sine wave and phase difference
    ######################################################################

    fracsctr_phase, fracsctr_phase_err, fracsctr_sinefit = get_phase(fracsctr,\
            fracsctr_err, num_spectra)

    print "simpl fracsctr phase: %.4e +- %.4e" % (fracsctr_phase, \
            fracsctr_phase_err)

    if var_bb_t:
        bb_t_phase, bb_t_phase_err, bb_t_sinefit = get_phase(bb_t, bb_t_err, \
                num_spectra)
        phase_2 = bb_t_phase
        print "diskBB Tin phase: %.4e +- %.4e" % (bb_t_phase, bb_t_phase_err)

    if var_bb_norm:
        print "\tIN HERE"
        bb_norm_phase, bb_norm_phase_err, bb_norm_sinefit = get_phase(bb_norm, \
                bb_norm_err, num_spectra)
        phase_2 = bb_norm_phase
        print "diskBB norm. phase: %.4e +- %.4e" % (bb_norm_phase, \
                bb_norm_phase_err)

    phase_diff = fracsctr_phase - phase_2
    if phase_diff < 0:
        phase_diff += 1.0
    elif phase_diff > 1.0:
        phase_diff -= 1.0

    print "Phase diff: %.4f" % phase_diff

    ###################
    ## Make the plots!
    ###################

    plot_name = log_file.replace('.log', '.eps')

    if var_fracsctr and var_gamma and var_bb_t:
        make_plots_var3(plot_name, num_spectra, fracsctr, fracsctr_err, gamma, \
                gamma_err, bb_t, bb_t_err, fracsctr_sinefit, bb_t_sinefit, \
                r"diskBB: T$_{in}$ (kev)")
    elif var_fracsctr and var_gamma and var_bb_norm:
        make_plots_var3(plot_name, num_spectra, fracsctr, fracsctr_err, gamma, \
                gamma_err, bb_norm, bb_norm_err, fracsctr_sinefit, \
                bb_norm_sinefit, "diskBB: norm.")
    elif var_fracsctr and var_gamma:
        make_plots_var2(plot_name, num_spectra, fracsctr, fracsctr_err, gamma, \
                gamma_err, r"simpl: $\Gamma$")
    elif var_fracsctr and var_bb_t:
        make_plots_var2(plot_name, num_spectra, fracsctr, fracsctr_err, bb_t, \
                bb_t_err, r"diskBB: T$_{in}$ (kev)")
    elif var_fracsctr and var_bb_norm:
        make_plots_var2(plot_name, num_spectra, fracsctr, fracsctr_err, bb_norm,
                bb_norm_err, "diskBB: norm.")


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
