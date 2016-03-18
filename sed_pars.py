#!/usr/bin/env python
"""
Classes used in the 'spectral' part of spectral-timing analysis.
"""
import numpy as np

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2015-2016"

class Parameter(object):
    def __init__(self, mod_name=None, label=None, par_name=None):
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
        self.phi_max = None
        self.phi_max_err = None
        self.phi_min = None
        self.phi_min_err = None

    def __str__(self):
        return "%s" % (self.label)

    def __repr__(self):
        return self.__str__()


class Phabs(object):
    def __init__(self):
        self.mod_name = "phabs"
        self.nH = Parameter(self.mod_name,
                label=r"phabs $N_H$ ($\times 10^{22}$)", par_name="nH")

    def __str__(self):
        return "%s" % (self.mod_name)

    def __repr__(self):
        return self.__str__()


class Simpl(object):
    def __init__(self):
        self.mod_name = "simpl "
        self.Gamma = Parameter(self.mod_name, label="PL photon index",
                par_name="Gamma")
        self.FracSctr = Parameter(self.mod_name, label="PL $F_{scatt}$ (norm)",
                par_name="FracSctr")
        self.UpScOnly = Parameter(self.mod_name, label="simpl UpScOnly",
                par_name="UpScOnly")

    def __str__(self):
        return "%s" % (self.mod_name)

    def __repr__(self):
        return self.__str__()


class Simpler(object):
    def __init__(self):
        self.mod_name = "simpler"
        self.UpScOnly = Parameter(self.mod_name, label="simpler UpScOnly",
                par_name="UpScOnly")
        self.Gamma = Parameter(self.mod_name, label="PL photon index",
                par_name="Gamma")
        self.FracSctr = Parameter(self.mod_name, label="PL scattering fraction",
                par_name="FracSctr")


    def __str__(self):
        return "%s" % (self.mod_name)

    def __repr__(self):
        return self.__str__()


class Nthcomp(object):
    def __init__(self):
        self.mod_name = "nthComp"
        self.Gamma = Parameter(self.mod_name, label="nthComp photon index",
                par_name="Gamma")
        self.kT_e = Parameter(self.mod_name,
                label=r"nthComp temp (keV)", par_name="kT_e")
        self.kT_bb = Parameter(self.mod_name,
                label=r"nthComp: kT$_{bb}$ (keV)", par_name="kT_bb")
        self.inp_type = Parameter(self.mod_name, label="nthComp inp type",
                par_name="inp_type")
        self.Redshift = Parameter(self.mod_name, label="nthComp Redshift",
                par_name="Redshift")
        self.norm = Parameter(self.mod_name, label="nthComp norm",
                par_name="norm")

    def __str__(self):
        return "%s" % (self.mod_name)

    def __repr__(self):
        return self.__str__()


class Diskbb(object):
    def __init__(self):
        self.mod_name = "diskbb"
        self.Tin = Parameter(self.mod_name, label="Disk temp (keV)", \
                par_name="Tin")
        self.norm = Parameter(self.mod_name, label="Disk norm",
                par_name="norm")

    def __str__(self):
        return "%s" % (self.mod_name)

    def __repr__(self):
        return self.__str__()


class Diskpbb(object):
    def __init__(self):
        self.mod_name = "diskpbb"
        self.Tin = Parameter(self.mod_name, label=r"DiskP temp (keV)",
                par_name="Tin")
        self.p = Parameter(self.mod_name, label="DiskP exponent", \
                par_name=" p ")
        self.norm = Parameter(self.mod_name, label="DiskP norm",
                par_name="norm")

    def __str__(self):
        return "%s" % (self.mod_name)

    def __repr__(self):
        return self.__str__()


class Bbodyrad(object):
    def __init__(self):
        self.mod_name = "bbodyrad"
        self.kT = Parameter(self.mod_name, label="BB temp (keV)",
                par_name="kT")
        self.norm = Parameter(self.mod_name, label="BB norm",
                par_name="norm")

    def __str__(self):
        return "%s" % (self.mod_name)

    def __repr__(self):
        return self.__str__()


class Gaussian(object):
    def __init__(self):
        self.mod_name = "gaussian"
        self.LineE = Parameter(self.mod_name, label="Line energy (keV)",
                par_name="LineE")
        self.Sigma = Parameter(self.mod_name, label="Equivalent width (keV)",
                par_name="Sigma")
        self.norm = Parameter(self.mod_name, label="Line norm",
                par_name="norm")

    def __str__(self):
        return "%s" % (self.mod_name)

    def __repr__(self):
        return self.__str__()


class Diskline(object):
    def __init__(self):
        self.mod_name = "diskline"
        self.LineE = Parameter(self.mod_name, label=r"Diskline energy (keV)",
                par_name="LineE")
        self.norm = Parameter(self.mod_name, label=r"Diskline norm",
                par_name="norm")

    def __str__(self):
        return "%s" % (self.mod_name)

    def __repr__(self):
        return self.__str__()


class Cutoffpl(object):
    def __init__(self):
        self.mod_name = "cutoffpl"
        self.PhoIndex = Parameter(self.mod_name, label=r"Cut-off PL photon index",
                par_name="PhoIndex")
        self.norm = Parameter(self.mod_name, label=r"Cut-off PL norm",
                par_name="norm")

    def __str__(self):
        return "%s" % (self.mod_name)

    def __repr__(self):
        return self.__str__()