import sys
import json
import copy

import ROOT
import xboa.common as common
from xboa.hit import Hit
from xboa.bunch import Bunch
import Configuration
import tracking
import maus_cpp.field
import maus_cpp.polynomial_map
import numpy

#TODO: fit at several differenet values of pz
#TODO: worried that this pushes x', y' to the limit... maybe we are not sensitive to x', y'?

class TMFitter(object):
    """
    Fit a transfer matrix to a transfer matrix
    """
    def __init__(self, config):
        """
        Initialises the fitter
        """
        self.config = config
        # object that handles the MAUS interface (geometry and tracking stuff)
        self.tracking = tracking.Tracking(config, self._get_configuration())
        # the measured transfer matrices, errors and tof12 (set on call to
        # fitter, list of dictionaries as per self.fit(...) docs)
        self.measured_tm_list = []
        # seed for the misalignment and scale factor (starting point for the
        # iterations)
        self.seed = {
          "dx":0.0,
          "dy":0.0,
          "dz":0.0,
          "dxp":0.0,
          "dyp":0.0,
          "scale_factor":self.config.fc_current_seed
        }
        # minuit is only allowed to vary parameters within $seed +- max_delta$
        self.max_delta = {
          "dx":20.0,
          "dy":20.0,
          "dz":50.0,
          "dxp":0.02,
          "dyp":0.02,
          "scale_factor":10.0
        }

        # best guess i.e. output from latest fit iteration
        self.best_guess = copy.deepcopy(self.seed)
        # gets filled with the estimated error on each parameter after fitting
        self.estimated_error = {}
        # maximum number of iterations for each step of the fitter
        self.max_iterations = self.config.fit_max_iterations
        # fitter will attempt to minimise sum chi2 within +- 0.1
        self.resolution = self.config.fit_resolution
        # gets filled with the calculated transfer matrices on each iteration of
        # the fitter (one for each tof)
        self.calculated_tm_list = []
        # index of the upstream virtual plane
        self.plane_us = self.config.fit_plane_us
        # index of the downstream virtual plane
        self.plane_ds = self.config.fit_plane_ds
        # iteration counter
        self._iteration = 0
        # tof12 distance
        # 8222.85 - straight tracks from run 7417
        # 8224.8 - fc data from run 7541 
        self.measured_data = None

    def get_best_guess(self):
        return self.best_guess

    def set_measured_data(self, measured_data):
        self.measured_data = measured_data
        self.measured_tm_list = []
        for data_set in self.measured_data:
            self.measured_tm_list.append(data_set["tm_entries"])

    def fit(self):
        """
        Attempt to fit the transfer matrix, given statistical errors on each
        term in the transfer matrix.
        - transfer_matrices: list of dictionaries. Each dict should have a 
                structure like:
                           {
                              'transfer_matrix':<matrix>, 
                              'matrix_errors':<errors>,
                              'tof12':<tof12>
                           }
                <matrix> should be a list of lists 4x5 containing the transfer
                         matrix.
                <errors> should be a list of lists 4x5 containing the error
                         of the transfer matrix.
                <tof12> should be a float corresponding to the tof12 time for
                        the matrix/errors
        The fit proceeds in a number of steps
        1. Attempt to fit scale factor; looks at the upper quadrant of the
           focussing part of the matrix i.e. terms (0,1; 0,2; 1,2, 1,2,) in the
           transfer matrix.
        2. Attempt to fit x and x'; looks at the horizontal kick terms i.e.
           terms (0,0; 1,0) in the transfer matrix.
        3. Attempt to fit y and y'; looks at the vertical kick terms i.e.
           terms (2,0; 3,0) in the transfer matrix.
        """
        self.best_guess = copy.deepcopy(self.seed)
        self.estimated_error = {}
        self.calculated_tm_list = []
        self._iteration = 0
        TMFitter.fitter = self
        print "Running minuit to find scale_factor"
        scale_factor_terms = []
        for row in range(0, 2):
            for col in range(1,3):
                scale_factor_terms.append((row, col))
        self._do_one_fit(["scale_factor"], scale_factor_terms)
        #print "Running minuit to find x, x', y, y'"
        #self._do_one_fit(["dx", "dxp", "dy", "dyp"], [(0, 0), (1, 0), (2, 0), (3, 0)])
        for i, data_set in enumerate(self.measured_data):
            self.measured_data[i]["fitted_tms"] = self.calculated_tm_list[i]
            self.measured_data[i]["alignment"] = self.best_guess
        print "Fields:"
        print maus_cpp.field.str(True)


    def save_fit(self, suffix):
        for i, data in enumerate(self.measured_data):
            fname = self.config.data_dir+"/tm_"+suffix+"_"+str(i)+".json"
            print "Saving to", fname
            fout = open(fname, "w")
            print >> fout, json.dumps(data)

    def get_station(self, tracks_out, z_position):
        hit_z = [abs(hit['position']['z']-z_position) for hit in tracks_out[0]]
        station = hit_z.index(min(hit_z))
        return station

    def get_tm(self, setup, momentum):
        """
        Calculate the transfer matrix, based on calculating Jacobian numerically
        using particle tracking
        """
        #beta = self.z_tof12/tof12/common.constants['c_light']
        #if beta > 1.:
        #    raise ValueError("tof12 is faster than speed of light")
        #self.tracking.pz = common.pdg_pid_to_mass[13] * beta * 1/(1-beta**2)**0.5
        self.tracking.pz = momentum
        self.tracking.reset_modules()
        for module in ["FCoil_0",  "FCoil_1"]:
            self.tracking.scale_mice_module(module, setup["scale_factor"])
        self.tracking.misalign_mice_module(setup["dx"], setup["dxp"], setup["dy"], setup["dyp"], setup["dz"], "AFC")
        hits_out = self.tracking.do_tracking(self.config.fit_delta_x, self.config.fit_delta_x, self.config.fit_delta_px, self.config.fit_delta_px)
        s_u = self.get_station(hits_out, self.config.z_tku)
        s_d = self.get_station(hits_out, self.config.z_tkd)
        points = [[track[s_u]['position']['x'], track[s_u]['momentum']['x']/track[s_u]['momentum']['z'],
                   track[s_u]['position']['y'], track[s_u]['momentum']['y']/track[s_u]['momentum']['z']] for track in hits_out]
        values = [[track[s_d]['position']['x'], track[s_d]['momentum']['x']/track[s_d]['momentum']['z'],
                   track[s_d]['position']['y'], track[s_d]['momentum']['y']/track[s_d]['momentum']['z']] for track in hits_out]
        print "Tracked from", hits_out[0][s_u]['position']['z'], "to", hits_out[0][s_d]['position']['z'], "(stations", s_u, "to", s_d, ")",
        print "p in:", hits_out[0][s_u]['momentum']['z'], "p out:", hits_out[0][s_d]['momentum']['z'], 
        matrix = maus_cpp.polynomial_map.PolynomialMap.least_squares_fit(points, values, 1).get_coefficients_as_matrix()
        return matrix

    def test_run(self):
        # test run
        TMFitter.fitter = self
        self.minuit = ROOT.TMinuit(0)
        self.current_fit = []
        self.score_function_terms = []
        self.run_once(True)
        for i, data_set in enumerate(self.measured_data):
            self.measured_data[i]["not_fitted_tms"] = self.calculated_tm_list[i]
            self.measured_data[i]["alignment"] = self.best_guess

    def _do_one_fit(self, fit_vars, score_function_terms):
        """
        Run a single iteration of the fitting code
        - fit_vars: list of vars from self.seed to fit for
        - score_function_terms: terms in the transfer matrix that will be used
          for calculating sum(chi2).
        """
        ierr = ROOT.Long()
        self.minuit = ROOT.TMinuit(len(fit_vars))
        self.current_fit = fit_vars
        # tell minuit the parameters that we vary
        for i, var in enumerate(fit_vars):
            j_min, j_ref, j_max = self._get_min_max(var)
            self.minuit.mnparm(i, var, j_ref, 1., j_min, j_max, ierr)
            self.minuit.SetPrintLevel(2) # -1 to silence
            self.minuit.SetFCN(self._score_function)
        args = numpy.array([self.max_iterations, self.resolution],
                            dtype=numpy.float64)
        # the score_function terms are the terms in transfer_matrix used in 
        # score_function to calculate sum(chi^2)
        self.score_function_terms = score_function_terms
        # run the minimiser
        self.minuit.mnexcm("SIMPLEX", args, 2, ierr)
        # update the best guess from minuit and errors and print them        
        self._update_from_minuit()
        print "Fits and estimated errors"
        for key in sorted(self.estimated_error.keys()):
            print key+":", round(self.best_guess[key], 3), \
                  "+/-", round(self.estimated_error[key], 3)

    def _get_configuration(self):
        """
        Generate a MAUS configuration for tracking
        """
        config_str = Configuration.Configuration().\
                                           getConfigJSON(command_line_args=True)
        config_json = json.loads(config_str)
        config_json['maximum_number_of_steps'] = 2000

        return json.dumps(config_json)

    def _get_min_max(self, var):
        """
        Get the minimum and maximum for a minuit parameter definition
        """
        d_ref = self.seed[var]
        d_min = self.seed[var] - self.max_delta[var]
        d_max = self.seed[var] + self.max_delta[var]
        return d_min, d_ref, d_max

    def _update_from_minuit(self):
        """
        update the best guess and estimated error based on minuit parameters
        """
        for i, var in enumerate(self.current_fit):
            current, error = ROOT.Double(), ROOT.Double()
            self.minuit.GetParameter(i, current, error)
            self.best_guess[var] = current
            self.estimated_error[var] = error

    @staticmethod
    def run_once(print_field = False):
        self = TMFitter.fitter
        # update the best_guess
        self._update_from_minuit()
        self.calculated_tm_list = []
        for key in self.current_fit:
            print key+":", self.best_guess[key],
        print

        for tm_list in self.measured_tm_list:
            self.calculated_tm_list.append([])
            for tm in tm_list:
                measured_tm = tm['transfer_matrix']
                measured_errors = tm['error_matrix']
                p_tot = tm['p']
                # calculate the transfer matrix based on best_guess parameters
                if self.config.fit_test_mode:
                    calculated_tm = measured_tm
                    for i in range(4):
                        for j in range(5):
                            calculated_tm['transfer_matrix'][i][j] *= 1.1
                else:
                    calculated_tm = self.get_tm(self.best_guess, p_tot)
                # calculate the sum(chi^2)
                score = 0.
                chi2 = [[0. for i in range(5)] for j in range(4)]
                for row, column in self.score_function_terms:
                    if measured_errors == None:
                        chi2[row][column] = -1.
                    else:
                        chi2[row][column] = (calculated_tm[row][column] - measured_tm[row][column])**2/measured_errors[row][column]**2
                        score += chi2[row][column]
                self.calculated_tm_list[-1].append({
                      'transfer_matrix':calculated_tm,
                      'normalised_residual':chi2,
                      'p':p_tot
                })
                # print some diagnostics
                print "Pz:", self.tracking.pz
                print "Calculated Matrix:".ljust(54), "  * ", "Measured Matrix:".ljust(54), "  * Normalised residuals Per Cell:"
                for row, cells in enumerate(chi2):
                    for column, element in enumerate(cells):
                        print str(round(calculated_tm[row][column], 5)).rjust(10),
                    print "  * ",
                    for column, element in enumerate(cells):
                        print str(round(measured_tm[row][column], 5)).rjust(10),
                    print "  * ",
                    for column, element in enumerate(cells):
                        print str(round(element, 5)).rjust(10),
                    print
        return score


    @staticmethod
    def _score_function(n_pars, pars, score, buff, err):
        """
        Calculate the transfer matrix and compare with the measured tranfer 
        matrix; use to calculate a chi2
        """
        score[0] = TMFitter.fitter.run_once()
        self._iteration += 1
        print "Iteration", self._iteration, "Score", score[0], "="*100
        print

