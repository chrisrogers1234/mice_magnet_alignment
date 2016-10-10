import sys
import os
import site
import json
import shutil

import numpy
import ROOT
ROOT.gROOT.SetBatch(True)

import xboa.common
import Configuration
import maus_cpp.globals
import maus_cpp.field
import maus_cpp.polynomial_map
import libxml2

from scripts.tm_calculator import TMCalculator
from scripts.data_plotter import DataPlotter
from scripts.tm_calculator import TOF12Predicate
import scripts.data_loader
import scripts.config_reco_fc_pos
import scripts.config_reco_fc_neg
import scripts.config_reco_fc_test
import scripts.config_reco_straight
import scripts.config_reco_8155
import scripts.config_reco_614

def maus_globals(config):
    str_conf = Configuration.Configuration().\
                                      getConfigJSON(command_line_args=False)
    json_conf = json.loads(str_conf)
    json_conf["simulation_geometry_filename"] = config.geometry
    json_conf["verbose_level"] = 1
    maus_conf = json.dumps(json_conf)
    maus_cpp.globals.birth(maus_conf)
    print maus_cpp.field.str(True)

def do_analysis(config, analysis_index):
    config_anal = config.analyses[analysis_index]
    all_event_count = 0

    for max_spill in [config.number_of_spills]:
        if max_spill == 0:
            continue
        data_loader = scripts.data_loader.DataLoader(config, analysis_index)
        data_loader.load_data(0, max_spill)
    sigma_x = 0.7
    sigma_xp = 5e-3
    print "Applying recon error with est x error", sigma_x, "est x' error", sigma_xp
    recon_error = [[0 for i in range(5)] for j in range(5)]
    recon_error[1][1] = sigma_x**2
    recon_error[2][2] = sigma_xp**2
    recon_error[3][3] = sigma_x**2
    recon_error[4][4] = sigma_xp**2

    tm_p_list = []
    print "Clearing old data"
    try:
        shutil.rmtree(config_anal["plot_dir"])
    except OSError:
        pass
    os.makedirs(config_anal["plot_dir"])
    print "Using p bins", config_anal["p_bins"]
    for p_low, p_high in config_anal["p_bins"]:
        continue
        calculator = TMCalculator(config, analysis_index, p_low, p_high)
        calculator.append_data(data_loader.events)
        tm, tm_error = calculator.calculate_fit_tm(recon_error)
        calculator.scraping_cut(config_anal["sigma_cut"], 150., config_anal["chi2_cut"], 5, recon_error) #n_sigma, max_radius, amp_cut, n_steps, use_recon
        for cut in config.lattice_scraping_cuts:
            p_tot = (p_high + p_low)/2.
            calculator.lattice_scraping_cut(cut, p_tot)
        tm, tm_error = calculator.calculate_fit_tm(recon_error)
        if tm == None or tm_error == None:
            continue
        calculator.print_tm()
        tm_p_list.append({"transfer_matrix":tm.get_coefficients_as_matrix(),
                          "error_matrix":tm_error,
                          "p_high":p_high,
                          "p_low":p_low,
                          "p":(p_high+p_low)/2.,
                          })
    xboa.common.clear_root()
    plot_predicate = TOF12Predicate(config_anal["tof12_cut_low"], config_anal["tof12_cut_high"])
    plotter = DataPlotter(config, analysis_index, data_loader.events, plot_predicate.test)
    plotter.print_cuts_summary()
    cov_us, cov_ds = plotter.print_covariance_matrix()
    sys.stdout.flush()
    plotter.plot_tof12()
    plotter.plot_var_2d(0, "upstream", 1, "upstream")
    plotter.plot_var_2d(2, "upstream", 3, "upstream")
    plotter.plot_var_2d(0, "downstream", 1, "downstream")
    plotter.plot_var_2d(2, "downstream", 3, "downstream")
    plotter.plot_var_2d(0, "downstream", 0, "residual")
    plotter.plot_var_2d(1, "downstream", 1, "residual")
    plotter.plot_var_2d(2, "downstream", 2, "residual")
    plotter.plot_var_2d(3, "downstream", 3, "residual")
    plotter.plot_var_1d(0, "residual")
    plotter.plot_var_1d(1, "residual")
    plotter.plot_var_1d(2, "residual")
    plotter.plot_var_1d(3, "residual")
    plotter.plot_var_1d(0, "upstream")
    plotter.plot_var_1d(1, "upstream")
    plotter.plot_var_1d(2, "upstream")
    plotter.plot_var_1d(3, "upstream")
    plotter.plot_var_1d(0, "downstream")
    plotter.plot_var_1d(1, "downstream")
    plotter.plot_var_1d(2, "downstream")
    plotter.plot_var_1d(3, "downstream")
    plotter.plot_r_ds()
    #plotter.plot_chi2()
    #plotter.plot_path_length()
    #plotter.plot_miss_distance()
    plotter.plot_p_tot()
    plotter.plot_p_tot_res()
    plotter.plot_p_tot_vs_tof()
    plotter.plot_p_tot_us_vs_ds()
    #for col in range(5):
    #    for row in range(4):
    #        plotter.plot_tm_element(row, col)
    accepted, accepted_data = plotter.get_cuts_summary()
    json_doc = {
        "covariance_matrix_us":[cov_us.tolist()], # should be done in p bins
        "covariance_matrix_ds":[cov_ds.tolist()], # should be done in p bins
        "tm_entries":tm_p_list,
        "accepted":accepted, # accepted by all cuts
        "accepted_data":accepted_data, # accepted after data cut
        "cuts_active":config.cuts_active,
    }
    fname = config.data_dir+"/tm_"+str(analysis_index)+".json"
    fout = open(fname, "w")
    print >> fout, json.dumps(json_doc, indent = 2)

def main():
    config = scripts.config_reco_614.Config()
    maus_globals(config)

    for i, anal in enumerate(config.analyses):
        do_analysis(config, i)

if __name__ == "__main__":
    main()
    print "Done - press <CR> to finish"
    #raw_input()

