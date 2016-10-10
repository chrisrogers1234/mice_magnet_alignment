import json
import os
import shutil

import scripts.fit_transfer_matrix.tm_fitter as tm_fitter
import scripts.tm_plotter as tm_plotter

import scripts.config_reco_straight
import scripts.config_reco_fc_pos
import scripts.config_reco_fc_neg
import scripts.config_reco_fc_test

def do_plot(config):
    try:
        shutil.rmtree(config.tm_plot_dir)
    except OSError:
        pass
    os.makedirs(config.tm_plot_dir)

    plotter = tm_plotter.TMPlotter(config)

def do_fit(config):
    data_list = []
    for analysis_index, analysis in enumerate(config.analyses):
        tm_file = config.data_dir+"/tm_"+str(analysis_index)+".json"
        fin = open(tm_file)
        data_list.append(json.loads(fin.read()))

    fitter = tm_fitter.TMFitter(config)
    fitter.set_measured_data(data_list)
    if config.will_do_not_fit:
        fitter.test_run()
        fitter.save_fit("fitted")
    if config.will_do_fit:
        fitter.fit()
        fitter.save_fit("fitted")

def print_cuts(config):
    data_list = []
    print "Cut".ljust(30),
    for analysis_index, analysis in enumerate(config.analyses):
        print analysis["name"].ljust(21),
        tm_file = config.data_dir+"/tm_"+str(analysis_index)+".json"
        fin = open(tm_file)
        data_list.append(json.loads(fin.read()))
    print
    print ",".ljust(30),

    for analysis_index, analysis in enumerate(config.analyses):
        print "All data".ljust(10), "After quality cuts".ljust(10),
    print

    for key in sorted(data_list[0]["cuts_active"]):
        if data_list[0]["cuts_active"][key] == False:
            continue
        print str(key).ljust(30),
        for data in data_list:
            print str(data["accepted"][key]).ljust(10),
            print str(data["accepted_data"][key]).ljust(10),
        print

  
def main():
    config = scripts.config_reco_fc_neg.Config()
    print_cuts(config)
    #do_fit(config)
    #do_plot(config)

if __name__ == "__main__":
    main()
    print "Finished - press <CR> to end"
    raw_input()

