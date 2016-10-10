import xboa.common
import json
import copy
import numpy

class DataPlotter(object):
    def __init__(self, config, analysis_index, events, will_cut):
        self.events = events
        self.will_cut = will_cut
        self.config = config
        self.config_analysis = config.analyses[analysis_index]
        self.plot_dir = config.analyses[analysis_index]["plot_dir"]

    def choose_data(self, choice, cuts):
        key_map = {
            "upstream":"data_in",
            "upstream normalised":"data_in",
            "downstream":"data_out",
            "downstream normalised":"data_out",
            "projected":"data_proj",
            "residual":"data_res",
        }
        key = key_map[choice]
        normalised = "normalised" in choice

        if cuts:
            predicate = lambda event: not self.will_cut(event) and event[key] != None and event['p_tot'] != None
        else:
            predicate = lambda event: event[key] != None and event['p_tot'] != None
        if normalised:
            normalise = lambda event: self.normalise(event, key)
        else:
            normalise = lambda event: event[key]
        data = [normalise(event) for event in self.events if predicate(event)]
        return data

    def normalise(self, event, key):
        a_dat = copy.deepcopy(event[key])
        if "data_out" in key:
            ptot = event['p_tot_ds']
        else:
            ptot = event['p_tot_us']
        # ptot^2 = px^2 + py^2 + pz^2 => pz = ptot/(x'^2+y'^2+1)**0.5
        pz = ptot/(a_dat[1]**2+a_dat[3]**2+1)**0.5
        a_dat[1] *= pz
        a_dat[3] *= pz
        return a_dat

    def print_optics(self, covariance_matrix):
        emit_perp = numpy.linalg.det(covariance_matrix)**0.25/xboa.common.pdg_pid_to_mass[13]
        print "emit 4D ", emit_perp

    def print_covariance_matrix(self):
        data_us = self.choose_data("upstream normalised", cuts = True)
        covariance_matrix_us = numpy.cov(data_us, rowvar = False)
        print "covariances upstream:"
        print covariance_matrix_us
        self.print_optics(covariance_matrix_us)
        print
        data_ds = self.choose_data("downstream normalised", cuts = True)
        covariance_matrix_ds = numpy.cov(data_ds, rowvar = False)
        print "covariances downstream:"
        print covariance_matrix_ds
        self.print_optics(covariance_matrix_ds)
        print
        return covariance_matrix_us, covariance_matrix_ds

    def get_cuts_summary(self):
        # there are a set of cuts that are to do with basic "data quality" e.g.
        # exactly one hit in TOF, etc and there are a set of cuts that are to do
        # with "bias correction"/etc.
        # I report the number of events that pass each data quality cut;
        # and the number of events that pass each "bias correction" cut AND the data quality cuts
        accepted_dict = {"any_cut":0, "data_cut":None, "all_events":0}
        accepted_data_dict = {"any_cut":None, "data_cut":0, "all_events":0}
        for key in self.events[0]["will_cut"]:
            accepted_dict[key] = 0
            accepted_data_dict[key] = 0
        for event in self.events:
            accepted_dict["all_events"]  += 1 # total number of events considered
            accepted_data_dict["all_events"]  += 1 # total number of events considered
            will_cut = event["will_cut"]
            for key in will_cut:
                if not will_cut[key]:
                    accepted_dict[key] += 1
                    if event["data_cut"]:
                        continue
                    accepted_data_dict[key] += 1
            if not event["any_cut"]:
                accepted_dict["any_cut"]  += 1 # total number accepted by all cuts
            if not event["data_cut"]:
                accepted_data_dict["data_cut"]  += 1 # total number accepted by data cuts
        return accepted_dict, accepted_data_dict

    def print_cuts_summary(self):
        cut_dict, cut_data_dict = self.get_cuts_summary()
        for key in sorted(cut_dict.keys()):
            if key in self.config.cuts_active:
                is_active = self.config.cuts_active[key]
            else:
                is_active = None
            print "   ", key.ljust(25), str(is_active).ljust(5), str(cut_dict[key]).ljust(8), str(cut_data_dict[key]).ljust(8)


    def plot_tof12(self):
        tof12_no_cut = [event["tof12"] for event in self.events if not self.will_cut(event)]
        tof12_all = [event["tof12"] for event in self.events]
        tof12_all = [tof12 for tof12 in tof12_all if tof12 != None]
        tof12_canvas = xboa.common.make_root_canvas("tof12_canvas")
        tof12_canvas.SetLogy()
        xmin = min(25., self.config_analysis["tof12_cut_low"])
        xmax = max(45., self.config_analysis["tof12_cut_high"])
        hist = xboa.common.make_root_histogram("tof12", tof12_all, "tof2 - tof1 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("tof12 in cut", tof12_no_cut, "tof2 - tof1 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetLineColor(4)
        hist.Draw("SAME")
        tof12_canvas.Update()
        for format in ["png", "root", "pdf"]:
            tof12_canvas.Print(self.plot_dir+"tof12."+format)

    def plot_chi2(self):
        residual_amplitudes_no_cut = [event["residual_amplitude"] for event in self.events if not self.will_cut(event)]
        residual_amplitudes_all = [event["residual_amplitude"] for event in self.events]
        residual_amplitudes_all = [res for res in residual_amplitudes_all if res != None]
        for event in self.events:
            if event["residual_amplitude"] == None and not self.will_cut(event):
                print json.dumps(event, indent=2)

        res_canvas = xboa.common.make_root_canvas("chi2_cut_canvas")
        res_canvas.SetLogy()
        hist = xboa.common.make_root_histogram("chi2", residual_amplitudes_all, "#chi^{2}", 100, xmin = 0., xmax = 50.)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("Included in analysis", residual_amplitudes_no_cut, "#chi^{2}", 100, xmin = 0., xmax = 50.)
        hist.SetLineColor(4)
        hist.Draw("SAME")
        res_canvas.Update()
        for format in ["png", "root", "pdf"]:
            res_canvas.Print(self.plot_dir+"chi2."+format)

    def plot_path_length(self):
        canvas = xboa.common.make_root_canvas("path length factor")
        z_0 = self.config.z_tof2 - self.config.z_tof1
        path_lengths_no_cut = [event["path_length"]/z_0 for event in self.events if not self.will_cut(event)]
        path_lengths_all = [event["path_length"] for event in self.events]
        path_lengths_all = [length/z_0 for length in path_lengths_all if length != None]

        hist = xboa.common.make_root_histogram("All events", path_lengths_all, "(est. path length)/(tof12 separation)", 100, xmin=1.0, xmax=1.01)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("Events in cut", path_lengths_no_cut, "(est. path length)/(tof12 separation)", 100, xmin=1.0, xmax=1.01)
        hist.SetLineColor(4)
        hist.Draw("SAME")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"path_length."+format)

    def plot_miss_distance(self):
        canvas = xboa.common.make_root_canvas("miss distance")
        miss_no_cut = [event["miss_distance"] for event in self.events if not self.will_cut(event)]
        miss_all = [event["miss_distance"] for event in self.events]
        miss_all = [length for length in miss_all if length != None]

        hist = xboa.common.make_root_histogram("All events", miss_all, "Miss distance [mm]", 100)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("Events in cut", miss_no_cut, "Miss distance [mm]", 100)
        hist.SetLineColor(4)
        hist.Draw("SAME")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"miss_distance."+format)

    def plot_p_tot(self):
        canvas = xboa.common.make_root_canvas("p_tot")
        p_tot_no_cut = [event["p_tot"] for event in self.events if not self.will_cut(event)]
        p_tot_all = [event["p_tot"] for event in self.events]
        p_tot_all = [p_tot for p_tot in p_tot_all if p_tot != None]

        hist = xboa.common.make_root_histogram("All events", p_tot_all, "p_{tot} [MeV/c]", 100, xmin=0., xmax=400.)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("Events in cut", p_tot_no_cut, "p_{tot} [MeV/c]", 100, xmin=0., xmax=400.)
        hist.SetLineColor(4)
        hist.Draw("SAME")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"p_tot."+format)

    def plot_p_tot_vs_tof(self):
        canvas = xboa.common.make_root_canvas("p_tot vs tof no cuts")
        tof_no_cut = [event["tof12"] for event in self.events if event["tof12"] != None and event["p_tot"] != None ]
        p_tot_no_cut = [event["p_tot"] for event in self.events if event["tof12"] != None and event["p_tot"] != None ]

        hist = xboa.common.make_root_histogram("Events in cut", tof_no_cut, "tof2 - tof1 [ns]", 100, p_tot_no_cut, "p_{tot} [MeV/c]", 100, ymin=0., ymax=300., xmin=25., xmax=45.)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"p_tot_vs_tof12_all."+format)

        canvas = xboa.common.make_root_canvas("p_tot vs tof cuts")
        tof_no_cut = [event["tof12"] for event in self.events if not self.will_cut(event)]
        p_tot_no_cut = [event["p_tot"] for event in self.events if not self.will_cut(event)]

        hist = xboa.common.make_root_histogram("Events in cut", tof_no_cut, "tof2 - tof1 [ns]", 100, p_tot_no_cut, "p_{tot} [MeV/c]", 100)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"p_tot_vs_tof12_cuts."+format)


    def plot_p_tot_res(self):
        canvas = xboa.common.make_root_canvas("p_tot us vs ds")
        cut_pred = lambda event: not event["p_tot_us"] == None and not event["p_tot_ds"] == None
        p_tot_us_no_cut = [event["p_tot_us"] for event in self.events if cut_pred(event)]
        p_tot_ds_no_cut = [event["p_tot_ds"] for event in self.events if cut_pred(event)]
        p_tot_res_no_cut = [p_tot - p_tot_ds_no_cut[i] for i, p_tot in enumerate(p_tot_us_no_cut)]

        p_tot_us_cut = [event["p_tot_us"] for event in self.events if not self.will_cut(event)]
        p_tot_ds_cut = [event["p_tot_ds"] for event in self.events if not self.will_cut(event)]
        p_tot_res_cut = [p_tot - p_tot_ds_cut[i] for i, p_tot in enumerate(p_tot_us_cut)]
        mean = numpy.mean(p_tot_res_cut)
        std = numpy.std(p_tot_res_cut)

        hist = xboa.common.make_root_histogram("All events", p_tot_res_no_cut, "p_{tot} [MeV/c]", 100, xmin=-50, xmax=50)
        hist.SetTitle(self.config_analysis['name'])
        hist.SetTitle("mean: "+str(mean)+" std: "+str(std))
        hist.Draw()
        hist = xboa.common.make_root_histogram("Events in cut", p_tot_res_cut, "p_{tot} [MeV/c]", 100, xmin=-50, xmax=50)
        hist.SetTitle("mean: "+str(mean)+" std: "+str(std))
        hist.SetLineColor(4)
        hist.Draw("SAME")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"p_tot_res."+format)

    def plot_p_tot_us_vs_ds(self):
        canvas = xboa.common.make_root_canvas("p_tot us vs ds")
        p_tot_us = [event["p_tot_us"] for event in self.events if not self.will_cut(event)]
        p_tot_ds = [event["p_tot_ds"] for event in self.events if not self.will_cut(event)]
        p_tot_res = [p_tot - p_tot_ds[i] for i, p_tot in enumerate(p_tot_us)]

        hist = xboa.common.make_root_histogram("Events in cut", p_tot_us, "p_{tot} upstream [MeV/c]", 100,
                                                                p_tot_ds, "p_{tot} downstream [MeV/c]", 100)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"p_tot_us_vs_ds."+format)

        hist = xboa.common.make_root_histogram("Events in cut", p_tot_us, "p_{tot} upstream [MeV/c]", 100,
                                                                p_tot_res, "p_{tot} resolution [MeV/c]", 100)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"p_tot_us_vs_res."+format)

    def plot_r_ds(self):
        canvas = xboa.common.make_root_canvas("radius downstream")
        nbins = 101
        data_out = self.choose_data("downstream", False)
        r_ds = [(u[0]**2+u[2]**2)**0.5 for u in data_out]
        hist = xboa.common.make_root_histogram("All events", r_ds, "r [mm]", 100, xmin = 0., xmax = 250.)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        data_out = self.choose_data("downstream", True)
        r_ds = [(u[0]**2+u[2]**2)**0.5 for u in data_out]
        hist = xboa.common.make_root_histogram("All events", r_ds, "r [mm]", 100, xmin = 0., xmax = 250.)
        hist.SetLineColor(4)
        hist.Draw("SAME")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"measured_radius."+format)

    def index_to_label(self, index):
        return ["x [mm]", "p_{x} [MeV/c]", "y [mm]", "p_{y} [MeV/c]"][index]

    def index_to_name(self, index):
        return ["x", "px", "y", "py"][index]

    def plot_var_2d(self, index_1, us_ds_1, index_2, us_ds_2,): # MAKE THIS A CONT PLOT AND SUPERIMPOSE EACH TOF BIN
        name =  us_ds_1+"_"+self.index_to_name(index_1)+"_"+ us_ds_2+"_"+self.index_to_name(index_2)
        canvas = xboa.common.make_root_canvas(name)
        data_1 = self.choose_data(us_ds_1, True)
        data_2 = self.choose_data(us_ds_2, True)
        var_1 = [u_in[index_1] for i, u_in in enumerate(data_1)]
        var_2 = [u_out[index_2] for i, u_out in enumerate(data_2)]
        lab_1 = us_ds_1+" "+self.index_to_label(index_1)
        lab_2 = us_ds_2+" "+self.index_to_label(index_2)
        print "Plot var", lab_1, "vs", lab_2, "for", len(var_1), "events"
        hist = xboa.common.make_root_histogram(name,
                                               var_1, lab_1, 100,
                                               var_2, lab_2, 100)
        hist.Draw("COL")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+name+"."+format)

    def plot_var_1d(self, index_1, us_ds_1):
        name =  us_ds_1+"_"+self.index_to_name(index_1)
        canvas = xboa.common.make_root_canvas(name)
        data_1 = self.choose_data(us_ds_1, False)
        var_1 = [u_in[index_1] for i, u_in in enumerate(data_1)]
        lab_1 = us_ds_1+" "+self.index_to_label(index_1)
        print "Plot var", lab_1, "for", len(var_1), "events"
        hist = xboa.common.make_root_histogram(name, var_1, lab_1, 100)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        data_1 = self.choose_data(us_ds_1, True)
        var_1 = [u_in[index_1] for i, u_in in enumerate(data_1)]
        hist = xboa.common.make_root_histogram(name, var_1, lab_1, 100)
        hist.SetLineColor(4)
        hist.Draw("SAME")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+name+"."+format)


