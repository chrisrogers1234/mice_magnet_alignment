import numpy
import ROOT
import xboa.common
import maus_cpp.polynomial_map

import scripts.lattice

class TMCalculator():
    def __init__(self, config, analysis_index, p_low, p_high):
        self.config = config
        self.plot_dir = config.analyses[analysis_index]["plot_dir"]
        self.pid = config.analyses[analysis_index]["pid"]
        self.p_low = p_low
        self.p_high = p_high
        self.will_cut = MomPredicate(p_low, p_high).test
        self.events = []
        self.fit = []
        self.tm = None
        self.tm_error = None

    def append_data(self, data):
        self.events += data
        print self.will_cut
        n_events = len([event for event in self.events if not self.will_cut(event)])
        print "TMCalculator has", len(self.events), \
              "events with p range", (self.p_low, self.p_high), \
              "of which", n_events, "pass cut"

    def calculate_lattice_tm(self, z_out, momentum):
        z_in = self.config.fit_z_us
        #beta = self.config.fit_z_tof12/tof12/xboa.common.constants['c_light']
        #gamma = 1/(1-beta**2)**0.5
        #mass = xboa.common.pdg_pid_to_mass[abs(self.pid)]
        #momentum = beta*gamma*mass
        currents = {
                "FocusCoil_US":self.config.fc_current_seed,
                "FocusCoil_DS":self.config.fc_current_seed,
        }
        lattice = scripts.lattice.Lattice(momentum, z_in, currents, self.config.geometry)
        lattice.calculate_transfer_matrices()
        tm = lattice.get_tm(z_out)
        return tm

    def calculate_fit_tm(self, recon_error_matrix):
        data_in = self.choose_data("upstream")
        data_out = self.choose_data("downstream")
        for data in data_in[0:2]:
            print "    data in", data
        for data in data_out[0:2]:
            print "    data out", data
        print "In:", len(data_in), "Out:", len(data_out), "events in tm calculation"
        if len(data_in) == 0 or len(data_in) != len(data_out):
            print "Bad data - not evaluating TM"
            return None, None
        self.tm = maus_cpp.polynomial_map.PolynomialMap.least_squares_fit(data_in, data_out, 1, recon_error_matrix)
        self.tm_error = self.calculate_tm_error(recon_error_matrix, data_in, data_out)
        return self.tm, self.tm_error

    def calculate_tm_error(self, recon_error_matrix, data_in, data_out):
        data_length = len(data_in)/self.config.n_subsamples

        print "Calculating errors in p bin:", self.p_low, self.p_high, \
              "all data length:", len(data_in), "subsample length:", data_length
        if data_length < self.config.min_subsample_size:
            print "Insufficient data in bin"
            return None
        tm_error = [[[] for col in range(5)] for row in range(4)]
        for i in range(self.config.n_subsamples):
            a_data_in = data_in[data_length*i:data_length*(i+1)]
            a_data_out = data_out[data_length*i:data_length*(i+1)]
            tm = maus_cpp.polynomial_map.PolynomialMap.least_squares_fit(a_data_in, a_data_out, 1, recon_error_matrix)
            tm = tm.get_coefficients_as_matrix()
            for row in range(len(tm)):
                for col in range(len(tm[row])):
                    tm_error[row][col].append(tm[row][col])
        err_coeff = self.config.n_subsamples**0.5
        for row in range(4):
            for col in range(5):
                tm_error[row][col] = numpy.std(tm_error[row][col])/err_coeff
        return tm_error

    def print_tm(self):
        print "TM:"
        for row in self.tm.get_coefficients_as_matrix():
            for cell in row:
                print str(round(cell, 4)).rjust(10),
            print
        if self.tm_error == None:
            print "No errors calculated"
            return
        print "Errors:"
        for row in self.tm_error:
            for cell in row:
                print str(round(cell, 4)).rjust(10),
            print
        print

    def get_res_fit(self, axis):
        data_res = self.choose_data("residual")
        res = [res[axis] for res in data_res]
        fit = ROOT.TF1("fit", "gaus")
        hist = xboa.common.make_root_histogram("residual fit", res, "Residual - "+str(axis), 1000, xmin = -500., xmax = 500.)
        hist.Draw()
        hist.Fit(fit, "Q")
        self.fit.append(fit)
        return fit.GetParameter(2)

    def residual_amplitude_cut(self, amplitude_cut):
        data_res_list = self.choose_data("residual")
        cov = numpy.cov(data_res_list, rowvar=0)
        cov_inv = numpy.linalg.inv(cov)
        for event in self.events:
            if self.will_cut(event):
                continue
            if event["data_res"] == None:
                event["residual_amplitude"] = None
                event["will_cut"]["residual_cut"] = True
                event["any_cut"] = True
                continue
            u = numpy.array(event["data_res"])
            u_t = numpy.transpose(u)
            amp = numpy.dot(numpy.dot(u_t, cov_inv), u)
            event["residual_amplitude"] = amp
            if amp > amplitude_cut:
                event["will_cut"]["residual_cut"] = True
                if self.config.cuts_active["residual_cut"]:
                    event["any_cut"] = True
            else:
                event["will_cut"]["residual_cut"] = False

    def lattice_scraping_cut(self, cut, p_tot):
        """
        cut on ((x_proj + stay_clear)**2 + (y_proj + stay_clear)**2)**0.5 < max_radius
        """
        max_radius = cut["max_radius"]
        stay_clear = cut["stay_clear"]
        z_pos = cut["z_pos"]
        cut_name = cut["cut_name"]
        tm = self.calculate_lattice_tm(z_pos, p_tot)
        for event in self.events:
                if self.will_cut(event):
                    continue
                proj = tm.evaluate(event["data_in"])
                radius = ((abs(proj[0])+stay_clear)**2 + (abs(proj[2])+stay_clear)**2)**0.5
                if radius > max_radius:
                    event["will_cut"][cut_name] = True
                    if self.config.cuts_active[cut_name]:
                        event["any_cut"] = True
                else:
                    event["will_cut"][cut_name] = False
                #print proj[0], proj[2], event["data_out"][0], event["data_out"][2], self.events

    def scraping_cut(self, n_sigma, max_radius, amp_cut, n_steps, recon_error):
        sigma_steps = [n_sigma*1.*i/n_steps for i in range(1, n_steps+1)]
        print "Scraping cut sigma steps", sigma_steps
        canvas = ROOT.TCanvas("scraping_cut_canvas", "scraping_cut_canvas")
        canvas.Draw()
        canvas.Divide(2, 1)
        for this_n_sigma in sigma_steps:
            # big memory usage start
            tm, tm_error = self.calculate_fit_tm(recon_error)
            # big memory usage end
            for event in self.events:
                if self.will_cut(event):
                    continue
                try:
                    data_proj = tm.evaluate(event["data_in"])
                    data_res = [data_proj[j] - event["data_out"][j] for j in range(4)]
                    event["data_proj"] = data_proj
                    event["data_res"] = data_res
                except TypeError:
                    event["data_res"] = None
                    event["data_proj"] = None
            canvas.cd(1)
            sigma_x = self.get_res_fit(0)
            canvas.cd(2)
            sigma_y = self.get_res_fit(2)
            canvas.Update()
            print "Found sigma x", sigma_x, "sigma y", sigma_y, "with", this_n_sigma
            for event in self.events:
                if self.will_cut(event):
                    continue
                proj = event["data_proj"]
                radius = ((sigma_x*this_n_sigma+abs(proj[0]))**2 + (sigma_y*this_n_sigma+abs(proj[2]))**2)**0.5
                if radius > max_radius:
                    event["will_cut"]["scraping_cut"] = True
                    if self.config.cuts_active["scraping_cut"]:
                        event["any_cut"] = True
                else:
                    event["will_cut"]["scraping_cut"] = False
                #print proj[0], proj[2], event["data_out"][0], event["data_out"][2], self.events
            residual_amplitudes = self.residual_amplitude_cut(amp_cut)
        for format in ["png", "root", "pdf"]:
            p_str = str(self.p_low)+"_"+str(self.p_high)
            canvas.Print(self.plot_dir+"scraping-cut_"+p_str+"."+format)

    def choose_data(self, choice):
        key_map = {
            "upstream":"data_in",
            "downstream":"data_out",
            "projected":"data_proj",
            "residual":"data_res",
        }
        key = key_map[choice]
        data = [event[key] for event in self.events if not self.will_cut(event)]
        return data

    def json_repr(self):
        pass

class TOF12Predicate(object):
    def __init__(self, tof12_min, tof12_max):
        self.tof12_min = tof12_min
        self.tof12_max = tof12_max

    def test(self, event):
        tof12 = event["tof12"]
        pred = event["any_cut"] or tof12 < self.tof12_min or tof12 > self.tof12_max
        return pred
        

class MomPredicate(object):
    def __init__(self, p_min, p_max):
        self.p_min = p_min
        self.p_max = p_max

    def test(self, event):
        ptot = event["p_tot"]
        pred = event["any_cut"] or ptot < self.p_min or ptot > self.p_max
        return pred
      

