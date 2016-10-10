import json

import ROOT

import xboa.common

class TMPlotter(object):
    def __init__(self, config):
        self.config = config
        self.data = []
        if config.plot_incl_fit:
            file_list = [config.data_dir+"/tm_fitted_"+str(analysis_index)+".json" \
                         for analysis_index, analysis in enumerate(config.analyses)]
        else:
            file_list = [config.data_dir+"/tm_"+str(analysis_index)+".json" \
                         for analysis_index, analysis in enumerate(config.analyses)]
        self.append_data(file_list)
        for row in range(4):
            for col in range(5):
                canvas = self.plot_tm_element(row, col)

    def append_data(self, file_list):
        for tm_file in file_list:
            print "Opening", tm_file
            fin = open(tm_file)
            self.data.append(json.loads(fin.read()))

    def get_one_data(self, data_key, a_data, row, column, tm_list, err_list, p_list):
        if data_key not in a_data:
            print data_key, "not in data with keys", a_data.keys()
            return False
        one_data = a_data[data_key]
        tm_list.append([tm["transfer_matrix"][row][column] for tm in one_data])
        err_list_this = []
        for tm in one_data:
            if "error_matrix" not in tm or tm["error_matrix"] == None:
                err_list_this.append(None)
            else:
                err_list_this.append(tm["error_matrix"][row][column])
        err_list.append(err_list_this)
        p_list.append([tm["p"] for tm in one_data])
        print data_key, tm_list[-1], err_list[-1]
        return True

    def get_data(self, row, column):
        name_list = []
        tm_list = []
        err_list = []
        p_list = []
        for key, name in [("tm_entries", "Measured M - "),
                          ("fitted_tms", "Fitted M - "),
                          ("not_fitted_tms", "Not Fitted M - ")]:
            for i, a_data in enumerate(self.data):
                okay = self.get_one_data(key, a_data, row, column, tm_list, err_list, p_list)
                if okay:
                    name_list.append(name+self.config.analyses[i]["name"])
        return name_list, p_list, tm_list, err_list

    def plot_tm_element(self, row, column):
        name_lookup = ["x", "x'", "y", "y'"]
        name_list, p_list, tm_list, err_list = self.get_data(row, column)
        name = "M_"+str(row)+str(column)
        canvas = xboa.common.make_root_canvas(name)
        canvas.SetWindowSize(1200, 800)
        canvas.SetCanvasSize(900, 700)
        canvas.SetBorderSize(1)
        canvas.SetBorderMode(1)
        canvas.Draw()
        canvas.Divide(2, 1)
        canvas.cd(1)
        canvas.GetPad(1).SetCanvasSize(1700, 700)
        canvas.GetPad(1).SetBorderSize(1)
        canvas.GetPad(1).SetBorderMode(1)
        #canvas.GetPad(2).SetCanvasSize(100, 700)
        canvas.GetPad(2).SetBorderSize(1)
        axis = "M_{"+str(row)+str(column)+"} aka M("
        axis += name_lookup[row]+"_{o}"
        if column > 0:
            axis += name_lookup[column-1]+"_{i}"
        axis += ")"
        if row % 2 == 0 and column % 2 == 0:
            axis += " [mm]"
        elif row % 2 == 1 and column % 2 == 1:
            axis += " [mm^{-1}]"
        all_tms, all_p = [], []
        for tm in tm_list:
            all_tms += tm
        for p_tot in p_list:
            all_p += p_tot
        hist, graph = xboa.common.make_root_graph(name, all_p, "p [MeV/c]", all_tms, axis)
        hist.GetYaxis().SetTitleOffset(1.4)
        hist.Draw()
        legend_list = []
        for i in range(len(tm_list)):
            anal_index = i % len(self.config.analyses)
            fit_index = i / len(self.config.analyses)
            a_p_list = p_list[i]
            a_tm_list = tm_list[i]
            a_err_list = err_list[i]
            this_name = name_list[i]
            p_delta = (a_p_list[1]-a_p_list[0])/2.
            n_points = len(a_p_list)
            graph = ROOT.TGraphErrors()
            graph.SetName(this_name)
            print this_name
            for i in range(n_points):
                if a_err_list[i] == None:
                    if self.config.plot_when_no_errors:
                        a_err_list[i] = 0.
                    else:
                        continue
                graph.SetPoint(i, a_p_list[i], a_tm_list[i])
                graph.SetPointError(i, p_delta, a_err_list[i])
                #print i, a_tof_list[i], a_tm_list[i], tof_delta, a_err_list[i]
            graph.SetMarkerStyle(23+fit_index)
            graph.SetMarkerColor(self.config.analyses[anal_index]["color"])
            graph.SetFillColor(0)
            #graph.SetLineColor(0)
            graph.Draw("SAMEP")
            self.root_objects.append(graph)
            legend_list.append(graph)
            canvas.Update()
        leg = xboa.common.make_root_legend(canvas, legend_list)
        canvas.cd(2)
        leg.Draw()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.config.tm_plot_dir+name+"."+format)
        return canvas

    root_objects = []

