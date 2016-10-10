import sys
import glob
import itertools
import math

import ROOT
import libMausCpp

import pz_calculator

class DataLoader(object):
    def __init__(self, config, analysis_index):
        self.config = config
        self.config_anal = config.analyses[analysis_index]
        self.tof_cut_count = 0
        self.tof12_cut_count = 0
        self.scifi_cluster_cut_count = 0
        self.scifi_space_point_cut_count = 0
        self.scifi_track_cut_count = 0
        self.scifi_track_point_cut_count = 0
        self.detector_count = {}
        self.event_count = 0
        self.accepted_count = 0
        self.pz_calc = pz_calculator.PzCalculator(config)
        self.events = []
        self.maus_version = ""
        self.cuts = {}

    def load_data(self, min_spill, max_spill):
        if min_spill >= max_spill:
            return
        file_name_list = []
        for fname in self.config_anal["reco_files"]:
            file_name_list += glob.glob(fname)
        file_name_list = sorted(file_name_list)
        print "Found", file_name_list, "files"
        if len(file_name_list) == 0:
            raise IOError("Couldn't find files for reco_files "+str(self.config_anal["reco_files"]))
        for file_name in file_name_list:
            print "Loading ROOT file", file_name
            sys.stdout.flush()
            root_file = ROOT.TFile(file_name, "READ") # pylint: disable = E1101

            data = ROOT.MAUS.Data() # pylint: disable = E1101
            tree = root_file.Get("Spill")
            try:
                tree.SetBranchAddress("data", data)
            except AttributeError:
                continue

            for i in range(min_spill, min(tree.GetEntries(), max_spill)):
                tree.GetEntry(i)
                spill = data.GetSpill()
                if spill.GetDaqEventType() == "physics_event":
                    for ev_number, reco_event in enumerate(spill.GetReconEvents()):
                        self.event_count += 1
                        try:
                            event = self.load_reco_event(reco_event)
                        except ValueError:
                            #sys.excepthook(*sys.exc_info())
                            print "spill", spill.GetSpillNumber(), "particle_number", reco_event.GetPartEventNumber()
                        except ZeroDivisionError:
                            pass
                        if event == None: # missing TOF1 - not considered further
                            continue 
                        event["spill"] = spill.GetSpillNumber()
                        self.events.append(event)
                        event["event_number"] = ev_number
            print "...loaded", tree.GetEntries(), "spills and", self.event_count, "events"
            sys.stdout.flush()

    def load_tof_event(self, tof_event):
        space_points = tof_event.GetTOFEventSpacePoint()
        tof_sp = []
        #if len(space_points.GetTOF2SpacePointArray()) > 0:
        #    print len(space_points.GetTOF0SpacePointArray()), len(space_points.GetTOF1SpacePointArray()), len(space_points.GetTOF2SpacePointArray())
        for tof0_sp in space_points.GetTOF0SpacePointArray():
            cov = [[0. for i in range(6)] for j in range(6)]
            cov[0][0] = 0.05
            cov[1][1] = tof0_sp.GetGlobalPosXErr()
            cov[2][2] = tof0_sp.GetGlobalPosXErr()
            tof0 = {
                "x":tof0_sp.GetGlobalPosX(),
                "y":tof0_sp.GetGlobalPosY(),
                "z":tof0_sp.GetGlobalPosZ(),
                "t":tof0_sp.GetTime(),
                "covariance":cov,
                "detector":"tof0"
            }
            tof_sp.append(tof0)
        for tof1_sp in space_points.GetTOF1SpacePointArray():
            cov = [[0. for i in range(6)] for j in range(6)]
            cov[0][0] = 0.05
            cov[1][1] = tof1_sp.GetGlobalPosXErr()
            cov[2][2] = tof1_sp.GetGlobalPosXErr()
            tof1 = {
                "x":tof1_sp.GetGlobalPosX(),
                "y":tof1_sp.GetGlobalPosY(),
                "z":tof1_sp.GetGlobalPosZ(),
                "t":tof1_sp.GetTime(),
                "covariance":cov,
                "detector":"tof1"
            }
            tof_sp.append(tof1)
        for tof2_sp in space_points.GetTOF2SpacePointArray():
            cov = [[0. for i in range(6)] for j in range(6)]
            cov[0][0] = 0.05
            cov[1][1] = tof2_sp.GetGlobalPosXErr()
            cov[2][2] = tof2_sp.GetGlobalPosXErr()           
            tof2 = {
                "x":tof2_sp.GetGlobalPosX(),
                "y":tof2_sp.GetGlobalPosY(),
                "z":tof2_sp.GetGlobalPosZ(),
                "t":tof2_sp.GetTime(),
                "covariance":cov,
                "detector":"tof2"
            }
            tof_sp.append(tof2)
        detectors = [x["detector"] for x in tof_sp]
        tof12_cut = False
        if "tof2" in detectors and "tof1" in detectors:
            tof12 = tof_sp[detectors.index("tof2")]["t"] - tof_sp[detectors.index("tof1")]["t"]
            if tof12 > self.config_anal["tof12_cut_high"] or \
               tof12 < self.config_anal["tof12_cut_low"]:
                self.tof12_cut_count += 1
                tof12_cut = True
        else:
            tof12_cut = True
        return (tof_sp, {
            "tof_0_sp":space_points.GetTOF0SpacePointArray().size() != 1,
            "tof_1_sp":space_points.GetTOF1SpacePointArray().size() != 1,
            "tof_2_sp":space_points.GetTOF2SpacePointArray().size() != 1,
            "tof12":tof12_cut,
          })

    def load_space_points(self, scifi_event):
        space_points_out = []
        for space_point in scifi_event.spacepoints():
            # if we require triplets and space point is not a triplet, skip it
            if self.config.will_require_triplets and \
                space_point.get_channels().GetEntries() != 3:
                continue
            position = space_point.get_global_position()
            if space_point.get_tracker() == 0:
                space_points_out.append({
                  "x":position.x(),
                  "y":position.y(),
                  "z":position.z(),
                  "t":None,
                  "detector":"tku"
                })
            else:
                space_points_out.append({
                  "x":position.x(),
                  "y":position.y(),
                  "z":position.z(),
                  "t":None,
                  "detector":"tkd"
                })
        space_points_out = sorted(space_points_out, key = lambda sp: sp["z"])
        return space_points_out

    def load_track_points(self, scifi_event):
        track_points_out = []
        for track in scifi_event.scifitracks():
            if track.tracker() == 0:
                detector = "tku"
            else:
                detector = "tkd"
            for track_point in track.scifitrackpoints():
                if track_point.station() != 1:
                    continue
                position = track_point.pos() # maybe global coordinates?
                momentum = track_point.mom() # maybe global coordinates?
                index = 0
                cov = [[0. for i in range(6)] for j in range(6)]
                index_mapping = [1, 4, 2, 5, 3]
                if track_point.covariance().size() == 25:
                    for i in range(5):
                        k = index_mapping[i]
                        for j in range(5):
                            l = index_mapping[j]
                            cov[k][l] = track_point.covariance()[index]
                            index += 1
                track_points_out.append({
                  "x":position.x(),
                  "y":position.y(),
                  "z":position.z(),
                  "px":momentum.x(),
                  "py":momentum.y(),
                  "pz":momentum.z(),
                  "t":None,
                  "detector":detector,
                  "covariance":cov,
                })
        track_points_out = sorted(track_points_out, key = lambda tp: tp["z"])
        return track_points_out

    def load_scifi_event(self, scifi_event):
        points = []
        will_cut_on_scifi_cluster = self.will_cut_on_scifi_clusters(scifi_event)
        will_cut_on_scifi_space_points = self.will_cut_on_scifi_space_points(scifi_event)
        will_cut_on_scifi_tracks = self.will_cut_on_scifi_tracks(scifi_event)
        will_cut_on_scifi_track_points = self.will_cut_on_scifi_track_points(scifi_event)

        if self.config.will_load_tk_space_points:
            points += self.load_space_points(scifi_event)
        if self.config.will_load_tk_track_points:
            points += self.load_track_points(scifi_event)
        return [
            points, {
                "scifi_space_points":will_cut_on_scifi_space_points,
                "scifi_space_clusters":will_cut_on_scifi_cluster, 
                "scifi_tracks_us":will_cut_on_scifi_tracks[0],
                "scifi_tracks_ds":will_cut_on_scifi_tracks[1],
                "scifi_track_points_us":will_cut_on_scifi_track_points[0],
                "scifi_track_points_ds":will_cut_on_scifi_track_points[1],
            }
        ]

    def will_cut_on_scifi_clusters(self, scifi_event):
        """Require exactly one cluster in each plane"""
        n_clusters = [[[0 for k in range(3)] for j in range(5)] for i in range(2)]
        for cluster in scifi_event.clusters():
            tracker = cluster.get_tracker()
            station = cluster.get_station()-1
            plane = cluster.get_plane()
            n_clusters[tracker][station][plane] += 1
        for i, tracker_cluster in enumerate(n_clusters):
            for station_cluster in tracker_cluster:
                if station_cluster != [1, 1, 1]:
                    self.scifi_cluster_cut_count += 1
                    return True
        return False

    def will_cut_on_scifi_space_points(self, scifi_event):
        """Require exactly one space point in each station"""
        n_space_points = [[0 for j in range(5)] for i in range(2)]
        for space_point in scifi_event.spacepoints():
            tracker = space_point.get_tracker()
            station = space_point.get_station()-1
            if not self.config.will_require_triplets:
                n_space_points[tracker][station] += 1
            # require triplet space point
            elif space_point.get_channels().GetEntries() == 3:
                n_space_points[tracker][station] += 1
        #print n_space_points
        for i, tracker_space_point in enumerate(n_space_points):
            if i in self.config.required_trackers:
                if tracker_space_point != [1, 1, 1, 1, 1]:
                    self.scifi_space_point_cut_count += 1
                    return True
        return False

    def will_cut_on_scifi_track_points(self, scifi_event):
        """Require minimum number of track points in at least one track in each tracker; require track point in station 1"""
        okay = [True, True]
        for track in scifi_event.scifitracks():
            points = track.scifitrackpoints()
            if len(points) >= self.config.required_number_of_track_points:
                for track_point in points:
                    if track_point.station() == 1:
                        okay[track.tracker()] = False
        return okay

        #if okay == [False, False]:
        #    return False # will NOT cut
        #else:
        #    return True

    def will_cut_on_scifi_tracks(self, scifi_event):
        """Require exactly one track in each tracker"""
        n_tracks = [0, 0]
        for track in scifi_event.scifitracks():
            n_tracks[track.tracker()] += 1
        n_tracks = [i != 1 for i in n_tracks]
        return n_tracks

    def load_reco_event(self, reco_event):
        tof_event = reco_event.GetTOFEvent()
        global_event = reco_event.GetGlobalEvent()
        scifi_event = reco_event.GetSciFiEvent()
        tof_loaded = self.load_tof_event(tof_event)
        if self.config.will_require_tof1:
            tofs = [ev["detector"] for ev in tof_loaded[0]]
            if "tof1" not in tofs:
                return None
        scifi_loaded = self.load_scifi_event(scifi_event)
        event = {"data":sorted(scifi_loaded[0]+tof_loaded[0], key = lambda tp: tp['z'])}
        event["particle_number"] = reco_event.GetPartEventNumber()
        
        cuts_chain = itertools.chain(tof_loaded[1].iteritems(), scifi_loaded[1].iteritems())
        cuts = (elem for elem in cuts_chain)
        event["will_cut"] = dict(cuts) # True means cut the thing
        event["any_cut"] = False
        event["data_cut"] = False
        event = self.tm_data(event)
        for key, value in event["will_cut"].iteritems():
            if value and self.config.cuts_active[key]:
                event["any_cut"] = True
                event["data_cut"] = True
        return event

    def will_do_p_cut(self, event):
        p_bins = self.config_anal["p_bins"]
        p_low = min(min(p_bins))
        p_high = max(max(p_bins))
        return event["p_tot"] < p_low or event["p_tot"] > p_high

    def will_do_p_cut_ds(self, event):
        return event["p_tot_ds"] == None or \
               event["p_tot_ds"] < self.config_anal["p_tot_ds_low"] or \
               event["p_tot_ds"] > self.config_anal["p_tot_ds_high"]

    def nan_check(self, float_list):
        if float_list == None:
            return
        for x in float_list:
            if math.isnan(x) or math.isinf(x) or x != x:
                raise ZeroDivisionError("Nan found in tracker: "+str(float_list))
            if x > 150.:
                raise ValueError("Tracker recon out of range: "+str(float_list))

    def tm_data(self, event):
        tku = None
        tkd = None
        tof1 = None
        tof2 = None
        event["p_tot"] = None
        event["p_tot_us"] = None
        event["p_tot_ds"] = None
        for point in event["data"]:
            if point["detector"] == "tku":
                tku = [point["x"], point["px"]/point["pz"], point["y"], point["py"]/point["pz"],]
                if self.config.momentum_from_tracker:
                    event["p_tot"] = (point["px"]**2+point["py"]**2+point["pz"]**2)**0.5 # TODO: get rid of p_tot
                    event["p_tot_us"] = (point["px"]**2+point["py"]**2+point["pz"]**2)**0.5
            elif point["detector"] == "tkd":
                tkd = [point["x"], point["px"]/point["pz"], point["y"], point["py"]/point["pz"],]
                if self.config.momentum_from_tracker:
                    event["p_tot_ds"] = (point["px"]**2+point["py"]**2+point["pz"]**2)**0.5
            elif point["detector"] == "tof1":
                tof1 = point["t"]
            elif point["detector"] == "tof2":
                tof2 = point["t"]
        try:
            event["tof12"] = tof2 - tof1
        except TypeError:
            event["tof12"] = None
        event["data_in"] = tku
        self.nan_check(tku)
        self.nan_check(tkd)
        event["data_out"] = tkd
        event["data_res"] = None
        event["data_proj"] = None
        event["residual_amplitude"] = None
        p_tot, path_length, miss_distance = self.pz_calc.get_p_tot(event)
        if not self.config.momentum_from_tracker:
            event["p_tot"] = p_tot # TODO: get rid of p_tot
            event["p_tot_us"] = p_tot
            event["p_tot_ds"] = p_tot
        event["path_length"] = path_length
        event["miss_distance"] = miss_distance
        event["will_cut"]["p_tot"] = self.will_do_p_cut(event)
        event["will_cut"]["p_tot_ds"] = self.will_do_p_cut_ds(event)
        event["will_cut"]["scraping_cut"] = False
        event["will_cut"]["residual_cut"] = False
        for cut in self.config.lattice_scraping_cuts:
            cut_name = cut["cut_name"]
            event["will_cut"][cut_name] = False
        return event


    cuts = None


