def file_names(run_number_list):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/"+run+"/"+run+"_recon.root"
        file_list.append(a_file)
    return file_list

class Config(object):
    geometry = "geometry_08020/ParentGeometryFile.dat"
    data_dir = "calculated_tms/fc_pos/"
    # prerequisite for data to load
    will_require_tof1 = True
    # prerequisite for space point cut
    will_require_triplets = False #True # require triplet space points
    cuts_active = {
          "any_cut":None,
          "scifi_space_clusters":False,
          "scraping_cut":True,
          "lattice_scraping_cut_z0":True,
          "lattice_scraping_cut_z400":True,
          "scifi_space_points":False,
          "scifi_tracks":True,
          "scifi_track_points":True,
          "tof12":True,
          "p_tot":True,
          "tof_0_sp":True,
          "tof_1_sp":True,
          "tof_2_sp":True,
          "residual_cut":True
    }
    lattice_scraping_cuts = [
        {"z_pos":0., "max_radius":160., "stay_clear":30., "cut_name":"lattice_scraping_cut_z0"},
        {"z_pos":400., "max_radius":160., "stay_clear":30., "cut_name":"lattice_scraping_cut_z400"},
    ]

    analyses = [
        {
            "plot_dir":data_dir+"plots_200MeVc/",
            "tof12_cut_low":30.5,
            "tof12_cut_high":32.0,
            "p_bins":[[175.+i*5., 175.+(i+1)*5.] for i in range(6)],
            "reco_files":file_names([8020, 8021]),
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"D2: 94.9 A",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_300MeVc/",
            "tof12_cut_low":28.8,
            "tof12_cut_high":29.8,
            "p_bins":[[260.+i*4., 260.+(i+1)*4.] for i in range(10)],
            "reco_files":file_names([8023, 8024]),
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"D2: 104.8 A",
            "color":6,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_400MeVc/",
            "tof12_cut_low":28.2,
            "tof12_cut_high":28.7,
            "p_bins":[[350.+i*10., 350.+(i+1)*10.] for i in range(5)],
            "reco_files":file_names([8025]),
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"D2: 170.0 A",
            "color":8,
            "pid":-13, # tof cut -> muon peak
        }
    ]
    required_trackers = [0, 1]
    required_number_of_track_points = 12
    eloss_forwards_model = "bethe_bloch_forwards"
    eloss_backwards_model = "bethe_bloch_backwards"
    will_load_tk_space_points = False
    will_load_tk_track_points = True
    number_of_spills = 50000

    z_tof1 = 12929.7
    z_tof2 = 21137.5
    z_tku = 15062.0
    z_tkd = 18848.5
    z_fc = (z_tku+z_tkd)/2.

    n_subsamples = 10
    min_subsample_size = 5

    plot_when_no_errors = True
    plot_incl_fit = True
    tm_plot_dir = data_dir+"/tm_plots/"
    will_do_fit = False
    will_do_not_fit = True
    fit_test_mode = False
    fit_z_tof12 = z_tof2 - z_tof1
    fit_z_us = -1818.5-1e-6
    fit_plane_us = 0
    fit_plane_ds = 38
    fit_resolution = 0.1
    fit_max_iterations = 20
    fit_pid = -13
    fit_delta_x = 0.1
    fit_delta_px = 0.1

    fc_current_seed = +50.

