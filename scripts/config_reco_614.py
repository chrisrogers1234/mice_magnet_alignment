def file_names(run_number_list):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/home/cr67/MAUS/work/reco/MAUS-v2.6.1/"+run+"/"+run+"_recon.root"
        file_list.append(a_file)
    return file_list

class Config(object):
    geometry = "geometry_08020/ParentGeometryFile.dat" # "scripts/fit_transfer_matrix/Step_IV_Coils.dat" #
    data_dir = "calculated_tms/6_1_4/"
    # prerequisite for data to load
    will_require_tof1 = True
    # prerequisite for space point cut
    will_require_triplets = False #True # require triplet space points
    cuts_active = {
          "any_cut":None,
          "scifi_space_clusters":False,
          "scraping_cut":False,
          "lattice_scraping_cut_z0":False,
          "lattice_scraping_cut_z400":False,
          "scifi_space_points":False,
          "scifi_tracks_us":True,
          "scifi_track_points_us":True,
          "scifi_tracks_ds":True,
          "scifi_track_points_ds":True,
          "tof12":True,
          "p_tot":True,
          "p_tot_ds":True,
          "tof_0_sp":True,
          "tof_1_sp":True,
          "tof_2_sp":True,
          "residual_cut":False
    }
    lattice_scraping_cuts = []
    #    {"z_pos":0., "max_radius":160., "stay_clear":30., "cut_name":"lattice_scraping_cut_z0"},
    #    {"z_pos":400., "max_radius":160., "stay_clear":30., "cut_name":"lattice_scraping_cut_z400"},
    #]

    analyses = [
        {
            "plot_dir":data_dir+"plots_140MeVc/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "p_bins":[[120.+i*25., 120.+(i+1)*25.] for i in range(1)],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8282, 8276, 8265, 8363, 8372, 8378]),
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"D2: 69.98 A",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_170MeVc/",
            "tof12_cut_low":30.5,
            "tof12_cut_high":34.5,
            "p_bins":[[140.+i*50., 140.+(i+1)*50.] for i in range(1)],
            "p_tot_ds_low":100.,
            "p_tot_ds_high":220.,
            "reco_files":file_names([8280, 8277, 8364, 8365, 8373, 8376]),
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"D2: 82.94 A",
            "color":6,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_200MeVc/",
            "tof12_cut_low":30.,
            "tof12_cut_high":33.2,
            "p_bins":[[170.+i*40., 170.+(i+1)*40.] for i in range(1)],
            "reco_files":file_names([8281, 8366, 8368, 8369, 8374, 8375]),
            "p_tot_ds_low":130.,
            "p_tot_ds_high":250.,
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"D2: 94.92 A",
            "color":8,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_240MeVc/",
            "tof12_cut_low":29.,
            "tof12_cut_high":31.8,
            "p_bins":[[210.+i*40., 210.+(i+1)*40.] for i in range(1)],
            "reco_files":file_names([8367, 8370, 8377]), #8341, 8342, 8344, 
            "p_tot_ds_low":180.,
            "p_tot_ds_high":300.,
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"D2: 110.82 A",
            "color":8,
            "pid":-13, # tof cut -> muon peak
        },
    ]
    required_trackers = [0, 1]
    required_number_of_track_points = 12
    eloss_forwards_model = "bethe_bloch_forwards"
    eloss_backwards_model = "bethe_bloch_backwards"
    will_load_tk_space_points = False
    will_load_tk_track_points = True
    number_of_spills = 50000
    momentum_from_tracker = True

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

    fc_current_seed = +72.00
    e2_us_current_seed = +183.82
    cc_us_current_seed = +205.53
    e1_us_current_seed = +173.66
    m2_us_current_seed = +247.81
    m1_us_current_seed = +99.83
    e1_ds_current_seed = +206.89
    cc_ds_current_seed = +206.89
    e2_ds_current_seed = +206.89

