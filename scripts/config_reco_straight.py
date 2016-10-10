class Config(object):
    geometry = "scripts/fit_transfer_matrix/Step_IV_Coils.dat" #"geometry_07534/ParentGeometryFile.dat"
    # prerequisite for data to load
    will_require_tof1 = True
    # prerequisite for space point cut
    will_require_triplets = False #True # require triplet space points
    cuts_active = {
          "any_cut":None,
          "scifi_space_clusters":False,
          "scraping_cut":True,
          "scifi_space_points":False,
          "scifi_tracks":True,
          "scifi_track_points":True,
          "tof12":True,
          "tof_0_sp":True,
          "tof_1_sp":True,
          "tof_2_sp":True,
          "residual_cut":True
    }
    data_dir = "calculated_tms/straight/"

    analyses = [
        {
            "plot_dir":data_dir+"plots_3a/",
            "tof12_cut_low":28.9,
            "tof12_cut_high":29.7,
            "tof_bins":[[28.9+i*0.2, 28.9+(i+1)*0.2] for i in range(4)],
            "reco_files":[
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08030/08030_recon.root",
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08030/08031_recon.root",
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08030/08032_recon.root",
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08030/08033_recon.root",
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08030/08034_recon.root",
            ],
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"3a",
            "color":4,
        },
        {
            "plot_dir":data_dir+"plots_3b/",
            "tof12_cut_low":30.9,
            "tof12_cut_high":31.5,
            "tof_bins":[[30.9+i*0.2, 30.9+(i+1)*0.2] for i in range(3)],
            "reco_files":[
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08035/08035_recon.root",
            ],
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"3b",
            "color":6,
        },
        {
            "plot_dir":data_dir+"plots_3c/",
            "tof12_cut_low":28.2,
            "tof12_cut_high":28.6,
            "tof_bins":[[28.2+i*0.2, 28.2+(i+1)*0.2] for i in range(2)],
            "reco_files":[
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08041/08041_recon.root",
            ],
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"3c",
            "color":8,
        },
        {
            "plot_dir":data_dir+"plots_3d/",
            "tof12_cut_low":29.6,
            "tof12_cut_high":32.0,
            "tof_bins":[[29.6+i*0.4, 29.6+(i+1)*0.4] for i in range(6)],
            "reco_files":[
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08036/08036_recon.root",
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08036/08037_recon.root",
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08036/08038_recon.root",
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08036/08039_recon.root",
              "/home/cr67/MAUS/work/reco/MAUS-v2.5.0/08036/08040_recon.root",
            ],
            "sigma_cut":1.0,
            "chi2_cut":20.,
            "name":"3d",
            "color":2,
        },
    ]
    required_trackers = [0, 1]
    required_number_of_track_points = 5
    number_of_events = 3
    eloss_forwards_model = "bethe_bloch_forwards"
    eloss_backwards_model = "bethe_bloch_backwards"
    will_load_tk_space_points = False
    will_load_tk_track_points = True
    number_of_spills = 50000

    n_subsamples = 10
    min_subsample_size = 5

    tm_plot_dir = data_dir+"/tm_plots/"
    will_do_fit = False
    fit_test_mode = False
    fit_z_tof12 = 8224.8
    fit_plane_us = 0
    fit_plane_ds = 38
    fit_resolution = 0.1
    fit_max_iterations = 20

