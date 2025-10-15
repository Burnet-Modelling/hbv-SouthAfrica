import pathlib
import hbv.constants as const
import hbv.utils 
import hbv.utils_plotting
import hbv.sensitivity


# Set Own Save Directory
# save_path  = 

# Creates directory (if not already existing)
save_dir = pathlib.Path(save_path)
save_dir.mkdir(parents = True, exist_ok = True)

# Run Main Analyses (number of runs can be modified in constants.py script)
# Note: This will return a pickle (.pkl) of extracted model results for post processing.

hbv.utils.hepbd_scenarios_ZAF(const.db_path, const.calib_path, const.calibration, save_dir)

# Data Analysis (Data for Plots and Tables Returned)

hbv.utils.zaf_epi_outcomes(res_save_dir) # Epidemiological Outcomes
hbv.zaf_economic_analysis(res_save_dir, h_disc = 0.03, c_disc = 0.03) # Economic outcomes, returns .pkl for cost-effectiveness fronteirs 

# Generate Outcome Plots
hbv.utils_plotting.cost_eff_frontier(res_save_dir) # Figure 1, Figure 2







country = "ZAF" 





db_path = _get_gitlab_folder()+"databooks/"
res_name = "20250918_add_negvax"
calib_path = _get_gitlab_folder()+"calibrations/revisions/"
calibration = "ZAF_calibration_vcont_test.xlsx" # denotes calibration to be used for this analysis

res_save_dir = pathlib.Path(_get_sharepoint_folder()+f"Applications/South Africa/results/{res_name}") #folder where results get saved
res_save_dir.mkdir(parents=True, exist_ok = True)

# Generate Databook
#hbv.generate_databook.generate_databook(country, savedir=db_path)

# Run central analyses only
## res_store = hbv.utils.central_run_only(db_path, calib_path, calibration)
#Save central run results
## test_save_dir = pathlib.Path(_get_sharepoint_folder()+f"Applications/South Africa/results/tests/")
## test_save_dir.mkdir(parents=True, exist_ok=True)
## save_name = "MTCT_test"
## hbv.utils_plotting.testing_central_plots(res_store, test_save_dir, save_name)

# Run scenarios and simulations
hbv.utils.hepbd_scenarios_ZAF(db_path, calib_path, calibration, res_save_dir)



#Calibration Plots
hbv.utils_plotting.calibration_plots(res_save_dir, calib_path, calibration, db_path)
hbv.utils_plotting.transmission_plots(res_save_dir)
hbv.utils_plotting.valid_plots(res_save_dir, calib_path, calibration, db_path)

#hbv.utils_plotting.nnv_plot(res_save_dir)
#hbv.utils_plotting.cost_eff_frontier(res_save_dir)

owsa_runs = hbv.sensitivity.owsa_mtct_run(db_path, calib_path, calibration)
hbv.sensitivity.owsa_mtct_outcomes(owsa_runs, res_save_dir)