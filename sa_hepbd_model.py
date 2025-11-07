import pathlib
import hbv.constants as const
import hbv.utils
import hbv.utils_plotting
import hbv.sensitivity


# Set Own Save Directory
#save_path  =

# Creates directory (if not already existing)
save_dir = pathlib.Path(save_path)
save_dir.mkdir(parents=True, exist_ok=True)

# Run Main Analyses (number of runs can be modified in constants.py script)
# Note: This will return a pickle (.pkl) of extracted model results for post processing.

hbv.utils.hepbd_scenarios_ZAF(save_dir)

# Data Analysis (Data for Plots and Tables Returned)

hbv.utils.zaf_epi_outcomes(save_dir) # Epidemiological Outcomes
hbv.utils.zaf_economic_analysis(save_dir, h_disc=0.03, c_disc=0.03) # Economic outcomes, returns .pkl for cost-effectiveness fronteirs


# Generate Outcome Plots
hbv.utils_plotting.cost_eff_frontier(save_dir) # Cost Effectiveness Frontiers
hbv.utils_plotting.transmission_plots(save_dir)  # Transmission modality in baseline scenario

# One Way Sensitivity Analysis (main analysis Figure 3)
owsa_runs = hbv.sensitivity.owsa_mtct_run() #run OWSA scenarios
hbv.sensitivity.owsa_mtct_outcomes(owsa_runs, save_dir) # calculate and plot outcomes







