from hbv import constants as const
from gitlab_utils import _get_sharepoint_folder, _get_gitlab_folder
import atomica as at
import pathlib
import pandas as pd
import numpy as np
import at_tools
import sciris as sc

# TODO: For next VIMC run (or any multisetting analysis, update calibration function)

def calibrate_single_country(country , savedir = _get_gitlab_folder()+"calibrations/revisions/"):

    # Make save directory if not already made (use revisions folder as it means nothing accidentally over written)
    savedir = pathlib.Path(savedir)
    savedir.mkdir(parents = True, exist_ok = True)

    # Import pre-requisites for calibration
    F = at.ProjectFramework(const.framework_path)
    P = at.Project(framework= F, databook =_get_gitlab_folder()+"databooks/"+f"{country}_databook.xlsx")
    P.settings.update_time_vector(start=const.sim_start, end = const.sim_end, dt = const.sim_dt)
    cal = P.make_parset()

    # Run calibrations using YAML files (make calib_yaml file and update structuring)
    pop_yaml = pathlib.Path(_get_gitlab_folder()+"yaml_calibrations/hbv_population_calibrate.yaml")
    cal = P.calibrate(cal, yaml = pop_yaml, savedir = savedir, log_output = True, save_intermediate = False)
    cal.save_calibration(savedir/f"{country}_calibration_pop.xlsx")

    infx_yaml = pathlib.Path(_get_gitlab_folder()+"yaml_calibrations/hbv_prevalence_calibrate.yaml")
    cal = P.calibrate(cal, yaml = infx_yaml, savedir = savedir, log_output = True, save_intermediate = False)
    cal.save_calibration(savedir/f"{country}_calibration_pop_prev.xlsx")


    burden_yaml = pathlib.Path(_get_gitlab_folder()+"yaml_calibrations/hbv_burden_calibrate.yaml")
    cal = P.calibrate(cal, yaml = burden_yaml, savedir = savedir, log_output = True, save_intermediate = False)
    cal.save_calibration(savedir/f"{country}_calibration.xlsx")


def central_run_only(db_path, calib_path, calibration):
    """
    Returns results items in a dictionary for plotting/datachecking - useful when wanting to try things, re-calibrate etc
    (no need to do 200 runs per thing).
    """

    F = at.ProjectFramework(const.framework_path)
    db_path = pathlib.Path(db_path+"ZAF_databook.xlsx")
    cal_dir = pathlib.Path(calib_path + calibration)

    res_store = {}
    for scen in const.zaf_scenarios:
        res_store[scen]={}
        res_store[scen]["P"] = []
        res_store[scen]["res"] = []

    for scen in const.zaf_scenarios:
        # Set up each modelled scenario
        if scen == "baseline":
            D = at.ProjectData.from_spreadsheet(db_path, framework= F) # no changes required
        # Selective HepB-BD, Current Guidelines
        if scen == "pessimistic_selective_SA":
            D = at.ProjectData.from_spreadsheet(db_path, framework= F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Targeted (selective) HepB-BD coverage to 40%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
        if scen == "optimistic_selective_SA":
            D = at.ProjectData.from_spreadsheet(db_path, framework= F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Targeted (selective) HepB-BD coverage to 40%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
        # Selective HepB-BD, WHO Guidelines
        if scen == "pessimistic_selective_WHO":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to one (treat all HBsAg+ diagnosed mothers)
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 1
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 1
            # Targeted (selective) HepB-BD coverage to 40%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
        if scen == "optimistic_selective_WHO":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to one (treat all HBsAg+ diagnosed mothers)
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 1
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 1
            # Targeted (selective) HepB-BD coverage to 80%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
        # Universal HepB-BD
        if scen == "pessimistic_universal":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Set universal HepB-BD coverage to 40%
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
        if scen == "optimistic_universal":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Set Universal HepB-BD coverage to 80%
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
        # Selective + Universal, Current Guidelines
        if scen == "pessimistic_combined_SA":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Set Targeted and Universal HepB-BD coverage to 40%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
        if scen == "optimistic_combined_SA":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Set Targeted and Universal HepB-BD coverage to 80%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
        # Selective + Universal, WHO Guidelines
        if scen == "pessimistic_combined_WHO":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 1
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 1
            # Set Targeted and Universal HepB-BD coverage to 40%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
        if scen == "optimistic_combined_WHO":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 1
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 1
            # Set Targeted and Universal HepB-BD coverage to 80%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")

        P = at.Project(framework=F, databook=D, do_run=False, sim_start=const.sim_start,
                       sim_end=const.sim_end, sim_dt=const.sim_dt)
        cal = P.parsets[0].load_calibration(cal_dir)
        res = P.run_sim(parset=cal, result_name=f"{scen}")

        res_store[scen]["P"].append(P)
        res_store[scen]["res"].append(res)
        
    return res_store


def hepbd_scenarios_ZAF(db_path, calib_path, calibration, res_save_dir):

    """
    This function carries out main analyses of HepB-BD South Africa modelling study. Not aimed to be a transferrable
    function.
    """
    # Sample transition and effectiveness parameters (all triangular)
    np.random.seed(const.seed)  # sets seed for random sampling
    trans_pars = pd.read_excel(_get_sharepoint_folder()+"Data\\hbv_flat datasheet.xlsx", sheet_name= "Global Pars")
    tpar_names = list(pd.unique(trans_pars.par))

    sampled_pars = pd.DataFrame(columns = ["pars", "central"] + [f"run_{i}" for i in range(1, const.runs+1, 1)])
    sampled_pars.pars = trans_pars.par  # parameter names
    sampled_pars.central = trans_pars.val  # central (non-sampled) values

    # If parameter has uncertainty will draw from triangular distribution, otherwise just repeats the central value for
    # all runs. Maps directly to trans_par sheet with Y/N for populations.
    for par in range(len(trans_pars)):
        if trans_pars.iloc[par,3] != trans_pars.iloc[par,4]:
            sampled_pars.iloc[par, 2:] = np.random.triangular(trans_pars.iloc[par,3], trans_pars.iloc[par,2], trans_pars.iloc[par,4], size= const.runs)
        else:
            sampled_pars.iloc[par, 2:const.runs+2] = trans_pars.iloc[par,2]

    # Framework, Databook and Calibration Paths for Running Model
    F = at.ProjectFramework(const.framework_path)
    db_path = pathlib.Path(db_path+"ZAF_databook.xlsx")
    cal_dir = pathlib.Path(calib_path + calibration)

    # Generate storage dictionary to save model runs
    res_store = {}
    for scen in const.zaf_scenarios:
        res_store[scen] = {}
        for store in list(const.zaf_data_extract.keys()):
            res_store[scen][store] = pd.DataFrame(columns = ["year", "central"] + [f"run_{run}" for run in range(1, const.runs +1, 1)])
            res_store[scen][store].year = np.arange(1990.5, 2101.5, 1)

    # Run the model (central and simulations)
    for scen in const.zaf_scenarios:
        # Set up each modelled scenario
        if scen == "baseline":
            D = at.ProjectData.from_spreadsheet(db_path, framework= F) # no changes required
        # Selective HepB-BD, Current Guidelines
        if scen == "pessimistic_selective_SA":
            D = at.ProjectData.from_spreadsheet(db_path, framework= F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Targeted (selective) HepB-BD coverage to 40%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
        if scen == "optimistic_selective_SA":
            D = at.ProjectData.from_spreadsheet(db_path, framework= F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Targeted (selective) HepB-BD coverage to 40%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
        # Selective HepB-BD, WHO Guidelines
        if scen == "pessimistic_selective_WHO":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to one (treat all HBsAg+ diagnosed mothers)
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 1
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 1
            # Targeted (selective) HepB-BD coverage to 40%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
        if scen == "optimistic_selective_WHO":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to one (treat all HBsAg+ diagnosed mothers)
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 1
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 1
            # Targeted (selective) HepB-BD coverage to 80%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
        # Universal HepB-BD
        if scen == "pessimistic_universal":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Set universal HepB-BD coverage to 40%
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
        if scen == "optimistic_universal":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Set Universal HepB-BD coverage to 80%
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
        # Selective + Universal, Current Guidelines
        if scen == "pessimistic_combined_SA":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Set Targeted and Universal HepB-BD coverage to 40%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
        if scen == "optimistic_combined_SA":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
            # Set Targeted and Universal HepB-BD coverage to 80%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
        # Selective + Universal, WHO Guidelines
        if scen == "pessimistic_combined_WHO":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 1
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 1
            # Set Targeted and Universal HepB-BD coverage to 40%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
        if scen == "optimistic_combined_WHO":
            D = at.ProjectData.from_spreadsheet(db_path, framework=F)
            # Set mav_trtall to zero
            D.tdve["mav_trtall"].ts["0-0M"].assumption = 1
            D.tdve["mav_trtall"].ts["0-0F"].assumption = 1
            # Set Targeted and Universal HepB-BD coverage to 80%
            D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")

        # Run each model simulation per scenario
        for run in range(const.runs + 1):
            for name in tpar_names:
                # Isolate each paramaeter
                refs_par = trans_pars[trans_pars.par == name]
                samp_par = sampled_pars[sampled_pars.pars == name]
                for bin in list(const.pop_bins.keys()):
                    for i in range(len(refs_par)):
                        if refs_par[f"{bin}M"].iloc[i] == "Y":
                            D.tdve[name].ts[f"{bin}M"].assumption = samp_par.iloc[i, run+1]
                        if refs_par[f"{bin}F"].iloc[i] == "Y":
                            D.tdve[name].ts[f"{bin}F"].assumption = samp_par.iloc[i, run + 1]

            P = at.Project(framework = F, databook = D, do_run = False, sim_start = const.sim_start,
                               sim_end = const.sim_end, sim_dt = const.sim_dt)
            cal = P.parsets[0].load_calibration(cal_dir)
            res = P.run_sim(parset=cal, result_name=f"{scen}_{run}")

            for store, vals in const.zaf_data_extract.items():
                res_store[scen][store].iloc[:, run+1] =  at.PlotData(res, vals[0], pops=vals[1], pop_aggregation=vals[2], t_bins=1).series[0].vals

            print(f"Run {run} of scenario {scen} completed")

    filename = sc.makefilepath(filename=f"ZAF_{const.runs}_sims.pkl", folder = res_save_dir/'result_pkls', makedirs= True)
    sc.save(filename=filename, obj=res_store)

def zaf_epi_outcomes(res_save_dir):
    scenarios = const.zaf_scenarios
    comparators = scenarios[1:] # to get differences
    data = sc.load(res_save_dir/'result_pkls'/f"ZAF_{const.runs}_sims.pkl")

    # Make Save Directory
    fig_tab_dir = res_save_dir/'tab_figs'
    fig_tab_dir.mkdir(parents=True, exist_ok=True)

    # Data for Table 2
    epi_measures = ["mtct_infx", "tot_inc", "hcc_incidence", "hbv_deaths"]

    epi_tot_short = pd.DataFrame(columns = ["scenarios"] + [f"{epi}_{bd}" for epi in epi_measures for bd in ["cent", "lb", "ub"]])
    epi_tot_short.scenarios = scenarios
    epi_tot_short = epi_tot_short.set_index("scenarios", drop = True)

    epi_comp_short = pd.DataFrame(columns = ["scenarios"] + [f"{epi}_{bd}" for epi in epi_measures for bd in ["cent", "lb", "ub"]])
    epi_comp_short.scenarios = scenarios
    epi_comp_short = epi_comp_short.set_index("scenarios", drop = True)

    for scen in scenarios:
        for meas in epi_measures:
            epi_tot_short.at[scen, f"{meas}_cent"] = np.sum(data[scen][meas].iloc[34:61, 1])
            epi_tot_short.at[scen, f"{meas}_lb"] = np.percentile(np.sum(data[scen][meas].iloc[34:61, 2:]), 2.5)
            epi_tot_short.at[scen, f"{meas}_ub"] = np.percentile(np.sum(data[scen][meas].iloc[34:61, 2:]), 97.5)

    for comp in comparators:
        for meas in epi_measures:
            epi_comp_short.at[comp, f"{meas}_cent"] = np.sum(data["baseline"][meas].iloc[34:61, 1]) -  np.sum(data[comp][meas].iloc[34:61, 1])
            epi_comp_short.at[comp, f"{meas}_lb"] = np.percentile(np.sum(data["baseline"][meas].iloc[34:61, 2:]) -  np.sum(data[comp][meas].iloc[34:61, 2:]), 2.5)
            epi_comp_short.at[comp, f"{meas}_ub"] = np.percentile(np.sum(data["baseline"][meas].iloc[34:61, 2:]) -  np.sum(data[comp][meas].iloc[34:61, 2:]), 97.5)

    epi_tot_long = pd.DataFrame(columns = ["scenarios"] + [f"{epi}_{bd}" for epi in epi_measures for bd in ["cent", "lb", "ub"]])
    epi_tot_long.scenarios = scenarios
    epi_tot_long = epi_tot_long.set_index("scenarios", drop = True)

    epi_comp_long = pd.DataFrame(columns = ["scenarios"] + [f"{epi}_{bd}" for epi in epi_measures for bd in ["cent", "lb", "ub"]])
    epi_comp_long.scenarios = scenarios
    epi_comp_long = epi_comp_long.set_index("scenarios", drop = True)

    for scen in scenarios:
        for meas in epi_measures:
            epi_tot_long.at[scen, f"{meas}_cent"] = np.sum(data[scen][meas].iloc[34:, 1])
            epi_tot_long.at[scen, f"{meas}_lb"] = np.percentile(np.sum(data[scen][meas].iloc[34:, 2:]), 2.5)
            epi_tot_long.at[scen, f"{meas}_ub"] = np.percentile(np.sum(data[scen][meas].iloc[34:, 2:]), 97.5)

    for comp in comparators:
        for meas in epi_measures:
            epi_comp_long.at[comp, f"{meas}_cent"] = np.sum(data["baseline"][meas].iloc[34:, 1]) -  np.sum(data[comp][meas].iloc[34:, 1])
            epi_comp_long.at[comp, f"{meas}_lb"] = np.percentile(np.sum(data["baseline"][meas].iloc[34:, 2:]) -  np.sum(data[comp][meas].iloc[34:, 2:]), 2.5)
            epi_comp_long.at[comp, f"{meas}_ub"] = np.percentile(np.sum(data["baseline"][meas].iloc[34:, 2:]) -  np.sum(data[comp][meas].iloc[34:, 2:]), 97.5)

    # Export to a consolidated excel spreadsheet
    with pd.ExcelWriter(res_save_dir/'tab_figs'/'table2epi.xlsx', engine='xlsxwriter') as writer:
        epi_tot_short.to_excel(writer, sheet_name = "Sum Res_2050", index=True)
        epi_tot_long.to_excel(writer, sheet_name = "Sum Res_2100", index=True)
        epi_comp_short.to_excel(writer, sheet_name = "Dif Res_2050", index=True)
        epi_comp_long.to_excel(writer, sheet_name = "Dif Res_2100", index=True)

    # TODO: Make a dictionary and save as an xlsx with sheet for each scenario
    # Data for Figure 2 and associated discussion (save as csv so plots can be modified without needing to re-run code)
    fig2_data = {}
    for scen in const.zaf_scenarios:
        fig2_data[scen] =  pd.DataFrame(columns = ["year", "mtct", "echt", "oth", "tot_sum", "validation"])
        fig2_data[scen].year = data[scen]["mtct_infx"].year
        fig2_data[scen].mtct = data[scen]["mtct_infx"].central
        fig2_data[scen].echt = data[scen]["echt_infx"].central
        fig2_data[scen].oth = data[scen]["oth_infx"].central
        fig2_data[scen].tot_sum = fig2_data[scen].mtct + fig2_data[scen].echt + fig2_data[scen].oth
        fig2_data[scen].validation = fig2_data[scen].tot_sum - data[scen]["tot_inc"].central

    excel_file_path = res_save_dir/'tab_figs'/'fig2data.xlsx'
    writer = pd.ExcelWriter(excel_file_path, engine='xlsxwriter')
    for sheet_name, df in fig2_data.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.close()

    # Data used in results paragraph
    context_data = np.zeros((3,3))
    # CHB infection rate among all births
    context_data[0,0] = np.sum(data["baseline"]["mtct_infx"].iloc[34:, 1])/np.sum(data["baseline"]["births"].iloc[34:, 1])
    context_data[0,1] = np.percentile(np.sum(data["baseline"]["mtct_infx"].iloc[34:, 2:])/np.sum(data["baseline"]["births"].iloc[34:, 2:]), 2.5)
    context_data[0,2] = np.percentile(np.sum(data["baseline"]["mtct_infx"].iloc[34:, 2:])/np.sum(data["baseline"]["births"].iloc[34:, 2:]), 97.5)
    # CHB infection rate among HBsAg+ mothers
    context_data[1,0] = np.sum(data["baseline"]["mtct_infx"].iloc[34:, 1])/np.sum(data["baseline"]["expos_births"].iloc[34:, 1])
    context_data[1,1] = np.percentile(np.sum(data["baseline"]["mtct_infx"].iloc[34:, 2:])/np.sum(data["baseline"]["expos_births"].iloc[34:, 2:]), 2.5)
    context_data[1,2] = np.percentile(np.sum(data["baseline"]["mtct_infx"].iloc[34:, 2:])/np.sum(data["baseline"]["expos_births"].iloc[34:, 2:]), 97.5)
    # MTCT proportion among all infections
    context_data[2,0] = np.sum(data["baseline"]["mtct_infx"].iloc[34:, 1])/np.sum(data["baseline"]["tot_inc"].iloc[34:, 1])
    context_data[2,1] = np.percentile(np.sum(data["baseline"]["mtct_infx"].iloc[34:, 2:])/np.sum(data["baseline"]["tot_inc"].iloc[34:, 2:]), 2.5)
    context_data[2,2] = np.percentile(np.sum(data["baseline"]["mtct_infx"].iloc[34:, 2:])/np.sum(data["baseline"]["tot_inc"].iloc[34:, 2:]), 97.5)

    column_headers = "Central, Lower, Upper"
    description = "CHB infection rate among all births (Row 1), CHB infection from HBsAg+ mothers (Row 2), MTCT among all incident infections (Row 3), from 2025-2100"

    np.savetxt(res_save_dir/'tab_figs'/'context_data.csv', context_data, delimiter=",", header = column_headers,
               comments = f'# {description} \n')

    # Data for Figure 3 and associated discussion (number needed to vaccinate)
    nnv_plot = pd.DataFrame(columns = ["scenario"] + [f"{nnv}_{bd}" for nnv in epi_measures for bd in ["cent", "lb", "ub"]])
    nnv_plot.scenario = scenarios
    nnv_plot = nnv_plot.set_index("scenario", drop=True)

    for comp in comparators:
        for meas in epi_measures:
            nnv_plot.at[comp, f"{meas}_cent"] = (np.sum(data[comp]["total_hepbd"].iloc[34:, 1])-np.sum(data["baseline"]["total_hepbd"].iloc[34:, 1]))/ (np.sum(data["baseline"][meas].iloc[34:, 1])-np.sum(data[comp][meas].iloc[34:, 1]))
            nnv_plot.at[comp, f"{meas}_lb"] = np.percentile((np.sum(data[comp]["total_hepbd"].iloc[34:, 2:])-np.sum(data["baseline"]["total_hepbd"].iloc[34:, 2:]))/ (np.sum(data["baseline"][meas].iloc[34:, 2:])-np.sum(data[comp][meas].iloc[34:, 2:])), 2.5)
            nnv_plot.at[comp, f"{meas}_ub"] = np.percentile((np.sum(data[comp]["total_hepbd"].iloc[34:, 2:])-np.sum(data["baseline"]["total_hepbd"].iloc[34:, 2:]))/ (np.sum(data["baseline"][meas].iloc[34:, 2:])-np.sum(data[comp][meas].iloc[34:, 2:])), 97.5)

    nnv_plot.to_csv(res_save_dir/'tab_figs'/'fig3data.csv', index=True)

    # Data for alternate Figure 3 (univ HepB-BD applied to HBsAg negative screened births)
    nnv_plot_add =  pd.DataFrame(columns = ["scenario"] + [f"{nnv}_{bd}" for nnv in epi_measures for bd in ["cent", "lb", "ub"]])
    nnv_plot_add.scenario = scenarios
    nnv_plot_add = nnv_plot_add.set_index("scenario", drop=True)

    for comp in comparators:
        for meas in epi_measures:
            nnv_plot_add.at[comp, f"{meas}_cent"] = ((np.sum(data[comp]["total_hepbd"].iloc[34:, 1])+np.sum(data[comp]["select_hepbd_univ"].iloc[34:, 1]))-np.sum(data["baseline"]["total_hepbd"].iloc[34:, 1]))/ (np.sum(data["baseline"][meas].iloc[34:, 1])-np.sum(data[comp][meas].iloc[34:, 1]))
            nnv_plot_add.at[comp, f"{meas}_lb"] = np.percentile(((np.sum(data[comp]["total_hepbd"].iloc[34:, 2:])+np.sum(data[comp]["select_hepbd_univ"].iloc[34:, 2:]))-np.sum(data["baseline"]["total_hepbd"].iloc[34:, 2:]))/ (np.sum(data["baseline"][meas].iloc[34:, 2:])-np.sum(data[comp][meas].iloc[34:, 2:])), 2.5)
            nnv_plot_add.at[comp, f"{meas}_ub"] = np.percentile(((np.sum(data[comp]["total_hepbd"].iloc[34:, 2:])+np.sum(data[comp]["select_hepbd_univ"].iloc[34:, 2:]))-np.sum(data["baseline"]["total_hepbd"].iloc[34:, 2:]))/ (np.sum(data["baseline"][meas].iloc[34:, 2:])-np.sum(data[comp][meas].iloc[34:, 2:])), 97.5)

    nnv_plot_add.to_csv(res_save_dir/'tab_figs'/'fig3data_add.csv', index=True)

    # Supplemental Plot: WHO Targets (mortality and incidence per 100k) from each scenario
    who_plot = {}
    who_measures = ["who_inc_u5", "who_inc_pop", "who_mort"]

    for whom in who_measures:
        who_plot[whom] = pd.DataFrame(columns = ["year"]+[f"{scen}_{bound}" for scen in const.zaf_scenarios for bound in ["cent", "lb", "ub"]])
        who_plot[whom].year = np.arange(1990.5, 2101.5, 1)

        for scen in const.zaf_scenarios:
            who_plot[whom][f"{scen}_cent"] = data[scen][whom].iloc[:, 1]
            who_plot[whom][f"{scen}_lb"] = np.percentile(data[scen][whom].iloc[:, 2:], 2.5, axis = 1)
            who_plot[whom][f"{scen}_ub"] = np.percentile(data[scen][whom].iloc[:, 2:], 97.5, axis = 1)

    excel_file_path = res_save_dir/'tab_figs'/'whotarget_data.xlsx'
    writer = pd.ExcelWriter(excel_file_path, engine='xlsxwriter')
    for sheet_name, df in who_plot.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.close()

def zaf_economic_analysis(res_save_dir, h_disc = 0.03, c_disc = 0.03):

    # Extract costs from flat datasheet
    flat_data = _get_sharepoint_folder() + "Data\\hbv_flat datasheet.xlsx"
    costs = pd.read_excel(flat_data, sheet_name="HBV Costs")

    # Conversion to USD
    usd_conversion = 1/costs[costs.cost_comp=="currency_conv"].iloc[0,4]

    np.random.seed(const.seed) # set seed for reproducibility
    # Sample DALY weights (triangular)
    daly_weights = costs[(costs.cost_comp=="dc_daly")|(costs.cost_comp=="hcc_daly")]
    daly_sample = pd.DataFrame(columns=["par", "central"] + [f"run_{run}" for run in range(1, const.runs + 1, 1)])
    daly_sample.par = daly_weights.cost_comp
    for i, daly in enumerate(list(pd.unique(daly_sample.par))):
        daly_sample.iloc[i, 1] = daly_weights.iloc[i, 3]
        daly_sample.iloc[i, 2:] = np.random.triangular(daly_weights.iloc[i, 4], daly_weights.iloc[i, 3],
                                                   daly_weights.iloc[i, 5], size=const.runs)

    # Sample Costs (uniform)
    unit_costs = costs[(costs.unit=="ZAR") & (costs.cost_comp!="currency_conv")]
    cost_sample = pd.DataFrame(columns=["par", "central"] + [f"run_{run}" for run in range(1, const.runs + 1, 1)])
    cost_sample.par = unit_costs.cost_comp
    for i, cost in enumerate(list(pd.unique(cost_sample.par))):
        cost_sample.iloc[i, 1] = unit_costs.iloc[i, 3]
        cost_sample.iloc[i, 2:] = np.random.uniform(unit_costs.iloc[i, 4], unit_costs.iloc[i, 5], size=const.runs)
    cost_sample = cost_sample.set_index('par', drop=False)

    # Discount Arrays
    data = sc.load(res_save_dir / 'result_pkls' / f"ZAF_{const.runs}_sims.pkl")
    disc_array = pd.DataFrame(columns=["year", "consumables", "health"])
    disc_array.year = data["baseline"]["total_hepbd"].year

    for i in range(len(disc_array)):
        if disc_array.iloc[i, 0] < 2024.5:
            disc_array.iloc[i, 1] = 0
            disc_array.iloc[i, 2] = 0
        else:
            disc_array.iloc[i, 1] = (1 - c_disc) ** (disc_array.iloc[i, 0] - 2024.5)
            disc_array.iloc[i, 2] = (1 - h_disc) ** (disc_array.iloc[i, 0] - 2024.5)

    # Estimate scenario costs and DALYs from model data
    scenarios = const.zaf_scenarios
    comparators = const.zaf_scenarios[1:]  # removes baseline when running comparisons, needed for ICERs

    # Cost Estimates
    cost_comps = ["hepbd", "hepbd_add", "mat_scr", "mat_av", "hepb3", "diag", "treat",
                  "ncirr", "cirr", "decomp", "canc", "bd_sum", "bd_sum_add", "hbv_sum",
                  "total", "total_add"]
    costs = {}
    for scen in scenarios:
        costs[scen] = {}
        for ccomp in cost_comps:
            costs[scen][ccomp] = pd.DataFrame(columns = ["year", "central"] + [f"run_{run}" for run in range(1, const.runs+1, 1)])
            costs[scen][ccomp].year = data["baseline"]["total_hepbd"].year

    # HepB-BD costs (HepB-BD vaccines administered at birth)
    for scen in scenarios:
        if (scen == "baseline"
            or scen == "pessimistic_universal"
            or scen == "optimistic_universal"
            or scen == "pessimistic_combined_SA"
            or scen == "optimistic_combined_SA"
            or scen == "pessimistic_combined_WHO"
            or scen == "optimistic_combined_WHO"):
            for run in range(const.runs+1):
                costs[scen]["hepbd"].iloc[:, run+1] = data[scen]["total_hepbd"].iloc[:, run+1]*cost_sample.loc["hepb_bd_univ"][run+1]*disc_array.iloc[:, 1]
        elif (scen == "pessimistic_selective_SA"
            or scen == "optimistic_selective_SA"
            or scen == "pessimistic_selective_WHO"
            or scen == "optimistic_selective_WHO"):
            for run in range(const.runs+1):
                costs[scen]["hepbd"].iloc[:, run+1] = data[scen]["total_hepbd"].iloc[:, run+1]*cost_sample.loc["hepb_bd_tgt"][run+1]*disc_array.iloc[:, 1]

    # HepB-BD costs (additional if universal vax also applies to test negative)
    for scen in scenarios:
        if (scen == "baseline"
            or scen == "pessimistic_universal"
            or scen == "optimistic_universal"):
            for run in range(const.runs+1):
                costs[scen]["hepbd_add"].iloc[:, run+1] = data[scen]["total_hepbd"].iloc[:, run+1]*cost_sample.loc["hepb_bd_univ"][run+1]*disc_array.iloc[:, 1]
        elif (scen == "pessimistic_selective_SA"
            or scen == "optimistic_selective_SA"
            or scen == "pessimistic_selective_WHO"
            or scen == "optimistic_selective_WHO"):
            for run in range(const.runs+1):
                costs[scen]["hepbd_add"].iloc[:, run+1] = data[scen]["total_hepbd"].iloc[:, run+1]*cost_sample.loc["hepb_bd_tgt"][run+1]*disc_array.iloc[:, 1]
        elif (scen == "pessimistic_combined_SA"
            or scen == "optimistic_combined_SA"
            or scen == "pessimistic_combined_WHO"
            or scen == "optimistic_combined_WHO"):
            for run in range(const.runs+1):
                costs[scen]["hepbd_add"].iloc[:, run+1] = (data[scen]["total_hepbd"].iloc[:, run+1]+data[scen]["select_hepbd_univ"].iloc[:, run+1])*cost_sample.loc["hepb_bd_univ"][run+1]*disc_array.iloc[:, 1]

    # HepB3 costs
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["hepb3"].iloc[:, run+1] = data[scen]["total_hepb3"].iloc[:, run+1]*cost_sample.loc["hepb3"][run+1]*disc_array.iloc[:, 1]

    # Maternal Screening Costs (assumed rapid test with or without HBV DNA)
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["mat_scr"].iloc[:, run+1] = (data[scen]["anc_screening"].iloc[:, run+1]*cost_sample.loc["diag_rapid"][run+1] + data[scen]["dna_screening"].iloc[:, run+1]*cost_sample.loc["hbv_dna"][run+1])*disc_array.iloc[:, 1]

    # Maternal Antiviral Costs
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["mat_av"].iloc[:, run+1] = data[scen]["mav_hepbd_cov"].iloc[:, run+1]*cost_sample.loc["tdf_pap"][run+1]*disc_array.iloc[:, 1]

    # Diagnosis costs (assumed rapid test)
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["diag"].iloc[:, run+1] = data[scen]["diag_cost"].iloc[:, run+1]*cost_sample.loc["diag_rapid"][run+1]*disc_array.iloc[:, 1]

    # Treatment costs
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["treat"].iloc[:, run+1] = data[scen]["treat_cost"].iloc[:, run+1]*cost_sample.loc["tdf_treat"][run+1]*disc_array.iloc[:, 1]

    # Non cirrhotic costs
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["ncirr"].iloc[:, run+1] = data[scen]["noncirr_cost"].iloc[:, run+1]*cost_sample.loc["chb_cost"][run+1]*disc_array.iloc[:, 1]

    # Cirrhotic costs
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["cirr"].iloc[:, run+1] = data[scen]["compcirr_cost"].iloc[:, run+1]*cost_sample.loc["cc_cost"][run+1]*disc_array.iloc[:, 1]

    # Decompensated costs
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["decomp"].iloc[:, run+1] = data[scen]["decomp_cost"].iloc[:, run+1]*cost_sample.loc["dc_cost"][run+1]*disc_array.iloc[:, 1]

    # HCC costs
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["canc"].iloc[:, run+1] = data[scen]["cancer_cost"].iloc[:, run+1]*cost_sample.loc["hcc_cost"][run+1]*disc_array.iloc[:, 1]

    # HepB-BD associated costs
    for scen in scenarios:
        for run in range(const.runs+1):
            costs[scen]["bd_sum"].iloc[:, run+1] = costs[scen]["hepbd"].iloc[:, run+1] + costs[scen]["mat_scr"].iloc[:, run+1] + costs[scen]["mat_av"].iloc[:, run+1]

    # Hep-BD associated costs (add - univ HBsAg vaccine)
    for scen in scenarios:
        for run in range(const.runs+1):
            costs[scen]["bd_sum_add"].iloc[:, run+1] = costs[scen]["hepbd_add"].iloc[:, run+1] + costs[scen]["mat_scr"].iloc[:, run+1] + costs[scen]["mat_av"].iloc[:, run+1]

    # Other hepatitis costs
    for scen in scenarios:
        for run in range(const.runs+1):
            costs[scen]["hbv_sum"].iloc[:, run+1] = costs[scen]["hepb3"].iloc[:, run+1] + costs[scen]["diag"].iloc[:, run+1] + costs[scen]["treat"].iloc[:, run+1] + costs[scen]["cirr"].iloc[:, run+1] + costs[scen]["decomp"].iloc[:, run+1] + costs[scen]["canc"].iloc[:, run+1]+ costs[scen]["ncirr"].iloc[:, run+1]

    # Total Costs
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["total"].iloc[:, run+1] = costs[scen]["bd_sum"].iloc[:, run+1] + costs[scen]["hbv_sum"].iloc[:, run+1]

    # Total Costs (add - univ HBsAg vaccine)
    for scen in scenarios:
        for run in range(const.runs + 1):
            costs[scen]["total_add"].iloc[:, run+1] = costs[scen]["bd_sum_add"].iloc[:, run+1] + costs[scen]["hbv_sum"].iloc[:, run+1]
    # DALYs
    dalys = {}
    daly_comps = ["dc_yld", "hcc_yld", "yll", "dalys"]
    for scen in scenarios:
        dalys[scen] = {}
        for dcomp in daly_comps:
            dalys[scen][dcomp] = pd.DataFrame(columns = ["year", "central"] + [f"run_{run}" for run in range(1,const.runs+1,1)])
            dalys[scen][dcomp].year = data["baseline"]["total_hepbd"].year

    for scen in scenarios:
        # ICER under assumption HBsAg- screened do not recieve HepB-BD under universal HepB-BD
        dalys[scen]["icer"] = pd.DataFrame(columns = ["time_horizon", "central"] + [f"run_{run}" for run in range(1,const.runs+1,1)])
        dalys[scen]["icer"].time_horizon = ["2025-2050", "2025-2100"]
        #ICER under assumption HBsAg- screened do recieve HepB-BD under universal HepB-BD
        dalys[scen]["icer_add"] = pd.DataFrame(columns = ["time_horizon", "central"] + [f"run_{run}" for run in range(1,const.runs+1,1)])
        dalys[scen]["icer_add"].time_horizon = ["2025-2050", "2025-2100"]



    # YLD and YLL (discounted)
    for scen in scenarios:
        for run in range(const.runs+1):
            dalys[scen]["dc_yld"].iloc[:, run+1] = data[scen]["decomp_daly"].iloc[:, run+1]*daly_sample.iloc[0, run+1]*disc_array.iloc[:,2]
            dalys[scen]["hcc_yld"].iloc[:, run+1] = data[scen]["cancer_daly"].iloc[:, run+1]*daly_sample.iloc[1, run+1]*disc_array.iloc[:,2]
            dalys[scen]["yll"].iloc[:, run+1] = data[scen]["yll_daly"].iloc[:, run+1]*disc_array.iloc[:,2]

    # DALYS
    for scen in scenarios:
        for run in range(const.runs+1):
            dalys[scen]["dalys"].iloc[:, run+1] =  dalys[scen]["dc_yld"].iloc[:, run+1] + dalys[scen]["hcc_yld"].iloc[:, run+1] + dalys[scen]["yll"].iloc[:, run+1]

    ## Outcomes and Outputs
    # TODO: Add in outcomes where HepB-BD vaccine under combination strategies is provided to HBsAg- births

    # ICERs (cost per DALY averted)
    for comps in comparators:
        for run in range(const.runs+1):
            #ICER under assumption HBsAg- screened do not recieve HepB-BD under universal HepB-BD
            dalys[comps]["icer"].iloc[0, run+1] = (np.sum(costs[comps]["total"].iloc[34:61, run+1])-np.sum(costs["baseline"]["total"].iloc[34:61, run+1]))/(np.sum(dalys["baseline"]["dalys"].iloc[34:61, run+1])-np.sum(dalys[comps]["dalys"].iloc[34:61, run+1]))
            dalys[comps]["icer"].iloc[1, run+1] = (np.sum(costs[comps]["total"].iloc[34:, run+1])-np.sum(costs["baseline"]["total"].iloc[34:, run+1]))/(np.sum(dalys["baseline"]["dalys"].iloc[34:, run+1])-np.sum(dalys[comps]["dalys"].iloc[34:, run+1]))
            #ICER under assumption HBsAg- screened do recieve HepB-BD under universal HepB-BD
            dalys[comps]["icer_add"].iloc[0, run+1] = (np.sum(costs[comps]["total_add"].iloc[34:61, run+1])-np.sum(costs["baseline"]["total_add"].iloc[34:61, run+1]))/(np.sum(dalys["baseline"]["dalys"].iloc[34:61, run+1])-np.sum(dalys[comps]["dalys"].iloc[34:61, run+1]))
            dalys[comps]["icer_add"].iloc[1, run+1] = (np.sum(costs[comps]["total_add"].iloc[34:, run+1])-np.sum(costs["baseline"]["total_add"].iloc[34:, run+1]))/(np.sum(dalys["baseline"]["dalys"].iloc[34:, run+1])-np.sum(dalys[comps]["dalys"].iloc[34:, run+1]))

    # Outputs for Table 2
    cost_measure = ["bd_sum", "bd_sum_add", "hbv_sum", "total", "total_add"]
    util_measure = ["dalys"]
    eff_measure = ["icer", "icer_add"]

    tot_short = pd.DataFrame(columns = ["scenario"] + [f"{meas}_{bnd}" for meas in cost_measure for bnd in ["cent", "lb", "ub"]]  + [f"{meas}_{bnd}" for meas in util_measure for bnd in ["cent", "lb", "ub"]])
    tot_short.scenario = const.zaf_scenarios
    tot_short = tot_short.set_index("scenario", drop=True)

    tot_long = pd.DataFrame(columns = ["scenario"] + [f"{meas}_{bnd}" for meas in cost_measure for bnd in ["cent", "lb", "ub"]] + [f"{meas}_{bnd}" for meas in util_measure for bnd in ["cent", "lb", "ub"]])
    tot_long.scenario = const.zaf_scenarios
    tot_long = tot_long.set_index("scenario", drop=True)

    for scen in const.zaf_scenarios:
        for meas in cost_measure:
            #2025-2050
            tot_short.at[scen, f"{meas}_cent"] = np.sum(costs[scen][meas].iloc[34:61, 1])
            tot_short.at[scen, f"{meas}_lb"] = np.percentile(np.sum(costs[scen][meas].iloc[34:61, 2:]), 2.5)
            tot_short.at[scen, f"{meas}_ub"] = np.percentile(np.sum(costs[scen][meas].iloc[34:61, 2:]), 97.5)
            #2025-2100
            tot_long.at[scen, f"{meas}_cent"] = np.sum(costs[scen][meas].iloc[34:, 1])
            tot_long.at[scen, f"{meas}_lb"] = np.percentile(np.sum(costs[scen][meas].iloc[34:, 2:]), 2.5)
            tot_long.at[scen, f"{meas}_ub"] = np.percentile(np.sum(costs[scen][meas].iloc[34:, 2:]), 97.5)
        for util in util_measure:
            #2025-2050
            tot_short.at[scen, f"{util}_cent"] = np.sum(dalys[scen][util].iloc[34:61, 1])
            tot_short.at[scen, f"{util}_lb"] = np.percentile(np.sum(dalys[scen][util].iloc[34:61, 2:]), 2.5)
            tot_short.at[scen, f"{util}_ub"] = np.percentile(np.sum(dalys[scen][util].iloc[34:61, 2:]), 97.5)
            # 2025-2100
            tot_long.at[scen, f"{util}_cent"] = np.sum(dalys[scen][util].iloc[34:, 1])
            tot_long.at[scen, f"{util}_lb"] = np.percentile(np.sum(dalys[scen][util].iloc[34:, 2:]), 2.5)
            tot_long.at[scen, f"{util}_ub"] = np.percentile(np.sum(dalys[scen][util].iloc[34:, 2:]), 97.5)

    comp_short = pd.DataFrame(columns = ["scenario"] + [f"{meas}_{bnd}" for meas in cost_measure for bnd in ["cent", "lb", "ub"]] + [f"{meas}_{bnd}" for meas in util_measure for bnd in ["cent", "lb", "ub"]] + [f"{meas}_{bnd}" for meas in eff_measure for bnd in ["cent", "lb", "ub"]])
    comp_short.scenario = comparators
    comp_short = comp_short.set_index("scenario", drop = True)

    comp_long = pd.DataFrame(columns = ["scenario"] + [f"{meas}_{bnd}" for meas in cost_measure for bnd in ["cent", "lb", "ub"]] + [f"{meas}_{bnd}" for meas in util_measure for bnd in ["cent", "lb", "ub"]] + [f"{meas}_{bnd}" for meas in eff_measure for bnd in ["cent", "lb", "ub"]])
    comp_long.scenario = comparators
    comp_long = comp_long.set_index("scenario", drop=True)

    for comp in comparators:
        for meas in cost_measure:
            #2025-2050
            comp_short.at[comp, f"{meas}_cent"] = np.sum(costs[comp][meas].iloc[34:61, 1]) - np.sum(costs["baseline"][meas].iloc[34:61, 1])
            comp_short.at[comp, f"{meas}_lb"] = np.percentile(np.sum(costs[comp][meas].iloc[34:61, 2:]) - np.sum(costs["baseline"][meas].iloc[34:61, 2:]), 2.5)
            comp_short.at[comp, f"{meas}_ub"] = np.percentile(np.sum(costs[comp][meas].iloc[34:61, 2:]) - np.sum(costs["baseline"][meas].iloc[34:61, 2:]),97.5)
            #2025-2100
            comp_long.at[comp, f"{meas}_cent"] = np.sum(costs[comp][meas].iloc[34:, 1]) - np.sum(costs["baseline"][meas].iloc[34:, 1])
            comp_long.at[comp, f"{meas}_lb"] = np.percentile(np.sum(costs[comp][meas].iloc[34:, 2:]) - np.sum(costs["baseline"][meas].iloc[34:, 2:]), 2.5)
            comp_long.at[comp, f"{meas}_ub"] = np.percentile(np.sum(costs[comp][meas].iloc[34:, 2:]) - np.sum(costs["baseline"][meas].iloc[34:, 2:]),97)

        for util in util_measure:
            #2025-2050
            comp_short.at[comp, f"{util}_cent"] = np.sum(np.sum(dalys["baseline"][util].iloc[34:61, 1] - dalys[comp][util].iloc[34:61, 1]))
            comp_short.at[comp, f"{util}_lb"] = np.percentile(np.sum(dalys["baseline"][util].iloc[34:61, 2:]) - np.sum(dalys[comp][util].iloc[34:61, 2:]), 2.5)
            comp_short.at[comp, f"{util}_ub"] = np.percentile(np.sum(dalys["baseline"][util].iloc[34:61, 2:]) - np.sum(dalys[comp][util].iloc[34:61, 2:]), 97.5)
            #2025-2100
            comp_long.at[comp, f"{util}_cent"] = np.sum(np.sum(dalys["baseline"][util].iloc[34:, 1] - dalys[comp][util].iloc[34:, 1]))
            comp_long.at[comp, f"{util}_lb"] = np.percentile(np.sum(dalys["baseline"][util].iloc[34:, 2:]) - np.sum(dalys[comp][util].iloc[34:, 2:]), 2.5)
            comp_long.at[comp, f"{util}_ub"] = np.percentile(np.sum(dalys["baseline"][util].iloc[34:, 2:]) - np.sum(dalys[comp][util].iloc[34:, 2:]), 97.5)

        for eff in eff_measure:
            # 2025-2050
            comp_short.at[comp, f"{eff}_cent"] = dalys[comp][eff].iloc[0,1]
            comp_short.at[comp, f"{eff}_lb"] = np.percentile(dalys[comp][eff].iloc[0,2:], 2.5)
            comp_short.at[comp, f"{eff}_ub"] = np.percentile(dalys[comp][eff].iloc[0,2:], 97.5)
            #2025-2100
            comp_long.at[comp, f"{eff}_cent"] = dalys[comp][eff].iloc[1,1]
            comp_long.at[comp, f"{eff}_lb"] = np.percentile(dalys[comp][eff].iloc[1,2:], 2.5)
            comp_long.at[comp, f"{eff}_ub"] = np.percentile(dalys[comp][eff].iloc[1,2:], 97.5)

    # Export to a consolidated excel spreadsheet
    with pd.ExcelWriter(res_save_dir/'tab_figs'/'table2econ.xlsx', engine='xlsxwriter') as writer:
        tot_short.to_excel(writer, sheet_name="Sum Res_2050", index=True)
        tot_long.to_excel(writer, sheet_name="Sum Res_2100", index=True)
        comp_short.to_excel(writer, sheet_name="Dif Res_2050", index=True)
        comp_long.to_excel(writer, sheet_name="Dif Res_2100", index=True)


    # Output for cost-effectiveness frontier plot
    wtp_daly = const.wtp_daly

    # Generate dictionary to iterate results to
    cef_data_short = {} # analysis limited to 2050
    cef_data_long = {} # analysis from 2025-2100

    cef_data_short_add = {} # limited to 2050, add HBsAg- univ vax costs
    cef_data_long_add = {} # to 2100, add HBsAg- univ vax costs

    for wtp in wtp_daly:
        cef_data_short[f"{wtp}_thresh"] = pd.DataFrame(columns = ["run"] + const.zaf_scenarios)
        cef_data_short[f"{wtp}_thresh"].run = [f"run_{run}" for run in range(1, const.runs+1,1)]
        cef_data_short[f"{wtp}_thresh"] = cef_data_short[f"{wtp}_thresh"].set_index("run", drop=True)

        cef_data_short_add[f"{wtp}_thresh"] = pd.DataFrame(columns = ["run"] + const.zaf_scenarios)
        cef_data_short_add[f"{wtp}_thresh"].run = [f"run_{run}" for run in range(1, const.runs+1,1)]
        cef_data_short_add[f"{wtp}_thresh"] = cef_data_short_add[f"{wtp}_thresh"].set_index("run", drop=True)

        cef_data_long[f"{wtp}_thresh"] = pd.DataFrame(columns = ["run"] + const.zaf_scenarios)
        cef_data_long[f"{wtp}_thresh"].run = [f"run_{run}" for run in range(1, const.runs+1,1)]
        cef_data_long[f"{wtp}_thresh"] = cef_data_long[f"{wtp}_thresh"].set_index("run", drop=True)

        cef_data_long_add[f"{wtp}_thresh"] = pd.DataFrame(columns = ["run"] + const.zaf_scenarios)
        cef_data_long_add[f"{wtp}_thresh"].run = [f"run_{run}" for run in range(1, const.runs+1,1)]
        cef_data_long_add[f"{wtp}_thresh"] = cef_data_long_add[f"{wtp}_thresh"].set_index("run", drop=True)

    for wtp in wtp_daly:
        for run in range(1, const.runs+1, 1):
            for scen in const.zaf_scenarios:
                cef_data_short[f"{wtp}_thresh"].at[f"run_{run}", scen] = (np.sum(dalys[scen]["dalys"].iloc[34:61, run+1])*wtp)+(np.sum(costs[scen]["total"].iloc[34:61, run+1]))
                cef_data_short_add[f"{wtp}_thresh"].at[f"run_{run}", scen] = (np.sum(dalys[scen]["dalys"].iloc[34:61, run+1])*wtp)+(np.sum(costs[scen]["total_add"].iloc[34:61, run+1]))

                cef_data_long[f"{wtp}_thresh"].at[f"run_{run}", scen] = (np.sum(dalys[scen]["dalys"].iloc[34:, run+1])*wtp)+(np.sum(costs[scen]["total"].iloc[34:, run+1]))
                cef_data_long_add[f"{wtp}_thresh"].at[f"run_{run}", scen] = (np.sum(dalys[scen]["dalys"].iloc[34:, run+1])*wtp)+(np.sum(costs[scen]["total_add"].iloc[34:, run+1]))


    filename_short = sc.makefilepath(filename=f"ZAF_{const.runs}_cef_short_data.pkl", folder=res_save_dir/'result_pkls', makedirs=True)
    sc.save(filename=filename_short, obj=cef_data_short)

    filename_short = sc.makefilepath(filename=f"ZAF_{const.runs}_cef_short_add_data.pkl", folder=res_save_dir/'result_pkls', makedirs=True)
    sc.save(filename=filename_short, obj=cef_data_short_add)

    filename_long = sc.makefilepath(filename=f"ZAF_{const.runs}_cef_long_data.pkl", folder=res_save_dir/'result_pkls', makedirs=True)
    sc.save(filename=filename_long, obj=cef_data_long)

    filename_long = sc.makefilepath(filename=f"ZAF_{const.runs}_cef_long_data_add.pkl", folder=res_save_dir/'result_pkls', makedirs=True)
    sc.save(filename=filename_long, obj=cef_data_long_add)

    # Export HepB-BD cost components for each scenario (central only)
    vax_comps = ["hepbd", "hepbd_add", "mat_scr", "mat_av"]
    vax_costs = {}

    for vcomp in vax_comps:
        vax_costs[vcomp] = pd.DataFrame(columns = ["year"] + const.zaf_scenarios)
        vax_costs[vcomp].year = np.arange(1990.5, 2101.5, 1)

        for scen in const.zaf_scenarios:
            vax_costs[vcomp][scen] = costs[scen][vcomp].iloc[:, 1] # Column 1 is the central estimate

    excel_file_path = res_save_dir/'tab_figs'/'vcost_data.xlsx'
    writer = pd.ExcelWriter(excel_file_path, engine='xlsxwriter')
    for sheet_name, df in vax_costs.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.close()











