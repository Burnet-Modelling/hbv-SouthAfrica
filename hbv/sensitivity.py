import hbv.constants as const
import numpy as np
import pandas as pd
import atomica as at
import pathlib
from gitlab_utils import _get_gitlab_folder, _get_sharepoint_folder
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 14}) # Sets the default font size to 14


## OWSA looking at rate of MTCT (low vs high) on health outcomes and costs (generate plot)
def owsa_mtct_run(db_path, calib_path, calibration):

    # Model Requisites
    F = at.ProjectFramework(const.framework_path)
    db_path = pathlib.Path(db_path + "ZAF_databook.xlsx")
    cal_dir = pathlib.Path(calib_path + calibration)

    # Rates of MTCT
    owsa_ins = {"higher":[0.85, 0.15], "lower":[0.07, 0.001], "main": [0.383, 0.048]} #higher risk, lower risk [hvl, lvl]

    # Initialize dictionary
    owsa_runs =  {}

    for (key,val) in owsa_ins.items():
        owsa_runs[key] = {}
        for scen in const.zaf_scenarios:
            owsa_runs[key][scen]= {}
            owsa_runs[key][scen]["P"] = []
            owsa_runs[key][scen]["result"] = []

    # Run the model and store results for each assumption and scenario pair
    for (key, val) in owsa_ins.items():
        if key == "higher":
            hvl_trisk = val[0]
            lvl_trisk = val[1]
        elif key == "lower":
            hvl_trisk = val[0]
            lvl_trisk = val[1]
        elif key == "main":
            hvl_trisk = val[0]
            lvl_trisk = val[1]

        for scen in const.zaf_scenarios:
        # Set up each modelled scenario
            if scen == "baseline":
                D = at.ProjectData.from_spreadsheet(db_path, framework= F) # no changes required
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
        # Selective HepB-BD, Current Guidelines
            if scen == "pessimistic_selective_SA":
                D = at.ProjectData.from_spreadsheet(db_path, framework= F)
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
                # Set mav_trtall to zero
                D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
                D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
                # Targeted (selective) HepB-BD coverage to 40%
                D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
                D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            if scen == "optimistic_selective_SA":
                D = at.ProjectData.from_spreadsheet(db_path, framework= F)
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
                # Set mav_trtall to zero
                D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
                D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
                # Targeted (selective) HepB-BD coverage to 40%
                D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
                D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            # Selective HepB-BD, WHO Guidelines
            if scen == "pessimistic_selective_WHO":
                D = at.ProjectData.from_spreadsheet(db_path, framework=F)
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
                # Set mav_trtall to one (treat all HBsAg+ diagnosed mothers)
                D.tdve["mav_trtall"].ts["0-0M"].assumption = 1
                D.tdve["mav_trtall"].ts["0-0F"].assumption = 1
                # Targeted (selective) HepB-BD coverage to 40%
                D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
                D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            if scen == "optimistic_selective_WHO":
                D = at.ProjectData.from_spreadsheet(db_path, framework=F)
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
                # Set mav_trtall to one (treat all HBsAg+ diagnosed mothers)
                D.tdve["mav_trtall"].ts["0-0M"].assumption = 1
                D.tdve["mav_trtall"].ts["0-0F"].assumption = 1
                # Targeted (selective) HepB-BD coverage to 80%
                D.tdve["tgt_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
                D.tdve["tgt_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            # Universal HepB-BD
            if scen == "pessimistic_universal":
                D = at.ProjectData.from_spreadsheet(db_path, framework=F)
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
                # Set mav_trtall to zero
                D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
                D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
                # Set universal HepB-BD coverage to 40%
                D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
                D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.4, 0.4], units="Number")
            if scen == "optimistic_universal":
                D = at.ProjectData.from_spreadsheet(db_path, framework=F)
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
                # Set mav_trtall to zero
                D.tdve["mav_trtall"].ts["0-0M"].assumption = 0
                D.tdve["mav_trtall"].ts["0-0F"].assumption = 0
                # Set Universal HepB-BD coverage to 80%
                D.tdve["univ_bd_cov"].ts["0-0M"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
                D.tdve["univ_bd_cov"].ts["0-0F"] = at.TimeSeries([1990, 2024, 2030, 2050], [0, 0, 0.8, 0.8], units="Number")
            # Selective + Universal, Current Guidelines
            if scen == "pessimistic_combined_SA":
                D = at.ProjectData.from_spreadsheet(db_path, framework=F)
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
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
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
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
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
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
                D.tdve["hvl_trisk"].ts["0-0M"].assumption = hvl_trisk
                D.tdve["hvl_trisk"].ts["0-0F"].assumption = hvl_trisk
                D.tdve["lvl_trisk"].ts["0-0M"].assumption = lvl_trisk
                D.tdve["lvl_trisk"].ts["0-0F"].assumption = lvl_trisk
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

            owsa_runs[key][scen]["P"].append(P)
            owsa_runs[key][scen]["result"].append(res)

    return owsa_runs


def owsa_mtct_outcomes(owsa_runs, res_save_dir):

    # Given this is a constrained analysis, all plotting code etc. can be contained in here

    results = ["mtct", "chb_inc", "hcc_inc", "deaths"]
    outputs = ["b_chb", "all_trans", "hcc_inc", "deaths_chb"]
    pops = ["total", "total", "total", "total"]
    sens = ["higher", "lower", "main"]

    # Epidemiological Outcomes Data (MTCT infections, total infections, HCC incidence, deaths)

    # Generate results dictionary
    epi_res = {}
    for sen in sens:
        epi_res[sen] = {}
        for res in results:
            epi_res[sen][res]= pd.DataFrame(columns = ["year"] + const.zaf_scenarios)
            epi_res[sen][res].year = np.arange(1990.5, 2101.5, 1)

    # Fill results dictionary
    for sen in sens:
        for idx,res in enumerate(results):
            for scen in const.zaf_scenarios:
                epi_res[sen][res][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], outputs[idx], pops[idx], pop_aggregation="sum", t_bins=1).series[0].vals

    # Summarise difference to 2100 (as lead in is the same, can just be sum across entirity)
    for sen in sens:
        for res in results:
            epi_res[sen][res] =  epi_res[sen][res].sum(numeric_only=True).to_frame().T
            epi_res[sen][res] = epi_res[sen][res].drop("year", axis=1)

            baseline = epi_res[sen][res]["baseline"]
            for scen in const.zaf_scenarios:
                epi_res[sen][res][scen] = baseline -  epi_res[sen][res][scen]


    # NNV per MTCT averted (assume no HepB-BD for neg screen HBsAg, different sens analysis)

    # Generate Empty Dictionary
    nnv_res = {}
    for sen in sens:
        nnv_res[sen] = {}
        nnv_res[sen]= pd.DataFrame(columns = ["year"] + const.zaf_scenarios)
        nnv_res[sen].year = np.arange(1990.5, 2101.5, 1)

    # Fill dictionary
    for sen in sens:
        for scen in const.zaf_scenarios:
            nnv_res[sen][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "bd_vax_total", "total", pop_aggregation="sum", t_bins=1).series[0].vals

    # Calculate NNV
    for sen in sens:
        nnv_res[sen] = nnv_res[sen].sum(numeric_only=True).to_frame().T
        nnv_res[sen] = nnv_res[sen].drop("year", axis=1)

        for scen in const.zaf_scenarios:
            nnv_res[sen][scen] = np.nan_to_num(nnv_res[sen][scen]/epi_res[sen]["mtct"][scen], 0)

    # Costs and cost-effectiveness threshold (additional cost/DALYs averted)

    # Costs
    econ_comps = ["hepbd", "mav", "anc_scr", "hepb3", "diag", "treat", "noncir", "cirr", "dc", "hcc", "yld_dc", "yld_hcc",
                  "yll", "tot_cost", "dalys"]

    # Generate dictionary for storage of costs
    econ_dict = {}
    for sen in sens:
        econ_dict[sen] = {}
        for comp in econ_comps:
            econ_dict[sen][comp] = pd.DataFrame(columns = ["year"]+const.zaf_scenarios)
            econ_dict[sen][comp].year = np.arange(1990.5, 2101.5, 1)

    # Discounting Array (fixed at 3% per annum)
    disc_array = pd.DataFrame(columns=["year", "disc"])
    disc_array.year =np.arange(1990.5, 2101.5, 1)

    for i in range(len(disc_array)):
        if disc_array.iloc[i, 0] < 2024.5:
            disc_array.iloc[i, 1] = 0
        else:
            disc_array.iloc[i, 1] = (1 - 0.03) ** (disc_array.iloc[i, 0] - 2024.5)

    # Import Costs (only need central)
    flat_data = _get_sharepoint_folder() + "Data\\hbv_flat datasheet.xlsx"
    costs = pd.read_excel(flat_data, sheet_name="HBV Costs")

    # Extract costs
    hepbd_univ, hepbd_tgt, hepb3, diag_lab, diag_rapid, hbv_dna, tdf_treat, tdf_pap, chb_cost, cc_cost, dc_cost, hcc_cost, dc_daly, hcc_daly = costs["est"].loc[1:15]

    # HepB-BD
    for sen in sens:
        for scen in const.zaf_scenarios:
            # Costed as universal HepB-BD
            if (scen == "baseline"
            or scen == "pessimistic_universal"
            or scen == "optimistic_universal"
            or scen == "pessimistic_combined_SA"
            or scen == "optimistic_combined_SA"
            or scen == "pessimistic_combined_WHO"
            or scen == "optimistic_combined_WHO"):
                econ_dict[sen]["hepbd"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0],"bd_vax_total", "total", pop_aggregation="sum", t_bins=1).series[0].vals * \
                                                hepbd_univ * disc_array.iloc[:, 1]
            # Costed as selective HepB-BD
            elif (scen == "pessimistic_selective_SA"
            or scen == "optimistic_selective_SA"
            or scen == "pessimistic_selective_WHO"
            or scen == "optimistic_selective_WHO"):
                econ_dict[sen]["hepbd"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0],"bd_vax_total", "total", pop_aggregation="sum", t_bins=1).series[0].vals * \
                                                hepbd_tgt * disc_array.iloc[:, 1]
    # Maternal Screening
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["anc_scr"][scen] = (at.PlotData(owsa_runs[sen][scen]["result"][0], "n_anc_scr", pops="total", pop_aggregation = "sum", t_bins=1).series[0].vals * diag_rapid +\
                                         at.PlotData(owsa_runs[sen][scen]["result"][0], "n_hbv_dna", pops="total", pop_aggregation = "sum", t_bins=1).series[0].vals * hbv_dna) * \
                                         disc_array.iloc[:, 1]
    # Maternal antivirals
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["mav"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "n_mav_bd", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          tdf_pap * disc_array.iloc[:, 1]

    # HepB3
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["hepb3"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "hepb3_vax_total", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          hepb3 * disc_array.iloc[:, 1]

    # Diagnoses
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["diag"][scen]= at.PlotData(owsa_runs[sen][scen]["result"][0], "new_diag", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          diag_rapid * disc_array.iloc[:, 1]

    # Treatment
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["treat"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "treated", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          tdf_treat * disc_array.iloc[:, 1]

    # CHB
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["noncir"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "cost_noncir", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          chb_cost * disc_array.iloc[:, 1]

    # CC
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["cirr"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "cost_cir", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          cc_cost * disc_array.iloc[:, 1]

    # DC
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["dc"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "cost_decomp", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          dc_cost * disc_array.iloc[:, 1]

    # HCC
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["hcc"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "cost_canc", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          hcc_cost * disc_array.iloc[:, 1]

    # DC YLDs
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["yld_dc"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "decompensated", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          dc_daly * disc_array.iloc[:, 1]

    # HCC YLDs
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["yld_hcc"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "cancer", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          hcc_daly * disc_array.iloc[:, 1]

    # YLLs
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["yll"][scen] = at.PlotData(owsa_runs[sen][scen]["result"][0], "yll", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals *\
                                          disc_array.iloc[:, 1]

    # Total Costs (already discounted)
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["tot_cost"][scen] = econ_dict[sen]["hepbd"][scen] + econ_dict[sen]["anc_scr"][scen] + econ_dict[sen]["mav"][scen] + econ_dict[sen]["hepb3"][scen] + econ_dict[sen]["diag"][scen] + \
                                               econ_dict[sen]["treat"][scen] + econ_dict[sen]["noncir"][scen] + econ_dict[sen]["cirr"][scen] + econ_dict[sen]["dc"][scen] + econ_dict[sen]["hcc"][scen]

    # Total DALYs (already discounted)
    for sen in sens:
        for scen in const.zaf_scenarios:
            econ_dict[sen]["dalys"][scen] = econ_dict[sen]["yld_dc"][scen] + econ_dict[sen]["yld_hcc"][scen] + econ_dict[sen]["yll"][scen]

    # Costs and DALYs averted by scenario (vs baseline)
    summ_data = {}
    measures = ["tot_cost", "dalys"]

    # Generate empty dataframe
    for sen in sens:
        summ_data[sen] = pd.DataFrame(columns=["measure"]+const.zaf_scenarios)
        summ_data[sen].measure = measures
        summ_data[sen] = summ_data[sen].set_index("measure", drop = False)

    # Sum results
    for sen in sens:
        for scen in const.zaf_scenarios:
            for meas in measures:
                summ_data[sen].at[meas, scen] = np.sum(econ_dict[sen][meas][scen])

    # Calculate difference vs baseline
    for sen in sens:
        bl_cost = summ_data[sen].at["tot_cost", "baseline"]
        bl_daly = summ_data[sen].at["dalys", "baseline"]
        for scen in const.zaf_scenarios:
            summ_data[sen].at["tot_cost", scen] =  bl_cost - summ_data[sen].at["tot_cost", scen] # costs averted
            summ_data[sen].at["dalys", scen] = bl_daly - summ_data[sen].at["dalys", scen] # DALYs averted


    # Plot data
    plot_scens = ["optimistic_selective_WHO", "optimistic_universal",
                  "optimistic_combined_WHO"]
    scens_full = ["Selective PAP+\nSelective HepB-BD", "Universal HepB-BD", "Selective PAP+\nUniversal HepB-BD"]

    scens_leg = ["Selective PAP+\nSelective HepB-BD", "Universal HepB-BD", "Selective PAP+\nUniversal HepB-BD"]

    fig = plt.figure(figsize=(20,18))
    msize = 100

    mtct = fig.add_subplot(3,2,1)
    for i, pl_sc in enumerate(plot_scens):
        mtct.scatter(scens_full[i], epi_res["lower"]["mtct"][pl_sc], color="red", s=msize, marker="^", label = "Lower MTCT")
        mtct.scatter(scens_full[i], epi_res["main"]["mtct"][pl_sc], color="blue", s=msize, marker="^", label = "Modelled MTCT")
        mtct.scatter(scens_full[i], epi_res["higher"]["mtct"][pl_sc], color="green", s=msize, marker="^", label = "Higher MTCT")
        mtct.set_ylabel("MTCT attributable \nCHB infections averted")
        mtct.set_xticks([])


    chb = fig.add_subplot(3,2,2)
    for i, pl_sc in enumerate(plot_scens):
        chb.scatter(scens_full[i], epi_res["lower"]["chb_inc"][pl_sc], color="red", s=msize, marker="^", label = "Lower MTCT")
        chb.scatter(scens_full[i], epi_res["main"]["chb_inc"][pl_sc], color="blue", s=msize, marker="^", label = "Modelled MTCT")
        chb.scatter(scens_full[i], epi_res["higher"]["chb_inc"][pl_sc], color="green", s=msize, marker="^", label = "Higher MTCT")
        chb.set_ylabel("Total CHB infections averted")
        chb.set_xticks([])


    hcc = fig.add_subplot(3,2,3)
    for i, pl_sc in enumerate(plot_scens):
        hcc.scatter(scens_full[i], epi_res["lower"]["hcc_inc"][pl_sc], color="red", s=msize, marker="^", label = "Lower MTCT")
        hcc.scatter(scens_full[i], epi_res["main"]["hcc_inc"][pl_sc], color="blue", s=msize, marker="^", label = "Modelled MTCT")
        hcc.scatter(scens_full[i], epi_res["higher"]["hcc_inc"][pl_sc], color="green", s=msize, marker="^", label = "Higher MTCT")
        hcc.set_ylabel("Total HCC cases averted")
        hcc.set_xticks([])

    deaths = fig.add_subplot(3,2,4)
    for i, pl_sc in enumerate(plot_scens):
        deaths.scatter(scens_full[i], epi_res["lower"]["deaths"][pl_sc], color="red", s=msize, marker="^", label = "Lower MTCT")
        deaths.scatter(scens_full[i], epi_res["main"]["deaths"][pl_sc], color="blue", s=msize, marker="^", label = "Modelled MTCT")
        deaths.scatter(scens_full[i], epi_res["higher"]["deaths"][pl_sc], color="green", s=msize, marker="^", label = "Higher MTCT")
        deaths.set_ylabel("Total hepatitis B attributable \ndeaths averted")

    nnv = fig.add_subplot(3,2,5)
    for i, pl_sc in enumerate(plot_scens):
        nnv.scatter(scens_full[i], nnv_res["lower"][pl_sc], color="red", s=msize, marker="^", label = "Lower MTCT")
        nnv.scatter(scens_full[i], nnv_res["main"][pl_sc], color="blue", s=msize, marker="^", label = "Modelled MTCT")
        nnv.scatter(scens_full[i], nnv_res["higher"][pl_sc], color="green", s=msize, marker="^", label = "Higher MTCT")
        nnv.set_ylabel("Number of HepB-BD doses per \n MTCT infection averted")
        if i ==0:
            nnv.legend(loc="best")

    fig.tight_layout()
    fig.savefig(res_save_dir/'tab_figs'/"Supp OWSA plot.png", dpi=400)

    # Cost Effectiveness Threshold (seperate plot)
    scen_marker = ["o", "*", "P", "x", "D"]

    fig3 = plt.figure(figsize=(20,15))

    # Create arrays to plot ICER on OWSA chart
    half_gdp = np.zeros((100,2))
    half_gdp[:,0] = np.arange(0, 3e9,3e9/100)
    half_gdp[:,1] = half_gdp[:,0]/((0.5*6250)*18.33)

    full_gdp = np.zeros((100, 2))
    full_gdp[:, 0] = np.arange(0, 3e9, 3e9 / 100)
    full_gdp[:, 1] = half_gdp[:, 0] / ((1* 6250) * 18.33)

    two_gdp = np.zeros((100, 2))
    two_gdp[:, 0] = np.arange(0, 3e9, 3e9 / 100)
    two_gdp[:, 1] = half_gdp[:, 0] / ((2 * 6250) * 18.33)

    cet = fig3.add_subplot(1,1,1)
    for i, pl_sc in enumerate(plot_scens):
        cet.scatter(-summ_data["lower"].at["tot_cost", pl_sc]/1e9, summ_data["lower"].at["dalys", pl_sc], s=msize, marker = scen_marker[i], color="red", label=scens_leg[i])
        cet.scatter(-summ_data["main"].at["tot_cost", pl_sc]/1e9, summ_data["main"].at["dalys", pl_sc], s=msize, marker = scen_marker[i], color="blue")
        cet.scatter(-summ_data["higher"].at["tot_cost", pl_sc]/1e9, summ_data["higher"].at["dalys", pl_sc], s=msize, marker = scen_marker[i], color="green")
        cet.set_ylabel("DALYs averted")
        cet.set_xlabel("Additional Costs (ZAR, Billions)")
        cet.vlines(x=0, ymin = 0, ymax=7e5, color="black", linestyles="--")
        cet.set_ylim(0, 6.5e5)
        cet.set_xlim(-14, 2)
        cet.legend(loc="upper right")

    cet.plot(half_gdp[:, 0]/1e9, half_gdp[:, 1], color="orange", linestyle="--", alpha=0.4)
    cet.plot(full_gdp[:, 0]/1e9, full_gdp[:, 1], color="orange", linestyle="--", alpha=0.4)
    cet.plot(two_gdp[:, 0]/1e9, two_gdp[:, 1], color="orange", linestyle="--", alpha=0.4)

    fig3.tight_layout()
    fig3.savefig(res_save_dir/'tab_figs'/"Figure 3 OWSA.pdf", dpi=400)























#
#
# def data_extract(res):
#     """ Extract data from model runs when doing OWSA sensitivity analyses """
#     data = {}
#
#     # Transmission via MTCT and ECHT in each population group
#     pop_bins = list(const.pop_bins.keys())
#     data["trans_paths"] = pd.DataFrame(columns = ["year"] + [f"{pop}{sex}_{route}" for pop in pop_bins for sex in ["M", "F"] for route in ["mtct", "horiz", "total", "echt_rate", "prev"]])
#     data["trans_paths"].year = np.arange(1990.5, 2101.5, 1)
#
#     for bin in list(const.pop_bins.keys()):
#         # Male
#         data["trans_paths"][f"{bin}M_mtct"] = at.PlotData(res, "b_chb", f"{bin}M", t_bins=1).series[0].vals
#         data["trans_paths"][f"{bin}M_horiz"]  = at.PlotData(res, "horiz_trans", f"{bin}M", t_bins=1).series[0].vals
#         data["trans_paths"][f"{bin}M_total"] = at.PlotData(res, "all_trans", f"{bin}M", t_bins=1).series[0].vals
#         data["trans_paths"][f"{bin}M_echt_rate"] = at.PlotData(res, "echt", f"{bin}M", t_bins=1).series[0].vals
#         data["trans_paths"][f"{bin}M_prev"] = at.PlotData(res, "hbsag_prev", f"{bin}M", t_bins=1).series[0].vals
#         # Female
#         data["trans_paths"][f"{bin}F_mtct"] = at.PlotData(res, "b_chb", f"{bin}F", t_bins=1).series[0].vals
#         data["trans_paths"][f"{bin}F_horiz"] = at.PlotData(res, "horiz_trans", f"{bin}F", t_bins=1).series[0].vals
#         data["trans_paths"][f"{bin}F_total"] = at.PlotData(res, "all_trans", f"{bin}F", t_bins=1).series[0].vals
#         data["trans_paths"][f"{bin}F_echt_rate"] = at.PlotData(res, "echt", f"{bin}F", t_bins=1).series[0].vals
#         data["trans_paths"][f"{bin}F_prev"] = at.PlotData(res, "hbsag_prev", f"{bin}F", t_bins=1).series[0].vals
#
#
#
#     # Flow of births in MTCT cascade
#     data["mtct_flow"] = pd.DataFrame(columns = ["year"] + [f"0-0{sex}_{route}" for sex in ["M", "F"] for route in ["chb", "rec", "prot", "sus", "vax", "tot"]])
#     data["mtct_flow"].year = np.arange(1990.5, 2101.5, 1)
#
#     for sex in ["M", "F"]:
#         data["mtct_flow"][f"0-0{sex}_chb"] = at.PlotData(res, "b_chb", f"0-0{sex}", t_bins=1).series[0].vals
#         data["mtct_flow"][f"0-0{sex}_rec"] = at.PlotData(res, "b_rec", f"0-0{sex}", t_bins=1).series[0].vals
#         data["mtct_flow"][f"0-0{sex}_prot"]  = at.PlotData(res, "b_prot", f"0-0{sex}", t_bins=1).series[0].vals
#         data["mtct_flow"][f"0-0{sex}_sus"]  = at.PlotData(res, "b_sus", f"0-0{sex}", t_bins=1).series[0].vals
#         data["mtct_flow"][f"0-0{sex}_vax"]  = at.PlotData(res, "b_vax", f"0-0{sex}", t_bins=1).series[0].vals
#         data["mtct_flow"][f"0-0{sex}_tot"]  = at.PlotData(res, "b_rate", f"0-0{sex}", t_bins=1).series[0].vals
#
#     return data
#
#
#
#
# def mtct_rate_owsa(res_save_dir, calib_path, calibration, db_path):
#
#     # Increased (Asia-Pacific) MTCT risk by viral load
#     mtct_hvl = 0.85
#     mtct_lvl = 0.15
#
#     # Generate baseline values for comparison (modelled values for sub-saharan Africa)
#     F = at.ProjectFramework(const.framework_path)
#     P = at.Project(framework= F, databook=pathlib.Path(db_path+"ZAF_databook.xlsx"), sim_start=const.sim_start,
#                    sim_end = const.sim_end, sim_dt = const.sim_dt, do_run = False)
#     cal_parset = P.make_parset()
#     cal_parset = cal_parset.load_calibration(pathlib.Path(calib_path+calibration))
#     res_bl = P.run_sim(parset=cal_parset, result_name="Baseline Comparison")
#
#     bl_data = data_extract(res_bl)
#
#     # Recalibrate model to match cumulative infections prior to HepB-BD (1990.5-2024.5) in 0-0M/0-0F
#     D_sens = at.ProjectData.from_spreadsheet(pathlib.Path(db_path+"ZAF_databook.xlsx"), framework= F)
#
#     # Replace relevant parameters
#     D_sens.tdve["hvl_trisk"].ts["0-0M"].assumption = mtct_hvl
#     D_sens.tdve["hvl_trisk"].ts["0-0F"].assumption = mtct_hvl
#     D_sens.tdve["lvl_trisk"].ts["0-0M"].assumption = mtct_lvl
#     D_sens.tdve["lvl_trisk"].ts["0-0F"].assumption = mtct_lvl
#
#     P_sens = at.Project(framework= F, databook=D_sens, sim_start=const.sim_start,
#                    sim_end = const.sim_end, sim_dt = const.sim_dt, do_run = False)
#     cal_parset_sens = P_sens.make_parset()
#     cal_parset_sens = cal_parset_sens.load_calibration(pathlib.Path(calib_path+calibration))
#     res_sens = P_sens.run_sim(parset=cal_parset_sens, result_name="Unadjusted MTCT increase")
#
#     sens_data =  data_extract(res_sens)
#
#
#
#     usens_data = pd.DataFrame(columns=["year", "0-0M_mtct", "0-0M_echt", "0-0M_total",
#                                       "0-0F_mtct", "0-0F_echt", "0-0F_total"])
#
#     usens_data.year = at.PlotData(res_sens, "chb_mtct", "0-0M", t_bins=1).series[0].tvec
#
#     # Mother to child transmission
#     usens_data["0-0M_mtct"] = at.PlotData(res_sens, "b_chb", "0-0M", t_bins=1).series[0].vals
#     usens_data["0-0F_mtct"] = at.PlotData(res_sens, "b_chb", "0-0F", t_bins=1).series[0].vals
#
#     # Early Childhood Horizontal transmission
#     usens_data["0-0M_echt"] = at.PlotData(res_sens, "horiz_trans", "0-0M", t_bins=1).series[0].vals
#     usens_data["0-0F_echt"] = at.PlotData(res_sens, "horiz_trans", "0-0F", t_bins=1).series[0].vals
#
#     # Total new CHB infections in 0-0y
#     usens_data["0-0M_total"] = usens_data["0-0M_mtct"] + usens_data["0-0M_echt"]
#     usens_data["0-0F_total"] = usens_data["0-0F_mtct"] + usens_data["0-0F_echt"]
#
#     # Test approach one: ratio of horizontal from run to run
#     male_ratio = (comp_data["0-0M_echt"] /usens_data["0-0M_echt"])*(usens_data["0-0M_mtct"]/comp_data["0-0M_total"])
#     female_ratio = (comp_data["0-0F_echt"] /usens_data["0-0F_echt"])*(usens_data["0-0F_mtct"]/comp_data["0-0F_total"])
#
#
#     # Re-run model calibration
#     D_sens.tdve["foi_cal"].ts["0-0M"] = at.TimeSeries(list(np.arange(1990, 2101,1)), list(male_ratio), units = "N.A.")
#     D_sens.tdve["foi_cal"].ts["0-0F"] = at.TimeSeries(list(np.arange(1990, 2101,1)), list(female_ratio), units = "N.A.")
#
#     P_sens = at.Project(framework=F, databook=D_sens, sim_start=const.sim_start,
#                         sim_end=const.sim_end, sim_dt=const.sim_dt, do_run=False)
#     res_sens = P_sens.run_sim(parset=cal_parset, result_name="Unadjusted MTCT increase")
#
#     usens_data = pd.DataFrame(columns=["year", "0-0M_mtct", "0-0M_echt", "0-0M_total",
#                                        "0-0F_mtct", "0-0F_echt", "0-0F_total"])
#
#     usens_data.year = at.PlotData(res_sens, "chb_mtct", "0-0M", t_bins=1).series[0].tvec
#
#     # Mother to child transmission
#     usens_data["0-0M_mtct"] = at.PlotData(res_sens, "b_chb", "0-0M", t_bins=1).series[0].vals
#     usens_data["0-0F_mtct"] = at.PlotData(res_sens, "b_chb", "0-0F", t_bins=1).series[0].vals
#
#     # Early Childhood Horizontal transmission
#     usens_data["0-0M_echt"] = at.PlotData(res_sens, "horiz_trans", "0-0M", t_bins=1).series[0].vals
#     usens_data["0-0F_echt"] = at.PlotData(res_sens, "horiz_trans", "0-0F", t_bins=1).series[0].vals
#
#     # Total new CHB infections in 0-0y
#     usens_data["0-0M_total"] = usens_data["0-0M_mtct"] + usens_data["0-0M_echt"]
#     usens_data["0-0F_total"] = usens_data["0-0F_mtct"] + usens_data["0-0F_echt"]
#
#
#
#
#
#
