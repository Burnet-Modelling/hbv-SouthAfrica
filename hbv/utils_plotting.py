from hbv import constants as const
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 16}) # Sets the default font size to 14
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import ticker
import seaborn as sns
import atomica as at
import pandas as pd
import numpy as np
from gitlab_utils import _get_sharepoint_folder, _get_gitlab_folder
import pathlib
import sciris as sc

def calibration_plots(res_save_dir, calib_path, calibration, db_path):

    F = at.ProjectFramework(const.framework_path)
    calib_dir = pathlib.Path(calib_path+calibration)
    db_dir = pathlib.Path(db_path+"ZAF_databook.xlsx")

    P = at.Project(framework= F, databook = db_dir, sim_start = const.sim_start, sim_end= const.sim_end, sim_dt = const.sim_dt)
    cal_parset = P.make_parset()
    cal_parset = cal_parset.load_calibration(calib_dir)

    res = P.run_sim(parset=cal_parset, result_name="YAML Calibrations")
    pop_bins = [f"{item}{suffix}" for item in list(const.pop_bins.keys()) for suffix in ["M", "F"]]

    save_path = res_save_dir/'tab_figs'
    pdf_path = save_path/'ZAF_calibration plots.pdf'

    with PdfPages(pdf_path, keep_empty=False) as pdf:
        # Population Size
        pop_cal = plt.figure(figsize=(20,15))
        for i, pop in enumerate(pop_bins):
            # Extract Data
            years = at.PlotData(res, outputs="alive", pops = [pop], t_bins=1).series[0].tvec
            model = at.PlotData(res, outputs="alive", pops = [pop], t_bins=1).series[0].vals
            data = np.array(P.parsets[0].get_par("alive").ts[pop].vals)

            ax = pop_cal.add_subplot(int(len(pop_bins)/4), 4, i+1)
            ax.plot(years, model/1e6, color = "orange", alpha = 0.7, label = "Model Output")
            ax.scatter(years, data/1e6, color = "darkgrey", alpha=0.4, label = "Calibration Data")
            ax.set_title(pop)
            ax.set_ylim(bottom = 0)

            if i == 0:
                ax.legend(loc="best")

        pop_cal.supylabel("Population (millions)")
        pop_cal.tight_layout()
        plt.savefig(save_path/"ZAF_pop cal.png", dpi=400)
        pdf.savefig(pop_cal)
        plt.close()

        # HBsAg prevalence
        prev_cal = plt.figure(figsize=(20,20))
        for i, pop in enumerate(pop_bins):
            model = at.PlotData(res, outputs="hbsag_prev", pops=[pop], t_bins=1).series[0].vals
            data = np.array(P.parsets[0].get_par("hbsag_prev").ts[pop].vals)
            data_years = P.parsets[0].get_par("hbsag_prev").ts[pop].t

            ax = prev_cal.add_subplot(int(len(pop_bins) / 4), 4, i + 1)
            ax.plot(years, model*1e2, color="orange", alpha=0.7, label="Model Output")
            ax.scatter(data_years, data*1e2, color="darkgrey", alpha=0.4, label="Calibration Data")
            ax.set_title(pop)
            ax.set_ylim(bottom=0)

            if i == 0:
                ax.legend(loc="best")

        prev_cal.supylabel("HBsAg seroprevalence (%)")
        prev_cal.tight_layout()
        plt.savefig(save_path/"ZAF_hbsag cal.png", dpi=400)
        pdf.savefig(prev_cal)
        plt.close()

        #HBeAg prevalence
        eag_cal = plt.figure(figsize=(20, 20))
        for i, pop in enumerate(pop_bins):
            model = at.PlotData(res, outputs="hbeag_prev", pops=[pop], t_bins=1).series[0].vals
            data = np.array(P.parsets[0].get_par("hbeag_prev").ts[pop].vals)
            data_years = P.parsets[0].get_par("hbeag_prev").ts[pop].t

            ax = eag_cal.add_subplot(int(len(pop_bins) / 4), 4, i + 1)
            ax.plot(years, model*1e2, color="orange", alpha=0.7, label="Model Output")
            ax.scatter(data_years, data*1e2, color="darkgrey", alpha=0.4, label="Calibration Data")
            ax.set_title(pop)
            ax.set_ylim(bottom=0)

            if i == 0:
                ax.legend(loc="best")

        eag_cal.supylabel("HBeAg seroprevalence among \nHBsAg positive (%)")
        eag_cal.tight_layout()
        plt.savefig(save_path/"ZAF_hbeag cal.png", dpi=400)
        pdf.savefig(eag_cal)
        plt.close()

        # HCC incidence
        hcc_cal = plt.figure(figsize=(20, 20))
        for i, pop in enumerate(pop_bins):
            model = at.PlotData(res, outputs="hcc_inc", pops=[pop], t_bins=1).series[0].vals
            data = P.parsets[0].get_par("hcc_inc").ts[pop].vals
            data_years = P.parsets[0].get_par("hcc_inc").ts[pop].t

            ax = hcc_cal.add_subplot(int(len(pop_bins) / 4), 4, i + 1)
            ax.plot(years, model, color="orange", alpha=0.7, label = "Model Output")
            ax.scatter(data_years, data, color="darkgrey", alpha=0.4, label=" Calibration Data")
            ax.set_title(pop)
            ax.set_ylim(bottom=0)

            if i == 0:
                ax.legend(loc="best")

        hcc_cal.supylabel("Annual incident HBV attributable \nHCC cases")
        hcc_cal.tight_layout()
        plt.savefig(save_path/"ZAF_hcc cal.png", dpi=400)
        pdf.savefig(hcc_cal)
        plt.close()

    # Treatment Coverage
        trt_cal = plt.figure(figsize=(20, 20))
        for i, pop in enumerate(pop_bins):
            model = at.PlotData(res, outputs="treated", pops=[pop], t_bins=1).series[0].vals
            data = P.parsets[0].get_par("treated").ts[pop].vals
            data_years = P.parsets[0].get_par("treated").ts[pop].t

            ax = trt_cal.add_subplot(int(len(pop_bins) / 4), 4, i + 1)
            ax.plot(years, model, color="orange", alpha=0.7)
            ax.scatter(data_years, data, color="darkgrey", alpha=0.4)
            ax.set_title(pop)
            ax.set_ylim(bottom=0)
        trt_cal.tight_layout()
        plt.savefig(save_path/"ZAF_treat cal.png")
        pdf.savefig(trt_cal)
        plt.close()

# Transmission Route Plot (Fig 2)
def transmission_plots(res_save_dir):

    # Figure 2 (total CHB by transmission route in baseline scenario to 2050)
    fig2_data = pd.read_excel(res_save_dir/'tab_figs'/'fig2data.xlsx', sheet_name="baseline")
    fig2 = plt.figure(figsize=(20,15))
    trans = fig2.add_subplot(1,1,1)
    trans.bar(fig2_data.year, fig2_data.mtct, color = "blue", alpha=0.6, label = "MTCT")
    trans.bar(fig2_data.year, fig2_data.echt, bottom = fig2_data.mtct, color="green", alpha = 0.6, label = "Early Childhood")
    trans.bar(fig2_data.year, fig2_data.oth, bottom = fig2_data.mtct + fig2_data.echt, color="orange", alpha = 0.6, label = "Other")
    trans.legend (loc="best")
    trans.set_xlim(1990, 2100)
    trans.set_ylabel("Total incident CHB infections")
    fig2.tight_layout()
    fig2.savefig(res_save_dir/'tab_figs'/'Figure 2.pdf', dpi=400)
    fig2.savefig(res_save_dir/'tab_figs'/'Figure 2.png', dpi=400)
    plt.close()

    # Sense Check - transmission by route by scenario comparisons
    ifx_routes = ["mtct", "echt", "oth", "tot_sum"]
    titles = ["Mother to Child", "Early Childhood", "Other Horizontal", "Total"]

    fig_ch = plt.figure(figsize=(20,20))
    for i,ifx in enumerate(ifx_routes):
        ax = fig_ch.add_subplot(2,2, i+1)
        for scen in const.zaf_scenarios:
            plt_data = pd.read_excel(res_save_dir/'tab_figs'/'fig2data.xlsx', sheet_name = scen)
            ax.plot(plt_data.year, plt_data[ifx], label = scen)
            ax.set_title(titles[i])
            if i ==3:
                ax.legend(loc="best")
    fig_ch.tight_layout()
    fig_ch.savefig(res_save_dir/'tab_figs'/'Infx check.pdf', dpi=400)
    plt.close()

# NNV Plot (Fig 3)
def nnv_plot (res_save_dir):

    nnv_data = pd.read_csv(res_save_dir/'tab_figs'/'fig3data.csv')
    nnv_data.scenario =  const.zaf_scenarios_full
    nnv_data = nnv_data[nnv_data.scenario != "Baseline \n (No HepB-BD)"]

    scen_rename = ["Selective HepB-BD \n (Current Guidelines)",
                   "Selective HepB-BD \n (Provisional WHO recommendations)",
                   "Universal HepB-BD",
                   "Combination HepB-BD \n (Current Guidelines)",
                   "Combination HepB-BD \n (Provisional WHO recommendations)"
                   ]
    nnv_data = nnv_data[~nnv_data.scenario.str.contains("Pessimistic")]
    nnv_data.scenario = scen_rename



    metrics = ["mtct_infx", "tot_inc", "hcc_incidence", "hbv_deaths"]
    full_name = ["CHB from MTCT infection", "Incident CHB infection", "Incident HCC case", "Hepatitis B associated death"]


    fig3 = plt.figure(figsize = (20,15))
    for i, metric in enumerate(metrics):
        error = [list(nnv_data[f"{metric}_cent"]-nnv_data[f"{metric}_lb"]), list(nnv_data[f"{metric}_ub"]-nnv_data[f"{metric}_cent"])]
        panel = fig3.add_subplot(2,2, i+1)
        panel.barh(nnv_data.scenario, nnv_data[f"{metric}_cent"], xerr=error, color="orange", alpha=0.6, ecolor = "grey", capsize=4)
        panel.set_title(full_name[i])

        if i == 2 or i == 3:
            panel.set_xlabel("Hepatitis B birth dose vaccines administered \n (per outcome averted)")

    fig3.tight_layout()
    fig3.savefig(res_save_dir/'tab_figs'/"Figure 3.pdf", dpi=400)
    fig3.savefig(res_save_dir/'tab_figs'/"Figure 3.png", dpi=400)

    plt.close()

# Cost Effectiveness Frontier (inbuilt probability function)
def cost_eff_frontier(res_save_dir):

    main_scens = ["selective_WHO", "universal", "combined_WHO"]
    sens_scens = ["selective_SA", "universal", "combined_SA"]

    half_gdp = (0.5*6250)*18.33
    full_gdp = (1*6250)*18.33
    two_gdp = (2*6250)*18.33
    x_upper_bound = (2.5*6250)*18.33

    # Figure 1: CEF of selective_WHO, universal, combined_WHO to 2100, no HepB-BD in HBsAg screened
    cef_data_f1 = sc.load(res_save_dir/'result_pkls'/f'ZAF_{const.runs}_cef_long_data.pkl')

    # Empty Dataframe Generate
    f1_data = pd.DataFrame(columns = ["wtp_thresh", "baseline"] + [f"optimistic_{scen}" for scen in main_scens])
    f1_data.wtp_thresh = const.wtp_daly
    f1_data = f1_data.set_index("wtp_thresh", drop = False)

    # Extract Data Required for CEF
    f1_keep_cols = ["baseline"] + [f"optimistic_{scen}" for scen in main_scens]
    f1_model_data = {k: v[f1_keep_cols] for k, v in cef_data_f1.items()}

    # Return proportion of runs with highest net benefit
    for wtp in const.wtp_daly:
        mins_f1 = f1_model_data[f"{wtp}_thresh"].eq(f1_model_data[f"{wtp}_thresh"].min(axis=1),axis=0).astype(int)
        mins_f1 = mins_f1.sum(numeric_only = True).to_frame().T/const.runs

        for scen in f1_keep_cols:
            f1_data.at[wtp, scen] = mins_f1[scen].values

    # Plot Data
    full_scen = ["Baseline\n(No HepB-BD)", "Selective PAP+\nSelective HepB-BD",
                 "Universal HepB-BD",
                 "Selective PAP+\nUniversal HepB-BD"]
    fig1 = plt.figure(figsize=(20,15))
    ax = fig1.add_subplot(1,1,1)
    for idx, scen in enumerate(["baseline"] + [f"optimistic_{scen}" for scen in main_scens]):
        ax.plot(f1_data.wtp_thresh, f1_data[scen]*1e2, label=full_scen[idx], lw=2)
        ax.legend(loc="center right")
        ax.set_ylabel (f"Proportion of {const.runs} simulations where \nscenario yielded highest net-benefit (%)")
        ax.set_xlabel("Willingness-to-pay per DALY averted (2024 ZAR)")
    ax.vlines(x=half_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    ax.vlines(x=full_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    ax.vlines(x=two_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    ax.text(x=half_gdp + 250, y=90, s="0.5x GDP \nper capita")
    ax.text(x=full_gdp + 250, y=90, s="1x GDP \nper capita")
    ax.text(x=two_gdp + 250, y=90, s="2x GDP \nper capita")
    ax.set_xlim(0, x_upper_bound)
    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    ax.set_ylim(0, 101)
    fig1.tight_layout()
    fig1.savefig(res_save_dir / 'tab_figs' / "Figure 1 CEF.pdf", dpi=400)
    plt.close()

    # Figure 2: CEF of selective_WHO, universal, combined_WHO to 2100, yes HepB-BD in HBsAg screened (sensitivity)
    cef_data_f2 = sc.load(res_save_dir/'result_pkls'/f'ZAF_{const.runs}_cef_long_data_add.pkl')

    f2_data = pd.DataFrame(columns=["wtp_thresh", "baseline"] + [f"optimistic_{scen}" for scen in main_scens])
    f2_data.wtp_thresh = const.wtp_daly
    f2_data = f2_data.set_index("wtp_thresh", drop=False)

    f2_keep_cols = ["baseline"] + [f"optimistic_{scen}" for scen in main_scens]
    f2_model_data = {k: v[f2_keep_cols] for k, v in cef_data_f2.items()}

    for wtp in const.wtp_daly:
        mins_f2 = f2_model_data[f"{wtp}_thresh"].eq(f2_model_data[f"{wtp}_thresh"].min(axis=1), axis=0).astype(int)
        mins_f2 = mins_f2.sum(numeric_only=True).to_frame().T / const.runs

        for scen in f2_keep_cols:
            f2_data.at[wtp, scen] = mins_f2[scen].values

    fig2 = plt.figure(figsize=(20,15))
    ax = fig2.add_subplot(1,1,1)
    for idx, scen in enumerate(["baseline"] + [f"optimistic_{scen}" for scen in main_scens]):
        ax.plot(f2_data.wtp_thresh, f2_data[scen]*1e2, label=full_scen[idx], lw=2)
        ax.legend(loc="center right")
        ax.set_ylabel (f"Proportion of {const.runs} simulations where \nscenario yielded highest net-benefit (%)")
        ax.set_xlabel("Willingness-to-pay per DALY averted (2024 ZAR)")
    ax.vlines(x=half_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    ax.vlines(x=full_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    ax.vlines(x=two_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    ax.text(x=half_gdp + 250, y=90, s="0.5x GDP \nper capita")
    ax.text(x=full_gdp + 250, y=90, s="1x GDP \nper capita")
    ax.text(x=two_gdp + 250, y=90, s="2x GDP \nper capita")
    ax.set_xlim(0, x_upper_bound)
    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    ax.set_ylim(0, 101)
    fig2.tight_layout()
    fig2.savefig(res_save_dir / 'tab_figs' / "Figure 2 CEF_sens.pdf", dpi=400)
    plt.close()

    # Supplement: CEF of selective_SA, universal, combined_SA to 2100, yes HepB-BD in HBsAg screened (sensitivity; uses cef_data_f1, cef_data_f2)

    s1_data = pd.DataFrame(columns=["wtp_thresh", "baseline"] + [f"optimistic_{scen}" for scen in sens_scens]) # top plot
    s1_data.wtp_thresh = const.wtp_daly
    s1_data = s1_data.set_index("wtp_thresh", drop=False)

    s2_data = pd.DataFrame(columns=["wtp_thresh", "baseline"] + [f"optimistic_{scen}" for scen in sens_scens])# bottom plot
    s2_data.wtp_thresh = const.wtp_daly
    s2_data = s2_data.set_index("wtp_thresh", drop=False)

    sup_keep_cols = ["baseline"] + [f"optimistic_{scen}" for scen in sens_scens]
    s1_model_data = {k: v[sup_keep_cols] for k, v in cef_data_f1.items()}
    s2_model_data = {k: v[sup_keep_cols] for k, v in cef_data_f2.items()}

    for wtp in const.wtp_daly:
        mins_s1 = s1_model_data[f"{wtp}_thresh"].eq(s1_model_data[f"{wtp}_thresh"].min(axis=1), axis=0).astype(int)
        mins_s1 = mins_s1.sum(numeric_only=True).to_frame().T / const.runs

        mins_s2 = s2_model_data[f"{wtp}_thresh"].eq(s2_model_data[f"{wtp}_thresh"].min(axis=1), axis=0).astype(int)
        mins_s2 = mins_s2.sum(numeric_only=True).to_frame().T / const.runs

        for scen in sup_keep_cols:
            s1_data.at[wtp, scen] = mins_s1[scen].values
            s2_data.at[wtp, scen] = mins_s2[scen].values

    sup_cef = plt.figure(figsize=(20,15))

    top = sup_cef.add_subplot(2,1,1) # HBsAg negative not costed
    for idx,scen in enumerate(["baseline"] + [f"optimistic_{scen}" for scen in sens_scens]):
        top.plot(s1_data.wtp_thresh, s1_data[scen] * 1e2, label=full_scen[idx], lw=2)
        #top.legend(loc="center right")
        #top.set_ylabel(f"Proportion of {const.runs} simulations where \nscenario yielded highest net-benefit (%)")
    top.vlines(x=half_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    top.vlines(x=full_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    top.vlines(x=two_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    top.text(x=half_gdp + 250, y=80, s="0.5x GDP \nper capita")
    top.text(x=full_gdp + 250, y=80, s="1x GDP \nper capita")
    top.text(x=two_gdp + 250, y=80, s="2x GDP \nper capita")
    top.set_xlim(0, x_upper_bound)
    top.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    top.set_ylim(0, 101)
    top.set_title("ANC screened HBsAg negative not vaccinated")

    bottom = sup_cef.add_subplot(2,1,2) # HBsAg negative not costed
    for idx,scen in enumerate(["baseline"] + [f"optimistic_{scen}" for scen in sens_scens]):
        bottom.plot(s2_data.wtp_thresh, s2_data[scen] * 1e2, label=full_scen[idx], lw=2)
        bottom.legend(loc="center right")
        #bottom.set_ylabel(f"Proportion of {const.runs} simulations where \nscenario yielded highest net-benefit (%)")
    bottom.vlines(x=half_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    bottom.vlines(x=full_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    bottom.vlines(x=two_gdp, ymin=0, ymax=105, color="black", linestyles="--", alpha=0.4)
    bottom.set_xlim(0, x_upper_bound)
    bottom.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    bottom.set_ylim(0, 101)
    bottom.set_title("ANC screened HBsAg negative HepB-BD eligible")
    bottom.set_xlabel("Willingness-to-pay per DALY averted (2024 ZAR)")

    sup_cef.supylabel(f"Proportion of {const.runs} simulations where \nscenario yielded highest net-benefit (%)")
    sup_cef.tight_layout()
    sup_cef.savefig(res_save_dir / 'tab_figs' / "Supp CEF_HBV DNA screen.png", dpi=400)
    plt.close()


def valid_plots(res_save_dir, calib_path, calibration, db_path):
    """
    Alternate calibration and mortality validation plots
    """

    # Total Mortality Estimate
    who_ghr_2022 =  25564

    # Modelled estimate in 2022
    data = sc.load(res_save_dir/'result_pkls'/f"ZAF_{const.runs}_sims.pkl")
    deaths = data["baseline"]["hbv_deaths"]
    lower = np.percentile(deaths.iloc[:, 2:], 2.5, axis=1).astype(float)
    upper= np.percentile(deaths.iloc[:, 2:], 97.5, axis=1).astype(float)

    # Distribution of mortality within the model (2022)
    F = at.ProjectFramework(const.framework_path)
    calib_dir = pathlib.Path(calib_path + calibration)
    db_dir = pathlib.Path(db_path + "ZAF_databook.xlsx")

    P = at.Project(framework=F, databook=db_dir, sim_start=const.sim_start, sim_end=const.sim_end, sim_dt=const.sim_dt)
    cal_parset = P.make_parset()
    cal_parset = cal_parset.load_calibration(calib_dir)

    res = P.run_sim(parset=cal_parset, result_name="YAML Calibrations")
    pop_bins = [f"{item}{suffix}" for item in list(const.pop_bins.keys()) for suffix in ["M", "F"]]

    death_dist = pd.DataFrame(columns=["year"]+pop_bins)
    death_dist.year = np.arange(1990.5, 2101.5, 1)

    for bin in const.pop_bins:
        death_dist[f"{bin}M"] = at.PlotData(res, outputs="deaths_chb", pops=f"{bin}M", t_bins=1).series[0].vals
        death_dist[f"{bin}F"] = at.PlotData(res, outputs="deaths_chb", pops=f"{bin}F", t_bins=1).series[0].vals

    # Extract Validation Data From Model
    total_deaths = at.PlotData(res, outputs="deaths_chb", pops="total", pop_aggregation="sum", t_bins=1).series[0].vals[32]
    pop_valids = ["0-4", "5-14", "15-29", "30-49", "50-59", "60-69", "70+"]
    male_valids = pd.DataFrame(columns = ["year"]+[f"{bin}M" for bin in pop_valids])
    male_valids.year = [2022]
    female_valids = pd.DataFrame(columns = ["year"]+[f"{bin}F" for bin in pop_valids])
    female_valids.year = [2022]

    #Male
    male_valids["0-4M"] = (death_dist["0-0M"].iloc[32]+death_dist["1-4M"].iloc[32])/total_deaths
    male_valids["5-14M"] = (death_dist["5-9M"].iloc[32]+0.5*death_dist["10-19M"].iloc[32])/total_deaths
    male_valids["15-29M"] = (0.5*death_dist["10-19M"].iloc[32] + death_dist["20-29M"].iloc[32])/total_deaths
    male_valids["30-49M"] = (death_dist["30-39M"].iloc[32] + death_dist["40-49M"].iloc[32])/total_deaths
    male_valids["50-59M"] = death_dist["50-59M"].iloc[32]/total_deaths
    male_valids["60-69M"] = death_dist["60-69M"].iloc[32]/total_deaths
    male_valids["70+M"] = (death_dist["70-79M"].iloc[32] + death_dist["80-89M"].iloc[32] + death_dist["90+M"].iloc[32])/total_deaths

    #Female
    female_valids["0-4F"] = (death_dist["0-0F"].iloc[32] + death_dist["1-4F"].iloc[32])/total_deaths
    female_valids["5-14F"] = (death_dist["5-9F"].iloc[32] + 0.5 * death_dist["10-19F"].iloc[32])/total_deaths
    female_valids["15-29F"] = (0.5 * death_dist["10-19F"].iloc[32] + death_dist["20-29F"].iloc[32])/total_deaths
    female_valids["30-49F"] = (death_dist["30-39F"].iloc[32] + death_dist["40-49F"].iloc[32])/total_deaths
    female_valids["50-59F"] = death_dist["50-59F"].iloc[32]/total_deaths
    female_valids["60-69F"] = death_dist["60-69F"].iloc[32]/total_deaths
    female_valids["70+F"] = (death_dist["70-79F"].iloc[32] + death_dist["80-89F"].iloc[32] + death_dist["90+F"].iloc[32])/total_deaths

    # Import Mortality Validation Data
    data_valid = pd.read_excel( _get_sharepoint_folder() + "Data\\hbv_flat datasheet.xlsx", sheet_name = "Validation_deaths")


    # Plot model projection vs available data
    fig = plt.figure(figsize=(20, 15))
    gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])

    mort_tot = fig.add_subplot(gs[0, :])
    mort_tot.plot(deaths.year, deaths.central, color="blue", label="Modelled Estimate (Calibration)")
    mort_tot.fill_between(deaths.year, lower, upper, alpha=0.2, color="blue")
    mort_tot.scatter([2022], who_ghr_2022, color="black", label = "WHO Global Hepatitis Report Estimate", marker="*", s=100)
    mort_tot.legend(loc="best")
    mort_tot.set_title("Total Hepatitis B Attributable Mortality")
    mort_tot.set_xlabel("Year")
    mort_tot.set_ylabel("Annual Deaths")
    mort_tot.set_xlim(2010, 2030)

    # Plot data (bars for model, scatter points for data sources)
    ax2 = fig.add_subplot(gs[1, 0])
    male_plot = male_valids.iloc[:, 1:].sum().plot(ax = ax2, kind="bar", color="blue", alpha=0.3)
    male_plot.scatter(range(len(data_valid.iloc[0, 2:9])), data_valid.iloc[0, 2:9], color="green", label="Global Burden of Disease")
    male_plot.scatter(range(len(data_valid.iloc[1, 2:9])), data_valid.iloc[1, 2:9], color="orange", label = "WHO Global Health Estimates")
    male_plot.set_ylabel("Proportion of hepatitis B \nattributable deaths")
    male_plot.legend(loc="best")

    ax3 = fig.add_subplot(gs[1, 1])
    female_plot = female_valids.iloc[:, 1:].sum().plot(ax = ax3, kind="bar", color="blue", alpha = 0.3)
    female_plot.scatter(range(len(data_valid.iloc[0, 9:16])), data_valid.iloc[0, 9:16], color="green", label="Global Burden of Disease")
    female_plot.scatter(range(len(data_valid.iloc[1, 9:16])), data_valid.iloc[1, 9:16], color="orange", label = "WHO Global Health Estimates")
    #female_plot.legend(loc="best")

    fig.tight_layout()
    fig.savefig(res_save_dir/'tab_figs'/"Validation Plot.png", dpi=400)







def testing_central_plots(res_store, test_save_dir, savename):

    pop_bins = [f"{item}{suffix}" for item in list(const.pop_bins.keys()) for suffix in ["M", "F"]]

    # Calibration plot (baseline scenario)
    pdf_path_cal = test_save_dir / f'{savename}_ZAF_calibration plots.pdf'
    with PdfPages(pdf_path_cal, keep_empty=False) as pdf:
        # Population
        pop_cal = plt.figure(figsize=(20, 20))
        for i, pop in enumerate(pop_bins):
            # Extract Data
            years = at.PlotData(res_store["baseline"]["res"][0], outputs="alive", pops=[pop], t_bins=1).series[0].tvec
            model = at.PlotData(res_store["baseline"]["res"][0], outputs="alive", pops=[pop], t_bins=1).series[0].vals
            data = res_store["baseline"]["P"][0].parsets[0].get_par("alive").ts[pop].vals

            ax = pop_cal.add_subplot(int(len(pop_bins) / 4), 4, i + 1)
            ax.plot(years, model, color="orange", alpha=0.7)
            ax.scatter(years, data, color="darkgrey", alpha=0.4)
            ax.set_title(pop)
            ax.set_ylim(bottom=0)
        pop_cal.tight_layout()
        pdf.savefig(pop_cal)
        plt.close()
        # HBsAg prevalence
        prev_cal = plt.figure(figsize=(20, 20))
        for i, pop in enumerate(pop_bins):
            # Extract Data
            years = at.PlotData(res_store["baseline"]["res"][0], outputs="hbsag_prev", pops=[pop], t_bins=1).series[0].tvec
            model = at.PlotData(res_store["baseline"]["res"][0], outputs="hbsag_prev", pops=[pop], t_bins=1).series[0].vals
            data = res_store["baseline"]["P"][0].parsets[0].get_par("hbsag_prev").ts[pop].vals
            data_years = res_store["baseline"]["P"][0].parsets[0].get_par("hbsag_prev").ts[pop].t

            ax = prev_cal.add_subplot(int(len(pop_bins) / 4), 4, i + 1)
            ax.plot(years, model, color="orange", alpha=0.7)
            ax.scatter(data_years, data, color="darkgrey", alpha=0.4)
            ax.set_title(pop)
            ax.set_ylim(bottom=0)
        prev_cal.tight_layout()
        pdf.savefig(prev_cal)
        plt.close()
        # HBeAg prevalence
        eag_cal = plt.figure(figsize=(20, 20))
        for i, pop in enumerate(pop_bins):
            # Extract Data
            years = at.PlotData(res_store["baseline"]["res"][0], outputs="hbeag_prev", pops=[pop], t_bins=1).series[0].tvec
            model = at.PlotData(res_store["baseline"]["res"][0], outputs="hbeag_prev", pops=[pop], t_bins=1).series[0].vals
            data = res_store["baseline"]["P"][0].parsets[0].get_par("hbeag_prev").ts[pop].vals
            data_years = res_store["baseline"]["P"][0].parsets[0].get_par("hbeag_prev").ts[pop].t

            ax = eag_cal.add_subplot(int(len(pop_bins) / 4), 4, i + 1)
            ax.plot(years, model, color="orange", alpha=0.7)
            ax.scatter(data_years, data, color="darkgrey", alpha=0.4)
            ax.set_title(pop)
            ax.set_ylim(bottom=0)
        eag_cal.tight_layout()
        pdf.savefig(eag_cal)
        plt.close()
        # HCC incidence
        hcc_cal = plt.figure(figsize=(20, 20))
        for i, pop in enumerate(pop_bins):
            # Extract Data
            years = at.PlotData(res_store["baseline"]["res"][0], outputs="hcc_inc", pops=[pop], t_bins=1).series[0].tvec
            model = at.PlotData(res_store["baseline"]["res"][0], outputs="hcc_inc", pops=[pop], t_bins=1).series[0].vals
            data = res_store["baseline"]["P"][0].parsets[0].get_par("hcc_inc").ts[pop].vals
            data_years = res_store["baseline"]["P"][0].parsets[0].get_par("hcc_inc").ts[pop].t

            ax = hcc_cal.add_subplot(int(len(pop_bins) / 4), 4, i + 1)
            ax.plot(years, model, color="orange", alpha=0.7)
            ax.scatter(data_years, data, color="darkgrey", alpha=0.4)
            ax.set_title(pop)
            ax.set_ylim(bottom=0)
        hcc_cal.tight_layout()
        pdf.savefig(hcc_cal)
        plt.close()
        # Treatment coverage
        treat_cal = plt.figure(figsize=(20, 20))
        for i, pop in enumerate(pop_bins):
            # Extract Data
            years = at.PlotData(res_store["baseline"]["res"][0], outputs="treated", pops=[pop], t_bins=1).series[0].tvec
            model = at.PlotData(res_store["baseline"]["res"][0], outputs="treated", pops=[pop], t_bins=1).series[0].vals
            data = res_store["baseline"]["P"][0].parsets[0].get_par("treated").ts[pop].vals
            data_years = res_store["baseline"]["P"][0].parsets[0].get_par("treated").ts[pop].t

            ax = treat_cal.add_subplot(int(len(pop_bins) / 4), 4, i + 1)
            ax.plot(years, model, color="orange", alpha=0.7)
            ax.scatter(data_years, data, color="darkgrey", alpha=0.4)
            ax.set_title(pop)
            ax.set_ylim(bottom=0)
        treat_cal.tight_layout()
        pdf.savefig(treat_cal)
        plt.close()

    # Outcomes projection plot
    scen_names = list(res_store.keys())
    fig1= plt.figure(figsize=(20, 20))
    mtct_plt = fig1.add_subplot(2,2,1)
    chb_plt = fig1.add_subplot(2,2,2)
    hcc_plt = fig1.add_subplot(2,2,3)
    dth_plt = fig1.add_subplot(2,2,4)
    for scen in scen_names:
        mtct_plt.plot(at.PlotData(res_store[scen]["res"][0], outputs="b_chb", pops="total", t_bins=1).series[0].tvec, at.PlotData(res_store[scen]["res"][0], outputs="b_chb", pops="total", t_bins=1).series[0].vals, label = scen)
        mtct_plt.set_title("Mother to Child Transmission")
        chb_plt.plot(at.PlotData(res_store[scen]["res"][0], outputs="all_trans", pops="total", t_bins=1).series[0].tvec, at.PlotData(res_store[scen]["res"][0], outputs="all_trans", pops="total", t_bins=1).series[0].vals)
        chb_plt.set_title("Total CHB incidence")
        hcc_plt.plot(at.PlotData(res_store[scen]["res"][0], outputs="hcc_inc", pops="total", t_bins=1).series[0].tvec, at.PlotData(res_store[scen]["res"][0], outputs="hcc_inc", pops="total", t_bins=1).series[0].vals)
        hcc_plt.set_title("HCC Incidence")
        dth_plt.plot(at.PlotData(res_store[scen]["res"][0], outputs="deaths_chb", pops="total", t_bins=1).series[0].tvec, at.PlotData(res_store[scen]["res"][0], outputs="deaths_chb", pops="total", t_bins=1).series[0].vals)
        dth_plt.set_title("HBV Attributable Deaths")
        mtct_plt.legend(loc="best")
    fig1.tight_layout()
    fig1.savefig(test_save_dir / f'{savename}_ZAF_projection plots.pdf', dpi=400)
    plt.close()














# Supplemental Figures (trends, costs, etc.)