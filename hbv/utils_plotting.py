from hbv import constants as const
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 16}) # Sets the default font size to 16
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import ticker
import atomica as at
import pandas as pd
import numpy as np
import pathlib
import sciris as sc



def cost_eff_frontier(res_save_dir):

    main_scens = ["sel pap + sel bd", "universal", "sel pap + univ bd"]
    sens_scens = ["sel pap dna + sel bd", "sel pap dna + univ bd"]

    half_gdp = (0.5*6250)*18.33
    full_gdp = (1*6250)*18.33
    two_gdp = (2*6250)*18.33
    x_upper_bound = (2.5*6250)*18.33

    # Make save folder for results presented in manuscript
    main_save_dir = pathlib.Path(res_save_dir/"Main Text Results/")
    main_save_dir = pathlib.Path(main_save_dir)
    main_save_dir.mkdir(parents=True, exist_ok=True)

    # Make save folder for results presented in supplement
    supp_save_dir = pathlib.Path(res_save_dir/"Supplemental Results/")
    supp_save_dir = pathlib.Path(supp_save_dir)
    supp_save_dir.mkdir(parents=True, exist_ok=True)

    # Figure 1: CEF of selective_WHO, universal, combined_WHO to 2100, no HepB-BD in HBsAg screened
    cef_data_f1 = sc.load(res_save_dir/'cost_pickles'/f"ZAF_{const.runs}_cef_2025_2100_main.pkl")

    # Empty Dataframe Generate
    f1_data = pd.DataFrame(columns = ["wtp_thresh", "baseline"] + main_scens)
    f1_data.wtp_thresh = const.wtp_daly
    f1_data = f1_data.set_index("wtp_thresh", drop = False)

    # Extract Data Required for CEF
    f1_keep_cols = ["baseline"] +  main_scens
    f1_model_data = {k: v[f1_keep_cols] for k, v in cef_data_f1.items()}

    # Return proportion of runs with highest net benefit
    for wtp in const.wtp_daly:
        mins_f1 = f1_model_data[f"{wtp}_thresh"].eq(f1_model_data[f"{wtp}_thresh"].min(axis=1),axis=0).astype(int)
        mins_f1 = mins_f1.sum(numeric_only = True).to_frame().T/const.runs

        for scen in f1_keep_cols:
            f1_data.at[wtp, scen] = mins_f1[scen].values

    # Plot Data
    full_scen = const.main_full
    fig1 = plt.figure(figsize=(20,15))
    ax = fig1.add_subplot(1,1,1)
    for idx, scen in enumerate(f1_keep_cols):
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
    fig1.savefig(main_save_dir / "Figure 1 CEF.pdf", dpi=400)
    plt.close()

    # Figure 2: CEF of selective_WHO, universal, combined_WHO to 2100, yes HepB-BD in HBsAg screened (sensitivity)
    cef_data_f2 = sc.load(res_save_dir/'cost_pickles'/f"ZAF_{const.runs}_cef_2025_2100_supp.pkl")

    f2_data = pd.DataFrame(columns=["wtp_thresh", "baseline"] +  main_scens)
    f2_data.wtp_thresh = const.wtp_daly
    f2_data = f2_data.set_index("wtp_thresh", drop=False)

    f2_keep_cols = ["baseline"] + main_scens
    f2_model_data = {k: v[f2_keep_cols] for k, v in cef_data_f2.items()}

    for wtp in const.wtp_daly:
        mins_f2 = f2_model_data[f"{wtp}_thresh"].eq(f2_model_data[f"{wtp}_thresh"].min(axis=1), axis=0).astype(int)
        mins_f2 = mins_f2.sum(numeric_only=True).to_frame().T / const.runs

        for scen in f2_keep_cols:
            f2_data.at[wtp, scen] = mins_f2[scen].values

    fig2 = plt.figure(figsize=(20,15))
    ax = fig2.add_subplot(1,1,1)
    for idx, scen in enumerate(f2_keep_cols):
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
    fig2.savefig(main_save_dir/ "Figure 2 CEF_sens.pdf", dpi=400)
    plt.close()

    # Supplement: CEF of selective_SA, universal, combined_SA to 2100, yes HepB-BD in HBsAg screened (sensitivity; uses cef_data_f1, cef_data_f2)

    s1_data = pd.DataFrame(columns=["wtp_thresh", "baseline"] + sens_scens) # top plot
    s1_data.wtp_thresh = const.wtp_daly
    s1_data = s1_data.set_index("wtp_thresh", drop=False)

    s2_data = pd.DataFrame(columns=["wtp_thresh", "baseline"] +  sens_scens)# bottom plot
    s2_data.wtp_thresh = const.wtp_daly
    s2_data = s2_data.set_index("wtp_thresh", drop=False)

    sup_keep_cols = ["baseline"] + sens_scens
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
    for idx,scen in enumerate(sup_keep_cols):
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
    for idx,scen in enumerate(sup_keep_cols):
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
    sup_cef.savefig(supp_save_dir/ "Supp CEF_HBV DNA screen.png", dpi=400)
    plt.close()

# Transmission Route Plot (Fig 2)
def transmission_plots(res_save_dir):

    # Supplemental Figure
    fig2_data = pd.read_excel(res_save_dir/'raw_tab_figs'/'CHB_transmission_source.xlsx', sheet_name="baseline")
    fig2 = plt.figure(figsize=(20,15))
    trans = fig2.add_subplot(1,1,1)
    trans.bar(fig2_data.year, fig2_data.mtct, color = "blue", alpha=0.6, label = "MTCT")
    trans.bar(fig2_data.year, fig2_data.echt, bottom = fig2_data.mtct, color="green", alpha = 0.6, label = "Early Childhood")
    trans.bar(fig2_data.year, fig2_data.oth, bottom = fig2_data.mtct + fig2_data.echt, color="orange", alpha = 0.6, label = "Other")
    trans.legend (loc="best")
    trans.set_xlim(1990, 2100)
    trans.set_ylabel("Total incident CHB infections")
    fig2.tight_layout()
    fig2.savefig(res_save_dir/"Supplemental Results/"/'CHB by source_baseline.png', dpi=400)
    plt.close()





