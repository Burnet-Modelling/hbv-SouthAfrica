import hbv
import numpy as np

sim_dt = 0.25 # simulation timestep
sim_start = 1990 # simulation start year
sim_end = 2101 # simulation end year
runs = 10  # number of samples for PSA
seed = 20250902 # set seed for numpy random sampling

framework_path = hbv.root/"framework"/"hbv_hepbd_zaf_fw.xlsx" # path to framework
db_path = hbv.root/"data"/"ZAF_databook.xlsx" # path to databook
calibration = hbv.root/"calibration"/"ZAF_calibration.xlsx" # path to calibration (y_) factors
samples_path = hbv.root/"data"/"ZAF_samples.xlsx" # path to model uncertainty bounds
costs_path = hbv.root/"data"/"ZAF_costs.xlsx" # path to cost inputs

# Population Bins
pop_bins = {"0-0":[0,0], "1-4":[1,4], "5-9":[5,9], "10-19":[10,19], "20-29":[20,29], "30-39":[30,39], "40-49":[40,49],
            "50-59":[50,59], "60-69":[60,69], "70-79":[70,79], "80-89":[80,89], "90+":[90,"100+"]} # All model age bins

birth_bins = {"10-19": [15, 19], "20-29": [20,29], "30-39": [30,39], "40-49": [40,49], "total": [15, 49]}  # WOCBA

trans_bins = {"0-0": [0], "1-4": [4], "5-9": [9], "10-19": [19], "20-29": [29], "30-39": [39], "40-49": [49],
              "50-59": [59], "60-69": [69], "70-79": [79], "80-89": [89]}  # Age Transfers

main_scenarios = ["baseline", "sel pap + sel bd", "universal", "sel pap + univ bd"]
add_scenarios = ["sel pap dna + sel bd", "sel pap dna + univ bd"]
supp_scenarios = ["baseline", "sel pap dna + sel bd", "universal", "sel pap dna + univ bd"]
zaf_scenarios = main_scenarios + add_scenarios

main_full = ["Baseline \n (No HepB-BD)", "Selective PAP \n+ Selective HepB-BD", "Universal \n HepB-BD",
             "Selective PAP \n+ Universal HepB-BD"]
add_full = ["Selective PAP (DNA) \n+ Selective HepB-BD", "Selective PAP (DNA) \n+ Universal HepB-BD"]
supp_full = ["Baseline \n (No HepB-BD)", "Selective PAP (DNA) \n+ Selective HepB-BD", "Universal \n HepB-BD", "Selective PAP (DNA) \n+ Universal HepB-BD"]
zaf_scenarios_full = main_full + add_full



zaf_scenarios_selective = ["pessimistic_selective_SA", "optimistic_selective_SA",
                           "pessimistic_selective_WHO", "optimistic_selective_WHO"] #hepbd cost hep_bd_tgt

zaf_scenarios_dna = ["pessimistic_selective_SA", "optimistic_selective_SA",
                     "pessimistic_combined_SA", "optimistic_combined_SA"] #maternal screening cost hbsag + hbv dna (expos_births)

zaf_data_extract = {# Canary in the coalmine
                "m_hcc": ["m_hcc", "total", "weighted"],
                #Transmission Outcomes
                "mtct_infx": ["b_chb", {"infants": ["0-0M", "0-0F"]}, "sum"],
                "echt_infx": ["horiz_trans", {"under 5": ["0-0M", "0-0F", "1-4M", "1-4F"]}, "sum"],
                "oth_infx": ["horiz_trans", {"over 5": ["5-9M", "5-9F", "10-19M", "10-19F", "20-29M", "20-29F", "30-39M", "30-39F",
                                                        "40-49M", "40-49F", "50-59M", "50-59F", "60-69M", "60-69F", "70-79M", "70-79F",
                                                        "80-89M", "80-89F", "90+M", "90+F"]}, "sum"],
                "tot_inc": ["all_trans", "total", "sum"],
                "expos_births": ["expos_births", {"infants": ["0-0M", "0-0F"]}, "sum"],
                "births": ["b_rate", {"infants": ["0-0M", "0-0F"]}, "sum" ],
                # Vaccination Coverage/Costs
                "anc_screening": ["n_anc_scr", {"infants": ["0-0M", "0-0F"]}, "sum"], # this is only HBsAg screening
                "dna_screening": ["n_hbv_dna",  {"infants": ["0-0M", "0-0F"]}, "sum"],
                "mav_hepbd_cov": ["n_mav_bd", {"infants": ["0-0M", "0-0F"]}, "sum" ],
                "selective_hepbd": ["n_select_bd",  {"infants": ["0-0M", "0-0F"]}, "sum"],
                "universal_hepbd": ["n_univ_bd",  {"infants": ["0-0M", "0-0F"]}, "sum"],
                "select_hepbd_univ": ["n_univ_bd_additional", {"infants": ["0-0M", "0-0F"]}, "sum"], #number of hepB-BD administered if univ included neg HBsAg tests
                "total_hepbd": ["bd_vax_total",  {"infants": ["0-0M", "0-0F"]}, "sum"],
                "total_hepb3": ["hepb3_vax_total",  {"infants": ["0-0M", "0-0F"]}, "sum"],
                # Disease Management Costs/Coverage
                "diag_cost":["new_diag", "total", "sum"],
                "treat_cost": ["treated", "total", "sum"],
                "noncirr_cost": ["cost_noncir", "total", "sum"],
                "compcirr_cost": ["cost_cir", "total", "sum"],
                "decomp_cost": ["cost_decomp", "total", "sum"],
                "cancer_cost": ["cost_canc", "total", "sum"],
                # Epidemiological Outcomes
                "hbv_deaths": ["deaths_chb", "total", "sum"],
                "hcc_incidence": ["hcc_inc", "total", "sum"],
                # Health Utility Inputs (DALYs)
                "decomp_daly":["decompensated", "total", "sum"],
                "cancer_daly":["cancer", "total", "sum"],
                "yll_daly": ["yll", "total", "sum"],
                # WHO Elimination Indicators
                "who_inc_u5": ["who_inc_meas", {"under 5": ["0-0M", "0-0F", "1-4M", "1-4F"]}, "weighted"],
                "who_inc_pop": ["who_inc_meas", "total", "weighted"],
                "who_mort": ["who_mort_meas", "total", "weighted"]
                }

wtp_daly = np.arange(0, 25e3*18.33, 0.5e2*18.33) #for cost-effectiveness threshold analysis