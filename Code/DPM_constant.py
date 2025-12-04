from DPM_lib import np
'''
This script defines the model parameter names, input parameter names, and the default values and ranges used in the DPM model.#
'''
'''
CHECKED
'''


# Generate the parameter names used for different number of drugs.#
def DPM_constant_parname(parname, celltype_name, drug_name):
    # Growth rate for each cell type.#
    g0_name = [f'g0_{i_celltype}' for i_celltype in celltype_name]

    # Initial cell count for each cell type.#
    X0_name = [f'{i_celltype}pop' for i_celltype in celltype_name]

    # 'Sa' matrix (number of cell types * number of drugs).#
    Sa_name = [f'Sa.{i_celltype}.{i_drugtype}.' for i_celltype in celltype_name for i_drugtype in drug_name]

    # 'T' matrix (number of cell types * Number of cell types). T.S..S. from S to S, T.S..R1. from R1 to S.#
    T_name = [f'T.{i_celltype}..{j_celltype}.' for i_celltype in celltype_name for j_celltype in celltype_name]

    return parname + X0_name + g0_name + Sa_name + T_name


PARFORMAT_LIST_0 = ['paramID', 'Num_drug', 'Num_cell_type']

# Total number of parameters in the 2-drug case.#
ALL_POSSIBLE_CELLTYPE_2DRUG = ['S', 'R1', 'R2', 'R12']
ALL_POSSIBLE_DRUG_2DRUG = ['D1', 'D2']
PARFORMAT_LIST_2DRUG = DPM_constant_parname(PARFORMAT_LIST_0, ALL_POSSIBLE_CELLTYPE_2DRUG, ALL_POSSIBLE_DRUG_2DRUG, )

# Total number of parameters in the 3-drug case.#
ALL_POSSIBLE_CELLTYPE_3DRUG = ['S', 'R1', 'R2', 'R3', 'R12', 'R13', 'R23', 'R123']
ALL_POSSIBLE_DRUG_3DRUG = ['D1', 'D2', 'D3']
PARFORMAT_LIST_3DRUG = DPM_constant_parname(PARFORMAT_LIST_0, ALL_POSSIBLE_CELLTYPE_3DRUG, ALL_POSSIBLE_DRUG_3DRUG)

# Total cell population of X0.#
X0TOTAL_DEFAULT_VAL = 5e9
X0TOTAL_UPPERTHRES = 1e13
X0T0TAL_LOWERTHRES = 1e5

# Number of drugs.#
NUM_DRUG_LIST = [2, 3]
NUM_DRUG_DEFAULT_VAL = 2
NUM_DRUG_UPPERTHRES = 3
NUM_DRUG_LOWERTHRES = 2

# Simulation duration (days).#
SIMDURATION_DEFAULT_VAL = 1800
SIMDURATION_UPPERTHRES = 7300
SIMDURATION_LOWERTHRES = 90

# Limit mortality (number of cells).#
LIMIT_MORTALITY_DEFAULT_VAL = 1e13
LIMIT_MORTALITY_UPPERTHRES = 1e14
LIMIT_MORTALITY_LOWERTHRES = 1e9

# Step size (days).#
STEPSIZE_DEFAULT_VAL = 45
STEPSIZE_UPPERTHRES = 200
STEPSIZE_LOWERTHRES = 15

# Simulation time step (days).#
SIMTIMESTEP_DEFAULT_VAL = 1
SIMTIMESTEP_UPPERTHRES = 2
SIMTIMESTEP_LOWERTHRES = 1

# Limit molecular detection (1/number of cells).#
LIMIT_MOLECULARDETECTION_DEFAULT_VAL = 1/1e4
LIMIT_MOLECULARDETECTION_UPPERTHRES = 1/1e1
LIMIT_MOLECULARDETECTION_LOWERTHRES = 1/1e9

# Limit radiologic detection (number of cells).#
LIMIT_RADIOLOGICDETECTION_DEFAULT_VAL = 1e9
LIMIT_RADIOLOGICDETECTION_UPPERTHRES = 1e9
LIMIT_RADIOLOGICDETECTION_LOWERTHRES = 1e8

# Parameter save block size: number of parameter sets saved per file.#
PAR_SAVE_BLOCK_SIZE_DEFAULT_VAL = 1e4
PAR_SAVE_BLOCK_SIZE_UPPERTHRES = 5e6
PAR_SAVE_BLOCK_SIZE_LOWERTHRES = 1e4

# Dose method: "continous" = "c" dosing; "discrete" = "d" dosing.#
DOSE_METHOD_LIST = ['continous', 'c', 'discrete', 'd']
# Convert the string to lowercase.#
DOSE_METHOD_LIST = [i.lower() for i in DOSE_METHOD_LIST]
DOSE_METHOD_DEFAULT_VAL = 'd'

# Dose interval: dose for a single drug is: (0:dose interval:1).#
DOSE_INTERVAL_LIST = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
DOSE_INTERVAL_DEFAULT_VAL = 0.1
DOSE_INTERVAL_UPPERTHRES = 0.5
DOSE_INTERVAL_LOWERTHRES = 0.01

# Threshold for the sum of dose combinations.#
# Example: drug1 + drug2 <= 2.#
# Example: drug1 + drug2 >= 1.#
DOSE_COMBINATION_SUM_UPPERTHRES = 2
DOSE_COMBINATION_SUM_LOWERTHRES = 1

# Dose thresholds for each drug within dose combination.#
# Example: dose combination consists of drug1 and drug2.#
# 0 <= drug1 <= 1.#
# 0 <= drug2 <= 1.#
DOSE_COMBINATION_ELEMENT_UPPERTHRES = 1
DOSE_COMBINATION_ELEMENT_LOWERTHRES = 0

# Number of steps that differ between CPM and DPM.#
NUM_STEPDIFF_DEFAULT_VAL = 2
NUM_STEPDIFF_UPPERTHRES = 5
NUM_STEPDIFF_LOWERTHRES = 0

# Sample size specified in the 'DPM_run_drugmis_analysis' function, this represents the number of patients required for #
# the different study arms.#
SAMPLESIZE = 125

# Number of bootstrap samples used to generate the confidence interval.#
NSAMPLE = 10000

# Specifies the number of iterations used to compute the confidence interval in the DPM_run_drugmis_ci function.#
N_ITERATIONS = 10000

# Alpha value used to compute the confidence interval in the DPM_run_drugmis_ci function.#
ALPHA_CI = 0.95

# The p-value used to determine the false-positive and false-negtive rates.#
PVALUE_FNFP = 0.025

# Strategy names.#
STRATEGY_LIST = ['CPM', 'DPM1', 'DPM2.1', 'DPM2.2', 'DPM3', 'DPM4']

# Parameter settings for the SciPy integrate function.#
USE_MATRIXEXP = False
SOLVER_ALL = {'ODEINT', 'SOLVE_IVP'}
SOLVER = 'SOLVE_IVP'  # 'ODEINT'  # 'SOLVE_IVP'
SCIPY_INTEGRATE_ODEINT = {'RTOL': 1e-3, 'ATOL': 1e-6, 'hmax': 1}
SCIPY_INTEGRATE_SOLVE_IVP = {"RTOL": 1e-3, "ATOL": 1e-6, "max_step": 1}

# Plotting parameters.#
PLT_INTERACTIVE = False
# Values of X smaller than MINX will be set to MINX.#
MINX = .1
XTOTAL_RATIO = 1.5
COLOR_X_2DRUG = ('b', 'g', 'c', 'm', 'r')
LABEL_X_2DRUG = ('S cells', 'R1 cells', 'R2 cells', 'R12 cells')
COLOR_DRUG_2DRUG = ('g', 'b')
LABEL_2DRUG = ('Drug 1', 'Drug 2')
LEGEND_COLNUM_2DRUG = 6
LEGEND_ORDER_2DRUG = np.array([0, 4, 1, 5, 2, 3])
FIG_WIDTH_2DRUG = 8.5
FIG_HEIGTH_2DRUG = 6.5

COLOR_X_3DRUG = ('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
LABEL_X_3DRUG = ('Total cells', 'S cells', 'R1 cells', 'R2 cells', 'R3 cells', 'R12 cells', 'R13 cells', 'R23 cells', 'R123 cells')
COLOR_DRUG_3DRUG = ('#1f77b4', '#ff7f0e', '#2ca02c')
LABEL_3DRUG = ('Drug 1', 'Drug 2', 'Drug 3')
LEGEND_COLNUM_3DRUG = 7
LEGEND_ORDER_3DRUG = np.reshape(np.array(range(0, 12)), (-1, 6)).flatten('F')
LEGEND_ORDER_3DRUG = np.concatenate((LEGEND_ORDER_3DRUG, np.array([12], int)))
FIG_WIDTH_3DRUG = 12
FIG_HEIGTH_3DRUG = 9

COLOR_STRATEGY = ('black', 'maroon', 'red', 'sienna', 'darkgoldenrod', 'darkolivegreen', 'forestgreen', 'darkcyan', 'dodgerblue',
                  'darkviolet', 'deeppink')
FIG_DPI = 300
LEGEND_BOXANCHOR = (0.5, 1.15)
YMINVAL = 1e-2
XMINVAL = -1

# Divide t to T_NORM to plot time in month instead of days.#
T_NORM = 1
LOGBASE = 10
YINCREASE = 2
XINCREASE = 2
TEXTFONT = 10
YMINVAL_XALL = 0
LEGEND_COLNUM_XALL = 6
FIG_WIDTH_XALL = 8.5

PLT_MARKERSIZE = 5
ALPHA = 0.5

# Header of the parameter input .csv file.#
HEADER_2DRUG_INPUT_PARAM_CSV = PARFORMAT_LIST_2DRUG[:1] + PARFORMAT_LIST_2DRUG[3:]
HEADER_3DRUG_INPUT_PARAM_CSV = PARFORMAT_LIST_3DRUG[:1] + PARFORMAT_LIST_3DRUG[3:]

# Header of the result parameter .csv file.#
HEADER_2DRUG_PARAM_CSV = HEADER_2DRUG_INPUT_PARAM_CSV
HEADER_3DRUG_PARAM_CSV = HEADER_3DRUG_INPUT_PARAM_CSV

# Header of the result stopt .csv file.#
HEADER_STOPT_CSV = ['paramID']
HEADER_STOPT_CSV.extend(STRATEGY_LIST[:-1])

# Header of the result dosage .csv file.#
HEADER_DOSAGE_CSV = ['paramID', 'Strategy name']

# Header of the result pop .csv file.#
HEADER_POP_CSV = ['paramID', 'Strategy name']

# Header of the result eachtimepoint .csv file.#
HEADER_EACHTIMEPOINT_CSV = ['paramID', 'Strategy name']

# Maximum file name length.#
MAXFILENAMELEN = 200  # os.pathconf('/', 'PC_NAME_MAX') - 50

# Input parameter keywords for the function 'DPM_run_par_csv_folder'.#
PARNAME_LIST_RUN_PAR_CSV_FOLDER = \
    ['pathload', 'pathsave', 'use_parallel', 'filename_pattern', 'misspecification_atdecision', 'misspecification_atsim',
     'Num_drug', 'Limit_moleculardetection', 'Limit_radiologicdetection', 'Limit_mortality', 'Simduration', 'Stepsize', 'dose_method',
     'dose_interval', 'dose_combination', 'use_input_dose_combination_only', 'PAR_criterisa', 'Strategy_name', 'run_sim',
     'save_filename_param', 'save_filename_stopt', 'save_filename_pop', 'save_filename_dosage', 'save_filename_eachtimepoint',
     'misspecification_sim_only', 'misspecification_fileload']

# Input parameter keywords for the function 'DPM_run_par_csv'.#
PARNAME_LIST_RUN_PAR_CSV = \
    ['filename_csv', 'pathsave', 'misspecification_filename_csv', 'misspecification_atdecision', 'misspecification_atsim',
     'Num_drug', 'Limit_moleculardetection', 'Limit_radiologicdetection', 'Limit_mortality', 'Simduration', 'Stepsize', 'dose_method',
     'dose_interval', 'dose_combination', 'use_input_dose_combination_only', 'PAR_criterisa', 'Strategy_name', 'run_sim',
     'save_filename_param', 'save_filename_stopt', 'save_filename_pop', 'save_filename_dosage', 'save_filename_eachtimepoint',
     'misspecification_sim_only', 'par_ind']

# Input parameter keywords for the function 'DPM_run_plot_1PAR'.#
PARNAME_LIST_RUN_PLOT_1PAR = \
    ['par', 'mis_par', 'Limit_radiologicdetection', 'Limit_moleculardetection', 'Limit_mortality', 'Simduration', 'Stepsize', 'dose_method',
     'dose_interval', 'dose_combination', 'use_input_dose_combination_only', 'PAR_criterisa', 'Strategy_name', 'save_filename_param',
     'save_filename_stopt', 'save_filename_pop', 'save_filename_dosage', 'save_filename_eachtimepoint', 'pathsave', 'plot', 'savename',
     'misspecification_atdecision', 'misspecification_atsim']

# Input parameter keywords for the function 'DPM_run_processing'.#
PARNAME_LIST_RUN_PROCESSING = ['Num_drug', 'Strategy_name', 'pathsave', 'pathload', 'Num_stepdiff', 'Simduration', 'Stepsize',
                               'misspecification_fileload', 'use_parallel']

# Total input parameter names.#
PARNAME_ALLFUN = [PARNAME_LIST_RUN_PAR_CSV_FOLDER, PARNAME_LIST_RUN_PAR_CSV, PARNAME_LIST_RUN_PLOT_1PAR, PARNAME_LIST_RUN_PROCESSING]
PARNAME_LIST = list({PAR for FUN_LIST in PARNAME_ALLFUN for PAR in FUN_LIST})


# Specify the colors used for printec text.#
class Colored:
    PURPLE: str = '\033[95m'
    CYAN: str = '\033[96m'
    DARKCYAN: str = '\033[36m'
    BLUE: str = '\033[94m'
    GREEN: str = '\033[92m'
    YELLOW: str = '\033[93m'
    RED: str = '\033[91m'
    BOLD: str = '\033[1m'
    UNDERLINE: str = '\033[4m'
    ITALICS: str = '\x1B[3m'
    END: str = '\033[0m'
