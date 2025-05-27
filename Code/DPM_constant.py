from DPM_lib import np
from scipy.stats import rv_continuous, loguniform

""" This script sets the model parameter names, input parameter names, parameter default values and ranges used in the
DPM (Dynamic Precision Medicine) model."""


# Generate the parameter names used in different drug numbers.
def DPM_constant_parname(parname, celltype_name, drug_name):
    # Growth rate for each cell type.
    g0_name = [f'g0_{i_celltype}' for i_celltype in celltype_name]

    # Initial cell number for each cell type.
    X0_name = [f'{i_celltype}pop' for i_celltype in celltype_name]

    # Sa matrix (Number of cell types * number of drugs).
    Sa_name = [f'Sa.{i_celltype}.{i_drugtype}.' for i_celltype in celltype_name for i_drugtype in drug_name]

    # T matrix (Number of cell types * Number of cell types). T.S..S. from S to S, T.S..R1. from R1 to S
    T_name = [f'T.{i_celltype}..{j_celltype}.' for i_celltype in celltype_name for j_celltype in celltype_name]

    return parname + X0_name + g0_name + Sa_name + T_name


PARFORMAT_LIST_0 = ['paramID', 'Num_drug', 'Num_cell_type']

# Total parameters in 4 drug case. These are the all the possiable parameters, include all parameter names used in 2 and 3 drug case.
# ALL_POSSIBLE_CELLTYPE = ['S', 'R1', 'R2', 'R3', 'R4', 'R12', 'R13', 'R14', 'R23', 'R24', 'R34', 'R123', 'R124', 'R134', 'R234', 'R1234']
# ALL_POSSIBLE_DRUG = ['D1', 'D2', 'D3', 'D4']
# PARFORMAT_LIST = DPM_constant_parname(PARFORMAT_LIST_0, ALL_POSSIBLE_CELLTYPE, ALL_POSSIBLE_DRUG)

# Total parameters in 2 drug case.
ALL_POSSIBLE_CELLTYPE_2DRUG = ['S', 'R1', 'R2', 'R12']
ALL_POSSIBLE_DRUG_2DRUG = ['D1', 'D2']
PARFORMAT_LIST_2DRUG = DPM_constant_parname(PARFORMAT_LIST_0, ALL_POSSIBLE_CELLTYPE_2DRUG, ALL_POSSIBLE_DRUG_2DRUG, )

# Total parameters in 3 drug case.
ALL_POSSIBLE_CELLTYPE_3DRUG = ['S', 'R1', 'R2', 'R3', 'R12', 'R13', 'R23', 'R123']
ALL_POSSIBLE_DRUG_3DRUG = ['D1', 'D2', 'D3']
PARFORMAT_LIST_3DRUG = DPM_constant_parname(PARFORMAT_LIST_0, ALL_POSSIBLE_CELLTYPE_3DRUG, ALL_POSSIBLE_DRUG_3DRUG)

# Total cell population of X0.
X0TOTAL_DEFAULT_VAL = 5e9
X0TOTAL_UPPERTHRES = 1e13
X0T0TAL_LOWERTHRES = 1e5

# Number of drug.
NUM_DRUG_LIST = [2, 3, 4]
NUM_DRUG_DEFAULT_VAL = 2
NUM_DRUG_UPPERTHRES = 4
NUM_DRUG_LOWERTHRES = 2

# Simulation duration: days.
SIMDURATION_DEFAULT_VAL = 1800
SIMDURATION_UPPERTHRES = 7300
SIMDURATION_LOWERTHRES = 90

# Limit mortality: number of cells.
LIMIT_MORTALITY_DEFAULT_VAL = 1e13
LIMIT_MORTALITY_UPPERTHRES = 1e14
LIMIT_MORTALITY_LOWERTHRES = 1e9

# Stepsize: days.
STEPSIZE_DEFAULT_VAL = 45
STEPSIZE_UPPERTHRES = 200
STEPSIZE_LOWERTHRES = 15

# Simulation timestep: days.
SIMTIMESTEP_DEFAULT_VAL = 1
SIMTIMESTEP_UPPERTHRES = 2
SIMTIMESTEP_LOWERTHRES = 1

# Limit molecular detection: 1/(number of cells).
LIMIT_MOLECULARDETECTION_DEFAULT_VAL = 1/1e4
LIMIT_MOLECULARDETECTION_UPPERTHRES = 1/1e1
LIMIT_MOLECULARDETECTION_LOWERTHRES = 1/1e9

# Limit radiologic detection: number of cells.
LIMIT_RADIOLOGICDETECTION_DEFAULT_VAL = 1e9
LIMIT_RADIOLOGICDETECTION_UPPERTHRES = 1e9
LIMIT_RADIOLOGICDETECTION_LOWERTHRES = 1e8

# Parameter save block size: every # of parameter saves to a file.
PAR_SAVE_BLOCK_SIZE_DEFAULT_VAL = 1e4
PAR_SAVE_BLOCK_SIZE_UPPERTHRES = 5e6
PAR_SAVE_BLOCK_SIZE_LOWERTHRES = 1e4

# Dose method: "continous" = "c" dosing; "discrete" = "d" dosing.
DOSE_METHOD_LIST = ['continous', 'c', 'discrete', 'd']
# change string to lower case
DOSE_METHOD_LIST = [i.lower() for i in DOSE_METHOD_LIST]
DOSE_METHOD_DEFAULT_VAL = 'd'

# Dose interval: dose for a single drug is 0: dose interval: 1.
DOSE_INTERVAL_LIST = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
DOSE_INTERVAL_DEFAULT_VAL = 0.1
DOSE_INTERVAL_UPPERTHRES = 0.5
DOSE_INTERVAL_LOWERTHRES = 0.01

# Dose combination sum together threshold.
# E.g. drug1 dose + drug2 dose + drug3 dose <= 5.
# E.g. drug1 dose + drug2 dose + drug3 dose >= 1.
DOSE_COMBINATION_SUM_UPPERTHRES = 5
DOSE_COMBINATION_SUM_LOWERTHRES = 1

# Dose threshold of each drug in dose combination.
# E.g. dose combination is drug1 + drug2 + drug3.
# 0 <= drug1 <= 2 in the combination.
# 0 <= drug2 <= 2 in the combination.
# 0 <= drug3 <= 2 in the combination.
DOSE_COMBINATION_ELEMENT_UPPERTHRES = 2
DOSE_COMBINATION_ELEMENT_LOWERTHRES = 0

# Look-ahead steps in multistep extension strategies (strategy 5, strategy 6, strategy 7, strategy 8).
LOOKAHEAD_STEP_DEFAULT_VAL = 5
LOOKAHEAD_STEP_UPPERTHRES = 5
LOOKAHEAD_STEP_LOWERTHRES = 2

# Maximum number of tracked subsequences in strategy 9 (ALTO: Adaptive Long-Term Optimization).
MAXNUM_SUBSEQ_DEFAULT_VAL_2DRUG = 500
MAXNUM_SUBSEQ_UPPERTHRES_2DRUG = 500
MAXNUM_SUBSEQ_LOWERTHRES_2DRUG = 100

MAXNUM_SUBSEQ_DEFAULT_VAL_3DRUG = 50
MAXNUM_SUBSEQ_UPPERTHRES_3DRUG = 50
MAXNUM_SUBSEQ_LOWERTHRES_3DRUG = 10

MAXNUM_SUBSEQ_DEFAULT_VAL_4DRUG = 20
MAXNUM_SUBSEQ_UPPERTHRES_4DRUG = 20
MAXNUM_SUBSEQ_LOWERTHRES_4DRUG = 5

# Subtree depth of strategy 9 (ALTO: Adaptive Long-Term Optimization).
SUBTREEDEPTH_DEFAULT_VAL = 5
SUBTREEDEPTH_UPPERTHRES = 5
SUBTREEDEPTH_LOWERTHRES = 2

# Number of steps that are different between CPM and DPM
NUM_STEPDIFF_DEFAULT_VAL = 2
NUM_STEPDIFF_UPPERTHRES = 5
NUM_STEPDIFF_LOWERTHRES = 0

# Subcolone limit of detection: fraction of subcolone on total cell number.
SUBCLONE_LOD_DEFAULT_VAL = 1e-4
SUBCLONE_LOD_UPPERTHRES = 5e-1
SUBCLONE_LOD_LOWERTHRES = 1e-8

# Mutation rate: per base pair per cell division
MUTATION_RATE_DEFAULT_VAL = 7.1e-7
MUTATION_RATE_UPPERTHRES = 1e-6
MUTATION_RATE_LOWERTHRES = 1e-7

# Liquid biopsy subcolone limit of detection ratio compared to tissue biopy
LIQUID_BIOPSY_LOD_RATIO_DEFAULT_VAl = 10
LIQUID_BIOPSY_LOD_RATIO_UPPERTHRES = 100
LIQUID_BIOPSY_LOD_RATIO_LOWERTHRES = 1


# Strategy names.
STRATEGY_LIST = ['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6', 'strategy7',
                 'strategy8', 'strategy9', 'optimal', 'all']
STRATEGY_LIST = [i.lower() for i in STRATEGY_LIST]

# LOD (limit of detection) values.
LOD_LIST = ['1e-01', '1e-02', '1e-03', '1e-04', '1e-05', '1e-06']
LOD_LIST_ALL = ['1e-06', '2e-06', '4e-06', '6e-06', '8e-06',
                '1e-05', '2e-05', '4e-05', '6e-05', '8e-05',
                '1e-04', '2e-04', '4e-04', '6e-04', '8e-04',
                '1e-03', '2e-03', '4e-03', '6e-03', '8e-03',
                '1e-02', '2e-02', '4e-02', '6e-02', '8e-02',
                '1e-01']

# misspecification_LOD values
MISSPECIFICATION_LOD = [0, 'pdf', 'loguni', 'max']
MISSPECIFICATION_LOD_STR = [str(i_val) for i_val in MISSPECIFICATION_LOD]

# Scipy integrate function parameter settings.
USE_MATRIXEXP = False
SOLVER_ALL = {'ODEINT', 'SOLVE_IVP'}
SOLVER = 'SOLVE_IVP'  # 'ODEINT'  # 'SOLVE_IVP'
SCIPY_INTEGRATE_ODEINT = {'RTOL': 1e-3, 'ATOL': 1e-6, 'hmax': 1}
SCIPY_INTEGRATE_SOLVE_IVP = {"RTOL": 1e-3, "ATOL": 1e-6, "max_step": 1}

# Plot parameters.
PLT_INTERACTIVE = True
# X smaller than MINX will set to MINX.
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

# Divided t to T_NORM, plot time in month not days.
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

# Heading of para input .csv file.
HEADING_2DRUG_INPUT_PARAM_CSV = PARFORMAT_LIST_2DRUG[:1] + PARFORMAT_LIST_2DRUG[3:]
HEADING_3DRUG_INPUT_PARAM_CSV = PARFORMAT_LIST_3DRUG[:1] + PARFORMAT_LIST_3DRUG[3:]

# Heading of result param .csv file.
HEADING_2DRUG_PARAM_CSV = HEADING_2DRUG_INPUT_PARAM_CSV
HEADING_3DRUG_PARAM_CSV = HEADING_3DRUG_INPUT_PARAM_CSV

# Heading of result stopt .csv file.
HEADING_STOPT_CSV = ['paramID']
HEADING_STOPT_CSV.extend(STRATEGY_LIST[:-1])

# Heading of result dosage .csv file.
HEADING_DOSAGE_CSV = ['paramID', 'Strategy name']

# Heading of result pop .csv file.
HEADING_POP_CSV = ['paramID', 'Strategy name']

# Heading of result eachtimepoint .csv file.
HEADING_EACHTIMEPOINT_CSV = ['paramID', 'Strategy name']

# Maximum file name length.
MAXFILENAMELEN = 200  # os.pathconf('/', 'PC_NAME_MAX') - 50

# Input parameter keywords for function DPM_run_par_csv_folder.
PARNAME_LIST_RUN_PAR_CSV_FOLDER = \
    ['pathload', 'pathsave', 'use_parallel', 'filename_pattern', 'misspecification_ofdecision', 'misspecification_atsim',
     'Num_drug', 'Limit_moleculardetection', 'Limit_radiologicdetection', 'Limit_mortality', 'Simduration', 'Stepsize', 'dose_method',
     'dose_interval', 'dose_combination', 'use_input_dose_combination_only', 'PAR_criterisa', 'Strategy_name', 'lookahead_step', 'Maxnum_subseq',
     'subtreedepth', 'run_sim', 'erase_preresult', 'fullinput', 'save_filename_param', 'save_filename_stopt', 'save_filename_pop',
     'save_filename_dosage', 'save_filename_eachtimepoint', 'misspecification_sim_only', 'subclone_LOD', 'misspecification_LOD', 'mutation_rate',
     'misspecification_fileload']

# Input parameter keywords for function DPM_run_par_csv.
PARNAME_LIST_RUN_PAR_CSV = \
    ['filename_csv', 'Num_drug', 'Limit_moleculardetection', 'Limit_radiologicdetection', 'Limit_mortality', 'Simduration', 'Stepsize',
     'dose_method', 'dose_interval', 'dose_combination', 'use_input_dose_combination_only', 'PAR_criterisa', 'Strategy_name',
     'lookahead_step', 'Maxnum_subseq', 'subtreedepth', 'run_sim', 'erase_preresult', 'fullinput', 'misspecification_filename_csv',
     'misspecification_ofdecision', 'misspecification_atsim', 'save_filename_param', 'save_filename_stopt', 'save_filename_pop',
     'misspecification_ofdecision', 'misspecification_atsim', 'save_filename_param', 'save_filename_stopt', 'save_filename_pop',
     'save_filename_dosage', 'save_filename_eachtimepoint', 'misspecification_sim_only', 'pathsave', 'par_ind', 'subclone_LOD',
     'misspecification_LOD', 'mutation_rate']

# Input parameter keywords function DPM_run_plot_1PAR.
PARNAME_LIST_RUN_PLOT_1PAR = \
    ['par', 'mis_par', 'Limit_radiologicdetection', 'Limit_mortality', 'Simduration', 'Stepsize', 'dose_method', 'dose_interval',
     'dose_combination', 'use_input_dose_combination_only', 'PAR_criterisa', 'Strategy_name', 'lookahead_step', 'Maxnum_subseq',
     'subtreedepth', 'erase_preresult', 'save_filename_param', 'save_filename_stopt', 'save_filename_pop', 'save_filename_dosage',
     'save_filename_eachtimepoint', 'pathsave', 'plot', 'savename', 'Limit_moleculardetection', 'misspecification_ofdecision',
     'misspecification_atsim', 'subclone_LOD']

# Input parameter keywords of function DPM_save_default_PARset.
PARNAME_LIST_SAVE_DEFAULT_PARSET = ['Num_drug', 'par_save_block_size', 'Limit_mortality', 'Simduration', 'pathsave']

# Input parameter keywords for function DPM_run_processing.
PARNAME_LIST_RUN_PROCESSING = ['Num_drug', 'Strategy_name', 'pathsave', 'pathload', 'Num_stepdiff', 'Simduration', 'Stepsize',
                               'pathloadmis', 'use_parallel']

# Input parameter keywords for function DPM_run_analysis.
PARNAME_LIST_RUN_ANALYSIS_LOD = ['Num_drug',
                                 'Simduration',
                                 'Strategy_name',
                                 'filename_pattern',
                                 'pathsave',
                                 'pathload',
                                 'pathloadmis',
                                 'LOD',
                                 'plot']

# Total input parameter names.
PARNAME_ALLFUN = [PARNAME_LIST_RUN_PAR_CSV_FOLDER, PARNAME_LIST_RUN_PAR_CSV, PARNAME_LIST_SAVE_DEFAULT_PARSET, PARNAME_LIST_RUN_PLOT_1PAR,
                  PARNAME_LIST_RUN_PROCESSING, PARNAME_LIST_RUN_ANALYSIS_LOD]
PARNAME_LIST = list({PAR for FUN_LIST in PARNAME_ALLFUN for PAR in FUN_LIST})

# FCS channels
FCSCHANNEL1 = ['FSC-A', 'FSC-H', 'FSC-W', 'SSC-A', 'SSC-H', 'SSC-W', 'GFP-A', 'mCherry-A', 'Helix NP NIR-A', 'Cell trace violet-A', 'Time']
FCSCHANNEL2 = ['FSC-A', 'FSC-H', 'FSC-W', 'SSC-A', 'SSC-H', 'SSC-W', 'GFP-A', 'PerCP-Cy5-5-A', 'APC-A', 'BV 421-A', 'Time']
GATE_FRACTION = 0.95

# Gaussian kernel to use when smoothing.
SIGMA_GKS = 10
HIST1D_BINS = 256
SCATTER2D_BINS = 1024

# FCS ODE model max generation
FCS_maxgen = 11
STEPSIZE_ODE_FCS = 0.1

# default mutation rate
MUTATION_RATE = 7.1 * 10 ** (-7)

# Define the color used in print text.
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
