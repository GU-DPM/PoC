from DPM_lib import argparse
from DPM_run import DPM_run_par_csv_folder, DPM_run_par_csv, DPM_run_processing
from DPM_constant import *
'''This script defines the argument format for command-line execution.'''
'''Command-line useage example:
Runing the function "run_par_csv_folder" for the misspecification case:
python3 DPM_cl.py run_par_csv_folder --num_drug 2 --pathload './pnas/' --pathsave './div10_atsim/' --strategy_name CPM DPM2.2 
                                     --misspecification_fileload './div10/' --filename_pattern '_para.csv' 
                                     --save_filename_param True --save_filename_stopt True --save_filename_dosage True 
                                     --save_filename_pop True --save_filename_eachtimepoint False --misspecification_sim_only True 
                                     --misspecification_atsim True
Runing the function "run_par_csv_folder" for the misspecification case:
python3 DPM_cl.py run_par_csv_folder --num_drug 2 --pathload './pnas/' --pathsave './pnas_sim/' --strategy_name CPM DPM2.2
                                     --filename_pattern '_para.csv' --save_filename_param True --save_filename_stopt True 
                                     --save_filename_dosage True --save_filename_pop True --save_filename_eachtimepoint False
                                     --misspecification_sim_only False --misspecification_atsim True

Runing the function "run_processing":
python3 DPM_cl.py run_processing --num_drug 2 --pathload ./div5_atsim/ --strategy_name CPM DPM2.2                                     
python3 DPM_cl.py run_processing --num_drug 2 --pathload ./div10_atsim/ --strategy_name CPM DPM2.2                                     
python3 DPM_cl.py run_processing --num_drug 2 --pathload ./div30_atsim/ --strategy_name CPM DPM2.2                                     
python3 DPM_cl.py run_processing --num_drug 2 --pathload ./5x_atsim/ --strategy_name CPM DPM2.2                                     
python3 DPM_cl.py run_processing --num_drug 2 --pathload ./10x_atsim/ --strategy_name CPM DPM2.2                                     
python3 DPM_cl.py run_processing --num_drug 2 --pathload ./30x_atsim/ --strategy_name CPM DPM2.2                                     
python3 DPM_cl.py run_processing --num_drug 2 --pathload ./pnas_sim/ --strategy_name CPM DPM2.2                                     
'''
'''
CHECKED
'''


# Define an action class that prevents duplicate argument inputs.#
class DPMclnorepeatinput(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: list, option_string: str = None)\
            -> None:
        if namespace.__getattribute__(self.dest) is not None:
            parser.error(f'Repeated input of --{self.dest}.')
        else:
            namespace.__setattr__(self.dest, values)
        return


# Define an action class that checks whether --misspecification_filename_csv input has the same length as --filename_csv input.#
class DPMclmisfilenamecsv(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: list, option_string: str = None)\
            -> None:
        if namespace.__getattribute__(self.dest) is not None:
            parser.error(f'Repeated input of --{self.dest}.')
        elif namespace.filename_csv is None:
            parser.error(f'Please provide the --filename_csv argument before specifying --{self.dest}.')
        elif len(namespace.filename_csv) != len(values):
            parser.error(f'The number of .csv files provided to --{self.dest} must match the number supplied to --filename_csv.')
        else:
            namespace.__setattr__(self.dest, values)
        return


# Define an action class that checks whether the number of --dose_combination inputs is an integer multiple of the drug count.#
class DPMcldosecombinationaction(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: list, option_string: str = None)\
            -> None:
        if namespace.__getattribute__(self.dest) is not None:
            parser.error(f'Repeated input of --{self.dest}.')
        elif namespace.num_drug is None:
            parser.error(f'Please input the --num_drug argument first if you want to specify this option --{self.dest}.')
        elif len(values) % namespace.num_drug[0] != 0:
            parser.error(f'The number of --dose_combination inputs is not a integer multiple of the drug number. --dose_combination contains '
                         f'{len(values)} values, which are {values}. But this count is not an interger multiple of the specified drug number '
                         f'{namespace.num_drug[0]}.')
        else:
            namespace.__setattr__(self.dest, values)
        return


# Define an action class that checks whether the parameter filter criterisa are specific to the drug number:#
# Criterisa 1-6 apply to the 2 drug case, and criterisa 7, 8 apply to the 3 drug case.#
class DPMclparcriterisa(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: list, option_string: str = None)\
            -> None:
        if namespace.__getattribute__(self.dest) is not None:
            parser.error(f'Repeated input of --{self.dest}.')
        elif namespace.num_drug is None:
            parser.error(f'Please input --num_drug argument first if you want to input this option --{self.dest}.')
        elif any({7, 8}.intersection(set(values))) and namespace.num_drug[0] == 2:
            parser.error(f'Parameter criterisa { {7, 8}.intersection(set(values))} not used for 2 durg case.')
        elif any(set(range(1, 7)).intersection(set(values))) and namespace.num_drug[0] == 3:
            parser.error(f'Parameter criterisa {set(range(1, 7)).intersection(set(values))} not used for 3 durg case.')
        else:
            namespace.__setattr__(self.dest, list(set(values)))
        return


# Define an action class for --use_input_dose_combination_only that checks whether the --dose_combination argument has been provided.#
class DPMCLinputdosecombinationonly(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values, option_string: str = None) \
            -> None:
        if namespace.__getattribute__(self.dest) is not False:
            parser.error(f'Repeated input of --{self.dest}.')
        elif namespace.dose_combination is None:
            parser.error('Please input --dose_combination argument first if you want to set --use_input_dose_combination_only to True.')
        else:
            namespace.__setattr__(self.dest, True)
        return


# Define an action class for --misspecification_sim_only that checks whether the --misspecification_filename_csv argument has been provided.#
class DPMCLmisspecificationsimonly(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: str, option_string: str = None) \
            -> None:
        values = DPM_cl_boolen(self, parser, namespace, values)
        if namespace.__contains__('misspecification_filename_csv') and namespace.misspecification_filename_csv is None and values is True:
            parser.error(f'Please input --misspecification_filename_csv first if you want to set --{self.dest} to True.')
        namespace.__setattr__(self.dest, values)
        return


# Define an action class for a boolean flag that takes no value and defaults to True.#
class DPMCLboolendefaulttrue(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: str, option_string: str = None) \
            -> None:
        if namespace.__getattribute__(self.dest) is not True:
            parser.error(f'Repeated input of --{self.dest}.')
        else:
            namespace.__setattr__(self.dest, False)
        return


# Define an action class for a boolean flag that takes no value and defaults to False.#
class DPMCLboolendefaultfalse(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: str, option_string: str = None) \
            -> None:
        if namespace.__getattribute__(self.dest) is not False:
            parser.error(f'Repeated input of --{self.dest}.')
        else:
            namespace.__setattr__(self.dest, True)
        return


# Define an action class for a boolean flag that requires an explicit input value.#
class DPMCLboolen(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: str, option_string: str = None) \
            -> None:
        value = DPM_cl_boolen(self, parser, namespace, values)
        namespace.__setattr__(self.dest, value)
        return


def DPM_cl_boolen(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: str) -> bool:
    if namespace.__getattribute__(self.dest) is not None:
        parser.error(f'Repeated input of --{self.dest}.')
    values_input = str(values[0]).upper()
    if 'TRUE'.startswith(values_input):
        values = True
    elif 'FALSE'.startswith(values_input):
        values = False
    else:
        parser.error(f'Please input --{self.dest} as a "True" or "False" Boolen value.')
    return values


# Return function handle for an ArgumentParser type-checking function that validates a numerical range:#
# minimum <= arg <= maximum. Define the function with default arguments.#
def DPM_cl_numerical_range(mini, maxi, typearg):
    def DPM_cl_range_checker(arg):
        # New argparse type function for a numerical argument within a predefined range.#
        try:
            f = float(arg)
        except ValueError:
            raise argparse.ArgumentTypeError('Input must be numerical.')
        if f < mini or f > maxi:
            if maxi < 1e4:
                raise argparse.ArgumentTypeError(f'Input must be in range [{mini:g}, {maxi:g}].')
            else:
                raise argparse.ArgumentTypeError(f'Input must be in range [{mini:.2e}, {maxi:.2e}].')
        return float(arg) if typearg == 'float' else int(arg) if typearg == 'int' else None
    # Return the function handle of the checking function.#
    return DPM_cl_range_checker


def DPM_cl_main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description='DPM (Dynamic Precision Medicine) model.\n'  #
        'S: cells sensitive to all drugs.\n'  #
        'Rx cells resistant to drug x but sensitive to other drugs.\n'  #
        'Rxy: cells resistant to drug x and y but sensitive to other drugs.\n'  #
        'N: total cells.\n'  #
        'CPM: (Current Personalized Medicine):'  #
        'Initially treats the dominate population with the most potent drug. The current treatment is maintained until either:'  #
        '(i) the total population reaches twices the nadir population (a nadir is a local minimum of the total population under a fixed '  #
        'treatment), or'  #
        '(ii) the total population reemerges from below the radiologic detection threshold a level.\n'  #
        'DPM1: Minimizes the predicted total cell population.\n'  #
        'DPM2: Mimimizes the risk of developing incurable cells unless there is an immediate threat of mortality.'  #
        'DPM2.1: immediate threat of mortality occurs at total cell numbers bigger than 1e9. '   #
        'DPM2.2: immediate threat of mortality occurs at total cell numbers bigger than 1e11.\n'   #
        'DPM3: Minimize the predicted total population unless the model predicts that the first incurable cell will form within next time step '
        '(stepsize).\n'  #
        'DPM4: Estimates the time to either incurability or death and responds to the most imminent threat, provided cure is still possible.\n'
        'For full details of these strategies, see'   #
        '(1) Impact of genetic dynamics and single-cell heterogeneity on development of nonstandard personalized medicine strategies for cancer.'
        'PNAS(2012). '   #
        '(2) Long range personalized cancer treatment strategies incorporating evolutionary dynamics. Biology Direct (2016).',   #
        formatter_class=argparse.RawTextHelpFormatter)
    subparsers: parser.add_subparsers = parser.add_subparsers(dest='function_name')
    subparsers.required = True

    # Specify the arguments used in the functions.#
    Run_par_csv_folder_parser: argparse.ArgumentParser = subparsers.add_parser('run_par_csv_folder')
    Run_par_csv_parser: argparse.ArgumentParser = subparsers.add_parser('run_par_csv')
    Run_processing_parser: argparse.ArgumentParser = subparsers.add_parser('run_processing')

    # Parameter 'num_drug'.#
    for parser_i in [Run_par_csv_folder_parser,
                     Run_par_csv_parser,
                     Run_processing_parser]:
        parser_i.add_argument('--num_drug', required=True, type=int, nargs=1, choices={2, 3}, action=DPMclnorepeatinput,
                              help='Number of drugs. Allowed values: 2 or 3.')

    # Parameter 'dose_method'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--dose_method', type=str, nargs=1, choices={'d', 'c'}, action=DPMclnorepeatinput,
                              help='Dose method. Allowed values: d (discrete) or c (continous). Default: discrete.')

    # Parameter 'par_ind'.#
    for parser_i in [Run_par_csv_parser]:
        parser_i.add_argument('--par_ind', type=int, nargs='*', action=DPMclnorepeatinput,
                              help='Integer indices of parameters. Default = None')

    # Parameter 'dose_interval'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--dose_interval', type=float, nargs=1, choices={0.01, 0.02, 0.05, 0.1, 0.2, 0.5}, action=DPMclnorepeatinput,
                              help='Dose interval for continous dosing. If --dose_method is set to "c" (continous), the d single durg is '   #
                                   'generated as 0:dose_interval:1. If --dose_method is set to "d" (discrete) is ignored. this argument. '   #
                                   'Default=0.1.')

    # Parameter 'use_input_dose_combination_only'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--use_input_dose_combination_only', action=DPMCLinputdosecombinationonly, default=False, nargs=0,
                              help='If this flag is set, only the input dose combinations will be used. Default = False.')

    # Parameter 'dose_combination'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--dose_combination', nargs='*', type=float, action=DPMcldosecombinationaction,
                              help='Dose combinations. The number of inputs must be an integer multiple of the drug number specified by '  #
                                   '--num_drug. For example, if the drug number is 2, then the number of --dose_combination inputs '  #
                                   'can be 1/2, 1/2; 1, 0; 2, 1. In each consecutive pair, the first value is the dose of drug 1 and the second '
                                   'is the dose of drug 2. Valid dose range: min = 0, max = 1. '  #
                                   'Default behavior:\n'  #
                                   '- For discrete dosing:\n'  #
                                   "In the 2 drug case, a single drug's default doses are: 0, 0.5, 1. \n"   #
                                   "In the 3 drug case, a single drug's default doses are: 0, 1/3, 1/2, 1. \n"  #
                                   'For continous dosing: \n'   #
                                   "Drug doses are generated as 0:dose_interval:1. "  #
                                   'The default dose combinations are all possiable drug dose combinations whose total dose equals 1.')

    # Parameter 'Simduration'.#
    for parser_i in [Run_par_csv_folder_parser,
                     Run_par_csv_parser,
                     Run_processing_parser]:
        parser_i.add_argument('--simduration', type=DPM_cl_numerical_range(SIMDURATION_LOWERTHRES, SIMDURATION_UPPERTHRES, 'int'), nargs=1,
                              action=DPMclnorepeatinput, help=f'Set the simulation duration (days). Value range: min: {SIMDURATION_LOWERTHRES},'
                                                              f'max: {SIMDURATION_UPPERTHRES}. Default = {SIMDURATION_DEFAULT_VAL}.')

    # Parameter 'Limit_mortality'. #
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--limit_mortality', type=DPM_cl_numerical_range(LIMIT_MORTALITY_LOWERTHRES, LIMIT_MORTALITY_UPPERTHRES, 'float'),
                              nargs=1, action=DPMclnorepeatinput,
                              help=f'Set the threshold of cell number that cause mortality. Value range: min: '
                                   f'{LIMIT_MORTALITY_LOWERTHRES:.2e}, max: {LIMIT_MORTALITY_UPPERTHRES:.2e}. '
                                   f'Default = {LIMIT_MORTALITY_DEFAULT_VAL:.2e}.')

    # Parameter 'step_size'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser, Run_processing_parser]:
        parser_i.add_argument('--stepsize', type=DPM_cl_numerical_range(STEPSIZE_LOWERTHRES, STEPSIZE_UPPERTHRES, 'float'), nargs=1,
                              action=DPMclnorepeatinput, help=f'Set the treatments adjusted duration (days). Value range: '
                                                              f'min: {STEPSIZE_LOWERTHRES}, max: {STEPSIZE_UPPERTHRES}. '
                                                              f'Default = {STEPSIZE_DEFAULT_VAL}.')

    # Parameter 'Limit_radiologicdetection'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--limit_radiologicdetection',
                              type=DPM_cl_numerical_range(LIMIT_RADIOLOGICDETECTION_LOWERTHRES, LIMIT_MORTALITY_UPPERTHRES, 'float'), nargs=1,
                              action=DPMclnorepeatinput,
                              help=f'Set the threshold of cell number that can be detected by radiology. '
                                   f' Value range: min: {LIMIT_RADIOLOGICDETECTION_LOWERTHRES:.2e}, max: '
                                   f'{LIMIT_RADIOLOGICDETECTION_UPPERTHRES:.2e}. Default = {LIMIT_RADIOLOGICDETECTION_DEFAULT_VAL:.2e}.')

    # Parmeter 'Limit_moleculardetection'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--limit_moleculardetection',
                              type=DPM_cl_numerical_range(LIMIT_MOLECULARDETECTION_LOWERTHRES, LIMIT_MOLECULARDETECTION_UPPERTHRES, 'float'),
                              nargs=1, action=DPMclnorepeatinput,
                              help=f'Set the ratio threshold of molecular cell types that can be detected. '  #
                                   f' Value range: min: {LIMIT_MOLECULARDETECTION_LOWERTHRES:.2e}, max: '
                                   f'{LIMIT_MOLECULARDETECTION_UPPERTHRES:.2e}. Default = {LIMIT_MOLECULARDETECTION_DEFAULT_VAL:.2e}.')

    # Parameter 'run_sim'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--run_sim', action=DPMCLboolen, nargs=1, help='Do simulation if True otherwise not. Default = True.')

    # Parameter 'misspecification_sim_only'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--misspecification_sim_only', action=DPMCLmisspecificationsimonly, nargs=1,
                              help='Only simulate the mis-specification parameters if True otherwise not. Default = False.')

    # Parameter number of step difference for processing.#
    for parser_i in [Run_processing_parser]:
        parser_i.add_argument('--num_stepdiff', type=DPM_cl_numerical_range(NUM_STEPDIFF_LOWERTHRES, NUM_STEPDIFF_UPPERTHRES, 'int'), nargs=1,
                              action=DPMclnorepeatinput, help=f'Set the number of step difference for processing. Value range: '
                                                              f'min: {NUM_STEPDIFF_LOWERTHRES}, max: {NUM_STEPDIFF_UPPERTHRES}. '
                                                              f'Default = {NUM_STEPDIFF_DEFAULT_VAL}.')

    # Parameter PAR_criterisa.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--par_criterisa', type=int, nargs='*', choices=set(range(1, 9)), action=DPMclparcriterisa,
                              help='Parameter criterisa for 2 drug case:\n'   #
                                   '0: Drug potency on S cells must satisfying: drug 1 >= drug2 (i.e., Sa_ratio_S_D2toS_D1 <= 1).\n'
                                   '1: Sensitive cells must be present in the initial cell population (R1(t=0)+R2(t=0)< N(t=0)).\n'
                                   '2: Drug 1 and Drug 2 kill their corresponding resistant populations more slowly than those populations grow '
                                   '(Sa(R1, D1) <= g0_R1 and Sa(R2, D2) <= g0_R2). \n'
                                   '3: Drug 1 and Drug 2 affect sensitive cells (Sa(S,D1) > 0 and Sa(S,D2) > 0).\n'
                                   '4: The S cell population can be eradicated by each drug (Sa(S, D1) > g0_S and Sa(S, D2) > g0_S).\n'
                                   '5: The multi-resistant population cannot be eradicated by either drug (Sa(R12, D1) < g0,Sa(R12, D2) < g0).\n'
                                   '6: Cure is achievable only under simultaneous full-dose administration of both drugs (which is toxic and '
                                   'thus not a valid option), while no valid static treatment can achieve cure.\n'
                                   'Parameter criteria for 3 drug case:\n'
                                   '7: Drug potency on S cells must follow the order drug 1 >= drug 2 >= drug 3.\n'
                                   '8: The initial size of a doubly resistant subpopulation must not exceed the initial sizes of either singly '
                                   'resistant subpopulation from which it arises. Default = no criterisa applied.')

    # Parameter 'Strategy_name'.#
    for parser_i in [Run_par_csv_folder_parser,
                     Run_par_csv_parser,
                     Run_processing_parser]:
        parser_i.add_argument('--strategy_name', type=str, nargs='*',
                              choices={'CPM', 'DPM1', 'DPM2.1', 'DPM2.2', 'DPM3', 'DPM4', 'all'}, action=DPMclnorepeatinput,
                              help='Select the strategy to be used. Default={CPM, DPM2.2}')

    # Parameter 'save_filename_param'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--save_filename_param', action=DPMCLboolen, nargs=1,
                              help='Save the parameter configuration to a .csv file if True, otherwise do not. Default = True.')

    # Parameter 'save_filename_stopt'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--save_filename_stopt', action=DPMCLboolen, nargs=1,
                              help='Save the survival times to a .csv file if True, otherwise do not. Default = True.')

    # Parameter 'save_filename_pop'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--save_filename_pop', action=DPMCLboolen, nargs=1,
                              help='Save the population composition dynamics at each step size to a .csv file if True, otherwise do not. '
                                   'Default = False.')

    # Parameter 'save_filename_dosage'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--save_filename_dosage', action=DPMCLboolen, nargs=1,
                              help='Save the dosage combination sequence to a .csv file if True, otherwise do not. Default = True.')

    # Parameter 'save_filename_eachtimepoint'.#
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--save_filename_eachtimepoint', action=DPMCLboolen, nargs=1,
                              help='Save the poulation composition dynamics at each time point to a .csv file if True, otherwise do not.'
                              'Default = False.')

    # Parameter 'pathsave'.#
    for i, parser_i in enumerate([Run_par_csv_folder_parser,
                                  Run_par_csv_parser,
                                  Run_processing_parser]):
        defaultpath = ''
        parser_i.add_argument('--pathsave', type=str, nargs=1, action=DPMclnorepeatinput,
                              help=f'Set the directory where the simulation results will be saved. Default = {defaultpath}.')

    # Parameter 'pathload'.#
    for parser_i in [Run_par_csv_folder_parser, Run_processing_parser]:
        parser_i.add_argument('--pathload', type=str, nargs=1, action=DPMclnorepeatinput,
                              help='Set the directory from which the data will be loaded. Default = ./pnas_result.')

    # Parameter 'misspecification_fileload'.#
    for parser_i in [Run_par_csv_folder_parser, Run_processing_parser]:
        parser_i.add_argument('--misspecification_fileload', type=str, nargs=1, action=DPMclnorepeatinput,
                              help='Set the loading directory of the misspecified parameters. Default = ''.')

    # Parameter 'filename_pattern'.#
    for parser_i in [Run_par_csv_folder_parser]:
        parser_i.add_argument('--filename_pattern', type=str, nargs=1, action=DPMclnorepeatinput,
                              help='Set the filename pattern for the input filename. Default = _para.csv.')

    # Parameter 'misspecification_atdecision'.#
    for parser_i in [Run_par_csv_folder_parser]:
        parser_i.add_argument('--misspecification_atdecision', action=DPMCLboolen, nargs=1,
                              help="Misspecification simulation in which the oncologist's decisions use parameter values from the pnas paper "
                                   "(not used). Default = False.")

    # Parameter 'misspecification_atsim'.#
    for parser_i in [Run_par_csv_folder_parser]:
        parser_i.add_argument('--misspecification_atsim', action=DPMCLboolen, nargs=1,
                              help="Misspecification simulation in which the simulation use parameter values from the pnas paper "
                                   "(not used). Default = False.")

    # Parameter 'use_parallel'.#
    for parser_i in [Run_par_csv_folder_parser, Run_processing_parser]:
        parser_i.add_argument('--use_parallel', action=DPMCLboolen, nargs=1, help='Use parallel if True otherwise do not. Default = False.')

    # Parameter 'filename_csv'.#
    Run_par_csv_parser.add_argument('--filename_csv', required=True, type=str, nargs='*', action=DPMclnorepeatinput,
                                    help='The .csv filename includes the parameters.')

    # Parameter 'misspecification_filename_csv'.#
    Run_par_csv_parser.add_argument('--misspecification_filename_csv', type=str, nargs='*', action=DPMclmisfilenamecsv,
                                    help='The .csv filename contains the mis-specified parameters. The number of parameters in this file must '
                                         'match exactly the number of parameters in the file specified by argument --filename.csv.')

    args: parser.parse_args = parser.parse_args()
    function_name = args.function_name
    # Ensure the input is unique and use the default value if no input is provided.#

    # Parameter 'Num_drug'.#
    if args.__contains__('num_drug') and args.num_drug is not None:
        Num_drug = args.num_drug[0]
    else:
        Num_drug = None

    # Parameter 'dose_method'.#
    if args.__contains__('dose_method') and args.dose_method is not None:
        dose_method = args.dose_method[0]
    else:
        dose_method = 'd'

    # Parameter 'dose_interval'.#
    if args.__contains__('dose_interval') and args.dose_interval is not None:
        dose_interval = args.dose_interval[0]
    else:
        dose_interval = 0.1

    # Parameter 'dose_combination'.#
    if args.__contains__('dose_combination') and args.dose_combination is not None:
        dose_combination = [tuple(args.dose_combination[i:i + Num_drug]) for i in range(0, len(args.dose_combination), Num_drug)]
    else:
        dose_combination = None

    # Parameter 'use_input_dose_combination_only'.#
    if args.__contains__('use_input_dose_combination_only'):
        use_input_dose_combination_only = args.use_input_dose_combination_only
    else:
        use_input_dose_combination_only = False

    # Parameter 'Simduration'.#
    if args.__contains__('simduration') and args.simduration is not None:
        Simduration = args.simduration[0]
    else:
        Simduration = SIMDURATION_DEFAULT_VAL

    # Parameter 'Limit_mortality'.#
    if args.__contains__('limit_mortality') and args.limit_mortality is not None:
        Limit_mortality = args.limit_mortality[0]
    else:
        Limit_mortality = LIMIT_MORTALITY_DEFAULT_VAL

    # Parameter 'Stepsize'.#
    if args.__contains__('stepszie') and args.stepsize is not None:
        Stepsize = args.stepsize[0]
    else:
        Stepsize = STEPSIZE_DEFAULT_VAL

    # Parameter 'Limit_radiologicdetection'.#
    if args.__contains__('limit_radiologicdetection') and args.limit_radiologicdetection is not None:
        Limit_radiologicdetection = args.limit_radiologicdetection[0]
    else:
        Limit_radiologicdetection = LIMIT_RADIOLOGICDETECTION_DEFAULT_VAL

    # Parameter 'Limit_moleculardetection'.#
    if args.__contains__('limit_moleculardetection') and args.limit_moleculardetection is not None:
        Limit_moleculardetection = args.limit_moleculardetection[0]
    else:
        Limit_moleculardetection = LIMIT_MOLECULARDETECTION_DEFAULT_VAL

    # Parameter 'run_sim'.#
    if args.__contains__('run_sim') and args.run_sim is not None:
        run_sim = args.run_sim
    else:
        run_sim = True

    # Parameter 'misspecification_sim_only'.#
    if args.__contains__('misspecification_sim_only') and args.misspecification_sim_only is not None:
        misspecification_sim_only = args.misspecification_sim_only
    else:
        misspecification_sim_only = False

    # Parameter specifying the number of step differences used for processing.#
    if args.__contains__('num_stepdiff') and args.num_stepdiff is not None:
        Num_stepdiff = args.num_stepdiff[0]
    else:
        Num_stepdiff = NUM_STEPDIFF_DEFAULT_VAL

    # Parameter 'PAR_criterisa'.#
    if args.__contains__('par_criterisa') and args.par_criterisa is not None:
        PAR_criterisa = args.par_criterisa
    else:
        PAR_criterisa = []

    # Parameter 'par_ind'.#
    if args.__contains__('par_ind') and args.par_ind is not None:
        par_ind = args.par_ind
    else:
        par_ind = None

    # Parameter 'Strategy_name'.#
    if args.__contains__('strategy_name') and args.strategy_name is not None:
        Strategy_name = args.strategy_name
    else:
        Strategy_name = ['CPM', 'DPM2.2']

    # Parameter 'save_filename_param'.#
    if args.__contains__('save_filename_param') and args.save_filename_param is not None:
        save_filename_param = args.save_filename_param
    else:
        save_filename_param = True

    # Parameter 'save_filename_stopt'.#
    if args.__contains__('save_filename_stopt') and args.save_filename_stopt is not None:
        save_filename_stopt = args.save_filename_stopt
    else:
        save_filename_stopt = True

    # Parameter 'save_filename_pop'.#
    if args.__contains__('save_filename_pop') and args.save_filename_pop is not None:
        save_filename_pop = args.save_filename_pop
    else:
        save_filename_pop = False

    # Parameter 'save_filename_dosage'.#
    if args.__contains__('save_filename_dosage') and args.save_filename_dosage is not None:
        save_filename_dosage = args.save_filename_dosage
    else:
        save_filename_dosage = True

    # Parameter 'save_filename_eachtimepoint'.#
    if args.__contains__('save_filename_eachtimepoint') and args.save_filename_eachtimepoint is not None:
        save_filename_eachtimepoint = args.save_filename_eachtimepoint
    else:
        save_filename_eachtimepoint = False

    # Parameter 'pathsave'.#
    if args.__contains__('pathsave') and args.pathsave is not None:
        pathsave = args.pathsave[0]
    else:
        pathsave = ''

    # Parameter 'pathload'.#
    if args.__contains__('pathload') and args.pathload is not None:
        pathload = args.pathload[0]
    else:
        pathload = ''

    # Parameter 'misspecification_fileload'.#
    if args.__contains__('misspecification_fileload') and args.misspecification_fileload is not None:
        misspecification_fileload = args.misspecification_fileload[0]
    else:
        misspecification_fileload = ''

    # Parameter 'filename_pattern'.#
    if args.__contains__('filename_pattern') and args.filename_pattern is not None:
        filename_pattern = args.filename_pattern[0]
    else:
        filename_pattern = ''

    # Parameter 'misspecification_atdecision'.#
    if args.__contains__('misspecification_atdecision') and args.misspecification_atdecision is not None:
        misspecification_atdecision = args.misspecification_atdecision
    else:
        misspecification_atdecision = False

    # Parameter 'misspecification_atsim'.#
    if args.__contains__('misspecification_atsim') and args.misspecification_atsim is not None:
        misspecification_atsim = args.misspecification_atsim
    else:
        misspecification_atsim = False

    # Parameter 'use_parallel'.#
    if args.__contains__('use_parallel') and args.use_parallel is not None:
        use_parallel = args.use_parallel
    else:
        use_parallel = False

    # Parameter 'filename_csv'.#
    if args.__contains__('filename_csv') and args.filename_csv is not None:
        used = set()
        filename_csv = args.filename_csv[0] if len(args.filename_csv) == 1 else \
            [x for x in args.filename_csv if x not in used and (used.add(x) or True)]
    else:
        filename_csv = None

    # Parameter 'misspecification_filename_csv'.#
    if args.__contains__('misspecification_filename_csv') and args.misspecification_filename_csv is not None:
        used = set()
        misspecification_filename_csv = args.misspecification_filename_csv[0] \
            if len(args.misspecification_filename_csv) == 1 else \
            [x for x in args.misspecification_filename_csv if x not in used and (used.add(x) or True)]
    else:
        misspecification_filename_csv = None

    # Run the function corresponding to the specified function name.#
    if function_name.lower() == 'run_par_csv_folder'.lower():
        parval = [pathload, pathsave, use_parallel, filename_pattern, misspecification_atdecision, misspecification_atsim, Num_drug,
                  Limit_moleculardetection, Limit_radiologicdetection, Limit_mortality, Simduration, Stepsize, dose_method, dose_interval,
                  dose_combination, use_input_dose_combination_only, PAR_criterisa, Strategy_name, run_sim, save_filename_param,
                  save_filename_stopt, save_filename_pop, save_filename_dosage, save_filename_eachtimepoint, misspecification_sim_only,
                  misspecification_fileload]

        parinput = dict(zip(PARNAME_LIST_RUN_PAR_CSV_FOLDER, parval)) if len(parval) == len(PARNAME_LIST_RUN_PAR_CSV_FOLDER) \
            else None
        DPM_run_par_csv_folder(**parinput)

    if function_name.lower() == 'run_par_csv'.lower():
        parval = [filename_csv, pathsave, misspecification_filename_csv, misspecification_atdecision, misspecification_atsim,
                  Num_drug, Limit_moleculardetection, Limit_radiologicdetection, Limit_mortality, Simduration, Stepsize, dose_method,
                  dose_interval, dose_combination, use_input_dose_combination_only, PAR_criterisa, Strategy_name, run_sim,
                  save_filename_param, save_filename_stopt, save_filename_pop, save_filename_dosage, save_filename_eachtimepoint,
                  misspecification_sim_only, par_ind]

        parinput = dict(zip(PARNAME_LIST_RUN_PAR_CSV, parval)) if len(parval) == len(PARNAME_LIST_RUN_PAR_CSV) else None
        DPM_run_par_csv(**parinput)

    if function_name.lower() == 'run_processing'.lower():
        parval = [Num_drug, Strategy_name, pathsave, pathload, Num_stepdiff, Simduration, Stepsize, misspecification_fileload, use_parallel]

        parinput = dict(zip(PARNAME_LIST_RUN_PROCESSING, parval)) if len(parval) == len(PARNAME_LIST_RUN_PROCESSING) else None
        DPM_run_processing(**parinput)
    return


if __name__ == '__main__':
    DPM_cl_main()
