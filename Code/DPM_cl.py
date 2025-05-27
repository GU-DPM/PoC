from DPM_lib import argparse
from DPM_run import DPM_run_par_csv_folder, DPM_run_par_csv, DPM_save_default_PARset, DPM_run_processing, DPM_run_analysis_LOD
from DPM_constant import *

""" This script sets arguments format for command line execution. """


# Define an action class that don't allow repeat argument input.
class DPMclnorepeatinput(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: list, option_string: str = None) \
            -> None:
        if namespace.__getattribute__(self.dest) is not None:
            parser.error(f'Repeated input of --{self.dest}.')
        else:
            namespace.__setattr__(self.dest, values)
        return


# Define an action class which check that --misspecification_filename_csv input should has the same length as --filename_csv input.
class DPMclmisfilenamecsv(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: list, option_string: str = None) \
            -> None:
        if namespace.__getattribute__(self.dest) is not None:
            parser.error(f'Repeated input of --{self.dest}.')
        elif namespace.filename_csv is None:
            parser.error(f'Please input --filename_csv argument first if you want to input --{self.dest}.')
        elif len(namespace.filename_csv) != len(values):
            parser.error(f'The number of .csv files input from --{self.dest} should be the same as the number of files input from '
                         f'--filename_csv.')
        else:
            namespace.__setattr__(self.dest, values)
        return


# Define an action class which check that the number of --dose_combination input should be a integer multiple of the drug number.
class DPMcldosecombinationaction(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: list, option_string: str = None) \
            -> None:
        if namespace.__getattribute__(self.dest) is not None:
            parser.error(f'Repeated input of --{self.dest}.')
        elif namespace.num_drug is None:
            parser.error(f'Please input --num_drug argument first if you want to input --{self.dest}.')
        elif len(values) % namespace.num_drug[0] != 0:
            parser.error(f'The number of --dose_combination input is not a integer multiple of the drug number. The --dose_combination has '
                         f'{len(values)} values, which are {values}. But its number is not a interger multiple of drug number which is set to '
                         f'{namespace.num_drug[0]}.')
        else:
            namespace.__setattr__(self.dest, values)
        return


# Define an action class which check that parameter filter criterisa is drug number specific. Criterisa 1-6 are used in 2 drug case and criterisa 7, 8
# are used in 3 drug case.
class DPMclparcriterisa(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: list, option_string: str = None) \
            -> None:
        if namespace.__getattribute__(self.dest) is not None:
            parser.error(f'Repeated input of --{self.dest}.')
        elif namespace.num_drug is None:
            parser.error(f'Please input --num_drug argument first if you want to input --{self.dest}.')
        elif any({7, 8}.intersection(set(values))) and namespace.num_drug[0] == 2:
            parser.error(f'Parameter criterisa { {7, 8}.intersection(set(values))} not used for 2 durg case.')
        elif any(set(range(1, 7)).intersection(set(values))) and namespace.num_drug[0] == 3:
            parser.error(f'Parameter criterisa {set(range(1, 7)).intersection(set(values))} not used for 3 durg case.')
        else:
            namespace.__setattr__(self.dest, list(set(values)))
        return


# Define an action class for --use_input_dose_combination_only parameter checking the argument --dose_combination is given.
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


# Define an action class for --misspecification_sim_only parameter checking the argument --misspecification_filename_csv is given.
class DPMCLmisspecificationsimonly(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: str, option_string: str = None) \
            -> None:
        values = DPM_cl_boolen(self, parser, namespace, values)
        if namespace.subcolone_LOD is None:
            if namespace.__contains__('misspecification_filename_csv') and namespace.misspecification_filename_csv is None and values is True:
                parser.error(f'Please input --misspecification_filename_csv first if you want to set --{self.dest} to True.')
        namespace.__setattr__(self.dest, values)
        return


# Define an action class for boolean input without input, default True.
class DPMCLboolendefaulttrue(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: str, option_string: str = None) \
            -> None:
        if namespace.__getattribute__(self.dest) is not True:
            parser.error(f'Repeated input of --{self.dest}.')
        else:
            namespace.__setattr__(self.dest, False)
        return


# Define an action class for boolean input without input, default False.
class DPMCLboolendefaultfalse(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: str, option_string: str = None) \
            -> None:
        if namespace.__getattribute__(self.dest) is not False:
            parser.error(f'Repeated input of --{self.dest}.')
        else:
            namespace.__setattr__(self.dest, True)
        return


# Define an action class for boolean input which needs one input.
class DPMCLboolen(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: str, option_string: str = None) \
            -> None:
        value = DPM_cl_boolen(self, parser, namespace, values)
        namespace.__setattr__(self.dest, value)
        return


class DPMclmaxnumsubseq(argparse.Action):
    def __call__(self: argparse.Action, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: list, option_string: str = None) \
            -> None:
        if namespace.__getattribute__(self.dest) is not None:
            parser.error(f'Repeated input of --{self.dest}.')
        elif namespace.num_drug is None:
            parser.error(f'Please input --num_drug argument first if you want to input --{self.dest}.')
        elif namespace.num_drug[0] == 2 and (values[0] < MAXNUM_SUBSEQ_LOWERTHRES_2DRUG or values[0] > MAXNUM_SUBSEQ_UPPERTHRES_2DRUG):
            parser.error(f'Input must be in range [{MAXNUM_SUBSEQ_LOWERTHRES_2DRUG:g}, {MAXNUM_SUBSEQ_UPPERTHRES_2DRUG:g}].')
        elif namespace.num_drug[0] == 3 and (values[0] < MAXNUM_SUBSEQ_LOWERTHRES_3DRUG or values[0] > MAXNUM_SUBSEQ_UPPERTHRES_3DRUG):
            parser.error(f'Input must be in range [{MAXNUM_SUBSEQ_LOWERTHRES_3DRUG:g}, {MAXNUM_SUBSEQ_UPPERTHRES_3DRUG:g}].')
        else:
            namespace.__setattr__(self.dest, values)
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


# Return function handle of an argument type function for ArgumentParser checking a numerical number range: minimum <= arg <= maximum.
# Define the function with default arguments.
def DPM_cl_numerical_range(mini, maxi, typearg):
    def DPM_cl_range_checker(arg):
        # New type function for argparse - a numerical arg within predefined range.
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
    # Return function handle of the checking function.
    return DPM_cl_range_checker


def DPM_cl_main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description='DPM (Dynamic Precision Medicine).\n'
        'S: cells sensitive to all drugs.\n'
        'Rx cells resistant to drug x but sensitive to other drugs.\n'
        'Rxy: cells resistant to drug x and y but sensitive to other drugs.\n'
        'Rxyz: cells resistant to drug x, y and z but sensitive to other drugs.\n'
        'N: total cells.\n'
        'Strategy 0: current personalized medicine strategy. Initially treat the dominate population using the most potent drug. Maintain the '
        'current treatment until either (i): the total population reaches twices the nadir population (A nadir is a local minimum of total '
        'population among the time-series profile where the current treatment is maintained.). (ii): the total population reemerges from a level '
        'below the radiologic detection threshold.\n'
        'Strategy 1: minimize the predicted total population.\n'
        'Strategy 2: mimimize the risk of incurable cells developing unless there is an immediate threat of mortality. Strategy 2.1: mortality '
        'happens at total cell numbers bigger than 1e9. Strategy 2.2: mortality happens at total cell numbers bigger than 1e11.\n'
        'Strategy 3: minimize the predicted total population unless there is a prediction that the first incurable cell will form within next '
        'stepsize.\n'
        'Strategy 4: estimate the time to either incurability or death, and react to the most proximal threat as long as there is a chance of cure.\n'
        'Strategy 5: multistep extension of strategy 1.\n'
        'Strategy 6: multistep extension of strategy 2.1.\n'
        'Strategy 7: multistep extension of strategy 2.2.\n'
        'Strategy 8: multistep extension of strategy 3.\n'
        'Strategy 9: ALTO(Adaptive Long-Term Optimization).\n'
        'For detail of the strategyies, see refs: (1) Impact of genetic dynamics and single-cell heterogeneity on development of nonstandard '
        'personalized medicine strategies for cancer. PNAS(2012). (2) Long range personalized cancer treatment strategies incorporating evolutionary '
        'dynamics. Biology Direct (2016).', formatter_class=argparse.RawTextHelpFormatter)
    subparsers: parser.add_subparsers = parser.add_subparsers(dest='function_name')
    subparsers.required = True

    # Specify the arguments used in "Run_par_csv" function.
    Run_par_csv_folder_parser: argparse.ArgumentParser = subparsers.add_parser('run_par_csv_folder')
    Run_par_csv_parser: argparse.ArgumentParser = subparsers.add_parser('run_par_csv')
    Save_default_PARset_parser: argparse.ArgumentParser = subparsers.add_parser('save_default_PARset')
    Run_processing_parser: argparse.ArgumentParser = subparsers.add_parser('run_processing')
    Run_analysis_LOD: argparse.ArgumentParser = subparsers.add_parser('run_analysis_LOD')

    # Parameter num_drug.
    for parser_i in [Run_par_csv_folder_parser,
                     Run_par_csv_parser,
                     Save_default_PARset_parser,
                     Run_processing_parser,
                     Run_analysis_LOD]:
        parser_i.add_argument('--num_drug', required=True, type=int, nargs=1, choices={2, 3}, action=DPMclnorepeatinput,
                              help='Number of drugs. Two drug number allowed: 2 or 3.')

    # Parameter dose_method.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--dose_method', type=str, nargs=1, choices={'d', 'c'}, action=DPMclnorepeatinput,
                              help='Dose method. Two dose methods allowed: d:discrete dose or c:continous dose. Default=discrete dose.')

    # Parameter par_ind.
    for parser_i in [Run_par_csv_parser]:
        parser_i.add_argument('--par_ind', type=int, nargs='*', action=DPMclnorepeatinput,
                              help='Integer index of parameters. E.g. 1 50. Default = None')

    # Parameter dose_interval.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--dose_interval', type=float, nargs=1, choices={0.01, 0.02, 0.05, 0.1, 0.2, 0.5}, action=DPMclnorepeatinput,
                              help='Dose interval for continous dose. If --dose_method is set to "c" which means continous dose, the dose for a '
                                   'single durg is 0:dose_interval:1. If --dose_method is set to "d" which means discrete dose, this argument is '
                                   'ignored. Default=0.1.')
    # Parameter use_input_dose_combination_only.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--use_input_dose_combination_only', action=DPMCLinputdosecombinationonly, default=False, nargs=0,
                              help='If flag this argument, will use the input dose combinations only. Default = False.')
    # Parameter dose_combination.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--dose_combination', nargs='*', type=float, action=DPMcldosecombinationaction,
                              help='Dose combinations. The number of input should be integer multiple of the drug number specified by --num_drug. '
                                   'E.g., if the drug number equals 2, the number of --dose_combination inputs can be 2, 4, 6, ..., 2*n. The first '
                                   'input in the two consecutive numbers is the dose of drug 1 and the second is the dose of drug 2. Value range: '
                                   'min: 0, max: 1. Default: if use discrete dose, the dose for a single drug will be 0, 0.5, 1 in 2 drug case and '
                                   '0, 1/3, 1/2, 1 in 3 drug case. If use continous dose, the dose for a single drug will be 0:dose_interval:1. '
                                   'The default dose combinations are all possiable drug dose combinations which gives total dose equal 1.')
    # Parameter Simduration.
    for parser_i in [Run_par_csv_folder_parser,
                     Run_par_csv_parser,
                     Save_default_PARset_parser,
                     Run_processing_parser,
                     Run_analysis_LOD]:
        parser_i.add_argument('--simduration', type=DPM_cl_numerical_range(SIMDURATION_LOWERTHRES, SIMDURATION_UPPERTHRES, 'int'), nargs=1,
                              action=DPMclnorepeatinput, help=f'Set the simulation duration (days). Value range: min: {SIMDURATION_LOWERTHRES}, max: '
                              f'{SIMDURATION_UPPERTHRES}. Default = {SIMDURATION_DEFAULT_VAL}.')
    # Parameter Limit_mortality.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser, Save_default_PARset_parser]:
        parser_i.add_argument('--limit_mortality', type=DPM_cl_numerical_range(LIMIT_MORTALITY_LOWERTHRES, LIMIT_MORTALITY_UPPERTHRES, 'float'),
                              nargs=1, action=DPMclnorepeatinput,
                              help=f'Set the threshold of cell number that cause mortality (#). Value range: min: {LIMIT_MORTALITY_LOWERTHRES:.2e}, '
                                   f'max: {LIMIT_MORTALITY_UPPERTHRES:.2e}. Default = {LIMIT_MORTALITY_DEFAULT_VAL:.2e}.')
    # Parameter step_size.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser, Run_processing_parser]:
        parser_i.add_argument('--stepsize', type=DPM_cl_numerical_range(STEPSIZE_LOWERTHRES, STEPSIZE_UPPERTHRES, 'float'), nargs=1,
                              action=DPMclnorepeatinput, help=f'Set the treatments adjusted duration (days). Value range: '
                                                              f'min: {STEPSIZE_LOWERTHRES}, max: {STEPSIZE_UPPERTHRES}. '
                                                              f'Default = {STEPSIZE_DEFAULT_VAL}.')
    # Parameter Limit_radiologicdetection.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--limit_radiologicdetection',
                              type=DPM_cl_numerical_range(LIMIT_RADIOLOGICDETECTION_LOWERTHRES, LIMIT_MORTALITY_UPPERTHRES, 'float'), nargs=1,
                              action=DPMclnorepeatinput,
                              help=f'Set the threshold of cell number that can be detected by radiology. '
                                   f' Value range: min: {LIMIT_RADIOLOGICDETECTION_LOWERTHRES:.2e}, max: '
                                   f'{LIMIT_RADIOLOGICDETECTION_UPPERTHRES:.2e}. Default = {LIMIT_RADIOLOGICDETECTION_DEFAULT_VAL:.2e}.')

    # Parmeter Limit_moleculardetection.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--limit_moleculardetection',
                              type=DPM_cl_numerical_range(LIMIT_MOLECULARDETECTION_LOWERTHRES, LIMIT_MOLECULARDETECTION_UPPERTHRES, 'float'), nargs=1,
                              action=DPMclnorepeatinput,
                              help=f'Set the ratio threshold of molecular cell type that can be detected. '
                                   f' Value range: min: {LIMIT_MOLECULARDETECTION_LOWERTHRES:.2e}, max: '
                                   f'{LIMIT_MOLECULARDETECTION_UPPERTHRES:.2e}. Default = {LIMIT_MOLECULARDETECTION_DEFAULT_VAL:.2e}.')

    # Parameter lookahead_step.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--lookahead_step', type=DPM_cl_numerical_range(LOOKAHEAD_STEP_LOWERTHRES, LOOKAHEAD_STEP_UPPERTHRES, 'int'),
                              nargs=1, action=DPMclnorepeatinput,
                              help=f'Set the look-ahead steps in multistep extension strategies (#). Value range: min: {LOOKAHEAD_STEP_LOWERTHRES}, '
                                   f'max: {LOOKAHEAD_STEP_UPPERTHRES}. Default = {LOOKAHEAD_STEP_DEFAULT_VAL}.')
    # Parameter Maxnum_subseq.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--maxnum_subseq', type=int, nargs=1, action=DPMclmaxnumsubseq,
                              help=f'Set the maximum number of tracked subsequences in strategy 9 (#). Value range: 2 drug case, '
                                   f'min: {MAXNUM_SUBSEQ_LOWERTHRES_2DRUG}, max: {MAXNUM_SUBSEQ_UPPERTHRES_2DRUG}, '
                                   f'default = {MAXNUM_SUBSEQ_DEFAULT_VAL_2DRUG}; 3 drug case, min: {MAXNUM_SUBSEQ_LOWERTHRES_3DRUG}, '
                                   f'max: {MAXNUM_SUBSEQ_UPPERTHRES_3DRUG}, default = {MAXNUM_SUBSEQ_DEFAULT_VAL_3DRUG}.')
    # Parameter subtreedepth.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--subtreedepth', type=DPM_cl_numerical_range(SUBTREEDEPTH_LOWERTHRES, SUBTREEDEPTH_UPPERTHRES, 'int'), nargs=1,
                              action=DPMclnorepeatinput, help=f'Set the subtree depth of strategy 9 (#). Value range: '
                                                              f'min: {SUBTREEDEPTH_LOWERTHRES}, max: {SUBTREEDEPTH_UPPERTHRES}. '
                                                              f'Default = {SUBTREEDEPTH_DEFAULT_VAL}.')
    # Parameter run_sim.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--run_sim', action=DPMCLboolen, nargs=1, help='Do simulation if True otherwise not. Default = True.')

    # Parameter misspecification_sim_only.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--misspecification_sim_only', action=DPMCLmisspecificationsimonly, nargs=1,
                              help='Only simulate the mis-specification parameters if True otherwise not. Default = False.')

    # Parameter subcolone limit of detection.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--subcolone_LOD',
                              type=DPM_cl_numerical_range(SUBCOLONE_LOD_LOWERTHRES, SUBCOLONE_LOD_UPPERTHRES, 'float'), nargs=1,
                              action=DPMclnorepeatinput,
                              help=f'Set the ratio threshold of subcolone that can be detected. '
                                   f' Value range: min: {SUBCOLONE_LOD_LOWERTHRES:.2e}, max: '
                                   f'{SUBCOLONE_LOD_UPPERTHRES:.2e}. Default = no limitation.')

    # Parameter LOD.
    for parser_i in [Run_analysis_LOD]:
        parser_i.add_argument('--LOD', type=str, nargs='*',
                              choices={'1e-01', '1e-02', '1e-03', '1e-04', '1e-05', '1e-06'}, action=DPMclnorepeatinput,
                              help='Set LOD (limit of detection). Default={1e-01, 1e-02, 1e-03, 1e-04, 1e-05, 1e-06}')

    # Parameter mutation rate.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--mutation_rate',
                              type=DPM_cl_numerical_range(MUTATION_RATE_LOWERTHRES, MUTATION_RATE_UPPERTHRES, 'float'), nargs=1,
                              action=DPMclnorepeatinput,
                              help=f'Set the mutation rate. '
                                   f' Value range: min: {MUTATION_RATE_LOWERTHRES:.2e}, max: '
                                   f'{MUTATION_RATE_UPPERTHRES:.2e}. Default = {MUTATION_RATE_DEFAULT_VAL:.2e}.')

    # Parameter misspecification method for subcolone under limit of detection.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--misspecification_LOD', type=str, nargs=1, action=DPMclnorepeatinput,
                              help='misspecification method for subcolone under limit of detection. 0 set subcolone to 0 all others sample from'
                                   ' a probability density function (pdf). Default=sample from pdf.')

    # Parameter number of step difference for processing
    for parser_i in [Run_processing_parser]:
        parser_i.add_argument('--num_stepdiff', type=DPM_cl_numerical_range(NUM_STEPDIFF_LOWERTHRES, NUM_STEPDIFF_UPPERTHRES, 'int'), nargs=1,
                              action=DPMclnorepeatinput, help=f'Set the number of step difference for processing. Value range: '
                                                              f'min: {NUM_STEPDIFF_LOWERTHRES}, max: {NUM_STEPDIFF_UPPERTHRES}. '
                                                              f'Default = {NUM_STEPDIFF_DEFAULT_VAL}.')

    # Parameter erase_preresult.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--erase_preresult', action=DPMCLboolen, nargs=1,
                              help='Erase the exist simulation result in the save files if True otherwise not. Default = True.')
    # Parameter PAR_criterisa.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--par_criterisa', type=int, nargs='*', choices=set(range(1, 9)), action=DPMclparcriterisa,
                              help='Parameter criterisa for 2 drug case, 0: The potency of drugs on S cell follows: drug 1 >= drug2 '
                                   '(Sa_ratio_S_D2toS_D1 <= 1). 1: There are sensitive cells at the initial cell population (R1(t=0) + R2(t=0) '
                                   '< N(t=0)). 2: Drug 1 and Drug 2\' rate of killing cells that resistant to drug 1 and drug 2 repectively is '
                                   'slower than its growth rate (Sa(R1, D1) <= g0_R1 and Sa(R2, D2) <= g0_R2). 3: Drug 1 and Drug 2 both have '
                                   'effect on sensitive cells (Sa(S,D1) > 0 and Sa(S,D2) > 0). 4: S cell type can be eradicated by each drug '
                                   '(Sa(S, D1) > g0_S and Sa(S, D2) > g0_S). 5: Multiply-resistant cell types cannot be eradicated by any drug '
                                   '(Sa(R12, D1) < g0 and Sa(R12, D2) < g0). 6: The patient can be cured by simultaneous full dosages of all drugs '
                                   'but cannot be cured by any valid static treatment. Parameter criterisa for 3 drug case, 7: The potency of drugs '
                                   'on S cell follows drug 1 >= drug 2 >= drug 3. 8: The initial size of a doubly-resistant subpopulation is no '
                                   'greater than the initial sizes of the two singly-resistant subpopulations from which it is derived. '
                                   'Default = no criterisa applied.')
    # Parameter Strategy_name.
    for parser_i in [Run_par_csv_folder_parser,
                     Run_par_csv_parser,
                     Run_processing_parser,
                     Run_analysis_LOD]:
        parser_i.add_argument('--strategy_name', type=str, nargs='*',
                              choices={'strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
                                       'strategy7', 'strategy8', 'strategy9', 'all'}, action=DPMclnorepeatinput,
                              help='Select the strategy to be used. Default={strategy0, strategy2.2}')
    # Parameter full_input.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--fullinput', action=DPMCLboolen, nargs=1,
                              help='Need to specify every parameters of the model in the input .csv file if True otherwise not. ' 'Default = False.')
    # Parameter par_save_block_size.
    Save_default_PARset_parser.add_argument('--par_save_block_size',
                                            type=DPM_cl_numerical_range(PAR_SAVE_BLOCK_SIZE_LOWERTHRES, PAR_SAVE_BLOCK_SIZE_UPPERTHRES, 'float'),
                                            nargs=1, action=DPMclnorepeatinput,
                                            help=f'Set the parameter save block size (#). Value range: min: {PAR_SAVE_BLOCK_SIZE_LOWERTHRES:.2e}, '
                                                 f'max: {PAR_SAVE_BLOCK_SIZE_UPPERTHRES:.2e}, Default = {PAR_SAVE_BLOCK_SIZE_DEFAULT_VAL:.2e}.')
    # Parameter plot.
    for parser_i in [Run_analysis_LOD]:
        parser_i.add_argument('--plot', action=DPMCLboolen, nargs=1,
                              help='Plot the result. Default = True.')
    # Parameter save_filename_param.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--save_filename_param', action=DPMCLboolen, nargs=1,
                              help='Save the parameter configuration to .txt file if True otherwise not. Default = True.')
    # Parameter save_filename_stopt.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--save_filename_stopt', action=DPMCLboolen, nargs=1,
                              help='Save the survival times to .txt file if True otherwise not. Default = True.')
    # Parameter save_filename_pop.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--save_filename_pop', action=DPMCLboolen, nargs=1,
                              help='Save the population composition dynamics at each stepsize to .txt file if True otherwise not. Default = False.')
    # Parameter save_filename_dosage.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--save_filename_dosage', action=DPMCLboolen, nargs=1,
                              help='Save the dosage combination sequence to .txt file if True otherwise not. Default = True.')
    # Parameter save_filename_eachtimepoint.
    for parser_i in [Run_par_csv_folder_parser, Run_par_csv_parser]:
        parser_i.add_argument('--save_filename_eachtimepoint', action=DPMCLboolen, nargs=1,
                              help='Save the poulation composition dynamics at each time points to .txt file if True otherwise not.'
                              'Default = False.')
    # Parameter pathsave.
    for i, parser_i in enumerate([Run_par_csv_folder_parser,
                                  Run_par_csv_parser,
                                  Save_default_PARset_parser,
                                  Run_processing_parser,
                                  Run_analysis_LOD]):
        if i in [0, 1]:
            defaultpath = "./csv_result"
        elif i == 2:
            defaultpath = "./parset_result"
        elif i == 3:
            defaultpath = "./csv_result/"
        else:
            defaultpath = ''
        parser_i.add_argument('--pathsave', type=str, nargs=1, action=DPMclnorepeatinput,
                              help=f'Set the saving directory of the simulation result. Default = {defaultpath}.')
    # Parameter pathload.
    for parser_i in [Run_par_csv_folder_parser,
                     Run_processing_parser,
                     Run_analysis_LOD]:
        parser_i.add_argument('--pathload', type=str, nargs=1, action=DPMclnorepeatinput,
                              help='Set the loading directory of the data. Default = ./parset_data.')

    # Parameter pathloadmis.
    for parser_i in [Run_processing_parser,
                     Run_analysis_LOD]:
        parser_i.add_argument('--pathloadmis', type=str, nargs=1, action=DPMclnorepeatinput,
                              help='Set the loading directory of the data. Default = ./parset_data.')

    # Parameter filename_pattern
    for parser_i in [Run_par_csv_folder_parser,
                     Run_analysis_LOD]:
        parser_i.add_argument('--filename_pattern', type=str, nargs=1, action=DPMclnorepeatinput,
                              help="Set the filename pattern for the input filename. Default = ''.")

    # Parameter misspecification.
    for parser_i in [Run_par_csv_folder_parser]:
        parser_i.add_argument('--misspecification', action=DPMCLboolen, nargs=1, help='Run misspecification .csv file in the folder. '
                                                                                      'Default = False.')
    # Parameter use_parallel.
    for parser_i in [Run_par_csv_folder_parser, Run_processing_parser]:
        parser_i.add_argument('--use_parallel', action=DPMCLboolen, nargs=1, help='Use parallel if True otherwise not. Default = False.')

    # Parameter filename_csv.
    Run_par_csv_parser.add_argument('--filename_csv', required=True, type=str, nargs='*', action=DPMclnorepeatinput,
                                    help='The .csv file name contains the parameters. E.g. example_DPM_parameter_2drug.csv.')
    # Parameter mis_filename_csv.
    Run_par_csv_parser.add_argument('--misspecification_filename_csv', type=str, nargs='*', action=DPMclmisfilenamecsv,
                                    help='The .csv file name contains the mis-specified parameters. The number of parameters in this file should be '
                                         'exactly the same as the number of parameters in the file specified by argument --filename.csv. E.g. '
                                         'example_DPM_parameter_2drug_mis. Default = None.')

    args: parser.parse_args = parser.parse_args()
    function_name = args.function_name
    # Make input unique and set to default value if no input.

    # Parameter Num_drug.
    if args.__contains__('num_drug') and args.num_drug is not None:
        Num_drug = args.num_drug[0]
    else:
        Num_drug = None

    # Parameter dose_method.
    if args.__contains__('dose_method') and args.dose_method is not None:
        dose_method = args.dose_method[0]
    else:
        dose_method = 'd'

    # Parameter dose_interval.
    if args.__contains__('dose_interval') and args.dose_interval is not None:
        dose_interval = args.dose_interval[0]
    else:
        dose_interval = 0.1

    # Parameter dose_combination.
    if args.__contains__('dose_combination') and args.dose_combination is not None:
        dose_combination = [tuple(args.dose_combination[i:i + Num_drug]) for i in range(0, len(args.dose_combination), Num_drug)]
    else:
        dose_combination = None

    # Parameter use_input_dose_combination_only.
    if args.__contains__('use_input_dose_combination_only'):
        use_input_dose_combination_only = args.use_input_dose_combination_only
    else:
        use_input_dose_combination_only = False

    # Parameter Simduration.
    if args.__contains__('simduration') and args.simduration is not None:
        Simduration = args.simduration[0]
    else:
        Simduration = SIMDURATION_DEFAULT_VAL

    # Parameter Limit_mortality.
    if args.__contains__('limit_mortality') and args.limit_mortality is not None:
        Limit_mortality = args.limit_mortality[0]
    else:
        Limit_mortality = LIMIT_MORTALITY_DEFAULT_VAL

    # Parameter Stepsize.
    if args.__contains__('stepszie') and args.stepsize is not None:
        Stepsize = args.stepsize[0]
    else:
        Stepsize = STEPSIZE_DEFAULT_VAL

    # Parameter Limit_radiologicdetection.
    if args.__contains__('limit_radiologicdetection') and args.limit_radiologicdetection is not None:
        Limit_radiologicdetection = args.limit_radiologicdetection[0]
    else:
        Limit_radiologicdetection = LIMIT_RADIOLOGICDETECTION_DEFAULT_VAL

    # Parameter Limit_moleculardetection
    if args.__contains__('limit_moleculardetection') and args.limit_moleculardetection is not None:
        Limit_moleculardetection = args.limit_moleculardetection[0]
    else:
        Limit_moleculardetection = LIMIT_MOLECULARDETECTION_DEFAULT_VAL

    # Parameter lookahead_step.
    if args.__contains__('lookahead_step') and args.lookahead_step is not None:
        lookahead_step = args.lookahead_step[0]
    else:
        lookahead_step = LOOKAHEAD_STEP_DEFAULT_VAL

    # Parameter Maxnum_subseq.
    if args.__contains__('maxnum_subseq') and args.maxnum_subseq is not None:
        Maxnum_subseq = args.maxnum_subseq[0]
    elif Num_drug == 2:
        Maxnum_subseq = MAXNUM_SUBSEQ_DEFAULT_VAL_2DRUG
    elif Num_drug == 3:
        Maxnum_subseq = MAXNUM_SUBSEQ_DEFAULT_VAL_3DRUG
    else:
        Maxnum_subseq = None

    # Parameter subtreedepth.
    if args.__contains__('subtreedepth') and args.subtreedepth is not None:
        subtreedepth = args.subtreedepth[0]
    else:
        subtreedepth = SUBTREEDEPTH_DEFAULT_VAL

    # Parameter run_sim.
    if args.__contains__('run_sim') and args.run_sim is not None:
        run_sim = args.run_sim
    else:
        run_sim = True

    # Parameter misspecification_sim_only.
    if args.__contains__('misspecification_sim_only') and args.misspecification_sim_only is not None:
        misspecification_sim_only = args.misspecification_sim_only
    else:
        misspecification_sim_only = False

    # Parameter subcolone_LOD.
    if args.__contains__('subcolone_LOD') and args.subcolone_LOD is not None:
        subcolone_LOD = args.subcolone_LOD[0]
    else:
        subcolone_LOD = ''

    # Parameter LOD
    if args.__contains__('LOD') and args.LOD is not None:
        LOD = list(set(args.LOD))
    else:
        LOD = ['1e-06', '1e-05', '1e-04', '1e-03', '1e-02', '1e-01']

    # Parameter mutation rate.
    if args.__contains__('mutation_rate') and args.mutation_rate is not None:
        mutation_rate = args.mutation_rate[0]
    else:
        mutation_rate = MUTATION_RATE_DEFAULT_VAL

    # Parameter misspecification_LOD
    if args.__contains__('misspecification_LOD') and args.misspecification_LOD is not None:
        misspecification_LOD = args.misspecification_LOD[0]
        try:
            misspecification_LOD = float(misspecification_LOD)
        except (ValueError, TypeError):
            misspecification_LOD = misspecification_LOD
    else:
        misspecification_LOD = 'pdf'

    # Parameter number of step difference for processing
    if args.__contains__('num_stepdiff') and args.num_stepdiff is not None:
        Num_stepdiff = args.num_stepdiff[0]
    else:
        Num_stepdiff = NUM_STEPDIFF_DEFAULT_VAL

    # Parameter erase_preresult.
    if args.__contains__('erase_preresult') and args.erase_preresult is not None:
        erase_preresult = args.erase_preresult
    else:
        erase_preresult = True

    # Parameter PAR_criterisa.
    if args.__contains__('par_criterisa') and args.par_criterisa is not None:
        PAR_criterisa = args.par_criterisa
    else:
        PAR_criterisa = []

    # Parameter par_ind.
    if args.__contains__('par_ind') and args.par_ind is not None:
        par_ind = args.par_ind
    else:
        par_ind = None

    # Parameter Strategy_name.
    if args.__contains__('strategy_name') and args.strategy_name is not None:
        Strategy_name = list(set(args.strategy_name))
        if 'all' in Strategy_name:
            Strategy_name = STRATEGY_LIST
    else:
        Strategy_name = ['strategy0', 'strategy2.2']

    # Parameter fullinput.
    if args.__contains__('fullinput') and args.fullinput is not None:
        fullinput = args.fullinput
    else:
        fullinput = False

    # Parameter par_save_block_size.
    if args.__contains__('par_save_block_size') and args.par_save_block_size is not None:
        par_save_block_size = args.par_save_block_size[0]
    else:
        par_save_block_size = PAR_SAVE_BLOCK_SIZE_DEFAULT_VAL

    # Parameter plot.
    if args.__contains__('plot') and args.plot is not None:
        plot = args.plot
    else:
        plot = True

    # Parameter save_filename_param.
    if args.__contains__('save_filename_param') and args.save_filename_param is not None:
        save_filename_param = args.save_filename_param
    else:
        save_filename_param = True

    # Parameter save_filename_stopt.
    if args.__contains__('save_filename_stopt') and args.save_filename_stopt is not None:
        save_filename_stopt = args.save_filename_stopt
    else:
        save_filename_stopt = True

    # Parameter save_filename_pop.
    if args.__contains__('save_filename_pop') and args.save_filename_pop is not None:
        save_filename_pop = args.save_filename_pop
    else:
        save_filename_pop = False

    # Parameter save_filename_dosage.
    if args.__contains__('save_filename_dosage') and args.save_filename_dosage is not None:
        save_filename_dosage = args.save_filename_dosage
    else:
        save_filename_dosage = True

    # Parameter save_filename_eachtimepoint.
    if args.__contains__('save_filename_eachtimepoint') and args.save_filename_eachtimepoint is not None:
        save_filename_eachtimepoint = args.save_filename_eachtimepoint
    else:
        save_filename_eachtimepoint = False

    # Parameter pathsave.
    if args.__contains__('pathsave') and args.pathsave is not None:
        pathsave = args.pathsave[0]
    elif function_name == 'run_par_csv':
        pathsave = './csv_result'
    elif function_name == 'run_default_PARset':
        pathsave = './parset_result'
    elif function_name == 'save_default_PARset':
        pathsave = './parset_data'
    else:
        pathsave = None

    # Parameter pathload.
    if args.__contains__('pathload') and args.pathload is not None:
        pathload = args.pathload[0]
    elif Num_drug == 2:
        pathload = './parset_data/2drug'
    elif Num_drug == 3:
        pathload = './parset_data/3drug'
    else:
        pathload = ''

    # Parameter pathloadmis.
    if args.__contains__('pathloadmis') and args.pathloadmis is not None:
        pathloadmis = args.pathloadmis[0]
    else:
        pathloadmis = ''

    # Parameter filename_pattern.
    if args.__contains__('filename_pattern') and args.filename_pattern is not None:
        filename_pattern = args.filename_pattern[0]
    else:
        filename_pattern = ''

    # Parameter misspecification.
    if args.__contains__('misspecification') and args.misspecification is not None:
        misspecification = args.misspecification
    else:
        misspecification = False

    # Parameter use_parallel.
    if args.__contains__('use_parallel') and args.use_parallel is not None:
        use_parallel = args.use_parallel
    else:
        use_parallel = False

    # Parameter filename_csv.
    if args.__contains__('filename_csv') and args.filename_csv is not None:
        used = set()
        filename_csv = args.filename_csv[0] if len(args.filename_csv) == 1 else \
            [x for x in args.filename_csv if x not in used and (used.add(x) or True)]
    else:
        filename_csv = None

    # Parameter misspecification_filename_csv.
    if args.__contains__('misspecification_filename_csv') and args.misspecification_filename_csv is not None:
        used = set()
        misspecification_filename_csv = args.misspecification_filename_csv[0] \
            if len(args.misspecification_filename_csv) == 1 else \
            [x for x in args.misspecification_filename_csv if x not in used and (used.add(x) or True)]
    else:
        misspecification_filename_csv = None

    # Run function according to different function name.
    if function_name.lower() == 'run_par_csv_folder'.lower():
        # use_parallel = True
        parval = [pathload, pathsave, use_parallel, filename_pattern, misspecification, Num_drug, Limit_moleculardetection,
                  Limit_radiologicdetection, Limit_mortality, Simduration, Stepsize, dose_method, dose_interval, dose_combination,
                  use_input_dose_combination_only, PAR_criterisa, Strategy_name, lookahead_step, Maxnum_subseq, subtreedepth, run_sim,
                  erase_preresult, fullinput, save_filename_param, save_filename_stopt, save_filename_pop, save_filename_dosage,
                  save_filename_eachtimepoint, misspecification_sim_only, subcolone_LOD, misspecification_LOD, mutation_rate]
        parinput = dict(zip(PARNAME_LIST_RUN_PAR_CSV_FOLDER, parval)) if len(parval) == len(PARNAME_LIST_RUN_PAR_CSV_FOLDER) \
            else None
        DPM_run_par_csv_folder(**parinput)

    if function_name.lower() == 'run_par_csv'.lower():
        parval = [filename_csv, Num_drug, Limit_moleculardetection, Limit_radiologicdetection, Limit_mortality, Simduration, Stepsize, dose_method,
                  dose_interval, dose_combination, use_input_dose_combination_only, PAR_criterisa, Strategy_name, lookahead_step, Maxnum_subseq,
                  subtreedepth, run_sim, erase_preresult, fullinput, misspecification_filename_csv, save_filename_param, save_filename_stopt,
                  save_filename_pop, save_filename_dosage, save_filename_eachtimepoint, misspecification_sim_only, pathsave, par_ind]
        parinput = dict(zip(PARNAME_LIST_RUN_PAR_CSV, parval)) if len(parval) == len(PARNAME_LIST_RUN_PAR_CSV) else None
        DPM_run_par_csv(**parinput)

    if function_name.lower() == 'run_processing'.lower():
        parval = [Num_drug, Strategy_name, pathsave, pathload, Num_stepdiff, Simduration, Stepsize, pathloadmis, use_parallel]
        parinput = dict(zip(PARNAME_LIST_RUN_PROCESSING, parval)) if len(parval) == len(PARNAME_LIST_RUN_PROCESSING) else None
        DPM_run_processing(**parinput)

    if function_name.lower() == 'save_default_PARset'.lower():
        parval = [Num_drug, par_save_block_size, Limit_mortality, Simduration, pathsave]
        parinput = dict(zip(PARNAME_LIST_SAVE_DEFAULT_PARSET, parval)) if len(parval) == len(PARNAME_LIST_SAVE_DEFAULT_PARSET) else None
        DPM_save_default_PARset(**parinput)

    if function_name.lower() == 'run_analysis_LOD'.lower():
        parval = [Num_drug, Simduration, Strategy_name, filename_pattern, pathsave, pathload, pathloadmis, LOD, plot]
        parinput = dict(zip(PARNAME_LIST_RUN_ANALYSIS_LOD, parval)) if len(parval) == len(PARNAME_LIST_RUN_ANALYSIS_LOD) else None
        DPM_run_analysis_LOD(**parinput)

    # pprint.pprint(vars(args))
    return


if __name__ == '__main__':
    DPM_cl_main()
