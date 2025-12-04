from DPM_lib import datetime, os
from DPM_generate import *
from DPM_model import *
from DPM_strategy import strategy_CPM
'''This script assigns and checks the arguments used in the DPM (Dynamic Precision Medicine) model.'''
'''
CHECKED
'''


# Check simulation parameter.#
# (i) The total cell number ('X0total') must be less than 'Limit_mortality' at t = 0; otherwise mortality occurs immediately at t = 0.#
# (ii) Simduration must be greater than stepsize; otherwise only a single step is simulated.#
# Return whether the parameter satisfy these conditions.#
def DPM_check_par(X0total, Limit_mortality, Simduration, Stepsize):
    Par_check_index = []
    # (1) If 'X0total' >= 'Limit_mortality', mortality occurs at t = 0.#
    if X0total >= Limit_mortality:
        Par_check_index.append('X0total >= Limit_mortality')
    # (2) 'Stepsize' > 'Simuduration': the 'stepsize' exceeds the entire simulation period, so even one full step cannot be completed.#
    if Stepsize > Simduration:
        Par_check_index.append('Stepsize > Simuduration')
    return Par_check_index


# Check whether the parameter is a dictionary.#
def DPM_check_dict_par(par):
    par_text = Colored.BOLD + Colored.PURPLE + 'par' + Colored.END
    print_quotation = True
    flag_inlist = True
    beginning_char = '\n'
    if type(par) is not dict:
        print('The type of the inputted ' + par_text + ' should be dictionary.')
        DPM_print_errorandclose()
        return
    keys = par.keys()
    keys_missing = set()
    if 'Num_drug' not in keys:
        Num_drug_text = Colored.BOLD + Colored.PURPLE + 'Num_drug' + Colored.END
        print(Num_drug_text + ' attribute is missing in the input ' + par_text + '.')
        DPM_print_errorandclose()
        return
    Num_drug = DPM_check_one_numerical_par(par['Num_drug'], 'Num_drug', '{:g}'.format, NUM_DRUG_UPPERTHRES, NUM_DRUG_LOWERTHRES, flag_inlist,
                                           NUM_DRUG_LIST, print_quotation, beginning_char)
    if Num_drug is None:
        return
    elif Num_drug == 2:
        keys_missing = set(PARFORMAT_LIST_2DRUG[1:2] + PARFORMAT_LIST_2DRUG[3:]).difference(set(keys))
    elif Num_drug == 3:
        keys_missing = set(PARFORMAT_LIST_3DRUG[1:2] + PARFORMAT_LIST_3DRUG[3:]).difference(set(keys))

    if not keys_missing == set():
        if len(keys_missing) == 1:
            print('The following attribute is missing in the input ' + par_text + ':')
        else:
            print('The following attributes are missing in the input ' + par_text + ':')
        DPM_print_list(list(keys_missing), print_quotation)
        DPM_print_errorandclose()
        return
    else:
        return par


# Check whether a numerical parameter is within its valid range or list.#
# Return the parameter if valid; otherwise, return None.#
def DPM_check_one_numerical_par(par, par_name, par_display, upper_limit, lower_limit, flag_inlist, lis=None, print_quotation=None,
                                beginning_char=''):
    if type(par) not in [int, float]:
        DPM_print_notnumerical(par_name)
        DPM_print_errorandclose()
        return
    if flag_inlist:
        par = DPM_check_in_list_par([par], par_name, lis, print_quotation, beginning_char)
        if par is None:
            return
        else:
            return par[0]
    else:
        if par > upper_limit:
            print('The inputted ' + Colored.BOLD + Colored.PURPLE + par_name + Colored.END + ' = ' + Colored.BOLD + Colored.PURPLE +
                  par_display(par) + Colored.END + ' is larger than the upper limit ' + Colored.BOLD + Colored.PURPLE +
                  par_display(upper_limit) + Colored.END + '.')
            DPM_print_errorandclose()
            return
        elif par < lower_limit:
            print('The inputted ' + Colored.BOLD + Colored.PURPLE + par_name + Colored.END + ' = ' + Colored.BOLD + Colored.PURPLE +
                  par_display(par) + Colored.END + ' is lower than the lower limit ' + Colored.BOLD + Colored.PURPLE +
                  par_display(lower_limit) + Colored.END + '.')
            DPM_print_errorandclose()
            return
        else:
            return par


# Check whether the parameter is a string or belongs to a specified list.#
# Return the parameter if valid; otherwise, return None.#
def DPM_check_one_string_par(par, par_name, flag_inlist, lis, print_quotation, beginning_char):
    if type(par) is not str:
        DPM_print_notstring(par_name)
        DPM_print_errorandclose()
        return
    if flag_inlist:
        par = DPM_check_in_list_par([par], par_name, lis, print_quotation, beginning_char)
        if par is None:
            return
        else:
            return par[0]
    else:
        return par


# Check whether the parameter is a boolean.#
# Return the parameter if it is, otherwise return None.#
def DPM_check_one_boolean_par(par, par_name):
    if type(par) is not bool:
        DPM_print_notboolean(par_name)
        DPM_print_errorandclose()
        return
    else:
        return par


# Check whether the input parameter is in the input list. Return it if present; otherwise return None.#
def DPM_check_in_list_par(par, par_name, lis, print_quotation, beginning_char=''):
    if not (all(elem in lis for elem in par)):
        print('The allowed input of ' + Colored.BOLD + Colored.PURPLE + par_name + Colored.END + ': ', end='')
        DPM_print_list(lis, print_quotation)
        element_not_in_list = [elem for elem in par if elem not in lis]
        print('\n', end='')
        print('The following input of ' + Colored.BOLD + Colored.PURPLE + par_name + Colored.END + ' not allowed: ', end='')
        DPM_print_list(element_not_in_list, print_quotation)
        DPM_print_errorandclose(beginning_char)
        return
    else:
        return par


# Check a list of numerical tuples used for dose combination input.#
# (i) The number of elements in each tuple must match the required length.#
# (ii) Each element in the tuple must be convertible to a numeric value.#
# Return the list of valid inputs if any exist; otherwise return None.#
def DPM_check_list_of_tuple_of_numerical_par(par, par_name, num_element_in_tuple, par_element_lowerthres, par_element_upperthres,
                                             par_sum_lowerthres, par_sum_upperthres, print_quotation):
    text_input = Colored.BOLD + Colored.PURPLE + par_name + Colored.END
    valid_input, invalid_input = [], []
    if par:
        for par_tuple in par:
            # The number of elements in the tuple matches the required number.#
            if type(par_tuple) is tuple and len(par_tuple) == num_element_in_tuple:
                try:
                    par_tuple = tuple(float(i) for i in par_tuple)
                    if all(par_element_lowerthres <= i <= par_element_upperthres for i in par_tuple) \
                            and par_sum_lowerthres <= sum(par_tuple) <= par_sum_upperthres:
                        valid_input.append(par_tuple)
                    else:
                        invalid_input.append(par_tuple)
                except ValueError:
                    invalid_input.append(par_tuple)
            else:
                invalid_input.append(par_tuple)

        if valid_input:
            if len(valid_input) > 1:
                print('The following ' + text_input + ' are valid inputs and accepted:')
            else:
                print('The following ' + text_input + ' is valid input and accepted:')
            DPM_print_list(valid_input, print_quotation)

        if invalid_input:
            if len(invalid_input) > 1:
                print('\nThe following ' + text_input + ' are invalid inputs and ignored:')
            else:
                print('\nThe following ' + text_input + ' is invalid input and ignored:')
            DPM_print_list(invalid_input, print_quotation)
            print('\n')
        return valid_input
    else:
        return


# Check the parameters in 'PARNAME_LIST'.#
def DPM_check_inputpar(kwargs, parname_list, caller=None):
    parname_list = [i_name.lower() for i_name in parname_list]

    print_quotation, flag_inlist, beginning_char = True, True, '\n'

    # Check that 'Num_drug' is numerical and within its valid range.#
    if 'Num_drug'.lower() in parname_list:
        if 'Num_drug'.lower() in kwargs.keys():
            Num_drug = kwargs['Num_drug'.lower()]
            Num_drug = DPM_check_one_numerical_par(Num_drug, 'Num_drug', '{:g}'.format, NUM_DRUG_UPPERTHRES, NUM_DRUG_LOWERTHRES, flag_inlist,
                                                   NUM_DRUG_LIST, print_quotation, beginning_char)
            if Num_drug is None:
                sys.exit()
        else:
            Num_drug = NUM_DRUG_DEFAULT_VAL
    else:
        Num_drug = ''

    # Check whether the parameter 'par' was provided.#
    if 'par'.lower() in parname_list:
        if 'par'.lower() in kwargs.keys():
            par = kwargs['par']
            par = DPM_check_dict_par(par)
            if par is None:
                sys.exit()
            Num_drug = par['Num_drug']
        else:
            par_text = Colored.BOLD + Colored.PURPLE + 'par' + Colored.END
            print('No ' + par_text + ' keyword argument was provided. The function requires ' + par_text + ' as input.')
            DPM_print_errorandclose()
            sys.exit()
    else:
        par = ''

    # Check whether the parameter 'par_mis' was provided.#
    if 'mis_par'.lower() in parname_list:
        if 'mis_par'.lower() in kwargs.keys():
            mis_par = kwargs['mis_par']
            mis_par = DPM_check_dict_par(mis_par)
            if mis_par is None:
                sys.exit()
        else:
            par_text = Colored.BOLD + Colored.PURPLE + 'mis_par' + Colored.END
            print('No ' + par_text + ' keyword argument was provided. The function requires ' + par_text + ' as input.')
            DPM_print_errorandclose()
            sys.exit()
    else:
        mis_par = ''

    # Check that 'dose_method' is a string and is either continous or discrete.#
    if 'dose_method'.lower() in parname_list:
        if 'dose_method'.lower() in kwargs.keys():
            dose_method = kwargs['dose_method'.lower()]
            dose_method = DPM_check_one_string_par(dose_method, 'dose_method', flag_inlist, DOSE_METHOD_LIST, print_quotation, beginning_char)
            if dose_method is None:
                sys.exit()
            elif dose_method in ['continous', 'c']:
                dose_method = 'c'
            else:
                dose_method = 'd'
        else:
            dose_method = DOSE_METHOD_DEFAULT_VAL
    else:
        dose_method = ''

    # Check that 'dose_interval' is numerical and within its valid range.#
    if 'dose_interval'.lower() in parname_list:
        if dose_method == 'c':
            if 'dose_interval'.lower() in kwargs.keys():
                dose_interval = kwargs['dose_interval'.lower()]
                dose_interval = DPM_check_one_numerical_par(dose_interval, 'dose_interval', '{:g}'.format, DOSE_INTERVAL_UPPERTHRES,
                                                            DOSE_INTERVAL_LOWERTHRES, flag_inlist, DOSE_INTERVAL_LIST, print_quotation,
                                                            beginning_char)
                if dose_interval is None:
                    sys.exit()
            else:
                dose_interval = DOSE_INTERVAL_DEFAULT_VAL
        else:
            dose_interval = ''
    else:
        dose_interval = ''

    # Generate the default 'dose_combination' parameter.#
    if dose_method == 'd' and Num_drug != '':
        dose_combination = DPM_generate_discrete_dose_combination(Num_drug)
    elif dose_method == 'c':
        dose_combination = DPM_generate_continuous_dose_combination(Num_drug, dose_interval)
    else:
        dose_combination = []

    # Check that 'use_input_dose_combination_only' is a boolean.#
    if 'use_input_dose_combination_only'.lower() in parname_list:
        if 'use_input_dose_combination_only'.lower() in kwargs.keys():
            use_input_dose_combination_only = kwargs['use_input_dose_combination_only'.lower()]
            use_input_dose_combination_only = DPM_check_one_boolean_par(use_input_dose_combination_only, 'use_input_dose_combination_only')
            if use_input_dose_combination_only is None:
                sys.exit()
        else:
            use_input_dose_combination_only = False
    else:
        use_input_dose_combination_only = ''

    # Check that dose_combination is a list of tuples and that each entry is within the valid range.#
    if 'dose_combination'.lower() in parname_list:
        if 'dose_combination'.lower() in kwargs.keys():
            dose_combination_input = kwargs['dose_combination'.lower()]
            if type(dose_combination_input) is not list:
                dose_combination_input = [dose_combination_input]
            dose_combination_input = DPM_check_list_of_tuple_of_numerical_par(
                dose_combination_input, 'dose_combination', Num_drug, DOSE_COMBINATION_ELEMENT_LOWERTHRES, DOSE_COMBINATION_ELEMENT_UPPERTHRES,
                DOSE_COMBINATION_SUM_LOWERTHRES, DOSE_COMBINATION_SUM_UPPERTHRES, print_quotation)
            if dose_combination_input:
                if not use_input_dose_combination_only:
                    dose_combination_input = [x for x in dose_combination_input if x not in dose_combination]
                    dose_combination.extend(dose_combination_input)
                elif use_input_dose_combination_only:
                    dose_combination = dose_combination_input
    else:
        dose_combination = ''

    flag_inlist = False

    # Check that 'Simduration' is numerical and within its valid range.#
    if 'Simduration'.lower() in parname_list:
        if 'Simduration'.lower() in kwargs.keys():
            Simduration = kwargs['Simduration'.lower()]
            Simduration = DPM_check_one_numerical_par(Simduration, 'Simduration', '{:g}'.format, SIMDURATION_UPPERTHRES,
                                                      SIMDURATION_LOWERTHRES, flag_inlist)
            if Simduration is None:
                sys.exit()
        else:
            Simduration = SIMDURATION_DEFAULT_VAL
    else:
        Simduration = ''

    # Check that 'Limit_mortality' is numerical and within its valid range.#
    if 'Limit_mortality'.lower() in parname_list:
        if 'Limit_mortality'.lower() in kwargs.keys():
            Limit_mortality = kwargs['Limit_mortality'.lower()]
            Limit_mortality = DPM_check_one_numerical_par(Limit_mortality, 'Limit_mortality', '{:.2e}'.format,
                                                          LIMIT_MORTALITY_UPPERTHRES, LIMIT_MORTALITY_LOWERTHRES, flag_inlist)
            if Limit_mortality is None:
                sys.exit()
        else:
            Limit_mortality = LIMIT_MORTALITY_DEFAULT_VAL
    else:
        Limit_mortality = ''

    # Check that 'Stepsize' is numerical and within its valid range.#
    if 'Stepsize'.lower() in parname_list:
        if 'Stepsize'.lower() in kwargs.keys():
            Stepsize = kwargs['Stepsize'.lower()]
            Stepsize = DPM_check_one_numerical_par(Stepsize, 'Stepsize', '{:g}'.format, STEPSIZE_UPPERTHRES,
                                                   STEPSIZE_LOWERTHRES, flag_inlist)
            if Stepsize is None:
                sys.exit()
        else:
            Stepsize = STEPSIZE_DEFAULT_VAL
    else:
        Stepsize = ''

    # Check that 'Limit_radiologicdetection' is numerical and within its valid range.#
    if 'Limit_radiologicdetection'.lower() in parname_list:
        if 'Limit_radiologicdetection'.lower() in kwargs.keys():
            Limit_radiologicdetection = kwargs['Limit_radiologicdetection'.lower()]
            Limit_radiologicdetection = DPM_check_one_numerical_par(Limit_radiologicdetection, 'Limit_radiologicdetection',
                                                                    '{:.2e}'.format, LIMIT_RADIOLOGICDETECTION_UPPERTHRES,
                                                                    LIMIT_RADIOLOGICDETECTION_LOWERTHRES, flag_inlist)
            if Limit_radiologicdetection is None:
                sys.exit()
        else:
            Limit_radiologicdetection = LIMIT_RADIOLOGICDETECTION_DEFAULT_VAL
    else:
        Limit_radiologicdetection = ''

    # Check that 'Limit_moleculardetection' is numerical and within its valid range.#
    if 'Limit_moleculardetection'.lower() in parname_list:
        if 'Limit_moleculardetection'.lower() in kwargs.keys():
            Limit_moleculardetection = kwargs['Limit_moleculardetection'.lower()]
            Limit_moleculardetection = DPM_check_one_numerical_par(Limit_moleculardetection, 'Limit_moleculardetection',
                                                                   '{:.2e}'.format, LIMIT_MOLECULARDETECTION_UPPERTHRES,
                                                                   LIMIT_MOLECULARDETECTION_LOWERTHRES, flag_inlist)
            if Limit_moleculardetection is None:
                sys.exit()
        else:
            Limit_moleculardetection = LIMIT_MOLECULARDETECTION_DEFAULT_VAL
    else:
        Limit_moleculardetection = ''

    # Check that 'run_sim' is a boolean.#
    if 'run_sim'.lower() in parname_list:
        if 'run_sim'.lower() in kwargs.keys():
            run_sim = kwargs['run_sim'.lower()]
            run_sim = DPM_check_one_boolean_par(run_sim, 'run_sim')
            if run_sim is None:
                sys.exit()
        else:
            run_sim = True
    else:
        run_sim = ''

    # Check that 'use_parallel' is a boolean.#
    if 'use_parallel'.lower() in parname_list:
        if 'use_parallel'.lower() in kwargs.keys():
            use_parallel = kwargs['use_parallel'.lower()]
            use_parallel = DPM_check_one_boolean_par(use_parallel, 'use_parallel')
            if use_parallel is None:
                sys.exit()
        else:
            use_parallel = False
    else:
        use_parallel = ''

    # Generate parameter criteria based on the number of drugs.#
    if Num_drug in [2, 3]:
        PAR_criterisa_list = DPM_generate_PAR_criterisa_list(Num_drug)
    else:
        PAR_criterisa_list = []

    # Check that 'misspecification_sim_only' is a boolean.#
    if 'misspecification_sim_only'.lower() in parname_list:
        if 'misspecification_sim_only'.lower() in kwargs.keys():
            misspecification_sim_only = kwargs['misspecification_sim_only'.lower()]
            misspecification_sim_only = DPM_check_one_boolean_par(misspecification_sim_only, 'misspecification_sim_only')
            if misspecification_sim_only is None:
                sys.exit()
        else:
            misspecification_sim_only = False
    else:
        misspecification_sim_only = ''

    # Check that all elements of 'PAR_criterisa' are included in the allowed list.#
    if 'PAR_criterisa'.lower() in parname_list:
        if 'PAR_criterisa'.lower() in kwargs.keys():
            PAR_criterisa = kwargs['PAR_criterisa'.lower()]
            if type(PAR_criterisa) is not list:
                PAR_criterisa = [PAR_criterisa]
            PAR_criterisa = DPM_check_in_list_par(PAR_criterisa, 'PAR_criterisa', PAR_criterisa_list, print_quotation, beginning_char)
            if PAR_criterisa is None:
                sys.exit()
            else:
                # Make 'PAR_criterisa' unique.#
                PAR_criterisa = list(set(PAR_criterisa))
        else:
            PAR_criterisa = []
    else:
        PAR_criterisa = ''

    # Check that 'Strategy_name' is in the list.#
    if 'Strategy_name'.lower() in parname_list:
        if 'Strategy_name'.lower() in kwargs.keys():
            Strategy_name = kwargs['Strategy_name'.lower()]
            if type(Strategy_name) is not list:
                Strategy_name = list(Strategy_name)
            Strategy_name = DPM_check_in_list_par([item.upper() for item in Strategy_name], 'Strategy_name', STRATEGY_LIST, print_quotation,
                                                  beginning_char)
            if Strategy_name is None:
                sys.exit()
        else:
            Strategy_name = ['CPM', 'DPM2.2']
        if 'all' in Strategy_name:
            Strategy_name = STRATEGY_LIST
    else:
        Strategy_name = ''

    # Check whether 'filename_csv' was provided.#
    if 'filename_csv'.lower() in parname_list:
        if 'filename_csv'.lower() in kwargs.keys():
            filename_csv = kwargs['filename_csv']
        else:
            filename_csv_text: str = Colored.BOLD + Colored.PURPLE + 'filename_csv' + Colored.END
            print('No ' + filename_csv_text + ' keyword argument input. Need at least input the name of the ' +
                  Colored.BOLD + Colored.PURPLE + '.csv' + Colored.END + ' file containing the parameters.')
            DPM_print_errorandclose()
            sys.exit()
    else:
        filename_csv = ''

    # Check that 'misspecification_filename_csv' has the same type and length as filename_csv.#
    if 'misspecification_filename_csv'.lower() in parname_list:
        if 'misspecification_filename_csv' in kwargs.keys():
            misspecification_filename_csv = kwargs['misspecification_filename_csv']
            if type(misspecification_filename_csv) != type(filename_csv):
                print('The type of the input keyword argument ' + Colored.BOLD + Colored.PURPLE + 'filename_csv' + Colored.END +
                      ' is different from the type of the input keyword argument ' + Colored.BOLD + Colored.PURPLE +
                      'misspecification_filename_csv' + Colored.END + '.')
                DPM_print_errorandclose()
                sys.exit()
            elif type(misspecification_filename_csv) is list and len(misspecification_filename_csv) != len(filename_csv):
                print('The length of the input keyword argument ' + Colored.BOLD + Colored.PURPLE + 'filename_csv' + Colored.END +
                      ' is different from the length of the input keyword argument ' + Colored.BOLD + Colored.PURPLE +
                      'misspecification_filename_csv' + Colored.END + '.')
                DPM_print_errorandclose()
                sys.exit()
        else:
            misspecification_filename_csv = ''
    else:
        misspecification_filename_csv = ''

    # Check that 'par_save_block_size' is numerical and within its valid range.#
    if 'par_save_block_size'.lower() in parname_list:
        flag_inlist = False
        if 'par_save_block_size'.lower() in kwargs.keys():
            par_save_block_size = kwargs['par_save_block_size'.lower()]
            par_save_block_size = DPM_check_one_numerical_par(par_save_block_size, 'par_save_block_size', '{:.2e}'.format,
                                                              PAR_SAVE_BLOCK_SIZE_UPPERTHRES, PAR_SAVE_BLOCK_SIZE_LOWERTHRES,
                                                              flag_inlist)
            if par_save_block_size is None:
                sys.exit()
        else:
            par_save_block_size = PAR_SAVE_BLOCK_SIZE_DEFAULT_VAL
    else:
        par_save_block_size = ''

    # Check that 'save_filename_param' is a boolean.#
    if 'save_filename_param'.lower() in parname_list:
        if 'save_filename_param'.lower() in kwargs.keys():
            save_filename_param = kwargs['save_filename_param'.lower()]
            save_filename_param = DPM_check_one_boolean_par(save_filename_param, 'save_filename_param')
            if save_filename_param is None:
                sys.exit()
        else:
            save_filename_param = True
    else:
        save_filename_param = ''

    # Check that 'save_filename_stopt' is a boolean.#
    if 'save_filename_stopt'.lower() in parname_list:
        if 'save_filename_stopt'.lower() in kwargs.keys():
            save_filename_stopt = kwargs['save_filename_stopt'.lower()]
            save_filename_stopt = DPM_check_one_boolean_par(save_filename_stopt, 'save_filename_stopt')
            if save_filename_stopt is None:
                sys.exit()
        else:
            save_filename_stopt = True
    else:
        save_filename_stopt = ''

    # Check that 'save_filename_pop' is a boolean.#
    if 'save_filename_pop'.lower() in parname_list:
        if 'save_filename_pop'.lower() in kwargs.keys():
            save_filename_pop = kwargs['save_filename_pop'.lower()]
            save_filename_pop = DPM_check_one_boolean_par(save_filename_pop, 'save_filename_pop')
            if save_filename_pop is None:
                sys.exit()
        else:
            save_filename_pop = True
    else:
        save_filename_pop = ''

    # Check that 'save_filename_dosage' is a boolean.#
    if 'save_filename_dosage'.lower() in parname_list:
        if 'save_filename_dosage'.lower() in kwargs.keys():
            save_filename_dosage = kwargs['save_filename_dosage'.lower()]
            save_filename_dosage = DPM_check_one_boolean_par(save_filename_dosage, 'save_filename_dosage')
            if save_filename_dosage is None:
                sys.exit()
        else:
            save_filename_dosage = True
    else:
        save_filename_dosage = ''

    # Check that 'save_filename_eachtimepoint' is a boolean.#
    if 'save_filename_eachtimepoint'.lower() in parname_list:
        if 'save_filename_eachtimepoint'.lower() in kwargs.keys():
            save_filename_eachtimepoint = kwargs['save_filename_eachtimepoint'.lower()]
            save_filename_eachtimepoint = DPM_check_one_boolean_par(save_filename_eachtimepoint, 'save_filename_eachtimepoint')
            if save_filename_eachtimepoint is None:
                sys.exit()
        else:
            save_filename_eachtimepoint = False
    else:
        save_filename_eachtimepoint = ''

    flag_inlist = False
    path_current = os.getcwd()
    # Check that 'pathsave' is a vaild directory. If it does not exist, attempt to creat it; otherwise use the current working directory.#
    if 'pathsave'.lower() in parname_list:
        if 'pathsave'.lower() in kwargs.keys():
            pathsave = kwargs['pathsave'.lower()]
            pathsave = DPM_check_one_string_par(pathsave, 'pathsave', flag_inlist, [], print_quotation, beginning_char)
            if pathsave is None:
                sys.exit()
            # If the save path does not exist.#
            pathsave = os.path.abspath(pathsave)
            if not os.path.exists(pathsave):
                try:
                    os.makedirs(pathsave)
                except (OSError, FileExistsError):
                    pathsave_text = Colored.BOLD + Colored.PURPLE + pathsave + Colored.END
                    print('The specified ' + Colored.BOLD + Colored.PURPLE + 'pathsave' + Colored.END + ':')
                    print(pathsave_text + ' does not exist and could not be created. The current working directory will be used instead.')  #
                    pathsave = path_current
        else:
            pathsave = path_current
        print('The directory for saving results is set to:')
        print(Colored.BOLD + Colored.PURPLE + pathsave + Colored.END)
    else:
        pathsave = ''

    path_current = os.getcwd()
    # Check that 'pathload' is a valid directory. If it does not exist, attempt to load data from the current directory.#
    if 'pathload'.lower() in parname_list:
        if 'pathload'.lower() in kwargs.keys():
            pathload = kwargs['pathload'.lower()]
            pathload = DPM_check_one_string_par(pathload, 'pathload', flag_inlist, [], print_quotation, beginning_char)
            if pathload is None:
                sys.exit()
            # If the data load path does not exist.#
            pathload = os.path.abspath(pathload)
            if not os.path.exists(pathload):
                pathload_text = Colored.BOLD + Colored.PURPLE + pathload + Colored.END
                print('The specified ' + Colored.BOLD + Colored.PURPLE + 'pathload' + Colored.END + ':')
                print(pathload_text + ' does not exist. The current working directory will be used instead.')  #
                pathload = path_current
        else:
            pathload = path_current
        print('The directory for loading files is set to:')
        print(Colored.BOLD + Colored.PURPLE + pathload + Colored.END)
    else:
        pathload = ''

    # Check that 'pathload' is a vaild directory. If it does not exist, attempt to load data from the current directory.#
    if 'misspecification_fileload'.lower() in parname_list:
        if 'misspecification_fileload'.lower() in kwargs.keys():
            misspecification_fileload = kwargs['misspecification_fileload'.lower()]
            misspecification_fileload = DPM_check_one_string_par(misspecification_fileload, 'misspecification_fileload',
                                                                 flag_inlist, [], print_quotation, beginning_char)
            if misspecification_fileload is None:
                sys.exit()
            if not os.path.exists(misspecification_fileload) and misspecification_fileload != '':
                misspecification_fileload_text = Colored.BOLD + Colored.PURPLE + misspecification_fileload + Colored.END
                print('The specified ' + Colored.BOLD + Colored.PURPLE + 'misspecification_fileload' + Colored.END + ':')
                print(misspecification_fileload_text + ' is not exist.')
                DPM_print_errorandclose()
                sys.exit()
        else:
            misspecification_fileload = ''
    else:
        misspecification_fileload = ''

    # Check that 'pathloadmis' is a vaild directory. If it does not exist, attempt to load data from the current directory.#
    if 'pathloadmis'.lower() in parname_list:
        if 'pathloadmis'.lower() in kwargs.keys():
            pathloadmis = kwargs['pathloadmis'.lower()]
            pathloadmis = DPM_check_one_string_par(pathloadmis, 'pathloadmis', False, [], print_quotation, beginning_char)
            if pathloadmis is None:
                sys.exit()
            # If the load path does not exist.#
            pathloadmis = os.path.abspath(pathloadmis) if pathloadmis != '' else ''
            if (not os.path.exists(pathloadmis)) and pathloadmis != '':
                pathloadmis_text = Colored.BOLD + Colored.PURPLE + pathloadmis + Colored.END
                print('The inputted ' + Colored.BOLD + Colored.PURPLE + 'pathloadmis' + Colored.END + ':')
                print(pathloadmis_text + ' is not exist')
                pathloadmis = ''
        else:
            pathloadmis = ''
    else:
        pathloadmis = ''

    # Check that 'par_ind' is numerical or a list of numerical values.#
    if 'par_ind'.lower() in parname_list:
        if 'par_ind'.lower() in kwargs.keys():
            par_ind_accept = []
            par_ind = kwargs['par_ind'.lower()]
            # Convert 'par_ind' to a list if it is not already one.#
            par_ind = [par_ind] if type(par_ind) is not list else par_ind
            for i_par_ind in par_ind:
                par_ind = DPM_check_one_numerical_par(i_par_ind, 'par_ind', '{:g}'.format, np.inf, -np.inf, flag_inlist, [],
                                                      print_quotation, beginning_char)
                if par_ind is not None and par_ind % 1 == 0:
                    par_ind_accept.append(i_par_ind)
            if len(par_ind_accept) == 0:
                print('No valid index is avaiable in the inputted ' + Colored.BOLD + Colored.PURPLE + 'par_ind' + Colored.END)
                DPM_print_errorandclose()
                sys.exit()
            used: set = set()
            par_ind = [int(x) for x in par_ind_accept if int(x) not in used and (used.add(int(x)) or True)]
        else:
            par_ind = ''
    else:
        par_ind = ''

    # Check that 'savename' is a string.#
    if 'savename'.lower() in parname_list:
        if 'savename'.lower() in kwargs.keys():
            savename = kwargs['savename'.lower()]
            savename = DPM_check_one_string_par(savename, 'savename', flag_inlist, [], print_quotation, beginning_char)
            if savename is None:
                sys.exit()
        else:
            # Use the local time as the save name.#
            savename = datetime.now().strftime('%Y-%m-%d-%H-%M-%S-%f')
    else:
        savename = ''

    # Check that 'save_filename_eachtimepoint' is a boolean.#
    if 'plot'.lower() in parname_list:
        if 'plot'.lower() in kwargs.keys():
            plot = kwargs['plot'.lower()]
            plot = DPM_check_one_boolean_par(plot, 'plot')
            if plot is None:
                sys.exit()
        else:
            plot = True if caller == 'DPM_run_plot_1PAR' else False
    else:
        plot = ''

    # Check that 'Num_stepdiff' is numerical and within its valid range.#
    if 'Num_stepdiff'.lower() in parname_list:
        if 'Num_stepdiff'.lower() in kwargs.keys():
            Num_stepdiff = kwargs['Num_stepdiff'.lower()]
            Num_stepdiff = DPM_check_one_numerical_par(Num_stepdiff, 'Num_stepdiff', '{:g}'.format, NUM_STEPDIFF_UPPERTHRES,
                                                       NUM_STEPDIFF_LOWERTHRES, flag_inlist)
            if Num_stepdiff is None:
                sys.exit()
        else:
            Num_stepdiff = NUM_STEPDIFF_DEFAULT_VAL
    else:
        Num_stepdiff = ''

    # Check that 'filename_pattern' is string.#
    if 'filename_pattern'.lower() in parname_list:
        if 'filename_pattern'.lower() in kwargs.keys():
            filename_pattern = kwargs['filename_pattern'.lower()]
            filename_pattern = DPM_check_one_string_par(filename_pattern, 'filename_pattern', flag_inlist, [], print_quotation, beginning_char)
            if filename_pattern is None:
                sys.exit()
        else:
            filename_pattern = ''
    else:
        filename_pattern = ''

    # Check that 'misspecification_atdecision' is a boolean.#
    if 'misspecification_atdecision'.lower() in parname_list:
        if 'misspecification_atdecision'.lower() in kwargs.keys():
            misspecification_atdecision = kwargs['misspecification_atdecision'.lower()]
            misspecification_atdecision = DPM_check_one_boolean_par(misspecification_atdecision, 'misspecification_atdesiction')
            if misspecification_atdecision is None:
                sys.exit()
        else:
            misspecification_atdecision = False
    else:
        misspecification_atdecision = ''

    # Check that 'misspecification_atsim' is a boolean.#
    if 'misspecification_atsim'.lower() in parname_list:
        if 'misspecification_atsim'.lower() in kwargs.keys():
            misspecification_atsim = kwargs['misspecification_atsim'.lower()]
            misspecification_atsim = DPM_check_one_boolean_par(misspecification_atsim, 'misspecification_atsim')
            if misspecification_atsim is None:
                sys.exit()
        else:
            misspecification_atsim = False
    else:
        misspecification_atsim = ''

    PAR_dict = {'Num_drug': Num_drug, 'dose_method': dose_method, 'dose_interval': dose_interval, 'dose_combination': dose_combination,
                'Simduration': Simduration, 'use_input_dose_combination_only': use_input_dose_combination_only, 'Stepsize': Stepsize,
                'Limit_mortality': Limit_mortality, 'Limit_radiologicdetection': Limit_radiologicdetection, 'run_sim': run_sim,
                'PAR_criterisa': PAR_criterisa, 'Strategy_name': Strategy_name, 'par_save_block_size': par_save_block_size,
                'save_filename_param': save_filename_param, 'save_filename_stopt': save_filename_stopt, 'save_filename_pop': save_filename_pop,
                'save_filename_dosage': save_filename_dosage, 'pathsave': pathsave, 'pathload': pathload,
                'save_filename_eachtimepoint': save_filename_eachtimepoint, 'use_parallel': use_parallel, 'filename_csv': filename_csv,
                'misspecification_filename_csv': misspecification_filename_csv, 'Limit_moleculardetection': Limit_moleculardetection,
                'misspecification_sim_only': misspecification_sim_only, 'savename': savename, 'par': par, 'mis_par': mis_par,
                'plot': plot, 'par_ind': par_ind, 'Num_stepdiff': Num_stepdiff, 'filename_pattern': filename_pattern,
                'misspecification_atdecision': misspecification_atdecision, 'misspecification_atsim': misspecification_atsim,
                'pathloadmis': pathloadmis, 'misspecification_fileload': misspecification_fileload}

    return PAR_dict


# Assign parameters from CSV files.#
def DPM_assign_par_csvread(par_csv, Num_drug):
    PARtotal = []
    for i_par in par_csv:
        if type(i_par) is dict:
            Num_cell_type = 4
            i_par.update({'Num_drug': Num_drug, 'Num_cell_type': Num_cell_type})
            PARtotal.append(i_par)
    return PARtotal


# Assign parameters in the 2 drug case fomart.#
def DPM_assign_par_default_2drug(par, X0total):
    par_out = dict()

    g0_S, ratioR1toX0, ratioR2toX0, Sa_ratio_S_D1tog0_S, Sa_ratio_S_D2toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, T_StoR1, T_StoR2 = par
    ratioStoX0, ratioR12toX0 = 1-ratioR1toX0-ratioR2toX0, 0

    g0_R1 = g0_R2 = g0_R12 = g0_S
    Sa_ratio_R1_D2toS_D2 = Sa_ratio_R2_D1toS_D1 = 1.0
    Sa_ratio_R12_D1toS_D1, Sa_ratio_R12_D2toS_D2 = Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2

    T_R1toR12, T_R2toR12 = T_StoR2, T_StoR1
    # X0#
    par_out['Spop'], par_out['R1pop'], par_out['R2pop'], par_out['R12pop'] = \
        X0total * np.array([ratioStoX0, ratioR1toX0,  ratioR2toX0, ratioR12toX0]).ravel()
    # g0#
    par_out['g0_S'], par_out['g0_R1'], par_out['g0_R2'], par_out['g0_R12'] = np.array([g0_S, g0_R1, g0_R2, g0_R12]).ravel()
    # Sa#
    par_out['Sa.S.D1.'] = Sa_ratio_S_D1tog0_S * g0_S
    par_out['Sa.S.D2.'] = Sa_ratio_S_D2toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R1.D1.'] = Sa_ratio_R1_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R1.D2.'] = Sa_ratio_R1_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R2.D1.'] = Sa_ratio_R2_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R2.D2.'] = Sa_ratio_R2_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R12.D1.'] = Sa_ratio_R12_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R12.D2.'] = Sa_ratio_R12_D2toS_D2 * par_out['Sa.S.D2.']
    # T#
    par_out['T.R1..S.'] = T_StoR1
    par_out['T.R2..S.'] = T_StoR2
    par_out['T.R12..R1.'] = T_R1toR12
    par_out['T.R12..R2.'] = T_R2toR12

    T = [i for i in HEADER_2DRUG_INPUT_PARAM_CSV if i.find('T.') != -1]
    for i_T in T:
        if i_T not in par_out.keys():
            par_out[i_T] = 0.0
    return par_out


# Assign parameters in the 3 drug case fomart.#
def DPM_assign_par_default_3drug(par, X0total):
    par_out = dict()

    g0_S, ratioStoX0, ratioR1toX0, ratioR2toX0, ratioR3toX0, ratioR12toX0, ratioR23toX0, ratioR13toX0, Sa_ratio_S_D1tog0_S, \
        Sa_ratio_S_D2toS_D1, Sa_ratio_S_D3toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, Sa_ratio_R3_D3toS_D3, T_StoR1, \
        T_StoR2, T_StoR3 = par
    ratioR123toX0 = 0.0
    g0_R1 = g0_R2 = g0_R3 = g0_R12 = g0_R13 = g0_R23 = g0_R123 = g0_S

    Sa_ratio_R1_D2toS_D2 = Sa_ratio_R1_D3toS_D3 = Sa_ratio_R2_D1toS_D1 = Sa_ratio_R2_D3toS_D3 = Sa_ratio_R3_D1toS_D1 = Sa_ratio_R3_D2toS_D2 = \
        Sa_ratio_R12_D3toS_D3 = Sa_ratio_R13_D2toS_D2 = Sa_ratio_R23_D1toS_D1 = 1.0
    Sa_ratio_R12_D1toS_D1 = Sa_ratio_R13_D1toS_D1 = Sa_ratio_R123_D1toS_D1 = Sa_ratio_R1_D1toS_D1
    Sa_ratio_R12_D2toS_D2 = Sa_ratio_R23_D2toS_D2 = Sa_ratio_R123_D2toS_D2 = Sa_ratio_R2_D2toS_D2
    Sa_ratio_R13_D3toS_D3 = Sa_ratio_R23_D3toS_D3 = Sa_ratio_R123_D3toS_D3 = Sa_ratio_R3_D3toS_D3

    T_R2toR12 = T_R3toR13 = T_R23toR123 = T_StoR1
    T_R1toR12 = T_R3toR23 = T_R13toR123 = T_StoR2
    T_R1toR13 = T_R2toR23 = T_R12toR123 = T_StoR3
    # X0#
    par_out['Spop'], par_out['R1pop'], par_out['R2pop'], par_out['R3pop'], par_out['R12pop'], par_out['R13pop'], par_out['R23pop'], \
        par_out['R123pop'] = X0total * np.array([ratioStoX0, ratioR1toX0, ratioR2toX0, ratioR3toX0, ratioR12toX0, ratioR23toX0,
                                                 ratioR13toX0, ratioR123toX0]).ravel()
    # g0#
    par_out['g0_S'], par_out['g0_R1'], par_out['g0_R2'], par_out['g0_R3'], par_out['g0_R12'], par_out['g0_R13'], par_out['g0_R23'],\
        par_out['g0_R123'] = np.array([g0_S, g0_R1, g0_R2, g0_R3, g0_R12, g0_R13, g0_R23, g0_R123]).ravel()
    # Sa#
    par_out['Sa.S.D1.'] = Sa_ratio_S_D1tog0_S * g0_S
    par_out['Sa.S.D2.'] = Sa_ratio_S_D2toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.S.D3.'] = Sa_ratio_S_D3toS_D1 * par_out['Sa.S.D1.']

    par_out['Sa.R1.D1.'] = Sa_ratio_R1_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R1.D2.'] = Sa_ratio_R1_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R1.D3.'] = Sa_ratio_R1_D3toS_D3 * par_out['Sa.S.D3.']

    par_out['Sa.R2.D1.'] = Sa_ratio_R2_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R2.D2.'] = Sa_ratio_R2_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R2.D3.'] = Sa_ratio_R2_D3toS_D3 * par_out['Sa.S.D3.']

    par_out['Sa.R3.D1.'] = Sa_ratio_R3_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R3.D2.'] = Sa_ratio_R3_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R3.D3.'] = Sa_ratio_R3_D3toS_D3 * par_out['Sa.S.D3.']

    par_out['Sa.R12.D1.'] = Sa_ratio_R12_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R12.D2.'] = Sa_ratio_R12_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R12.D3.'] = Sa_ratio_R12_D3toS_D3 * par_out['Sa.S.D3.']

    par_out['Sa.R13.D1.'] = Sa_ratio_R13_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R13.D2.'] = Sa_ratio_R13_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R13.D3.'] = Sa_ratio_R13_D3toS_D3 * par_out['Sa.S.D3.']

    par_out['Sa.R23.D1.'] = Sa_ratio_R23_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R23.D2.'] = Sa_ratio_R23_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R23.D3.'] = Sa_ratio_R23_D3toS_D3 * par_out['Sa.S.D3.']

    par_out['Sa.R123.D1.'] = Sa_ratio_R123_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R123.D2.'] = Sa_ratio_R123_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R123.D3.'] = Sa_ratio_R123_D3toS_D3 * par_out['Sa.S.D3.']
    # T#
    par_out['T.R1..S.'] = T_StoR1
    par_out['T.R2..S.'] = T_StoR2
    par_out['T.R3..S.'] = T_StoR3
    par_out['T.R12..R1.'] = T_R1toR12
    par_out['T.R13..R1.'] = T_R1toR13
    par_out['T.R12..R2.'] = T_R2toR12
    par_out['T.R23..R2.'] = T_R2toR23
    par_out['T.R13..R3.'] = T_R3toR13
    par_out['T.R23..R3.'] = T_R3toR23
    par_out['T.R123..R12.'] = T_R12toR123
    par_out['T.R123..R13.'] = T_R13toR123
    par_out['T.R123..R23.'] = T_R23toR123

    T = [i for i in HEADER_3DRUG_INPUT_PARAM_CSV if i.find('T.') != -1]

    for i_T in T:
        if i_T not in par_out.keys():
            par_out[i_T] = 0.0

    return par_out


# Assign strategy results for saving to .csv files.#
def DPM_assign_saveresult_csv(Strategy, Strategy_text, Simduration, Stepsize, cured_survival, non_report, para_PAR_index, filename_stopt,
                              filename_dosage, filename_pop, filename_eachtimepoint, Limit_mortality):
    if Strategy is not None:
        t, X, d = Strategy
        t_atstepsize = np.arange(0, Simduration+Stepsize, Stepsize).tolist()

        if filename_dosage:
            # Find the index of t.#
            pos = [int(np.where(t == i_step)[0]) if np.isin(i_step, t, assume_unique=True) and i_step < max(t) else -1
                   for i_step in t_atstepsize[0:-1]]
            dose_Strategy = [tuple(d[:, x].T) if x != -1 else -1 for x in pos]
            # Round the tuple of dose values to two digits.#
            dose_Strategy = [tuple(round(y, 2) for y in x) if type(x) is tuple else -1 for x in dose_Strategy]
            dose_Strategy = para_PAR_index + Strategy_text + ', '.join(str(x).replace(' ', '') for x in dose_Strategy)
        else:
            dose_Strategy = ''

        if filename_pop:
            pos = [int(np.where(t == i_step)[0]) if np.isin(i_step, t, assume_unique=True) else -1 for i_step in t_atstepsize[1:]]
            pop_Strategy = [tuple(X[:, x].T) if x != -1 else -1 for x in pos]
            pop_Strategy = [tuple(round(y, 2) for y in x) if type(x) is tuple else -1 for x in pop_Strategy]
            pop_Strategy = para_PAR_index + Strategy_text + ', '.join('(' + ','.join('{:.2e}'.format(y) for y in x) + ')'
                                                                      if type(x) is tuple else str(x) for x in pop_Strategy)
        else:
            pop_Strategy = ''

        if filename_stopt:
            if all(X[:, -1] < 1):
                survival_Strategy = cured_survival
            else:
                survival_Strategy = min(t[-1], Simduration)
                if survival_Strategy == Simduration and sum(X[:, -1]) < Limit_mortality:
                    survival_Strategy += 1
        else:
            survival_Strategy = ''

        if filename_eachtimepoint:
            X = map(tuple, X.T)
            X = tuple(X)
            X = [tuple(round(y, 2) for y in x) if type(x) is tuple else -1 for x in X]
            d = map(tuple, d.T)
            d = tuple(d)
            d = d + (d[-1],)
            a = tuple(zip(d, X))
            eachtimepoint_Strategy = para_PAR_index + Strategy_text + \
                ', '.join(', '.join(str(y).replace(' ', '') for y in x) if type(x) is tuple else str(x) for x in a)
        else:
            eachtimepoint_Strategy = ''
    else:
        survival_Strategy = non_report
        dose_Strategy = para_PAR_index + Strategy_text + non_report
        pop_Strategy = para_PAR_index + Strategy_text + non_report
        eachtimepoint_Strategy = para_PAR_index + Strategy_text + non_report

    return {'survival': survival_Strategy, 'dose': dose_Strategy, 'pop': pop_Strategy, 'eachtimepoint':  eachtimepoint_Strategy}


# Check whether the parameter satisfies the criterisa in 2 drug case.#
def DPM_check_criterisa_2drug(PAR_criterisa, PAR, Simduration, Limit_mortality):
    Num_criterisa = 9
    criterisa = [True] * Num_criterisa

    Num_drug, g0_S, g0_R1, g0_R2, g0_R12 = [PAR.get(key) for key in ['Num_drug', 'g0_S', 'g0_R1', 'g0_R2', 'g0_R12']]

    # Criterion(0): Drug potency on S cells must satisfy drug 1 >= drug 2.#
    criterisa_0 = (0 not in PAR_criterisa) or (PAR['Sa.S.D1.'] >= PAR['Sa.S.D2.'])
    criterisa[0] = criterisa_0
    if not criterisa_0:
        return criterisa

    # Criterion PNAS(1): S, R1, R2, R12 cell populations must be >= 0.#
    criterisa_1 = (1 not in PAR_criterisa) or all(i >= 0 for i in [PAR['Spop'], PAR['R1pop'], PAR['R2pop'], PAR['R12pop']])
    criterisa[1] = criterisa_1
    if not criterisa_1:
        return criterisa

    # Criterion PNAS(2): Sa(R1,D1) <= g0_R1 and Sa(R2,D2) <= g0_R2.#
    criterisa_2 = (2 not in PAR_criterisa) or (PAR['Sa.R1.D1.'] <= PAR['g0_R1'] and PAR['Sa.R2.D2.'] <= PAR['g0_R2'])
    criterisa[2] = criterisa_2
    if not criterisa_2:
        return criterisa

    # Criterion PNAS(3): Sa(S,D1) > 0 and Sa(S,D2) > 0.#
    criterisa_3 = (3 not in PAR_criterisa) or (PAR['Sa.S.D1.'] > 0 and PAR['Sa.S.D2.'] > 0)
    criterisa[3] = criterisa_3
    if not criterisa_3:
        return criterisa

    # Criterion Biology Direct(4): Sa(S,D1) > g0_S and Sa(S,D2) > g0_S must both hold.#
    criterisa_4 = (4 not in PAR_criterisa) or (PAR['Sa.S.D1.'] > PAR['g0_S'] and PAR['Sa.S.D2.'] > PAR['g0_S'])
    criterisa[4] = criterisa_4
    if not criterisa_4:
        return criterisa

    # Criterion Biology Direct(5): Sa(R12,D1) < g0_R12 and Sa(R12,D2) < g0_R12 must both hold.#
    criterisa_5 = (5 not in PAR_criterisa) or (PAR['Sa.R12.D1.'] < PAR['g0_R12'] and PAR['Sa.R12.D2.'] < PAR['g0_R12'])
    criterisa[5] = criterisa_5
    if not criterisa_5:
        return criterisa

    # Criterion Biology Direct(6): The patient can be cured by simultaneous full doses of all drugs#
    # (an invalid option due to toxicity) but cannot be cured by any valid static treatment.#
    if 6 in PAR_criterisa:
        LSsim = True
        X0, g0, Sa, T = DPM_generate_X0_2drug(PAR), DPM_generate_g0_2drug(PAR), DPM_generate_Sa_2drug(PAR), DPM_generate_T_2drug(PAR)
        t = np.arange(0, Simduration+SIMTIMESTEP_DEFAULT_VAL, SIMTIMESTEP_DEFAULT_VAL)
        d_fulldose = np.full([Num_drug, t.shape[0]-1], 1, dtype=float)
        X_fulldose, t_fulldose, _ = DPM_sim_model2012(X0, T, g0, Sa, d_fulldose, t, Limit_mortality, LSsim)
        if not np.all(X_fulldose[:, -1] < 1):
            criterisa[6] = False
            return criterisa
        else:
            d_drug1fulldose = np.tile(np.array([[1, 0]]).T, (1, t.shape[0]-1))
            X_drug1fulldose, t_drug1fulldose, _ = DPM_sim_model2012(X0, T, g0, Sa, d_drug1fulldose, t, Limit_mortality, LSsim)
            if not np.any(X_drug1fulldose[:, -1] > 1):
                criterisa[6] = False
                return criterisa
            else:
                d_drug2fulldose = np.tile(np.array([[0, 1]]).T, (1, t.shape[0]-1))
                X_drug2fulldose, t_drug2fulldose, _ = DPM_sim_model2012(X0, T, g0, Sa, d_drug2fulldose, t, Limit_mortality, LSsim)
                if not np.any(X_drug2fulldose[:, -1] > 1):
                    criterisa[6] = False
                    return criterisa
                else:
                    d_halfdose = np.tile(np.array([[1/2, 1/2]]).T, (1, t.shape[0]-1))
                    X_halfdose, t_halfdose, _ = DPM_sim_model2012(X0, T, g0, Sa, d_halfdose, t, Limit_mortality, LSsim)
                    if not np.any(X_halfdose[:, -1] > 1):
                        criterisa[6] = False
                        return criterisa

    # Criterion PNAS(7): A full dose of either drug1 and drug2 must produce at least a 25% increase in progression-free survival compared to#
    # no treatment.#
    if 7 in PAR_criterisa:
        LSsim = True
        X0, g0, Sa, T = DPM_generate_X0_2drug(PAR), DPM_generate_g0_2drug(PAR), DPM_generate_Sa_2drug(PAR), DPM_generate_T_2drug(PAR)
        t = np.arange(0, Simduration+SIMTIMESTEP_DEFAULT_VAL, SIMTIMESTEP_DEFAULT_VAL)
        # No treatment.#
        d_zerodose = np.zeros([Num_drug, t.shape[0]-1])
        _, t_notreat, _ = DPM_sim_model(X0, T, g0, Sa, d_zerodose, t, Limit_mortality, LSsim)
        if t_notreat[-1] > Simduration/1.25:
            criterisa[7] = False
            return criterisa
        else:
            # Full dose of drug 1.#
            d_fulldosedrug1 = np.tile(np.array([[1, 0]]).T, (1, t.shape[0] - 1))
            _, t_fulldosedrug1, _ = DPM_sim_model(X0, T, g0, Sa, d_fulldosedrug1, t, Limit_mortality, LSsim)
            if t_fulldosedrug1[-1] < 1.25 * t_notreat[-1]:
                criterisa[7] = False
                return criterisa
            else:
                # Full dose of drug 2.#
                d_fulldosedrug2 = np.tile(np.array([[0, 1]]).T, (1, t.shape[0] - 1))
                _, t_fulldosedrug2, _ = DPM_sim_model(X0, T, g0, Sa, d_fulldosedrug2, t, Limit_mortality, LSsim)
                if t_fulldosedrug2[-1] < 1.25 * t_notreat[-1]:
                    criterisa[7] = False
                    return criterisa

    # Criterion PNAS(8): the CPM survival time should be less than 75% of the total simulation time to allow evaluation of a potential 25%#
    # relative improvement by other strategies.#
    if 8 in PAR_criterisa:
        LSsim = True
        t_CPM, X_CPM, d_CPM = strategy_CPM(PAR, Simduration, STEPSIZE_DEFAULT_VAL, SIMTIMESTEP_DEFAULT_VAL, Limit_mortality,
                                           LIMIT_RADIOLOGICDETECTION_DEFAULT_VAL, LSsim, False, False, None)
        if (t_CPM[-1] > Simduration/1.25) or (not all(X_CPM[:, -1] < 1)):
            criterisa[8] = False

    return criterisa


# Check whether the parameter satisfies the criteria for the 3 drug case.#
def DPM_check_criterisa_3drug(PAR_criterisa, PAR, Simduration, Limit_mortality):
    Num_criterisa = 4
    criterisa = [True] * Num_criterisa

    # Criterisa(0): Drug potency on S cells must satisfy drug 1 >= drug 2 >= drug 3.#
    criterisa_0 = (0 not in PAR_criterisa) or (PAR['Sa.S.D1.'] >= PAR['Sa.S.D2.'] >= PAR['Sa.S.D3.'])
    criterisa[0] = criterisa_0
    if not criterisa_0:
        return criterisa

    # Criterisa(1): S, R1, R2, R3, R12, R23, R12, R123 cell populations must be >= 0.#
    criterisa_1 = (1 not in PAR_criterisa) or all(i >= 0 for i in [PAR['Spop'], PAR['R1pop'], PAR['R2pop'], PAR['R3pop'], PAR['R12pop'],
                                                                   PAR['R13pop'], PAR['R23pop'], PAR['R123pop']])
    criterisa[1] = criterisa_1
    if not criterisa_1:
        return criterisa

    # Criterisa (2): a full dose of drug1, drug2 or drug3 must produce at least a 25% increase in progression-free survival compared with no
    # treatment.#
    if 2 in PAR_criterisa:
        LSsim = True
        X0, g0, Sa, T = DPM_generate_X0_3drug(PAR), DPM_generate_g0_3drug(PAR), DPM_generate_Sa_3drug(PAR), DPM_generate_T_3drug(PAR)
        t = np.arange(0, Simduration + SIMTIMESTEP_DEFAULT_VAL, SIMTIMESTEP_DEFAULT_VAL)
        # No treatment.#
        d_zerodose = np.zeros([PAR['Num_drug'], t.shape[0] - 1])
        _, t_notreat, _ = DPM_sim_model(X0, T, g0, Sa, d_zerodose, t, Limit_mortality, LSsim)
        if t_notreat[-1] > Simduration/1.25:
            criterisa_2 = False
            criterisa[2] = criterisa_2
            return criterisa
        else:
            # Full dose of drug 1.#
            d_fulldosedrug1 = np.tile(np.array([[1, 0, 0]]).T, (1, t.shape[0] - 1))
            _, t_fulldosedrug1, _ = DPM_sim_model(X0, T, g0, Sa, d_fulldosedrug1, t, Limit_mortality, LSsim)
            if t_fulldosedrug1[-1] < 1.25 * t_notreat[-1]:
                criterisa_2 = False
                criterisa[2] = criterisa_2
                return criterisa
            else:
                # Full dose of drug 2.#
                d_fulldosedrug2 = np.tile(np.array([[0, 1, 0]]).T, (1, t.shape[0] - 1))
                _, t_fulldosedrug2, _ = DPM_sim_model(X0, T, g0, Sa, d_fulldosedrug2, t, Limit_mortality, LSsim)
                if t_fulldosedrug2[-1] < 1.25 * t_notreat[-1]:
                    criterisa_2 = False
                    criterisa[2] = criterisa_2
                    return criterisa
                else:
                    # Full dose of drug 3.#
                    d_fulldosedrug3 = np.tile(np.array([[0, 0, 1]]).T, (1, t.shape[0] - 1))
                    _, t_fulldosedrug3, _ = DPM_sim_model(X0, T, g0, Sa, d_fulldosedrug3, t, Limit_mortality, LSsim)
                    if t_fulldosedrug3[-1] < 1.25 * t_notreat[-1]:
                        criterisa_2 = False
                        criterisa[2] = criterisa_2
                        return criterisa
    return criterisa
