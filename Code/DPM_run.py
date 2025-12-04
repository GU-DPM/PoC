from DPM_lib import inspect, Counter, time, ThreadPool, multiprocessing, Parallel, delayed, tqdm, plt, deepcopy, bz2, pickle
from DPM_read_save import *
from DPM_strategy import *
from DPM_plot import DPM_plot_all, DPM_plot_1strategy, DPM_plot_allstrategy
from DPM_miscellaneous import DPM_miscellaneous_allequal, DPM_miscellaneous_slicedosage
from DPM_assign_check import *
from DPM_analysis import *
'''
CHECKED
'''


def DPM_run_generate_misDrug2(pathload_='./pnas/', pathsave_='./pnas_div5/', ratio=1, filename_pattern_='_para.csv'):
    # Generate the misspecified potency parameters for Drug 2.#
    file_list = os.listdir(pathload_)
    file_list = [filename for filename in file_list if re.search(filename_pattern_, filename) is not None]
    try:
        file_list.sort(key=lambda i: int(i.split('_')[-3]))
    except (ValueError, IndexError):
        pass
    file_list = [os.path.join(pathload_, filename) for filename in file_list]

    with tqdm(total=len(file_list), ncols=150) as pbar:
        for i_filename in file_list:
            par_csv = pd.read_csv(i_filename)
            par_csv['Sa.S.D2.'] = par_csv['Sa.S.D2.'] * ratio
            par_csv['Sa.R1.D2.'] = par_csv['Sa.R1.D2.'] * ratio
            par_csv['Sa.R2.D2.'] = par_csv['Sa.R2.D2.'] * ratio
            par_csv['Sa.R12.D2.'] = par_csv['Sa.R12.D2.'] * ratio

            par_csv['paramID'] = par_csv['paramID'].astype('int')

            i_filename_save = i_filename[i_filename.rfind('/') + 1:]
            if ratio >= 1:
                i_filename_save = f'mis_{int(ratio)}x_' + i_filename_save
            else:
                i_filename_save = f'mis_div{int(1/ratio)}_' + i_filename_save
            pathsave_ = os.path.abspath(pathsave_)
            if not os.path.exists(pathsave_):
                os.makedirs(pathsave_)
            i_filename_save = os.path.join(pathsave_, i_filename_save)
            par_csv.to_csv(i_filename_save, index=False)
            pbar.update(1)
    return


def DPM_run_par_csv_folder(*argv, **kwargs):
    print_quotation = True
    parname_list = PARNAME_LIST_RUN_PAR_CSV_FOLDER

    # Passing a non-keyword argument will result in an error.#
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return
    kwargs = {k: v for k, v in kwargs.items() if v is not None}

    # If any input keywords are not in the list of allowed argument names, raise an error.#
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    # Delete any key in kwargs with a value of None.#
    kwargs = {k.lower(): v for k, v in kwargs.items() if v is not None}

    par = DPM_check_inputpar(kwargs, parname_list)
    pathload_ = par['pathload']
    pathsave_ = par['pathsave']
    filename_pattern_ = par['filename_pattern']
    Num_drug = par['Num_drug']
    dose_method = par['dose_method']
    dose_combination = par['dose_combination']
    dose_interval = par['dose_interval']
    use_input_dose_combination_only = par['use_input_dose_combination_only']
    Simduration = par['Simduration']
    Limit_moleculardetection = par['Limit_moleculardetection']
    Limit_mortality = par['Limit_mortality']
    Limit_radiologicdetection = par['Limit_radiologicdetection']
    Stepsize = par['Stepsize']
    run_sim = par['run_sim']
    PAR_criterisa = par['PAR_criterisa']
    Strategy_name = par['Strategy_name']
    save_filename_param = par['save_filename_param']
    save_filename_stopt = par['save_filename_stopt']
    save_filename_dosage = par['save_filename_dosage']
    save_filename_pop = par['save_filename_pop']
    save_filename_eachtimepoint = par['save_filename_eachtimepoint']
    misspecification_sim_only = par['misspecification_sim_only']
    misspecification_atdecision = par['misspecification_atdecision']
    misspecification_atsim = par['misspecification_atsim']
    misspecification_fileload = par['misspecification_fileload']
    use_parallel = par['use_parallel']

    file_list = os.listdir(pathload_)
    file_list = [filename for filename in file_list if re.search(filename_pattern_, filename) is not None]
    if not file_list:
        print('No default dataset found in the current directory:')  #
        print(Colored.BOLD + Colored.PURPLE + pathload_ + Colored.END)
        DPM_print_errorandclose()
        sys.exit()
    try:
        file_list.sort(key=lambda x: int(x.split('_')[-3]))
    except (ValueError, IndexError):
        pass

    if misspecification_fileload:
        misfile_list_sort = []
        misfile_list = os.listdir(misspecification_fileload)
        for i_file in file_list:
            i_file = '_' + i_file.split('_')[4] + '_'
            i_misfile_list = [filename for filename in misfile_list if re.search(i_file, filename) is not None]
            assert len(i_misfile_list) == 1
            misfile_list_sort.append(os.path.join(misspecification_fileload, i_misfile_list[0]))
    else:
        misfile_list_sort = [''] * len(file_list)

    file_list = [os.path.join(pathload_, filename) for filename in file_list]

    # file_list = file_list[212:]
    # misfile_list_sort = misfile_list_sort[212:]
    if use_parallel and len(file_list) > 1:
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(DPM_run_par_csv)
                                   (filename_csv=i_filename,
                                    misspecification_filename_csv=misfile_list_sort[i],
                                    misspecification_atdecision=misspecification_atdecision,
                                    misspecification_atsim=misspecification_atsim,
                                    Num_drug=Num_drug,
                                    dose_method=dose_method,
                                    dose_combination=dose_combination,
                                    dose_interval=dose_interval,
                                    use_input_dose_combination_only=use_input_dose_combination_only,
                                    Simduration=Simduration,
                                    Limit_moleculardetection=Limit_moleculardetection,
                                    Limit_mortality=Limit_mortality,
                                    Limit_radiologicdetection=Limit_radiologicdetection,
                                    Stepsize=Stepsize,
                                    run_sim=run_sim,
                                    PAR_criterisa=PAR_criterisa,
                                    Strategy_name=Strategy_name,
                                    save_filename_param=save_filename_param,
                                    save_filename_stopt=save_filename_stopt,
                                    save_filename_pop=save_filename_pop,
                                    save_filename_dosage=save_filename_dosage,
                                    save_filename_eachtimepoint=save_filename_eachtimepoint,
                                    misspecification_sim_only=misspecification_sim_only,
                                    pathsave=pathsave_)
                                   for i, i_filename in enumerate(file_list))
        print('Finished.')
    else:
        count_file = 1
        for i, i_filename in enumerate(file_list):
            print('Processing the ' + Colored.BOLD + Colored.PURPLE + str(count_file) + 'th' + Colored.END +
                  ' file of ' + Colored.BOLD + Colored.PURPLE + str(len(file_list)) + Colored.END + ' total files...')
            DPM_run_par_csv(filename_csv=i_filename,
                            misspecification_filename_csv=misfile_list_sort[i],
                            misspecification_atdecision=misspecification_atdecision,
                            misspecification_atsim=misspecification_atsim,
                            Num_drug=Num_drug,
                            dose_method=dose_method,
                            dose_combination=dose_combination,
                            dose_interval=dose_interval,
                            use_input_dose_combination_only=use_input_dose_combination_only,
                            Simduration=Simduration,
                            Limit_moleculardetection=Limit_moleculardetection,
                            Limit_mortality=Limit_mortality,
                            Limit_radiologicdetection=Limit_radiologicdetection,
                            Stepsize=Stepsize,
                            run_sim=run_sim,
                            PAR_criterisa=PAR_criterisa,
                            Strategy_name=Strategy_name,
                            save_filename_param=save_filename_param,
                            save_filename_stopt=save_filename_stopt,
                            save_filename_pop=save_filename_pop,
                            save_filename_dosage=save_filename_dosage,
                            save_filename_eachtimepoint=save_filename_eachtimepoint,
                            misspecification_sim_only=misspecification_sim_only,
                            pathsave=pathsave_)
            count_file += 1
        print('Finished.')
    return


# Run simulation using parameters from a CSV file.#
def DPM_run_par_csv(*argv, **kwargs):
    def DPM_run_par_csv_1(filename_):
        # Delete the string '.csv' from the filename.#
        filename_ = filename_.replace('.csv', '')
        if sys.platform in ['darwin', 'linux']:
            filename_ = filename_[filename_.rfind('/') + 1:]
        else:
            filename_ = filename_[filename_.rfind('\\') + 1:]
        filename_ = filename_[:MAXFILENAMELEN] if len(filename_) >= MAXFILENAMELEN else filename_
        return filename_

    print_quotation = True
    LSsim = True
    parname_list = PARNAME_LIST_RUN_PAR_CSV

    # Raise an error if a non-keyword argument is provided.#
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return
    kwargs = {k: v for k, v in kwargs.items() if v is not None}

    # If the input keywords are not in the list of allowed argument names, raise an error.#
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    # Remove any keys from 'kwargs' whose values are None.#
    kwargs = {k.lower(): v for k, v in kwargs.items() if v is not None}

    par = DPM_check_inputpar(kwargs, parname_list)
    Num_drug = par['Num_drug']
    par_ind = par['par_ind']
    dose_combination = par['dose_combination']
    Simduration = par['Simduration']
    Limit_mortality = par['Limit_mortality']
    Stepsize = par['Stepsize']
    Limit_radiologicdetection = par['Limit_radiologicdetection']
    run_sim = par['run_sim']
    Strategy_name = par['Strategy_name']
    save_filename_param = par['save_filename_param']
    save_filename_stopt = par['save_filename_stopt']
    save_filename_pop = par['save_filename_pop']
    save_filename_dosage = par['save_filename_dosage']
    save_filename_eachtimepoint = par['save_filename_eachtimepoint']
    pathsave = par['pathsave']
    filename_csv = par['filename_csv']
    misspecification_filename_csv = par['misspecification_filename_csv']
    misspecification_atsim = par['misspecification_atsim']
    misspecification_atdecision = par['misspecification_atdecision']
    misspecification_sim_only = par['misspecification_sim_only']

    print(filename_csv)
    par_csv = DPM_read_par_csv(filename_csv, 'filename_csv')
    par_csv = DPM_assign_par_csvread(par_csv, Num_drug)
    if not par_csv:
        filename_csv_text = Colored.BOLD + Colored.PURPLE + 'filename_csv' + Colored.END
        print('No avaiable parameter from keyword input ' + filename_csv_text + '.')
        DPM_print_errorandclose()
        return

    # Check for mis-specification files.#
    mis_par_csv = None
    if misspecification_filename_csv:
        mis_par_csv = DPM_read_par_csv(misspecification_filename_csv, 'misspecification_filename_csv')
        mis_par_csv = DPM_assign_par_csvread(mis_par_csv, Num_drug)
        if not mis_par_csv:
            misspecification_filename_csv_text = Colored.BOLD + Colored.PURPLE + 'misspecification_filename_csv' + Colored.END
            print('No avaiable parameter from keyword input ' + misspecification_filename_csv_text + '.')
            DPM_print_errorandclose()
            return

    # Verify that the number of misspecified parameters matches the number of true parameters.#
    if 'misspecification_filename_csv' in kwargs.keys() and kwargs['misspecification_filename_csv'] != '':
        if len(par_csv) != len(mis_par_csv):
            filename_csv_text = Colored.BOLD + Colored.PURPLE + filename_csv + Colored.END
            misspecification_filename_csv_text = Colored.BOLD + Colored.PURPLE + misspecification_filename_csv + Colored.END
            print('The number of parameters ' + filename_csv_text +
                  ' is not the same as the number of parameters ' + misspecification_filename_csv_text + '.')
            DPM_print_errorandclose()
            return

    mis_specification_csv = True if mis_par_csv else False
    if type(filename_csv) is list:
        filename_csv = [i_file.replace('.csv', '') for i_file in filename_csv]
        if sys.platform in ['darwin', 'linux']:
            filename_csv = [i_file[i_file.rfind('/')+1:] for i_file in filename_csv]
        else:
            filename_csv = [i_file[i_file.rfind('\\') + 1:] for i_file in filename_csv]
        filename = '__'.join(filename_csv)
    else:
        filename = DPM_run_par_csv_1(filename_csv)

    if len(par_ind) > 0 and par_ind != [-1]:
        par_ind_accept = [x for x in par_ind if 0 <= x < len(par_csv)]
        if len(par_ind_accept) == 0:
            print('No avaiable index in the inputted ' + Colored.BOLD + Colored.PURPLE + 'par_ind' + Colored.END + '.')
            DPM_print_errorandclose()
            sys.exit()
        par_csv = [par_csv[i] for i in par_ind_accept]
        mis_par_csv = [mis_par_csv[i] for i in par_ind_accept] if mis_specification_csv else None

    if len(par_ind) > 0:
        if len(filename) >= MAXFILENAMELEN:
            filename = filename[:MAXFILENAMELEN]

    if par_csv and run_sim:
        timenow = datetime.now()
        timenow = timenow.strftime('%Y%m%d')
        filename_param = filename + '_result_para_' + timenow + '.csv' if save_filename_param else ''
        filename_stopt = filename + '_result_stopt_' + timenow + '.csv' if save_filename_stopt else ''
        filename_dosage = filename + '_result_dosage_' + timenow + '.csv' if save_filename_dosage else ''
        filename_pop = filename + '_result_pop_' + timenow + '.csv' if save_filename_pop else ''
        filename_eachtimepoint = filename + '_result_eachtimepoint_' + timenow + '.csv' if save_filename_eachtimepoint else ''

        mis_filename_param = os.path.join(pathsave, 'mis_' + filename_param) if (mis_specification_csv and save_filename_param) else ''
        mis_filename_stopt = os.path.join(pathsave, 'mis_' + filename_stopt) if (mis_specification_csv and save_filename_stopt) else ''
        mis_filename_dosage = os.path.join(pathsave, 'mis_' + filename_dosage) if (mis_specification_csv and save_filename_dosage) else ''
        mis_filename_pop = os.path.join(pathsave, 'mis_' + filename_pop) if (mis_specification_csv and save_filename_pop) else ''
        mis_filename_eachtimepoint = os.path.join(pathsave, 'mis_' + filename_eachtimepoint) if \
            (mis_specification_csv and filename_eachtimepoint) else ''

        filename_param = os.path.join(pathsave, filename_param) if save_filename_param else ''
        filename_stopt = os.path.join(pathsave, filename_stopt) if save_filename_stopt else ''
        filename_dosage = os.path.join(pathsave, filename_dosage) if save_filename_dosage else ''
        filename_pop = os.path.join(pathsave, filename_pop) if save_filename_pop else ''
        filename_eachtimepoint = os.path.join(pathsave, filename_eachtimepoint) if save_filename_eachtimepoint else ''

        Header_param_csv, Header_stopt_csv, Header_dosage_csv, Header_pop_csv, Header_eachtimepoint_csv = \
            DPM_generate_header_csv(Num_drug, Simduration, Stepsize)

        filename_zip = zip((save_filename_param, save_filename_stopt, save_filename_dosage, save_filename_pop, save_filename_eachtimepoint),
                           (filename_param, filename_stopt, filename_dosage, filename_pop, filename_eachtimepoint),
                           (Header_param_csv, Header_stopt_csv, Header_dosage_csv, Header_pop_csv, Header_eachtimepoint_csv))
        mis_filename_zip = zip((save_filename_param, save_filename_stopt, save_filename_dosage, save_filename_pop, save_filename_eachtimepoint),
                               (mis_filename_param, mis_filename_stopt, mis_filename_dosage, mis_filename_pop, mis_filename_eachtimepoint),
                               (Header_param_csv, Header_stopt_csv, Header_dosage_csv, Header_pop_csv, Header_eachtimepoint_csv))
        # If running only the mis-specification file and clearing prior results, create empty replacement files.#
        not mis_specification_csv and DPM_read_erase(filename_zip)
        # If run mis-specification file, create the empty files.#
        mis_specification_csv and DPM_read_erase(mis_filename_zip)

        # par_csv = par_csv[180:]
        # mis_par_csv = mis_par_csv[180:]
        with tqdm(total=len(par_csv), ncols=150, desc='Runing simulation using parameters from the .csv input file.') as pbar:  #
            for i in range(len(par_csv)):
                i_par = par_csv[i]
                i_paramID = int(i_par['paramID'])

                if mis_specification_csv:
                    i_mis_par = mis_par_csv[i]
                    i_mis_paramID = int(i_mis_par['paramID'])
                    assert i_paramID == i_mis_paramID

                    i_mis_CPM, i_mis_DPM1, i_mis_DPM2_1, i_mis_DPM2_2, i_mis_DPM3, i_mis_DPM4 = \
                        DPM_run_simulation(i_par, Strategy_name, dose_combination, Simduration, Stepsize, Limit_mortality,
                                           Limit_radiologicdetection, LSsim, misspecification_atdecision, misspecification_atsim, i_mis_par)

                    DPM_save_result_csv(i_mis_par, Stepsize, Simduration, i_mis_CPM, i_mis_DPM1, i_mis_DPM2_1, i_mis_DPM2_2, i_mis_DPM3,
                                        i_mis_DPM4, i_mis_paramID, mis_filename_param, mis_filename_stopt, mis_filename_dosage, mis_filename_pop,
                                        mis_filename_eachtimepoint, Limit_mortality)

                if not misspecification_sim_only or not mis_specification_csv:
                    i_CPM, i_DPM1, i_DPM2_1, i_DPM2_2, i_DPM3, i_DPM4 = \
                        DPM_run_simulation(i_par, Strategy_name, dose_combination, Simduration, Stepsize, Limit_mortality,
                                           Limit_radiologicdetection, LSsim, False, False, None)

                    DPM_save_result_csv(i_par, Stepsize, Simduration, i_CPM, i_DPM1, i_DPM2_1, i_DPM2_2, i_DPM3, i_DPM4, i_paramID,
                                        filename_param, filename_stopt, filename_dosage, filename_pop, filename_eachtimepoint, Limit_mortality)
                pbar.update(1)
    return


# Run simulation for a single parameter and plot the result.#
def DPM_run_plot_1PAR(*argv, **kwargs):
    print_quotation = True
    LSsim = True
    parname_list = PARNAME_LIST_RUN_PLOT_1PAR
    # If a non-keyword argument is provided, raise an error.#
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return

    # If any provided keyword is not in the list of allowed argument names, raise an error.#
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    kwargs = {k.lower(): v for k, v in kwargs.items()}

    # The caller of the 'DPM_check_inputpar' function should be 'DPM_run_plot_1PAR', the current function name.#
    callername = inspect.stack()[0][3]
    par = DPM_check_inputpar(kwargs, parname_list, callername)
    Num_drug = par['Num_drug']
    i_par = par['par']
    Limit_radiologicdetection = par['Limit_radiologicdetection']
    Limit_moleculardetection = par['Limit_moleculardetection']
    dose_combination = par['dose_combination']
    Simduration = par['Simduration']
    Limit_mortality = par['Limit_mortality']
    Stepsize = par['Stepsize']
    PAR_criterisa = par['PAR_criterisa']
    Strategy_name = par['Strategy_name']
    save_filename_param = par['save_filename_param']
    save_filename_stopt = par['save_filename_param']
    save_filename_pop = par['save_filename_pop']
    save_filename_dosage = par['save_filename_dosage']
    save_filename_eachtimepoint = par['save_filename_eachtimepoint']
    pathsave = par['pathsave']
    savename = par['savename']

    criterisa = DPM_check_criterisa_2drug(PAR_criterisa, i_par, Simduration, Limit_mortality) if Num_drug == 2 else \
        DPM_check_criterisa_3drug(PAR_criterisa, i_par, Simduration, Limit_mortality) if Num_drug == 3 else False
    X0total = i_par['Spop'] + i_par['R1pop'] + i_par['R2pop'] + i_par['R12pop']
    par_check_index = DPM_check_par(X0total, Limit_mortality, Simduration, Stepsize)
    if not (all(criterisa) and not par_check_index):
        print('The input ' + Colored.BOLD + Colored.PURPLE + 'par' + Colored.END + " doesn't satisfy the required criteria.")
        DPM_print_errorandclose()
        sys.exit()

    i_par['Num_cell_type'] = 4
    i_par['paramID'] = 1 if 'paramID' not in i_par.keys() else i_par['paramID']

    filename_param = savename + '_result_para' + '.csv' if save_filename_param else ''
    filename_stopt = savename + '_result_stopt' + '.csv' if save_filename_stopt else ''
    filename_dosage = savename + '_result_dosage' + '.csv' if save_filename_dosage else ''
    filename_pop = savename + '_result_pop' + '.csv' if save_filename_stopt else ''
    filename_eachtimepoint = savename + '_result_eachtimepoint' + '.csv' if save_filename_eachtimepoint else ''

    filename_param = os.path.join(pathsave, filename_param) if save_filename_param else ''
    filename_stopt = os.path.join(pathsave, filename_stopt) if save_filename_stopt else ''
    filename_dosage = os.path.join(pathsave, filename_dosage) if save_filename_dosage else ''
    filename_pop = os.path.join(pathsave, filename_pop) if save_filename_pop else ''
    filename_eachtimepoint = os.path.join(pathsave, filename_eachtimepoint) if save_filename_eachtimepoint else ''

    Header_param_csv, Header_stopt_csv, Header_dosage_csv, Header_pop_csv, Header_eachtimepoint_csv = \
        DPM_generate_header_csv(Num_drug, Simduration, Stepsize)

    filename_zip = zip((save_filename_param, save_filename_stopt, save_filename_dosage, save_filename_pop, save_filename_eachtimepoint),
                       (filename_param, filename_stopt, filename_dosage, filename_pop, filename_eachtimepoint),
                       (Header_param_csv, Header_stopt_csv, Header_dosage_csv, Header_pop_csv, Header_eachtimepoint_csv))

    DPM_read_erase(filename_zip)

    print(Colored.BOLD + Colored.PURPLE + 'Running...' + Colored.END)
    i_CPM, i_DPM1, i_DPM2_1, i_DPM2_2, i_DPM3, i_DPM4 = DPM_run_simulation(i_par, Strategy_name, dose_combination, Simduration, Stepsize,
                                                                           Limit_mortality, Limit_radiologicdetection, LSsim, False, False, None)
    i_paramID = 1
    DPM_save_result_csv(i_par, Stepsize, Simduration, i_CPM, i_DPM1, i_DPM2_1, i_DPM2_2, i_DPM3, i_DPM4, i_paramID, filename_param,
                        filename_stopt, filename_dosage, filename_pop, filename_eachtimepoint, Limit_mortality)

    strategy = {'CPM': i_CPM, 'DPM1': i_DPM1, 'DPM2.1': i_DPM2_1, 'DPM2.2': i_DPM2_2, 'DPM3': i_DPM3, 'DPM4': i_DPM4}
    DPM_plot_all(Num_drug, strategy, Simduration, Limit_moleculardetection, Limit_radiologicdetection, Limit_mortality, pathsave)
    return


# Simulate the model for one parameter value and plot the drug misspecification case used in Deepak's paper.#
def DPM_run_plot_1PAR_drugmis(*argv, **kwargs):
    print_quotation = True
    LSsim = True
    parname_list = PARNAME_LIST_RUN_PLOT_1PAR
    # Raise an error is a non-keyword argument is provided.#
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return

    # Raise an error if any provided keyword is not in the list of allowed argument names.#
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    kwargs = {k.lower(): v for k, v in kwargs.items()}

    par = DPM_check_inputpar(kwargs, parname_list)

    i_par = par['par']
    i_mis_par = par['mis_par']
    Num_drug = par['Num_drug']
    Limit_radiologicdetection = par['Limit_radiologicdetection']
    Limit_moleculardetection = par['Limit_moleculardetection']
    dose_combination = par['dose_combination']
    Simduration = par['Simduration']
    Limit_mortality = par['Limit_mortality']
    Stepsize = par['Stepsize']
    Strategy_name = par['Strategy_name']
    pathsave = par['pathsave']
    misspecification_atdecision = par['misspecification_atdecision']
    misspecification_atsim = par['misspecification_atsim']

    i_par['Num_cell_type'] = 4
    i_mis_par['Num_cell_type'] = i_par['Num_cell_type']

    i_CPM, i_DPM1, i_DPM2_1, i_DPM2_2, i_DPM3, i_DPM4 = DPM_run_simulation(i_par, Strategy_name, dose_combination, Simduration, Stepsize,
                                                                           Limit_mortality, Limit_radiologicdetection, LSsim, False, False, None)

    i_mis_CPM, i_mis_DPM1, i_mis_DPM2_1, i_mis_DPM2_2, i_mis_DPM3, i_mis_DPM4 = \
        DPM_run_simulation(i_par, Strategy_name, dose_combination, Simduration, Stepsize, Limit_mortality, Limit_radiologicdetection, LSsim,
                           misspecification_atdecision, misspecification_atsim, i_mis_par)

    DPM_plot_1strategy(Num_drug, i_CPM, Simduration, Limit_moleculardetection, Limit_radiologicdetection, Limit_mortality, 'CPM', pathsave)
    plt.savefig(os.path.join(pathsave, 'CPM.pdf'), dpi=FIG_DPI)
    plt.close()

    DPM_plot_1strategy(Num_drug, i_DPM2_2, Simduration, Limit_moleculardetection, Limit_radiologicdetection, Limit_mortality, 'DPM2', pathsave)
    plt.savefig(os.path.join(pathsave, 'DPM2.pdf'), dpi=FIG_DPI)
    plt.close()

    DPM_plot_1strategy(Num_drug, i_mis_CPM, Simduration, Limit_moleculardetection, Limit_radiologicdetection, Limit_mortality, 'CPM mis',
                       pathsave)
    plt.savefig(os.path.join(pathsave, 'CPM mis.pdf'), dpi=FIG_DPI)
    plt.close()

    DPM_plot_1strategy(Num_drug, i_mis_DPM2_2, Simduration, Limit_moleculardetection, Limit_radiologicdetection,
                       Limit_mortality, 'DPM2 mis', pathsave)
    plt.savefig(os.path.join(pathsave, 'DPM2 mis.pdf'), dpi=FIG_DPI)
    plt.close()

    X_total = dict()
    X_total['CPM'] = (i_CPM[0], i_CPM[1].sum(axis=0))
    X_total['CPM mis'] = (i_mis_CPM[0], i_mis_CPM[1].sum(axis=0))
    X_total['DPM2'] = (i_DPM2_2[0], i_DPM2_2[1].sum(axis=0))
    X_total['DPM2 mis'] = (i_mis_DPM2_2[0], i_mis_DPM2_2[1].sum(axis=0))
    DPM_plot_allstrategy(X_total, Simduration, Limit_radiologicdetection, Limit_mortality, pathsave)
    plt.savefig(os.path.join(pathsave, 'total.pdf'), dpi=FIG_DPI)
    plt.close()

    X_R12 = dict()
    X_R12['CPM'] = (i_CPM[0], i_CPM[1][-1, :])
    X_R12['CPM mis'] = (i_mis_CPM[0], i_mis_CPM[1][-1, :])
    X_R12['DPM2'] = (i_DPM2_2[0], i_DPM2_2[1][-1, :])
    X_R12['DPM2 mis'] = (i_mis_DPM2_2[0], i_mis_DPM2_2[1][-1, :])
    DPM_plot_allstrategy(X_R12, Simduration, Limit_radiologicdetection, Limit_mortality, pathsave)
    plt.savefig(os.path.join(pathsave, 'R12.pdf'), dpi=FIG_DPI)
    plt.close()
    return


# Run different strategy simulations.#
def DPM_run_simulation(PAR, Strategy_name, dose_combination, Simduration, Stepsize, Limit_mortality, Limit_radiologicdetection, LSsim,
                       misspecification_atdecision, misspecification_atsim, mis_par):
    CPM, DPM1, DPM2_1, DPM2_2, DPM3, DPM4 = None, None, None, None, None, None

    if 'CPM' in Strategy_name:
        t_CPM, X_CPM, d_CPM = strategy_CPM(PAR, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL, Limit_mortality, Limit_radiologicdetection,
                                           LSsim, misspecification_atdecision, misspecification_atsim, mis_par)
        CPM = [t_CPM, X_CPM, d_CPM]

    if 'DPM1' in Strategy_name:
        t_DPM1, X_DPM1, d_DPM1 = strategy_DPM1(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL, Limit_mortality, LSsim,
                                               misspecification_atdecision, mis_par)
        DPM1 = [t_DPM1, X_DPM1, d_DPM1]

    if 'DPM2.1' in Strategy_name:
        DPM2threshold = 1e9
        t_DPM2_1, X_DPM2_1, d_DPM2_1 = strategy_DPM2(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL, Limit_mortality,
                                                     LSsim, DPM2threshold, misspecification_atdecision, misspecification_atsim, mis_par)
        DPM2_1 = [t_DPM2_1, X_DPM2_1, d_DPM2_1]

    if 'DPM2.2' in Strategy_name:
        DPM2threshold = 1e11
        t_DPM2_2, X_DPM2_2, d_DPM2_2 = strategy_DPM2(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL, Limit_mortality,
                                                     LSsim, DPM2threshold, misspecification_atdecision, misspecification_atsim, mis_par)
        DPM2_2 = [t_DPM2_2, X_DPM2_2, d_DPM2_2]

    if 'DPM3' in Strategy_name:
        t_DPM3, X_DPM3, d_DPM3 = strategy_DPM3(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL, Limit_mortality,
                                               LSsim, misspecification_atdecision, mis_par)
        DPM3 = [t_DPM3, X_DPM3, d_DPM3]

    if 'DPM4' in Strategy_name:
        t_DPM4, X_DPM4, d_DPM4 = strategy_DPM4(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL, Limit_mortality,
                                               LSsim, misspecification_atdecision, mis_par)
        DPM4 = [t_DPM4, X_DPM4, d_DPM4]

    return CPM, DPM1, DPM2_1, DPM2_2, DPM3, DPM4


def DPM_run_processing(*argv, **kwargs):
    def DPM_run_processing_1(i_stoptime_2):
        def DPM_run_processing_1_1(ind, pop_strategy_, strategy_):
            i_para_t0 = [i_para[i_key_][ind] for i_key_ in list(para.keys())[1:5]]
            i_para_t0_tuple = tuple(i_para_t0)
            if pop_dict[strategy_][0] is None:
                pop_dict[strategy_][0] = [i_para_t0_tuple]
            else:
                pop_dict[strategy_][0].append(i_para_t0_tuple)
            pop_strategy_ = pop_strategy_.split(';')
            for i_, i_step in enumerate(pop_strategy_):
                if i_step != '-1':
                    i_step = [float(val) for val in re.sub('[()]', '', i_step).split(',')]
                    if pop_dict[strategy_][i_ + 1] is None:
                        pop_dict[strategy_][i_ + 1] = [tuple(i_step)]
                    else:
                        pop_dict[strategy_][i_ + 1].append(tuple(i_step))
            return

        i_filename_pattern = i_stoptime_2.split('_stopt', 1)[0]
        i_file_stoptime_2 = os.path.join(pathload, i_stoptime_2)
        i_stoptime = DPM_read_stoptime_csv(i_file_stoptime_2, Strategy_name)
        i_file_para = i_filename_pattern + '_para'
        i_file_para = [i_file for i_file in file_para if i_file_para in i_file]

        i_para = pd.read_csv(os.path.join(pathload, i_file_para[0]), usecols=para_key).to_dict('list')
        for i_key in para.keys():
            para[i_key].extend(i_para[i_key])

        i_file_dosage, i_file_pop = i_filename_pattern + '_dosage', i_filename_pattern + '_pop'
        i_file_dosage = [i_file for i_file in file_dosage if i_file_dosage in i_file]
        i_file_pop = [i_file for i_file in file_pop if i_file_pop in i_file]
        assert len(i_file_dosage) == 1 and len(i_file_pop) == 1
        i_file_dosage, i_file_pop = os.path.join(pathload, i_file_dosage[0]), os.path.join(pathload, i_file_pop[0])

        i_dosage = DPM_read_dosage_pop_csv(i_file_dosage, Strategy_name)
        i_pop = DPM_read_dosage_pop_csv(i_file_pop, Strategy_name)

        # num_cores = multiprocessing.cpu_count()
        for i_strategy in Strategy_name:
            i_pop_strategy = i_pop[i_strategy]
            for ii, i_pop_strategy_para in enumerate(i_pop_strategy):
                DPM_run_processing_1_1(ii, i_pop_strategy_para, i_strategy)
            # Parallel(n_jobs=num_cores)(delayed(DPM_run_processing_2_1)(ii, i_pop_strategy_para, i_strategy)
            #                            for ii, i_pop_strategy_para in enumerate(i_pop_strategy))

        for i_key in dosage.keys():
            dosage[i_key].extend(i_dosage[i_key])
            stoptime[i_key].extend(i_stoptime[i_key])
            pop[i_key].extend(i_pop[i_key])

        i_para = pd.read_csv(os.path.join(pathload, i_file_para[0]), usecols=para_key)
        i_para['X'] = i_para.loc[:, para_key[1:5]].sum(axis=1)
        i_para['R1pop per'] = i_para['R1pop']/i_para['X']
        i_para['R2pop per'] = i_para['R2pop']/i_para['X']
        i_para['R12pop per'] = i_para['R12pop']/i_para['X']

        pbar.update(1)
        return

    def DPM_run_processing_2(para_i, dosage_i, stoptime_i, num_dosage):
        para_i_df = pd.DataFrame(para_i)
        para_i_df['Xtotal'] = para_i_df['Spop'] + para_i_df['R1pop'] + para_i_df['R2pop'] + para_i_df['R12pop']
        # (1) Subset of patients with R1 >= 7.1e-7 and R2 >= 7.1e-7. #
        bool_R1_R2_e_7 = (para_i_df['R1pop'] / para_i_df['Xtotal'] >= 7.1e-7) & (para_i_df['R2pop'] / para_i_df['Xtotal'] >= 7.1e-7)
        bool_R1_R2_e_7 = bool_R1_R2_e_7.values.tolist()
        # (2) Subset of patients with R1 >= 7.1e-5 and R2 >= 7.1e-5. #
        bool_R1_R2_e_5 = (para_i_df['R1pop'] / para_i_df['Xtotal'] >= 7.1e-5) & (para_i_df['R2pop'] / para_i_df['Xtotal'] >= 7.1e-5)
        bool_R1_R2_e_5 = bool_R1_R2_e_5.values.tolist()

        total = [True for _ in range(len(para_i_df['Xtotal']))]
        del para_i, para_i_df

        dosage_i_sliced = DPM_miscellaneous_slicedosage(dosage_i, num_dosage, Strategy_name)
        bool_diff_dosage = [x != y for x, y in zip(dosage_i_sliced['CPM'], dosage_i_sliced['DPM2.2'])]
        bool_same_dosage = [x == y for x, y in zip(dosage_i_sliced['CPM'], dosage_i_sliced['DPM2.2'])]

        setind_i = dict(zip(setname, [total, bool_diff_dosage, bool_same_dosage, bool_R1_R2_e_7, bool_R1_R2_e_5]))

        km_i = {i_setname: {**{i_strategy: None for i_strategy in Strategy_name}, **{'hz_ratio': [], 'p': [], 'num': 0}} for i_setname in
                setname}
        info_i = {i_setname: {**{i_strategy: None for i_strategy in Strategy_name}, 'paramID sigbetter': None} for i_setname in setname}
        for i_setname in setname:
            i_stoptime_set = {item: list(itertools.compress(value, setind_i[i_setname])) for (item, value) in stoptime_i.items()}
            i_dosage_set = {item: list(itertools.compress(value, setind_i[i_setname])) for (item, value) in dosage_i.items()}
            i_hz, i_p = DPM_analysis_hazard_ratio(i_stoptime_set, Strategy_name, Simduration)
            km_i[i_setname]['hz_ratio'], km_i[i_setname]['p'], km_i[i_setname]['num'] = i_hz, i_p, len(i_stoptime_set['paramID'])
            info_i[i_setname]['paramID sigbetter'] = DPM_analysis_sigbetter(i_stoptime_set, Strategy_name)
            if i_setname == 'total':
                print(f"DPM sigbetter:{len(info_i[i_setname]['paramID sigbetter'])}")
            for i_strategy in Strategy_name:
                i_stoptime_strategy = i_stoptime_set[i_strategy]
                i_dosage_set_strategy = i_dosage_set[i_strategy]
                km_i[i_setname][i_strategy] = DPM_analysis_KM(i_stoptime_strategy, Simduration)
                info_i[i_setname][i_strategy] = DPM_analysis_dose(i_dosage_set_strategy, i_strategy)
        return setind_i, km_i, info_i

    print_quotation = True
    parname_list = PARNAME_LIST_RUN_PROCESSING
    # If a non-keyword argument is provided as input, raise an error.#
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return

    # Raise an error if any sypplied keyword does not appear in the list of permitted argument names.#
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    # Remove any entries from 'kwargs' whose values are None.#
    kwargs = {k.lower(): v for k, v in kwargs.items() if v is not None}

    par = DPM_check_inputpar(kwargs, parname_list)
    Strategy_name = par['Strategy_name']
    pathload = par['pathload']
    Num_stepdiff = par['Num_stepdiff']
    Simduration = par['Simduration']
    Stepsize = par['Stepsize']
    use_parallel = par['use_parallel']

    file_format = '.csv'
    file_list = os.listdir(pathload)
    partial_filename_stoptime = ['result_stopt', file_format]
    file_stoptime = [filename for filename in file_list if all([x in filename for x in partial_filename_stoptime])]
    partial_filename_dosage = ['result_dosage', file_format]
    file_dosage = [filename for filename in file_list if all([x in filename for x in partial_filename_dosage])]
    partial_filename_pop = ['result_pop', file_format]
    file_pop = [filename for filename in file_list if all([x in filename for x in partial_filename_pop])]
    partial_filename_para = ['result_para', file_format]
    file_para = [filename for filename in file_list if all([x in filename for x in partial_filename_para])]

    # If no files are detected in the current directory, terminate the process.#
    if not file_stoptime or not file_dosage or not file_pop:
        print('No stopping time file was located in the current directory:')
        print(Colored.BOLD + Colored.PURPLE + pathload + Colored.END)
        DPM_print_errorandclose()
        sys.exit()
    assert len(file_dosage) == len(file_stoptime) == len(file_pop)

    try:
        filename_i = file_stoptime[0].split('_')
        for i, i_val in enumerate(filename_i):
            try:
                int(i_val)
                break
            except ValueError:
                pass
        file_stoptime.sort(key=lambda x: int(x.split('_')[i]))
        file_dosage.sort(key=lambda x: int(x.split('_')[i]))
        file_pop.sort(key=lambda x: int(x.split('_')[i]))
    except ValueError:
        pass

    # Load the stopping time data from the .csv file.#
    result = {**{'paramID': []}, **{i_strategy: [] for i_strategy in Strategy_name}}
    stoptime, dosage, pop = deepcopy(result), deepcopy(result), deepcopy(result)

    para_key = ['paramID', 'Spop', 'R1pop', 'R2pop', 'R12pop', 'g0_S', 'Sa.S.D1.', 'Sa.S.D2.', 'Sa.R1.D1.', 'Sa.R2.D2.', 'T.R1..S.', 'T.R2..S.']
    para = {i_key: [] for i_key in para_key}

    # result = pd.DataFrame(columns=list(para.keys())[1:5])
    pop_dict = dict(zip(Strategy_name, [None] * len(Strategy_name)))
    for i_strategy in Strategy_name:
        pop_dict[i_strategy] = [[] for _ in range((int(Simduration/Stepsize) + 1))]

    with tqdm(total=len(file_stoptime), ncols=150, desc='Processing results') as pbar:
        st = time()
        if use_parallel:
            with ThreadPool() as pool:
                pool.map(DPM_run_processing_1, file_stoptime)
        else:
            for i, i_file_stoptime in enumerate(file_stoptime):
                try:
                    DPM_run_processing_1(i_file_stoptime)
                except (AssertionError,  ValueError):
                    print()
                    i_file_text = Colored.BOLD + Colored.PURPLE + i_file_stoptime + Colored.END
                    print(i_file_text)
                    raise ValueError('Error.')
        elapsed_time = round((time() - st)/60, 1)
        print('Execution time:', elapsed_time, 'mins')

    para_check = {i_key: [] for i_key in para_key}
    for i, i_file_stoptime in enumerate(file_stoptime):
        i_file_para_c = i_file_stoptime.split('_stopt', 1)[0] + '_para'
        i_file_para_c = [i_file for i_file in file_para if i_file_para_c in i_file]
        i_para_c = pd.read_csv(os.path.join(pathload, i_file_para_c[0]), usecols=para_key).to_dict('list')
        for i_key in para_check.keys():
            para_check[i_key].extend(i_para_c[i_key])

    check_list = [Counter(dosage['paramID']), Counter(stoptime['paramID']), Counter(para['paramID']), Counter(para_check['paramID'])]
    if not DPM_miscellaneous_allequal(check_list):
        raise ValueError('paramID mismatch.')

    setname = ['total', 'first 2 diff', 'first 2 same', 'R1,R2>=1e-7', 'R1,R2>=1e-5']

    setind, km, info = DPM_run_processing_2(para, dosage, stoptime, Num_stepdiff)
    filename_para, filename_stoptime, filename_dosage, filename_pop, filename_popdict, filename_setind, filename_km, filename_info = (
        os.path.join(pathload, 'result_para.pckl'), os.path.join(pathload, 'result_stoptime.pckl'), os.path.join(pathload, 'result_dosage.pckl'),
        os.path.join(pathload, 'result_pop.pckl'), os.path.join(pathload, 'result_popdict.pckl'), os.path.join(pathload, 'result_setind.pckl'),
        os.path.join(pathload, 'result_km.pckl'), os.path.join(pathload, 'result_info.pckl'))
    with bz2.BZ2File(filename_para, 'wb') as f:
        pickle.dump(para, f)
    with bz2.BZ2File(filename_stoptime, 'wb') as f:
        pickle.dump(stoptime, f)
    with bz2.BZ2File(filename_dosage, 'wb') as f:
        pickle.dump(dosage, f)
    with bz2.BZ2File(filename_pop, 'wb') as f:
        pickle.dump(pop, f)
    with bz2.BZ2File(filename_popdict, 'wb') as f:
        pickle.dump(pop_dict, f)
    with bz2.BZ2File(filename_setind, 'wb') as f:
        pickle.dump(setind, f)
    with bz2.BZ2File(filename_km, 'wb') as f:
        pickle.dump(km, f)
    with bz2.BZ2File(filename_info, 'wb') as f:
        pickle.dump(info, f)
    return


if __name__ == "__main__":
    pass
