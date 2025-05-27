from DPM_lib import inspect, Counter, time, ThreadPool, re
from DPM_read_save import *
from DPM_strategy import *
from DPM_plot import DPM_plot_all, DPM_plot_1strategy, DPM_plot_allstrategy
from DPM_miscellaneous import DPM_miscellaneous_allequal, DPM_miscellaneous_slicedosage
from DPM_assign_check import *
from DPM_analysis import *


def DPM_run_par_csv_folder(*argv, **kwargs):
    print_quotation = True
    parname_list = PARNAME_LIST_RUN_PAR_CSV_FOLDER

    # If input non-keywrod argument, error.
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return
    kwargs = {k: v for k, v in kwargs.items() if v is not None}

    # If the inputted keywords not in the allowed arguments name list, error.
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    # Delete the keys in kwargs if its value is None.
    kwargs = {k.lower(): v for k, v in kwargs.items() if v is not None}

    par = DPM_check_inputpar(kwargs, parname_list)
    pathload = par['pathload']
    pathsave = par['pathsave']
    filename_pattern = par['filename_pattern']
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
    erase_preresult = par['erase_preresult']
    PAR_criterisa = par['PAR_criterisa']
    Strategy_name = par['Strategy_name']
    fullinput = par['fullinput']
    save_filename_param = par['save_filename_param']
    save_filename_stopt = par['save_filename_stopt']
    save_filename_dosage = par['save_filename_dosage']
    save_filename_pop = par['save_filename_pop']
    save_filename_eachtimepoint = par['save_filename_eachtimepoint']
    misspecification_sim_only = par['misspecification_sim_only']
    lookahead_step = par['lookahead_step']
    Maxnum_subseq = par['Maxnum_subseq']
    subtreedepth = par['subtreedepth']
    misspecification_ofdecision = par['misspecification_ofdecision']
    misspecification_atsim = par['misspecification_atsim']
    subclone_LOD = par['subclone_LOD']
    misspecification_LOD = par['misspecification_LOD']
    mutation_rate = par['mutation_rate']
    misspecification_fileload = par['misspecification_fileload']
    use_parallel = par['use_parallel']

    file_list = os.listdir(pathload)
    file_list = [filename for filename in file_list if re.search(filename_pattern, filename) is not None]
    if not file_list:
        print('Cannot find any default dataset in the current directory:')
        print(Colored.BOLD + Colored.PURPLE + pathload + Colored.END)
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

    file_list = [os.path.join(pathload, filename) for filename in file_list]

    file_list = file_list[:10]
    misfile_list_sort = misfile_list_sort[:10]
    if use_parallel and len(file_list) > 1:
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(DPM_run_par_csv)
                                   (filename_csv=i_filename,
                                    misspecification_filename_csv=misfile_list_sort[i],
                                    misspecification_ofdecision=misspecification_ofdecision,
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
                                    erase_preresult=erase_preresult,
                                    PAR_criterisa=PAR_criterisa,
                                    Strategy_name=Strategy_name,
                                    fullinput=fullinput,
                                    save_filename_param=save_filename_param,
                                    save_filename_stopt=save_filename_stopt,
                                    save_filename_pop=save_filename_pop,
                                    save_filename_dosage=save_filename_dosage,
                                    save_filename_eachtimepoint=save_filename_eachtimepoint,
                                    misspecification_sim_only=misspecification_sim_only,
                                    pathsave=pathsave,
                                    lookahead_step=lookahead_step,
                                    Maxnum_subseq=Maxnum_subseq,
                                    subtreedepth=subtreedepth,
                                    subclone_LOD=subclone_LOD,
                                    misspecification_LOD=misspecification_LOD,
                                    mutation_rate=mutation_rate)
                                   for i, i_filename in enumerate(file_list))
        print('Finished.')
    else:
        count_file = 1
        for i, i_filename in enumerate(file_list):
            print('Processing the ' + Colored.BOLD + Colored.PURPLE + str(count_file) + 'th' + Colored.END +
                  ' file of ' + Colored.BOLD + Colored.PURPLE + str(len(file_list)) + Colored.END + ' total files...')
            DPM_run_par_csv(filename_csv=i_filename,
                            misspecification_filename_csv=misfile_list_sort[i],
                            misspecification_ofdecision=misspecification_ofdecision,
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
                            erase_preresult=erase_preresult,
                            PAR_criterisa=PAR_criterisa,
                            Strategy_name=Strategy_name,
                            fullinput=fullinput,
                            save_filename_param=save_filename_param,
                            save_filename_stopt=save_filename_stopt,
                            save_filename_pop=save_filename_pop,
                            save_filename_dosage=save_filename_dosage,
                            save_filename_eachtimepoint=save_filename_eachtimepoint,
                            misspecification_sim_only=misspecification_sim_only,
                            pathsave=pathsave,
                            lookahead_step=lookahead_step,
                            Maxnum_subseq=Maxnum_subseq,
                            subtreedepth=subtreedepth,
                            subclone_LOD=subclone_LOD,
                            misspecification_LOD=misspecification_LOD,
                            mutation_rate=mutation_rate)
            count_file += 1
        print('Finished.')
    return


# Run simulation of parameters from csv file.
def DPM_run_par_csv(*argv, **kwargs):
    def DPM_run_par_csv_1(filename_):
        # Delete string '.csv' in filename.
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

    # If input non-keywrod argument, error.
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return
    kwargs = {k: v for k, v in kwargs.items() if v is not None}

    # If the inputted keywords not in the allowed arguments name list, error.
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    # Delete the keys in kwargs if its value is None.
    kwargs = {k.lower(): v for k, v in kwargs.items() if v is not None}

    par = DPM_check_inputpar(kwargs, parname_list)
    Num_drug = par['Num_drug']
    par_ind = par['par_ind']
    dose_combination = par['dose_combination']
    Simduration = par['Simduration']
    Limit_moleculardetection = par['Limit_moleculardetection']
    Limit_mortality = par['Limit_mortality']
    Stepsize = par['Stepsize']
    Limit_radiologicdetection = par['Limit_radiologicdetection']
    lookahead_step = par['lookahead_step']
    Maxnum_subseq = par['Maxnum_subseq']
    subtreedepth = par['subtreedepth']
    run_sim = par['run_sim']
    erase_preresult = par['erase_preresult']
    PAR_criterisa = par['PAR_criterisa']
    Strategy_name = par['Strategy_name']
    fullinput = par['fullinput']
    save_filename_param = par['save_filename_param']
    save_filename_stopt = par['save_filename_stopt']
    save_filename_pop = par['save_filename_pop']
    save_filename_dosage = par['save_filename_dosage']
    save_filename_eachtimepoint = par['save_filename_eachtimepoint']
    pathsave = par['pathsave']
    filename_csv = par['filename_csv']
    misspecification_filename_csv = par['misspecification_filename_csv']
    misspecification_ofdecision = par['misspecification_ofdecision']
    misspecification_atsim = par['misspecification_atsim']
    misspecification_sim_only = par['misspecification_sim_only']
    subclone_LOD = par['subclone_LOD']
    misspecification_LOD = par['misspecification_LOD']
    mutation_rate = par['mutation_rate']

    print(filename_csv)
    par_csv = DPM_read_par_csv(filename_csv, 'filename_csv')
    par_csv = DPM_assign_par_csvread(par_csv, Num_drug)
    if not par_csv:
        filename_csv_text = Colored.BOLD + Colored.PURPLE + 'filename_csv' + Colored.END
        print('No avaiable parameter from keyword input ' + filename_csv_text + '.')
        DPM_print_errorandclose()
        return

    # Check have mis_specification files.
    mis_par_csv = None
    if misspecification_filename_csv:
        mis_par_csv = DPM_read_par_csv(misspecification_filename_csv, 'misspecification_filename_csv')
        mis_par_csv = DPM_assign_par_csvread(mis_par_csv, Num_drug)
        if not mis_par_csv:
            misspecification_filename_csv_text = Colored.BOLD + Colored.PURPLE + 'misspecification_filename_csv' + Colored.END
            print('No avaiable parameter from keyword input ' + misspecification_filename_csv_text + '.')
            DPM_print_errorandclose()
            return

    # Verify if the number of misspecified parameters matches the number of true parameters.
    if 'misspecification_filename_csv' in kwargs.keys() and kwargs['misspecification_filename_csv'] != '':
        if len(par_csv) != len(mis_par_csv):
            filename_csv_text = Colored.BOLD + Colored.PURPLE + filename_csv + Colored.END
            misspecification_filename_csv_text = Colored.BOLD + Colored.PURPLE + misspecification_filename_csv + Colored.END
            print('The number of parameters ' + filename_csv_text +
                  ' is not the same as the number of parameters ' + misspecification_filename_csv_text + '.')
            DPM_print_errorandclose()
            return

    celltype = ['R1pop', 'R2pop', 'R12pop']
    result = dict(zip(celltype, [{'total': [], 'percent': [], 'est': []} for _ in range(len(MISSPECIFICATION_LOD_STR))]))
    keys = ['all']
    [keys.append(i_strategy) for i_strategy in Strategy_name]
    i_indsame, mis_specfiy_pop = None, dict(zip(keys, [deepcopy(result) for _ in range(len(keys))]))

    if subclone_LOD:
        celltype = ['R1pop', 'R2pop']
        if mis_par_csv is None:
            mis_par_csv, mis_specfiy_pop = \
                DPM_generate_misspecification_subcolone(deepcopy(par_csv), subclone_LOD, misspecification_LOD, mutation_rate, celltype,
                                                        mis_specfiy_pop, 'all')

            par_csv_df = pd.DataFrame.from_dict(par_csv)
            par_csv_df['X'] = par_csv_df.loc[:, ['Spop', 'R1pop', 'R2pop', 'R12pop']].sum(axis=1)
            par_csv_df['R1pop per'] = par_csv_df['R1pop']/par_csv_df['X']
            par_csv_df['R2pop per'] = par_csv_df['R2pop']/par_csv_df['X']
            par_csv_df['R12pop per'] = par_csv_df['R12pop']/par_csv_df['X']
            assert all(par_csv_df['R12pop per'] == 0)

            i_ind = (par_csv_df['R1pop per'] >= 2 * subclone_LOD) & (par_csv_df['R2pop per'] >= 2 * subclone_LOD)
            mis_par_csv_df = pd.DataFrame.from_dict(mis_par_csv)

            flag = pd.DataFrame(par_csv_df[i_ind.values][mis_par_csv_df.columns] == mis_par_csv_df[i_ind.values][mis_par_csv_df.columns]).\
                all('columns')
            assert flag.all('rows')
            i_indsame = i_ind[i_ind].index.values
        else:
            mis_par_csv = DPM_generate_misspecification_subcolone(deepcopy(mis_par_csv), subclone_LOD, misspecification_LOD, mutation_rate,
                                                                  celltype, mis_specfiy_pop)

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
        if mis_specification_csv:
            mis_filename_specfiy_pop = os.path.join(pathsave, filename + '_result_misspecfiy_pop_' + timenow + '.pckl')
        else:
            mis_filename_specfiy_pop = ''

        filename_param = os.path.join(pathsave, filename_param) if save_filename_param else ''
        filename_stopt = os.path.join(pathsave, filename_stopt) if save_filename_stopt else ''
        filename_dosage = os.path.join(pathsave, filename_dosage) if save_filename_dosage else ''
        filename_pop = os.path.join(pathsave, filename_pop) if save_filename_pop else ''
        filename_eachtimepoint = os.path.join(pathsave, filename_eachtimepoint) if save_filename_eachtimepoint else ''

        Heading_param_csv, Heading_stopt_csv, Heading_dosage_csv, Heading_pop_csv, Heading_eachtimepoint_csv = \
            DPM_generate_heading_csv(Num_drug, Simduration, Stepsize)

        filename_zip = zip((save_filename_param, save_filename_stopt, save_filename_dosage, save_filename_pop, save_filename_eachtimepoint),
                           (filename_param, filename_stopt, filename_dosage, filename_pop, filename_eachtimepoint),
                           (Heading_param_csv, Heading_stopt_csv, Heading_dosage_csv, Heading_pop_csv, Heading_eachtimepoint_csv))
        mis_filename_zip = zip((save_filename_param, save_filename_stopt, save_filename_dosage, save_filename_pop, save_filename_eachtimepoint),
                               (mis_filename_param, mis_filename_stopt, mis_filename_dosage, mis_filename_pop, mis_filename_eachtimepoint),
                               (Heading_param_csv, Heading_stopt_csv, Heading_dosage_csv, Heading_pop_csv, Heading_eachtimepoint_csv))
        # If not run the mis-specification file only and would like to erase the previous results, create the empty files.
        not mis_specification_csv and erase_preresult and DPM_read_erase(filename_zip)
        # If run mis-specification file and would like to erase the previous results, create the empty files.
        mis_specification_csv and erase_preresult and DPM_read_erase(mis_filename_zip)

        # par_csv = par_csv[180:]
        # mis_par_csv = mis_par_csv[180:]
        with tqdm(total=len(par_csv), ncols=150, desc='Runing simulation for the parameters from .csv file input') as pbar:
            for i in range(len(par_csv)):
                i_par = par_csv[i]
                i_paramID = int(i_par['paramID'])

                if mis_specification_csv:
                    i_mis_par = mis_par_csv[i]
                    i_mis_paramID = int(i_mis_par['paramID'])
                    assert i_paramID == i_mis_paramID

                    i_mis_strategy0, i_mis_strategy1, i_mis_strategy2_1, i_mis_strategy2_2, i_mis_strategy3, i_mis_strategy4, \
                        i_mis_strategy5, i_mis_strategy6, i_mis_strategy7, i_mis_strategy8, i_mis_strategy9, i_mis_strategyoptimal = \
                        DPM_run_simulation(i_par, Strategy_name, dose_combination, Simduration, Stepsize, Limit_mortality,
                                           Limit_radiologicdetection, lookahead_step, Maxnum_subseq, subtreedepth, LSsim,
                                           misspecification_ofdecision, misspecification_atsim, i_mis_par, subclone_LOD,
                                           misspecification_LOD, mutation_rate, mis_specfiy_pop)

                    i_mis_specify_pop_strategy0 = i_mis_strategy0[-1]
                    i_mis_specify_pop_strategy2_2 = i_mis_strategy2_2[-1]
                    if subclone_LOD != '':
                        assert mis_specfiy_pop['all'] == i_mis_specify_pop_strategy0['all'] == i_mis_specify_pop_strategy2_2['all']
                        mis_specfiy_pop['strategy0'] = i_mis_specify_pop_strategy0['strategy0']
                        mis_specfiy_pop['strategy2.2'] = i_mis_specify_pop_strategy2_2['strategy2.2']

                    i_mis_strategy0 = i_mis_strategy0[:-1]
                    i_mis_strategy2_2 = i_mis_strategy2_2[:-1]

                    i_mis_stopt = DPM_save_result_csv(i_mis_par, Stepsize, Simduration, i_mis_strategy0, i_mis_strategy1,
                                                      i_mis_strategy2_1, i_mis_strategy2_2, i_mis_strategy3, i_mis_strategy4,
                                                      i_mis_strategy5, i_mis_strategy6, i_mis_strategy7, i_mis_strategy8,
                                                      i_mis_strategy9, i_mis_paramID, mis_filename_param, mis_filename_stopt,
                                                      mis_filename_dosage, mis_filename_pop, mis_filename_eachtimepoint,
                                                      Limit_mortality)
                #
                i_strategy0, i_strategy1, i_strategy2_1, i_strategy2_2, i_strategy3, i_strategy4, i_strategy5, \
                    i_strategy6, i_strategy7, i_strategy8, i_strategy9, i_strategyoptimal = \
                    DPM_run_simulation(i_par, Strategy_name, dose_combination, Simduration, Stepsize, Limit_mortality,
                                       Limit_radiologicdetection, lookahead_step, Maxnum_subseq, subtreedepth, LSsim, False,
                                       False, None, '', '', '', mis_specfiy_pop)
                # Drug misspecification during simulation should not affect the initial drug decision.
                if mis_specification_csv and misspecification_atsim and subclone_LOD == '':
                    assert np.array_equal(i_strategy0[2][:, 0], i_mis_strategy0[2][:, 0])
                    assert np.array_equal(i_strategy2_2[2][:, 0], i_mis_strategy2_2[2][:, 0])

                # if not misspecification_sim_only:
                #     i_strategy0, i_strategy2_2 = i_strategy0[:-1], i_strategy2_2[:-1]
                #     i_stopt = DPM_save_result_csv(i_par, Stepsize, Simduration, i_strategy0, i_strategy1, i_strategy2_1, i_strategy2_2,
                #                                   i_strategy3, i_strategy4, i_strategy5, i_strategy6, i_strategy7, i_strategy8, i_strategy9,
                #                                   i_paramID, filename_param, filename_stopt, filename_dosage, filename_pop,
                #                                   filename_eachtimepoint, Limit_mortality, False)
                # if not mis_specification_csv:
                #     DPM_save_result_csv(i_par, Stepsize, Simduration, i_strategy0, i_strategy1, i_strategy2_1, i_strategy2_2,
                #                         i_strategy3, i_strategy4, i_strategy5, i_strategy6, i_strategy7, i_strategy8, i_strategy9,
                #                         i_paramID, filename_param, filename_stopt, filename_dosage, filename_pop, filename_eachtimepoint,
                #                         Limit_mortality)

                # if mis_specification_csv and i_indsame is not None and i in i_indsame:
                #     if 'strategy0' in Strategy_name:
                #         assert i_strategy0[1].shape[1] == i_mis_strategy0[1].shape[1]
                #         assert all(np.equal(i_strategy0[1][:, -1], i_mis_strategy0[1][:, -1]))
                #     if 'strategy2.2' in Strategy_name:
                #         assert i_strategy2_2[1].shape[1] == i_mis_strategy2_2[1].shape[1]
                #         assert all(np.equal(i_strategy2_2[1][:, -1], i_mis_strategy2_2[1][:, -1]))
                #     assert i_mis_stopt == i_stopt

                pbar.update(1)

        # if mis_specification_csv:
        #     with bz2.BZ2File(mis_filename_specfiy_pop, 'wb') as f:
        #         pickle.dump(par, f)
    return


# Run simulation for 1 parameter and plot the result.
def DPM_run_plot_1PAR(*argv, **kwargs):
    print_quotation = True
    LSsim = True
    parname_list = PARNAME_LIST_RUN_PLOT_1PAR
    # If input non-keyword argument, error.
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return

    # If the inputted keywords not in the allowed arguments name list, error.
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    kwargs = {k.lower(): v for k, v in kwargs.items()}

    # Caller of the DPM_check_inputpar function, this should be 'DPM_run_plot_1PAR', the current function name.
    callername = inspect.stack()[0][3]
    par = DPM_check_inputpar(kwargs, parname_list, callername)
    Num_drug, i_par, Limit_radiologicdetection, dose_combination, Simduration, Limit_mortality, Stepsize, PAR_criterisa, Strategy_name, \
        lookahead_step, Maxnum_subseq, subtreedepth, erase_preresult, save_filename_param, save_filename_stopt, save_filename_pop, \
        save_filename_dosage, save_filename_eachtimepoint, pathsave, plot, savename, Limit_moleculardetection \
        = par['Num_drug'], par['par'], par['Limit_radiologicdetection'], par['dose_combination'], par['Simduration'], par['Limit_mortality'], \
        par['Stepsize'], par['PAR_criterisa'], par['Strategy_name'], par['lookahead_step'], par['Maxnum_subseq'], par['subtreedepth'], \
        par['erase_preresult'], par['save_filename_param'], par['save_filename_stopt'], par['save_filename_pop'], par['save_filename_dosage'], \
        par['save_filename_eachtimepoint'], par['pathsave'], par['plot'], par['savename'], par['Limit_moleculardetection']

    criterisa = DPM_check_criterisa_2drug(PAR_criterisa, i_par, Simduration, Limit_mortality) if Num_drug == 2 else \
        DPM_check_criterisa_3drug(PAR_criterisa, i_par, Simduration, Limit_mortality) if Num_drug == 3 else False
    X0total = i_par['Spop'] + i_par['R1pop'] + i_par['R2pop'] + i_par['R12pop'] if Num_drug == 2 else \
        i_par['Spop'] + i_par['R1pop'] + i_par['R2pop'] + i_par['R3pop'] + i_par['R12pop'] + i_par['R13pop'] + i_par['R23pop'] + \
        i_par['R123pop'] if Num_drug == 3 else None
    par_check_index = DPM_check_par(X0total, Limit_mortality, Simduration, Stepsize)
    if not (all(criterisa) and not par_check_index):
        print('The input ' + Colored.BOLD + Colored.PURPLE + 'par' + Colored.END + " doesn't satisfy the required criterisa.")
        DPM_print_errorandclose()
        sys.exit()

    i_par['Num_cell_type'] = 4 if Num_drug == 2 else 8 if Num_drug == 3 else None
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

    Heading_param_csv, Heading_stopt_csv, Heading_dosage_csv, Heading_pop_csv, Heading_eachtimepoint_csv = \
        DPM_generate_heading_csv(Num_drug, Simduration, Stepsize)

    filename_zip = zip((save_filename_param, save_filename_stopt, save_filename_dosage, save_filename_pop, save_filename_eachtimepoint),
                       (filename_param, filename_stopt, filename_dosage, filename_pop, filename_eachtimepoint),
                       (Heading_param_csv, Heading_stopt_csv, Heading_dosage_csv, Heading_pop_csv, Heading_eachtimepoint_csv))

    # If not run the mis-specification file only and would like to erase the previous results, create the empty files.
    erase_preresult and DPM_read_erase(filename_zip)

    print(Colored.BOLD + Colored.PURPLE + 'Running...' + Colored.END)
    i_strategy0, i_strategy1, i_strategy2_1, i_strategy2_2, i_strategy3, i_strategy4, i_strategy5, i_strategy6, i_strategy7, i_strategy8, \
        i_strategy9, i_strategyoptimal = DPM_run_simulation(i_par, Strategy_name, dose_combination, Simduration, Stepsize, Limit_mortality,
                                                            Limit_radiologicdetection, lookahead_step, Maxnum_subseq, subtreedepth,
                                                            LSsim, False, None, None, None, None)

    DPM_save_result_csv(i_par, Stepsize, Simduration, i_strategy0, i_strategy1, i_strategy2_1, i_strategy2_2, i_strategy3, i_strategy4,
                        i_strategy5, i_strategy6, i_strategy7, i_strategy8, i_strategy9, 1, filename_param, filename_stopt, filename_dosage,
                        filename_pop, filename_eachtimepoint, Limit_mortality)
    strategy = {'strategy0': i_strategy0, 'strategy1': i_strategy1, 'strategy2.1': i_strategy2_1, 'strategy2.2': i_strategy2_2,
                'strategy3': i_strategy3, 'strategy4': i_strategy4, 'strategy5': i_strategy5, 'strategy6': i_strategy6, 'strategy7': i_strategy7,
                'strategy8': i_strategy8, 'strategy9': i_strategy9}
    DPM_plot_all(Num_drug, strategy, Simduration, Limit_moleculardetection, Limit_radiologicdetection, Limit_mortality, pathsave, savename)

    print(Colored.BOLD + Colored.PURPLE + 'Finished.' + Colored.END)
    return


# Run simulation for 1 parameter and plot the result.
def DPM_run_plot_1PAR_drugmis(*argv, **kwargs):
    print_quotation = True
    LSsim = True
    parname_list = PARNAME_LIST_RUN_PLOT_1PAR
    # If input non-keyword argument, error.
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return

    # If the inputted keywords not in the allowed arguments name list, error.
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    kwargs = {k.lower(): v for k, v in kwargs.items()}

    # Caller of the DPM_check_inputpar function, this should be 'DPM_run_plot_1PAR', the current function name.
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
    lookahead_step = par['lookahead_step']
    Maxnum_subseq = par['Maxnum_subseq']
    subtreedepth = par['subtreedepth']
    pathsave = par['pathsave']
    misspecification_ofdecision = par['misspecification_ofdecision']
    misspecification_atsim = par['misspecification_atsim']
    subclone_LOD = par['subclone_LOD']
    misspecification_LOD = par['misspecification_LOD']
    mutation_rate = par['mutation_rate']

    i_par['Num_cell_type'] = 4 if Num_drug == 2 else 8 if Num_drug == 3 else None
    i_mis_par['Num_cell_type'] = i_par['Num_cell_type']

    celltype = ['R1pop', 'R2pop', 'R12pop']
    result = dict(zip(celltype, [{'total': [], 'percent': [], 'est': []} for _ in range(len(MISSPECIFICATION_LOD_STR))]))
    keys = ['all']
    [keys.append(i_strategy) for i_strategy in Strategy_name]
    mis_specfiy_pop = dict(zip(keys, [deepcopy(result) for _ in range(len(keys))]))

    i_strategy0, i_strategy1, i_strategy2_1, i_strategy2_2, i_strategy3, i_strategy4, i_strategy5, \
        i_strategy6, i_strategy7, i_strategy8, i_strategy9, i_strategyoptimal = \
        DPM_run_simulation(i_par, Strategy_name, dose_combination, Simduration, Stepsize, Limit_mortality, Limit_radiologicdetection,
                           lookahead_step, Maxnum_subseq, subtreedepth, LSsim, False, False, None, '', '', '', mis_specfiy_pop)

    strategy = {'strategy0': i_strategy0[:3], 'strategy1': i_strategy1, 'strategy2.1': i_strategy2_1, 'strategy2.2': i_strategy2_2[:3],
                'strategy3': i_strategy3, 'strategy4': i_strategy4, 'strategy5': i_strategy5, 'strategy6': i_strategy6, 'strategy7': i_strategy7,
                'strategy8': i_strategy8, 'strategy9': i_strategy9}
    strategy = {key: strategy[key] for key in Strategy_name}

    i_mis_strategy0, i_mis_strategy1, i_mis_strategy2_1, i_mis_strategy2_2, i_mis_strategy3, i_mis_strategy4, i_mis_strategy5, \
        i_mis_strategy6, i_mis_strategy7, i_mis_strategy8, i_mis_strategy9, i_mis_strategyoptimal = \
        DPM_run_simulation(i_par, Strategy_name, dose_combination, Simduration, Stepsize, Limit_mortality,
                           Limit_radiologicdetection, lookahead_step, Maxnum_subseq, subtreedepth, LSsim,
                           misspecification_ofdecision, misspecification_atsim, i_mis_par, subclone_LOD,
                           misspecification_LOD, mutation_rate, mis_specfiy_pop)

    mis_strategy = {'strategy0': i_mis_strategy0[:3], 'strategy1': i_mis_strategy1, 'strategy2.1': i_mis_strategy2_1,
                    'strategy2.2': i_mis_strategy2_2[:3], 'strategy3': i_mis_strategy3, 'strategy4': i_mis_strategy4,
                    'strategy5': i_mis_strategy5, 'strategy6': i_mis_strategy6, 'strategy7': i_mis_strategy7,
                    'strategy8': i_mis_strategy8, 'strategy9': i_mis_strategy9}
    mis_strategy = {key: mis_strategy[key] for key in Strategy_name}

    index_mis_strategy_sim = [0, 3, 2]
    DPM_plot_1strategy(Num_drug, i_strategy0[:3], Simduration, Limit_moleculardetection, Limit_radiologicdetection,
                       Limit_mortality, 'strategy0', pathsave)
    plt.savefig(os.path.join(pathsave, 'strategy0.pdf'), dpi=FIG_DPI)
    plt.close()

    DPM_plot_1strategy(Num_drug, i_strategy2_2[:3], Simduration, Limit_moleculardetection, Limit_radiologicdetection,
                       Limit_mortality, 'strategy2', pathsave)
    plt.savefig(os.path.join(pathsave, 'strategy2.pdf'), dpi=FIG_DPI)
    plt.close()

    DPM_plot_1strategy(Num_drug, i_mis_strategy0[:3], Simduration, Limit_moleculardetection, Limit_radiologicdetection,
                       Limit_mortality, 'strategy0 mis', pathsave)
    plt.savefig(os.path.join(pathsave, 'strategy0 mis.pdf'), dpi=FIG_DPI)
    plt.close()

    DPM_plot_1strategy(Num_drug, i_mis_strategy2_2[:3], Simduration, Limit_moleculardetection, Limit_radiologicdetection,
                       Limit_mortality, 'strategy2 mis', pathsave)
    plt.savefig(os.path.join(pathsave, 'strategy2 mis.pdf'), dpi=FIG_DPI)
    plt.close()

    X_total = dict()
    X_total['strategy0'] = (i_strategy0[0], i_strategy0[1].sum(axis=0))
    X_total['strategy0 mis'] = (i_mis_strategy0[0], i_mis_strategy0[1].sum(axis=0))
    X_total['strategy2'] = (i_strategy2_2[0], i_strategy2_2[1].sum(axis=0))
    X_total['strategy2 mis'] = (i_mis_strategy2_2[0], i_mis_strategy2_2[1].sum(axis=0))
    DPM_plot_allstrategy(Num_drug, X_total, Simduration, Limit_radiologicdetection, Limit_mortality, pathsave)
    plt.savefig(os.path.join(pathsave, 'total.pdf'), dpi=FIG_DPI)
    plt.close()

    X_R12 = dict()
    X_R12['strategy0'] = (i_strategy0[0], i_strategy0[1][-1, :])
    X_R12['strategy0 ,is'] = (i_mis_strategy0[0], i_mis_strategy0[1][-1, :])
    X_R12['strategy2'] = (i_strategy2_2[0], i_strategy2_2[1][-1, :])
    X_R12['strategy2 mis'] = (i_mis_strategy2_2[0], i_mis_strategy2_2[1][-1, :])
    DPM_plot_allstrategy(Num_drug, X_R12, Simduration, Limit_radiologicdetection, Limit_mortality, pathsave)
    plt.savefig(os.path.join(pathsave, 'R12.pdf'), dpi=FIG_DPI)
    plt.close()
    return


# Run different strategy simulations.
def DPM_run_simulation(PAR, Strategy_name, dose_combination, Simduration, Stepsize, Limit_mortality, Limit_radiologicdetection, lookahead_step,
                       Maxnum_subseq, subtreedepth, LSsim, misspecification_ofdecision, misspecification_atsim, mis_par, subclone_LOD,
                       misspecification_LOD, mutation_rate, mis_specfiy_pop):
    Strategy0, Strategy1, Strategy2_1, Strategy2_2, Strategy3, Strategy4, Strategy5, Strategy6, Strategy7, Strategy8, Strategy9, \
        Strategy_optimal = None, None, None, None, None, None, None, None, None, None, None, None

    # Optimal control results from Dymos package
    if 'optimal'.lower() in Strategy_name:
        t_optimal, X_optimal, d_optimal = DPM_strategy_optimal(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                               Limit_mortality, LSsim, misspecification_ofdecision, mis_par)
        Strategy_optimal = [t_optimal, X_optimal, d_optimal]

    if 'strategy0'.lower() in Strategy_name:
        if not subclone_LOD:
            t_Strategy0, X_Strategy0, d_Strategy0 = DPM_strategy_0(PAR,
                                                                   Simduration,
                                                                   Stepsize,
                                                                   SIMTIMESTEP_DEFAULT_VAL,
                                                                   Limit_mortality,
                                                                   Limit_radiologicdetection,
                                                                   LSsim,
                                                                   misspecification_ofdecision,
                                                                   misspecification_atsim,
                                                                   mis_par)
            Strategy0 = [t_Strategy0, X_Strategy0, d_Strategy0, None]
        # Misspecification from the subclone LOD.
        else:
            t_Strategy0, X_Strategy0, d_Strategy0, mis_specfiy_pop = DPM_strategy_LOD_0(PAR,
                                                                                        Simduration,
                                                                                        Stepsize,
                                                                                        SIMTIMESTEP_DEFAULT_VAL,
                                                                                        Limit_mortality,
                                                                                        Limit_radiologicdetection,
                                                                                        mis_par,
                                                                                        subclone_LOD,
                                                                                        misspecification_LOD,
                                                                                        mutation_rate,
                                                                                        mis_specfiy_pop,
                                                                                        'strategy0')
            Strategy0 = [t_Strategy0, X_Strategy0, d_Strategy0, mis_specfiy_pop]

    if 'strategy1'.lower() in Strategy_name:
        t_Strategy1, X_Strategy1, d_Strategy1 = DPM_strategy_1(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                               Limit_mortality, LSsim, misspecification_ofdecision, mis_par)
        Strategy1 = [t_Strategy1, X_Strategy1, d_Strategy1]

    if 'strategy2.1'.lower() in Strategy_name:
        Strategy2threshold = 1e9
        t_Strategy2_1, X_Strategy2_1, d_Strategy2_1 = DPM_strategy_2(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                                     Limit_mortality, LSsim, Strategy2threshold, misspecification_ofdecision,
                                                                     mis_par)
        Strategy2_1 = [t_Strategy2_1, X_Strategy2_1, d_Strategy2_1]

    if 'strategy2.2'.lower() in Strategy_name:
        Strategy2threshold = 1e11
        if not subclone_LOD:
            t_Strategy2_2, X_Strategy2_2, d_Strategy2_2 = DPM_strategy_2(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                                         Limit_mortality, LSsim, Strategy2threshold, misspecification_ofdecision,
                                                                         misspecification_atsim, mis_par)
            Strategy2_2 = [t_Strategy2_2, X_Strategy2_2, d_Strategy2_2, None]
        # misspecification from subcolone LOD
        else:
            t_Strategy2_2, X_Strategy2_2, d_Strategy2_2, mis_specfiy_pop = \
                DPM_strategy_LOD_2(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL, Limit_mortality,
                                   Limit_radiologicdetection, Strategy2threshold, mis_par, subclone_LOD, misspecification_LOD,
                                   mutation_rate, mis_specfiy_pop, 'strategy2.2')
            Strategy2_2 = [t_Strategy2_2, X_Strategy2_2, d_Strategy2_2, mis_specfiy_pop]

    if 'strategy3'.lower() in Strategy_name:
        t_Strategy3, X_Strategy3, d_Strategy3 = DPM_strategy_3(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                               Limit_mortality, LSsim, misspecification_ofdecision, mis_par)
        Strategy3 = [t_Strategy3, X_Strategy3, d_Strategy3]

    if 'strategy4'.lower() in Strategy_name:
        t_Strategy4, X_Strategy4, d_Strategy4 = DPM_strategy_4(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                               Limit_mortality, LSsim, misspecification_ofdecision, mis_par)
        Strategy4 = [t_Strategy4, X_Strategy4, d_Strategy4]

    # Multistep extension of strategy 1.
    if 'strategy5'.lower() in Strategy_name:
        t_Strategy5, X_Strategy5, d_Strategy5 = DPM_strategy_5(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                               Limit_mortality, lookahead_step, LSsim, misspecification_ofdecision, mis_par)
        Strategy5 = [t_Strategy5, X_Strategy5, d_Strategy5]

    # Multistep extension of strategy 2.1.
    if 'strategy6'.lower() in Strategy_name:
        Strategy6threshold = 1e9
        t_Strategy6, X_Strategy6, d_Strategy6 = DPM_strategy_6(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                               Limit_mortality, lookahead_step, LSsim, Strategy6threshold,
                                                               misspecification_ofdecision, mis_par)
        Strategy6 = [t_Strategy6, X_Strategy6, d_Strategy6]

    # Multistep extension of strategy 2.2.
    if 'strategy7'.lower() in Strategy_name:
        Strategy7threshold = 1e11
        t_Strategy7, X_Strategy7, d_Strategy7 = DPM_strategy_6(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                               Limit_mortality, lookahead_step, LSsim, Strategy7threshold,
                                                               misspecification_ofdecision, mis_par)
        Strategy7 = [t_Strategy7, X_Strategy7, d_Strategy7]

    # Multistep extension of strategy 3.
    if 'strategy8'.lower() in Strategy_name:
        t_Strategy8, X_Strategy8, d_Strategy8 = DPM_strategy_8(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                               Limit_mortality, lookahead_step, LSsim, misspecification_ofdecision, mis_par)
        Strategy8 = [t_Strategy8, X_Strategy8, d_Strategy8]

    # ALTO (Adaptive Long-Term Optimization).
    if 'strategy9'.lower() in Strategy_name:
        t_Strategy9, X_Strategy9, d_Strategy9 = DPM_strategy_9(PAR, dose_combination, Simduration, Stepsize, SIMTIMESTEP_DEFAULT_VAL,
                                                               Limit_mortality, Maxnum_subseq, subtreedepth, LSsim, misspecification_ofdecision,
                                                               mis_par)
        Strategy9 = [t_Strategy9, X_Strategy9, d_Strategy9]

    return Strategy0, Strategy1, Strategy2_1, Strategy2_2, Strategy3, Strategy4, Strategy5, Strategy6, Strategy7, Strategy8, Strategy9, \
        Strategy_optimal


def DPM_run_processing(*argv, **kwargs):
    def DPM_run_processing_1(partial_filename):
        file_mis = dict(zip(MISSPECIFICATION_LOD_STR, [None] * len(MISSPECIFICATION_LOD_STR)))
        for i_key in MISSPECIFICATION_LOD_STR:
            file_mis[i_key] = [filename for filename in os.listdir(os.path.join(pathloadmis, i_key)) if partial_filename in filename]
        return file_mis

    def DPM_run_processing_2(i_stoptime_2):
        def DPM_run_processing_2_1(ind, pop_strategy_, strategy_):
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

        i_dosage = DPM_read_dosage_csv(i_file_dosage, Strategy_name, Num_stepdiff)
        i_pop = DPM_read_dosage_csv(i_file_pop, Strategy_name, Num_stepdiff)

        # num_cores = multiprocessing.cpu_count()
        for i_strategy in Strategy_name:
            i_pop_strategy = i_pop[i_strategy]
            for ii, i_pop_strategy_para in enumerate(i_pop_strategy):
                DPM_run_processing_2_1(ii, i_pop_strategy_para, i_strategy)
            # Parallel(n_jobs=num_cores)(delayed(DPM_run_processing_2_1)(ii, i_pop_strategy_para, i_strategy)
            #                            for ii, i_pop_strategy_para in enumerate(i_pop_strategy))

        for i_key in dosage.keys():
            dosage[i_key].extend(i_dosage[i_key])
            stoptime[i_key].extend(i_stoptime[i_key])
            pop[i_key].extend(i_pop[i_key])

        # load mis_specification data
        i_para_mis, i_stoptime_mis = None, None
        if pathloadmis:
            try:
                i_stoptime_mis, i_dosage_mis, i_para_mis = \
                    DPM_read_misspecify_LOD(pathloadmis, Strategy_name, Num_stepdiff, file_stoptime_mis, file_dosage_mis,
                                            file_para_mis, i_filename_pattern, para_key)

                for i_key in para_mis.keys():
                    for j_key in para.keys():
                        para_mis[i_key][j_key].extend(i_para_mis[i_key][j_key])
                for i_key in stoptime_mis.keys():
                    for j_key in stoptime.keys():
                        dosage_mis[i_key][j_key].extend(i_dosage_mis[i_key][j_key])
                        stoptime_mis[i_key][j_key].extend(i_stoptime_mis[i_key][j_key])

            except ValueError:
                print('\n', i_filename_pattern)

        i_para = pd.read_csv(os.path.join(pathload, i_file_para[0]), usecols=para_key)
        i_para['X'] = i_para.loc[:, para_key[1:5]].sum(axis=1)
        i_para['R1pop per'] = i_para['R1pop']/i_para['X']
        i_para['R2pop per'] = i_para['R2pop']/i_para['X']
        i_para['R12pop per'] = i_para['R12pop']/i_para['X']

        if 'pdf' in pathload:
            assert pd.DataFrame(i_para['R1pop per'] > MUTATION_RATE_DEFAULT_VAL).all('rows').values and \
                   pd.DataFrame(i_para['R2pop per'] > MUTATION_RATE_DEFAULT_VAL).all('rows').values and \
                   all(i_para['R12pop per'] == 0)

        if pathloadmis:
            i_LOD = float(pathloadmis.split('/')[-1])
            i_ind = (i_para['R1pop per'] >= 2*i_LOD) & (i_para['R2pop per'] >= 2*i_LOD)
            # i_ind_R1 = i_para['R1pop per'] < 2*i_LOD
            # i_ind_R2 = i_para['R2pop per'] < 2*i_LOD
            for key, values in i_para_mis.items():
                values = pd.DataFrame.from_dict(values)
                i_stoptime_mis_key = i_stoptime_mis[key]

                assert pd.DataFrame(i_para[i_ind.values]['Spop'] == values[i_ind.values]['Spop']).all('rows').values
                assert pd.DataFrame(i_para[i_ind.values]['R1pop'] == values[i_ind.values]['R1pop']).all('rows').values
                assert pd.DataFrame(i_para[i_ind.values]['R2pop'] == values[i_ind.values]['R2pop']).all('rows').values
                assert pd.DataFrame(i_para[i_ind.values]['R12pop'] == values[i_ind.values]['R12pop']).all('rows').values

                # idx_true, = np.where(i_ind.values)
                # f_true = i_para[i_ind.values]['Spop'] != values[i_ind.values]['Spop']
                # f_true[f_true].index

                # assert pd.DataFrame(i_para[i_ind_R1.values]['R1pop'] != values[i_ind_R1.values]['R1pop']).all('rows').values
                # assert pd.DataFrame(i_para[i_ind_R2.values]['R2pop'] != values[i_ind_R2.values]['R2pop']).all('rows').values

                # ii = i_para[i_ind_R1.values]['R1pop'] == values[i_ind_R1.values]['R1pop']
                # i_para[i_ind_R1.values][ii.values]['R1pop']
                # values[i_ind_R1.values][ii.values]['R1pop']

                for i_strategy in Strategy_name:
                    idx_true, = np.where(i_ind.values)
                    flag = [i_stoptime[i_strategy][ii] == i_stoptime_mis_key[i_strategy][ii] for ii in idx_true]

                    # j_stoptime = [i_stoptime[i_strategy][ii] for ii in idx_true]
                    # j_mis_stoptime = [i_stoptime_mis_key[i_strategy][ii] for ii in idx_true]
                    # flag = [i_stoptime[i_strategy][ii] != i_stoptime_mis_key[i_strategy][ii] for ii in idx_true]
                    # ind_flag = [ii for ii, x in enumerate(flag) if x]
                    # if ind_flag:
                    #     idx_true[ind_flag[0]]
                    #     j_stoptime[ind_flag[0]]
                    #     j_mis_stoptime[ind_flag[0]]

                    assert all(flag)
        pbar.update(1)
        return

    def DPM_run_processing_3(para_i, dosage_i, stoptime_i, num_dosage):
        para_i_df = pd.DataFrame(para_i)
        para_i_df['Xtotal'] = para_i_df['Spop'] + para_i_df['R1pop'] + para_i_df['R2pop'] + para_i_df['R12pop']
        # Select subdataset for mis_specification under LOD
        # (1) Subset of patients with R1 >= 7.1e-7 and R2 >= 7.1e-7
        bool_R1_R2_e_7 = (para_i_df['R1pop'] / para_i_df['Xtotal'] >= 7.1e-7) & (para_i_df['R2pop'] / para_i_df['Xtotal'] >= 7.1e-7)
        bool_R1_R2_e_7 = bool_R1_R2_e_7.values.tolist()
        # (2) Subset of patients with R1 >= 7.1e-5 and R2 >= 7.1e-5
        bool_R1_R2_e_5 = (para_i_df['R1pop'] / para_i_df['Xtotal'] >= 7.1e-5) & (para_i_df['R2pop'] / para_i_df['Xtotal'] >= 7.1e-5)
        bool_R1_R2_e_5 = bool_R1_R2_e_5.values.tolist()

        total = [True for _ in range(len(para_i_df['Xtotal']))]
        del para_i, para_i_df

        dosage_i_sliced = DPM_miscellaneous_slicedosage(dosage_i, num_dosage, Strategy_name)
        bool_diff_dosage = [x != y for x, y in zip(dosage_i_sliced['strategy0'], dosage_i_sliced['strategy2.2'])]
        bool_same_dosage = [x == y for x, y in zip(dosage_i_sliced['strategy0'], dosage_i_sliced['strategy2.2'])]

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
                print(f"strategy2 sigbetter:{len(info_i[i_setname]['paramID sigbetter'])}")
            for i_strategy in Strategy_name:
                i_stoptime_strategy = i_stoptime_set[i_strategy]
                i_dosage_set_strategy = i_dosage_set[i_strategy]
                km_i[i_setname][i_strategy] = DPM_analysis_KM(i_stoptime_strategy, Simduration)
                info_i[i_setname][i_strategy] = DPM_analysis_dose(i_dosage_set_strategy, i_strategy)

        return setind_i, km_i, info_i

    print_quotation = True
    parname_list = PARNAME_LIST_RUN_PROCESSING
    # If input non-keyword argument, error.
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return

    # If the inputted keywords not in the allowed arguments name list, error.
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    # Delete the keys in kwargs if its value is None.
    kwargs = {k.lower(): v for k, v in kwargs.items() if v is not None}

    par = DPM_check_inputpar(kwargs, parname_list)
    Num_drug, Strategy_name, pathsave, pathload, pathloadmis, Num_stepdiff, Simduration, Stepsize, use_parallel \
        = par['Num_drug'], par['Strategy_name'], par['pathsave'], par['pathload'], par['pathloadmis'], par['Num_stepdiff'], par['Simduration'], \
        par['Stepsize'], par['use_parallel']

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

    file_stoptime_mis, file_dosage_mis, file_para_mis, LOD = dict(), dict(), dict(), None
    if pathloadmis:
        LOD = pathloadmis[pathloadmis.find('e-0')-1:]
        LOD = LOD if LOD[-1] != '/' else LOD[:-1]
        file_stoptime_mis = DPM_run_processing_1('result_stopt')
        file_dosage_mis = DPM_run_processing_1('result_dosage')
        file_para_mis = DPM_run_processing_1('result_para')

    # If there are no files found in the current directory, exit.
    if not file_stoptime or not file_dosage or not file_pop:
        print('Cannot find any stopping time file in the current directory:')
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

    # Read csv stopping time file
    result = {**{'paramID': []}, **{i_strategy: [] for i_strategy in Strategy_name}}
    stoptime, dosage, pop = deepcopy(result), deepcopy(result), deepcopy(result)

    stoptime_mis, dosage_mis, para_mis, para_mis_0, para_mis_pdf = dict(), dict(), dict(), pd.DataFrame(), pd.DataFrame()
    para_key = ['paramID', 'Spop', 'R1pop', 'R2pop', 'R12pop', 'g0_S', 'Sa.S.D1.', 'Sa.S.D2.', 'Sa.R1.D1.', 'Sa.R2.D2.', 'T.R1..S.', 'T.R2..S.']
    para = {i_key: [] for i_key in para_key}
    if pathloadmis:
        stoptime_mis = dict(zip(MISSPECIFICATION_LOD_STR, [deepcopy(result) for _ in range(len(MISSPECIFICATION_LOD_STR))]))
        dosage_mis = dict(zip(MISSPECIFICATION_LOD_STR, [deepcopy(result) for _ in range(len(MISSPECIFICATION_LOD_STR))]))
        para_mis = dict(zip(MISSPECIFICATION_LOD_STR, [deepcopy(para) for _ in range(len(MISSPECIFICATION_LOD_STR))]))

    # result = pd.DataFrame(columns=list(para.keys())[1:5])
    pop_dict = dict(zip(Strategy_name, [None] * len(Strategy_name)))
    for i_strategy in Strategy_name:
        pop_dict[i_strategy] = [[] for _ in range((int(Simduration/Stepsize) + 1))]

    with tqdm(total=len(file_stoptime), ncols=150, desc='Runing simulation for the parameters from .csv file input') as pbar:
        st = time()
        if use_parallel:
            with ThreadPool() as pool:
                pool.map(DPM_run_processing_2, file_stoptime)
        else:
            for i, i_file_stoptime in enumerate(file_stoptime):
                try:
                    DPM_run_processing_2(i_file_stoptime)
                except (AssertionError,  ValueError):
                    print()
                    i_file_text = Colored.BOLD + Colored.PURPLE + i_file_stoptime + Colored.END
                    print(i_file_text)
                    pbar.update(1)
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

    if pathloadmis:
        check_list = []
        for i_key in MISSPECIFICATION_LOD_STR:
            check_list.append(Counter(stoptime_mis[i_key]['paramID']))
            check_list.append(Counter(dosage_mis[i_key]['paramID']))
            check_list.append(Counter(para_mis[i_key]['paramID']))
        check_list.append(Counter(para['paramID']))
        if not DPM_miscellaneous_allequal(check_list):
            raise ValueError('paramID mismatch.')

    setname = ['total', 'first 2 diff', 'first 2 same', 'R1,R2>=1e-7', 'R1,R2>=1e-5']

    if not pathloadmis:
        setind, km, info = DPM_run_processing_3(para, dosage, stoptime, Num_stepdiff)
        filename_para, filename_stoptime, filename_dosage, filename_pop, filename_popdict, filename_setind, filename_km, filename_info = \
            os.path.join(pathload, 'result_para.pckl'), \
            os.path.join(pathload, 'result_stoptime.pckl'), \
            os.path.join(pathload, 'result_dosage.pckl'), \
            os.path.join(pathload, 'result_pop.pckl'), \
            os.path.join(pathload, 'result_popdict.pckl'), \
            os.path.join(pathload, 'result_setind.pckl'), \
            os.path.join(pathload, 'result_km.pckl'), \
            os.path.join(pathload, 'result_info.pckl')
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

    if pathloadmis:
        setind_mis, km_mis, info_mis = [], [], []
        for i_key in MISSPECIFICATION_LOD_STR:
            setind_mis_i, km_mis_i, info_mis_i = DPM_run_processing_3(para, dosage_mis[i_key], stoptime_mis[i_key], Num_stepdiff)
            setind_mis.append(setind_mis_i)
            km_mis.append(km_mis_i)
            info_mis.append(info_mis_i)
        setind_mis = dict(zip(MISSPECIFICATION_LOD_STR, setind_mis))
        km_mis = dict(zip(MISSPECIFICATION_LOD_STR, km_mis))
        info_mis = dict(zip(MISSPECIFICATION_LOD_STR, info_mis))
        filename_para, filename_stoptime, filename_dosage, filename_setind, filename_km, filename_info = \
            os.path.join(pathloadmis, 'result_' + LOD + '_para.pckl'), \
            os.path.join(pathloadmis, 'result_' + LOD + '_stoptime.pckl'), \
            os.path.join(pathloadmis, 'result_' + LOD + '_dosage.pckl'), \
            os.path.join(pathloadmis, 'result_' + LOD + '_setind.pckl'), \
            os.path.join(pathloadmis, 'result_' + LOD + '_km.pckl'), \
            os.path.join(pathloadmis, 'result_' + LOD + '_info.pckl')
        with bz2.BZ2File(filename_para, 'wb') as f:
            pickle.dump(para_mis, f)
        with bz2.BZ2File(filename_stoptime, 'wb') as f:
            pickle.dump(stoptime_mis, f)
        with bz2.BZ2File(filename_dosage, 'wb') as f:
            pickle.dump(dosage_mis, f)
        with bz2.BZ2File(filename_setind, 'wb') as f:
            pickle.dump(setind_mis, f)
        with bz2.BZ2File(filename_km, 'wb') as f:
            pickle.dump(km_mis, f)
        with bz2.BZ2File(filename_info, 'wb') as f:
            pickle.dump(info_mis, f)

    if not pathloadmis:
        pass

        # DPM_analysis_pop(pop, Strategy_name)

        # bool_same_dosage = [x == y for x, y in zip(dosage['strategy0'], dosage['strategy2.2'])]

        # bool_same_dosage, bool_diff_dosage = pd.Series(bool_same_dosage), pd.Series(setind['first 2 diff'])
        # df_dosage = pd.DataFrame.from_dict(dosage)
        # df_dosage[bool_same_dosage]['strategy0'] == df_dosage[bool_same_dosage]['strategy2.2']
        # df_dosage[setind['first 2 diff']]['strategy0'] != df_dosage[setind['first 2 diff']]['strategy2.2']
        #
        # DPM_analysis_stat(stoptime, Strategy_name, Simduration, None, titlestr='Total')
        # DPM_analysis_stat(stoptime, Strategy_name, Simduration, bool_diff_dosage, titlestr='First two move diff')
        # DPM_analysis_stat(stoptime, Strategy_name, Simduration, bool_same_dosage, titlestr='First two move same')
    return


def DPM_run_analysis(pathload=None, pathsave=None, filename_pattern=None, Strategy_name=None,
                     Stepsize=STEPSIZE_DEFAULT_VAL, Simduration=SIMDURATION_DEFAULT_VAL):
    filename_para = os.path.join(pathload, filename_pattern) + '_para.pckl'
    filename_stopt = os.path.join(pathload, filename_pattern) + '_stopt.pckl'
    filename_pop = os.path.join(pathload, filename_pattern) + '_pop.pckl'
    with bz2.BZ2File(filename_para, 'rb') as f:
        para = pickle.load(f)
    with bz2.BZ2File(filename_pop, 'rb') as f:
        pop = pickle.load(f)

    # pop_dict = dict(zip(Strategy_name, [None]*len(Strategy_name)))
    # result = pd.DataFrame(columns=list(para.keys())[1:5])
    # for i_strategy in Strategy_name:
    #     pop_total = [result] * (int(Simduration/Stepsize)+1)
    #     i_pop = pop[i_strategy]
    #     for i, i_para in enumerate(i_pop):
    #         print(i)
    #         para_t0 = [para[key][i] for key in list(para.keys())[1:5]]
    #         para_t0_df = pd.DataFrame([para_t0], columns=list(para.keys())[1:5])
    #         pop_total[0] = pd.concat([pop_total[0], para_t0_df])
    #         i_para = i_para.split(';')
    #         for j, j_step in enumerate(i_para):
    #             if j_step != '-1':
    #                 j_step = [float(val) for val in re.sub('[()]', '', j_step).split(',')]
    #                 j_step_df = pd.DataFrame([j_step], columns=list(para.keys())[1:5])
    #                 pop_total[j+1] = pd.concat([pop_total[j+1], j_step_df])
    #     pop_dict[i_strategy] = pop_total


def DPM_run_analysis_LOD(*argv, **kwargs):
    print_quotation = True
    parname_list = PARNAME_LIST_RUN_ANALYSIS_LOD
    # If input non-keyword argument, error.
    if argv:
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return

    # If the inputted keywords not in the allowed arguments name list, error.
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    # Delete the keys in kwargs if its value is None.
    kwargs = {k.lower(): v for k, v in kwargs.items() if v is not None}

    par = DPM_check_inputpar(kwargs, parname_list)
    Num_drug, LOD, Simduration, Strategy_name, filename_pattern, pathload, pathsave, pathloadmis, plot \
        = par['Num_drug'], par['LOD'], par['Simduration'], par['Strategy_name'], par['filename_pattern'], par['pathload'], par['pathsave'], \
        par['pathloadmis'], par['plot']

    def DPM_run_analysis_LOD_1(pathloadin, misspecification_LOD_i=None):
        filename_para = pathloadin + '_para.pckl'
        # filename_stoptime = pathloadin + '_stoptime.pckl'
        filename_setind = pathloadin + '_setind.pckl'
        filename_km = pathloadin + '_km.pckl'
        filename_info = pathloadin + '_info.pckl'
        filename_popdict = pathloadin + '_popdict.pckl'

        if misspecification_LOD_i is None:
            with bz2.BZ2File(filename_para, 'rb') as f:
                para_i = pickle.load(f)
                para_i = para_i  # if misspecification_LOD_i is None else para_i[misspecification_LOD_i]]
        else:
            para_i = None

        # with bz2.BZ2File(filename_stoptime, 'rb') as f:
        #     stoptime_1 = pickle.load(f)
        #     stoptime_1 = stoptime_1 if misspecification_LOD_1 is None else stoptime_1[misspecification_LOD_1]

        if misspecification_LOD_i is None:
            with bz2.BZ2File(filename_setind, 'rb') as f:
                setind_i = pickle.load(f)
        else:
            setind_i = None

        with bz2.BZ2File(filename_km, 'rb') as f:
            km_i = pickle.load(f)
            km_i = km_i if misspecification_LOD_i is None else km_i[misspecification_LOD_i]

        with bz2.BZ2File(filename_info, 'rb') as f:
            info_i = pickle.load(f)
            info_i = info_i if misspecification_LOD_i is None else info_i[misspecification_LOD_i]

        if misspecification_LOD_i is None:
            with bz2.BZ2File(filename_popdict, 'rb') as f:
                pop_i = pickle.load(f)
        else:
            pop_i = None

        return km_i, info_i, para_i, setind_i, pop_i  # stoptime_1

    # load true simulation reuslts
    km, info, para, setind, pop = DPM_run_analysis_LOD_1(os.path.join(pathload, filename_pattern))
    setname = list(km.keys())
    km_mis, info_mis, para_mis = [[] for _ in range(len(MISSPECIFICATION_LOD_STR))], \
                                 [[] for _ in range(len(MISSPECIFICATION_LOD_STR))], \
                                 [[] for _ in range(len(MISSPECIFICATION_LOD_STR))]

    misspecification_LOD = MISSPECIFICATION_LOD_STR
    for i_LOD in LOD:
        for i, i_misspecification_LOD in enumerate(misspecification_LOD):
            # print(i_LOD, i_misspecification_LOD)
            i_km_LOD, i_info_LOD, i_para_LOD, *_ = \
                DPM_run_analysis_LOD_1(os.path.join(pathloadmis, i_LOD, filename_pattern + '_' + i_LOD), i_misspecification_LOD)
            km_mis[i].append(i_km_LOD)
            info_mis[i].append(i_info_LOD)
            para_mis[i].append(i_para_LOD)

    DPM_analysis_misspec_pop(para, pop, info, LOD, pathsave, pathload)
    DPM_analysis_misspec_par(para, info, info_mis, LOD, pathsave, pathload)
    DPM_analysis_misspec(km, info, km_mis, info_mis, LOD, setname, Strategy_name, Simduration, pathsave, plot)
    return


def DPM_run_generate_R1R2_pdftrue(pathload='./csv_data/2drug3million/', pathsave='./csv_data/2drug3million/pdftrue/',
                                  filename_pattern='_para', Num_drug=NUM_DRUG_DEFAULT_VAL, Simduration=SIMDURATION_DEFAULT_VAL,
                                  Limit_mortality=LIMIT_MORTALITY_DEFAULT_VAL, Stepsize=STEPSIZE_DEFAULT_VAL):
    PAR_criterisa = []
    celltype = ['R1pop', 'R2pop']
    mutation_rate = MUTATION_RATE_DEFAULT_VAL
    subclone_LOD = 0.5
    misspecification_LOD = 'pdf'
    file_list = os.listdir(pathload)
    file_list = [filename for filename in file_list if re.search(filename_pattern, filename) is not None]
    try:
        file_list.sort(key=lambda i: int(i.split('_')[-3]))
    except (ValueError, IndexError):
        pass
    file_list = [os.path.join(pathload, filename) for filename in file_list]

    x = np.linspace(mutation_rate, subclone_LOD, int(1e8))
    norm_pdf = simps(DPM_generate_pdf(x, mutation_rate, subclone_LOD), x)
    with tqdm(total=len(file_list), ncols=150) as pbar:
        for i_filename in file_list:
            par_csv = DPM_read_par_csv(i_filename, 'filename_csv')
            par_csv = DPM_assign_par_csvread(par_csv, Num_drug)
            par_csv = DPM_generate_misspecification_subcolone(par_csv, subclone_LOD, misspecification_LOD, mutation_rate, celltype,
                                                              norm_pdf=norm_pdf, Xtotal=5e9)
            par_csv = pd.DataFrame.from_records(par_csv)
            par_csv.drop(columns=['Num_drug', 'Num_cell_type'])
            par_csv['paramID'] = par_csv['paramID'].astype('int')

            i_filename_save = i_filename[i_filename.rfind('/') + 1:]
            pathsave = os.path.abspath(pathsave)
            if not os.path.exists(pathsave):
                os.makedirs(pathsave)
            i_filename_save = os.path.join(pathsave, i_filename_save)
            par_csv.to_csv(i_filename_save, index=False)
            pbar.update(1)
    return


def DPM_run_generate_misDrug2(pathload='./pnas/', pathsave='./pnas_div5/', ratio=1, filename_pattern='_para'):

    file_list = os.listdir(pathload)
    file_list = [filename for filename in file_list if re.search(filename_pattern, filename) is not None]
    try:
        file_list.sort(key=lambda i: int(i.split('_')[-3]))
    except (ValueError, IndexError):
        pass
    file_list = [os.path.join(pathload, filename) for filename in file_list]

    with tqdm(total=len(file_list), ncols=150) as pbar:
        for i_filename in file_list:
            par_csv = pd.read_csv(i_filename)
            par_csv['Sa.S.D2.'] = par_csv['Sa.S.D2.'] * ratio
            par_csv['Sa.R1.D2.'] = par_csv['Sa.R1.D2.'] * ratio
            par_csv['Sa.R2.D2.'] = par_csv['Sa.R2.D2.'] * ratio
            par_csv['Sa.R12.D2.'] = par_csv['Sa.R12.D2.'] * ratio

            par_csv = pd.DataFrame.from_records(par_csv)
            par_csv['paramID'] = par_csv['paramID'].astype('int')

            i_filename_save = i_filename[i_filename.rfind('/') + 1:]
            if ratio >= 1:
                i_filename_save = f'mis_{int(ratio)}x_' + i_filename_save
            else:
                i_filename_save = f'mis_div{int(1/ratio)}_' + i_filename_save
            pathsave = os.path.abspath(pathsave)
            if not os.path.exists(pathsave):
                os.makedirs(pathsave)
            i_filename_save = os.path.join(pathsave, i_filename_save)
            par_csv.to_csv(i_filename_save, index=False)
            pbar.update(1)
    return


if __name__ == "__main__":
    filename_pattern = '_para'
    subclone_LOD = [1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01]
    misspecification_LOD = ['pdf', 'loguni']
    pathload = './pdftrue'
    for i_subclone_LOD in subclone_LOD:
        for i_misspecification_LOD in misspecification_LOD:
            pathsave = f'./misspecification/pdftrue_mis_LOD/{i_subclone_LOD:.0e}/{i_misspecification_LOD}/'
            DPM_run_par_csv_folder(pathload=pathload,
                                   pathsave=pathsave,
                                   filename_pattern=filename_pattern,
                                   run_sim=True,
                                   subclone_LOD=i_subclone_LOD,
                                   misspecification_LOD=i_misspecification_LOD,
                                   Strategy_name=['strategy0', 'strategy2.2'],
                                   fullinput=True,
                                   save_filename_param=True,
                                   save_filename_stopt=True,
                                   save_filename_dosage=True,
                                   save_filename_pop=False,
                                   save_filename_eachtimepoint=False,
                                   misspecification_sim_only=True,
                                   use_parallel=False)

    # pathload = './misspecification/PNAS/'
    # pathsave = './misspecification/PNAS/figure'
    # pathload_misLOD = './misspecification/PNAS_mis_LOD/'
    #
    # LOD = ['1e-06', '1e-05', '1e-04', '1e-03', '1e-02', '1e-01']
    # filename_pattern = 'result'
    # DPM_run_analysis_LOD(Num_drug=2,
    #                      LOD=LOD,
    #                      filename_pattern=filename_pattern,
    #                      pathload=pathload,
    #                      pathsave=pathsave,
    #                      pathloadmis=pathload_misLOD,
    #                      Strategy_name=['strategy0', 'strategy2.2'],
    #                      plot=True)

    # pathload_ = './csv_result/2drug3million/pdftrue/'
    # pathsave_ = './csv_result/2drug3million/pdftrue/figure'
    # DPM_run_processing(Num_drug=2,
    #                    pathload=pathload_,
    #                    pathsave=pathsave_,
    #                    Num_stepdiff=2,
    #                    Strategy_name=['strategy0', 'strategy2.2'],
    #                    pathloadmis='',
    #                    use_parallel=False)

    # pathload_ = './csv_data/2drug3million/pnas/'
    # misspecification_fileload_ = './csv_data/2drug3million/pnas_div30/'
    # pathsave_ = './Peter Mon/pnas_div30/'
    # filename_pattern_ = '_para'
    # DPM_run_par_csv_folder(pathload=pathload_,
    #                        pathsave=pathsave_,
    #                        misspecification_fileload=misspecification_fileload_,
    #                        filename_pattern=filename_pattern_,
    #                        run_sim=True,
    #                        Strategy_name=['strategy0', 'strategy2.2'],
    #                        fullinput=True,
    #                        save_filename_param=True,
    #                        save_filename_stopt=True,
    #                        save_filename_dosage=True,
    #                        save_filename_pop=False,
    #                        save_filename_eachtimepoint=False,
    #                        misspecification_sim_only=True,
    #                        misspecification=False,
    #                        use_parallel=True)
