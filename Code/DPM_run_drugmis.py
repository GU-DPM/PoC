from DPM_lib import ThreadPool, resample
from DPM_read_save import *
from DPM_plot import DPM_plot_df2pdf, DPM_plot_multi, DPM_plot_KM, DPM_plot_KM_multi
from DPM_assign_check import *
from DPM_analysis import *
from DPM_miscellaneous import DPM_miscellaneous_fillful
from DPM_run import DPM_run_plot_1PAR_drugmis
from DPM_lib import deepcopy, plt, bz2, tabulate, pickle, tqdm, random
'''
CHECKED
'''


def DPM_run_drugmis_processing(pathload='', Strategy_name=None, Simduration=SIMDURATION_DEFAULT_VAL, use_parallel=False):
    def DPM_run_drugmis_processing_1(i_stoptime_):
        i_filename_pattern = i_stoptime_.split('_stopt', 1)[0]
        i_file_stoptime_ = os.path.join(pathload, i_stoptime_)
        i_stoptime = DPM_read_stoptime_csv(i_file_stoptime_, Strategy_name)
        i_file_para = i_filename_pattern + '_para'
        i_file_para = [i_file for i_file in file_para if i_file_para in i_file]

        i_para = pd.read_csv(os.path.join(pathload, i_file_para[0]), usecols=para_key).to_dict('list')
        for i_key in para.keys():
            para[i_key].extend(i_para[i_key])

        i_file_dosage = i_filename_pattern + '_dosage'
        i_file_dosage = [i_file for i_file in file_dosage if i_file_dosage in i_file]
        i_file_dosage = os.path.join(pathload, i_file_dosage[0])

        i_dosage = DPM_read_dosage_pop_csv(i_file_dosage, Strategy_name)

        for i_key in dosage.keys():
            dosage[i_key].extend(i_dosage[i_key])
            stoptime[i_key].extend(i_stoptime[i_key])

        i_para = pd.read_csv(os.path.join(pathload, i_file_para[0]), usecols=para_key)
        i_para['X'] = i_para.loc[:, para_key[1:5]].sum(axis=1)
        i_para['R1pop per'] = i_para['R1pop']/i_para['X']
        i_para['R2pop per'] = i_para['R2pop']/i_para['X']
        i_para['R12pop per'] = i_para['R12pop']/i_para['X']

        pbar.update(1)
        return

    def DPM_run_drugmis_processing_2(dosage_, stoptime_):
        km_ = {**{i_strategy: None for i_strategy in Strategy_name}, **{'hz_ratio': [], 'p': [], 'num': 0}}
        info_ = {**{i_strategy: None for i_strategy in Strategy_name}, 'paramID sigbetter': None}
        hz, p = DPM_analysis_hazard_ratio(stoptime_, Strategy_name, Simduration)
        km_['hz_ratio'], km_['p'], km_['num'] = hz, p, len(stoptime_['paramID'])
        info_['paramID sigbetter'] = DPM_analysis_sigbetter(stoptime_, Strategy_name)

        for i_strategy in Strategy_name:
            i_stoptime_strategy = stoptime_[i_strategy]
            i_dosage_set_strategy = dosage_[i_strategy]
            km_[i_strategy] = DPM_analysis_KM(i_stoptime_strategy, Simduration)
            info_[i_strategy] = DPM_analysis_dose(i_dosage_set_strategy, i_strategy)

        return km_, info_

    file_format = '.csv'
    file_list = os.listdir(pathload)
    partial_filename_stoptime = ['result_stopt', file_format]
    file_stoptime = [filename for filename in file_list if all([x in filename for x in partial_filename_stoptime])]
    partial_filename_dosage = ['result_dosage', file_format]
    file_dosage = [filename for filename in file_list if all([x in filename for x in partial_filename_dosage])]
    partial_filename_para = ['result_para', file_format]
    file_para = [filename for filename in file_list if all([x in filename for x in partial_filename_para])]

    assert len(file_stoptime) == len(file_dosage) and len(file_stoptime) > 0

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
    except ValueError:
        pass

    # Load the stopping-time data from the CSV file.#
    result = {**{'paramID': []}, **{i_strategy: [] for i_strategy in Strategy_name}}
    stoptime, dosage = deepcopy(result), deepcopy(result)

    para_key = ['paramID', 'Spop', 'R1pop', 'R2pop', 'R12pop', 'g0_S', 'Sa.S.D1.', 'Sa.S.D2.', 'Sa.R1.D1.', 'Sa.R2.D2.', 'T.R1..S.', 'T.R2..S.']
    para = {i_key: [] for i_key in para_key}

    with tqdm(total=len(file_stoptime), ncols=150, desc='processing...') as pbar:
        if use_parallel:
            with ThreadPool() as pool:
                pool.map(DPM_run_drugmis_processing_1, file_stoptime)
        else:
            for i, i_file_stoptime in enumerate(file_stoptime):
                try:
                    DPM_run_drugmis_processing_1(i_file_stoptime)
                except (AssertionError,  ValueError):
                    print()
                    i_file_text = Colored.BOLD + Colored.PURPLE + i_file_stoptime + Colored.END
                    print(i_file_text)
                    pbar.update(1)

    km, info = DPM_run_drugmis_processing_2(dosage, stoptime)
    filename_para, filename_stoptime, filename_dosage, filename_km, filename_info = \
        os.path.join(pathload, 'result_para.pckl'), \
        os.path.join(pathload, 'result_stoptime.pckl'), \
        os.path.join(pathload, 'result_dosage.pckl'), \
        os.path.join(pathload, 'result_km.pckl'), \
        os.path.join(pathload, 'result_info.pckl')
    with bz2.BZ2File(filename_para, 'wb') as f:
        pickle.dump(para, f)
    with bz2.BZ2File(filename_stoptime, 'wb') as f:
        pickle.dump(stoptime, f)
    with bz2.BZ2File(filename_dosage, 'wb') as f:
        pickle.dump(dosage, f)
    with bz2.BZ2File(filename_km, 'wb') as f:
        pickle.dump(km, f)
    with bz2.BZ2File(filename_info, 'wb') as f:
        pickle.dump(info, f)
    return


def DPM_run_drugmis_analysis(mis=None, Strategy_name=None, Simduration=SIMDURATION_DEFAULT_VAL, filename_pattern='result',
                             pathload='./pnas/',
                             pathsave='./mis_oftrue/',
                             pathloadmis='./',
                             samplesize=SAMPLESIZE,
                             nsample=NSAMPLE):

    def DPM_run_drugmis_analysis_1(pathload_, mis_flag=False):
        filename_para = pathload_ + '_para.pckl'
        filename_dosage = pathload_ + '_dosage.pckl'
        filename_stoptime = pathload_ + '_stoptime.pckl'
        filename_setind = pathload_ + '_setind.pckl'
        filename_km = pathload_ + '_km.pckl'
        # filename_info = pathload_ + '_info.pckl'

        if mis_flag is False:
            with bz2.BZ2File(filename_para, 'rb') as f_:
                para_ = pickle.load(f_)
        else:
            para_ = None

        with bz2.BZ2File(filename_dosage, 'rb') as f_:
            dosage_ = pickle.load(f_)

        with bz2.BZ2File(filename_stoptime, 'rb') as f_:
            stoptime_ = pickle.load(f_)

        if not mis_flag:
            with bz2.BZ2File(filename_setind, 'rb') as f_:
                setind_ = pickle.load(f_)
        else:
            setind_ = None

        with bz2.BZ2File(filename_km, 'rb') as f_:
            km_ = pickle.load(f_)

        # if not mis_flag:
        #     with bz2.BZ2File(filename_info, 'rb') as f_:
        #         info_ = pickle.load(f_)
        # else:
        #     info_ = None

        return para_, stoptime_, dosage_, setind_, km_  # info_

    def DPM_run_drugmis_analysis_2(stoptime_, dosage_, hz_ratio, p_inner, hz_inner, firstuse, average_dose_intensity, average_move_num):
        hz_ratio_, p_, _, hz_, firstuse_, average_dose_intensity_, average_move_num_ = \
            DPM_run_drugmis_info(stoptime_, dosage_, Strategy_name, Simduration)

        hz_ratio.append(hz_ratio_), p_inner.append(p_)
        for i_strategy in Strategy_name:
            hz_inner[i_strategy].append(hz_[i_strategy])
            firstuse[i_strategy].append(firstuse_[i_strategy])
            average_dose_intensity[i_strategy].append(average_dose_intensity_[i_strategy])
            average_move_num[i_strategy].append(average_move_num_[i_strategy])
        return hz_ratio, p_inner, hz_inner, firstuse, average_dose_intensity, average_move_num

    def DPM_run_drugmis_analysis_3(result_, f_, mis_=None):
        result_format = {i_strategy: None for i_strategy in Strategy_name}
        hz_ratio_, p_, hz_, firstuse_, average_dose_intensity_, average_move_num_ = result_
        if mis_ is not None:
            hz_ratio_, p_, hz_, firstuse_, average_dose_intensity_, average_move_num_ = \
                hz_ratio_[mis_], p_[mis_], hz_[mis_], firstuse_[mis_], average_dose_intensity_[mis_], average_move_num_[mis_]

        hz_ratio_mean_ = np.mean(hz_ratio_)
        hz_ratio_ci_ = DPM_run_drugmis_ci(hz_ratio_)

        hz_mean_, hz_ci_, firstuse_mean_, average_dose_intensity_mean_, average_move_num_mean_ = \
            deepcopy(result_format), deepcopy(result_format), deepcopy(result_format), deepcopy(result_format), deepcopy(result_format)
        for i_strategy in Strategy_name:
            hz_mean_[i_strategy] = np.mean(hz_[i_strategy])
            firstuse_mean_[i_strategy] = np.mean(firstuse_[i_strategy])
            average_dose_intensity_mean_[i_strategy] = np.mean(average_dose_intensity_[i_strategy])
            average_move_num_mean_[i_strategy] = np.mean(average_move_num_[i_strategy])
            hz_ci_[i_strategy] = DPM_run_drugmis_ci(hz_[i_strategy])

        frate = None
        if f_ == 'fn':
            frate = sum([i_p > PVALUE_FNFP for i_p in p_])/len(p_) * 100
        elif f_ == 'fp':
            frate = sum([i_p < PVALUE_FNFP for i_p in p_])/len(p_) * 100

        return hz_ratio_mean_, hz_ratio_ci_, frate, hz_mean_, hz_ci_, firstuse_mean_, average_dose_intensity_mean_, average_move_num_mean_

    def DPM_run_drugmis_analysis_4(result_Ben_, result_NoBen_, result_Ben_mis_, result_NoBen_mis_, mis_, colname_, filename_, color_, linestyle_,
                                   name_, fig_name_, format_sci=False, plot=True):
        data = list()
        if type(result_Ben_) is dict:
            val_BenDPM = [result_Ben_['DPM2.2']] + [i_val['DPM2.2'] for _, i_val in result_Ben_mis_.items()]
            if type(val_BenDPM[0]) is tuple:
                val_BenDPM = ['('+'{:.2e}'.format(i[0])+',\n'+'{:.2e}'.format(i[1])+')' for i in val_BenDPM]
            else:
                val_BenDPM = ['{:.2e}'.format(i) for i in val_BenDPM] if format_sci else [round(i, 2) for i in val_BenDPM]

            val_BenCPM = [result_Ben_['CPM']] + [i_val['CPM'] for _, i_val in result_Ben_mis_.items()]
            if type(val_BenCPM[0]) is tuple:
                val_BenCPM = ['('+'{:.2e}'.format(i[0])+',\n'+'{:.2e}'.format(i[1])+')' for i in val_BenCPM]
            else:
                val_BenCPM = ['{:.2e}'.format(i) for i in val_BenCPM] if format_sci else [round(i, 2) for i in val_BenCPM]

            val_NoBenDPM = [result_NoBen_['DPM2.2']] + [i_val['DPM2.2'] for _, i_val in result_NoBen_mis_.items()]
            if type(val_NoBenDPM[0]) is tuple:
                val_NoBenDPM = ['('+'{:.2e}'.format(i[0])+',\n'+'{:.2e}'.format(i[1])+')' for i in val_NoBenDPM]
            else:
                val_NoBenDPM = ['{:.2e}'.format(i) for i in val_NoBenDPM] if format_sci else [round(i, 2) for i in val_NoBenDPM]

            val_NoBenCPM = [result_NoBen_['CPM']] + [i_val['CPM'] for _, i_val in result_NoBen_mis_.items()]
            if type(val_NoBenCPM[0]) is tuple:
                val_NoBenCPM = ['('+'{:.2e}'.format(i[0])+',\n'+'{:.2e}'.format(i[1])+')' for i in val_NoBenCPM]
            else:
                val_NoBenCPM = ['{:.2e}'.format(i) for i in val_NoBenCPM] if format_sci else [round(i, 2) for i in val_NoBenCPM]

            data.append(['BenDPM'] + val_BenDPM)
            data.append(['BenCPM'] + val_BenCPM)
            data.append(['NoBenDPM'] + val_NoBenDPM)
            data.append(['NoBenCPM'] + val_NoBenCPM)
        else:
            val_Ben = [result_Ben_] + [i_val for _, i_val in result_Ben_mis_.items()]
            if type(val_Ben[0]) is tuple:
                if format_sci:
                    val_Ben = ['('+'{:.2e}'.format(i[0])+',\n'+'{:.2e}'.format(i[1])+')' for i in val_Ben]
                else:
                    val_Ben = ['(' + str(round(i[0], 3)) + ',\n' + str(round(i[1], 3)) + ')' for i in val_Ben]
            else:
                val_Ben = ['{:.2e}'.format(i) for i in val_Ben] if format_sci else [round(i, 2) for i in val_Ben]

            val_NoBen = [result_NoBen_] + [i_val for _, i_val in result_NoBen_mis_.items()]
            if type(val_NoBen[0]) is tuple:
                if format_sci:
                    val_NoBen = ['(' + '{:.2e}'.format(i[0]) + ',\n' + '{:.2e}'.format(i[1]) + ')' for i in val_NoBen]
                else:
                    val_NoBen = ['(' + str(round(i[0], 3)) + ',\n' + str(round(i[1], 3)) + ')' for i in val_NoBen]
            else:
                val_NoBen = ['{:.2e}'.format(i) for i in val_NoBen] if format_sci else [round(i, 2) for i in val_NoBen]

            data.append(['Ben'] + val_Ben)
            data.append(['NoBen'] + val_NoBen)

        col_names_ = [colname_]
        colname_ = ['control'] + [d.removeprefix('./').removesuffix('_atsim/') for d in mis_]
        col_names_.extend(colname_)
        print(tabulate(data, headers=col_names_, tablefmt='rst', numalign='center'))
        data_df_ = pd.DataFrame(data, columns=col_names_)
        DPM_plot_df2pdf(data_df_, filename_, (['white'], ['lightgray']))

        if plot:
            leg = data_df_.iloc[:, 0].to_list()
            x = data_df_.columns.tolist()[1:]
            data = data_df_.iloc[:, 1:].to_numpy(dtype=float)
            DPM_plot_multi(x, data, leg, name_, color_, linestyle_)
            plt.title(fig_name_ + name)
            pathsave_ = os.path.join(pathsave_fig, filename_[:filename_.find('/', filename_.rfind('/')+1)+1] + fig_name_ + name_ + '.pdf')
            plt.savefig(pathsave_, format='pdf')
            plt.close()
        return

    if not os.path.exists(pathsave):
        os.makedirs(pathsave)
    filename_save = os.path.join(pathsave, 'result.pckl')
    pathload = os.path.join(pathload, filename_pattern)
    # Path where anaylzed results saved.#
    file_exist = os.path.isfile(filename_save)
    # If the simulation results are not yet analyzed, run the analysis. If they are already processed, #
    # skip this step and move on to generating the figures for the paper. #
    if not file_exist:
        # Load simulation reuslts without missepcification.#
        para, stoptime, dosage, setind, km = DPM_run_drugmis_analysis_1(pathload)
        setname = list(setind.keys())[:3]

        Benidx = [i for i, x in enumerate(setind['first 2 diff']) if x]
        NoBenidx = [i for i, x in enumerate(setind['first 2 same']) if x]

        stoptime_mis = {i_mis: None for i_mis in mis}
        dosage_mis = {i_mis: None for i_mis in mis}

        km_mis = {i_set: {i_mis: {i_key: None for i_key in list(km['total'].keys())} for i_mis in mis} for i_set in setname}
        para_sel = {i_mis: {i_strategy: None for i_strategy in Strategy_name} for i_mis in mis}
        survival_numset = {i_survival_diffset: 0 for i_survival_diffset in ['misworse', 'misbetter', 'miseven']}
        survival_num = {i_mis: {i_strategy: deepcopy(survival_numset) for i_strategy in Strategy_name} for i_mis in mis}
        survival_diff = {i_mis: {i_strategy: [] for i_strategy in Strategy_name} for i_mis in mis}
        # Load simulation reuslts with missepcification.#
        for i_mis in mis:
            i_pathload = os.path.join(pathloadmis, i_mis, filename_pattern)
            _, i_stoptime_mis, i_dosage_mis, _, i_km_mis = DPM_run_drugmis_analysis_1(i_pathload, True)
            # Use only the 'total' dataset in the 'mis_km'.#
            i_km_mis = i_km_mis['total']
            stoptime_mis[i_mis] = i_stoptime_mis
            dosage_mis[i_mis] = i_dosage_mis
            for i_set in setname:
                if i_set == 'total':
                    km_mis[i_set][i_mis] = i_km_mis
                    for i_strategy in Strategy_name:
                        for i, i_paramID in enumerate(stoptime['paramID']):
                            assert stoptime['paramID'][i] == i_stoptime_mis['paramID'][i]
                            survival_diff[i_mis][i_strategy].append(stoptime[i_strategy][i] - i_stoptime_mis[i_strategy][i])
                            if stoptime[i_strategy][i] > i_stoptime_mis[i_strategy][i]:
                                survival_num[i_mis][i_strategy]['misworse'] += 1
                            elif stoptime[i_strategy][i] < i_stoptime_mis[i_strategy][i]:
                                survival_num[i_mis][i_strategy]['misbetter'] += 1
                            else:
                                survival_num[i_mis][i_strategy]['miseven'] += 1

                        if i_strategy == 'CPM':
                            i_ind, = np.where(np.array(stoptime[i_strategy]) < Simduration)
                        elif i_strategy == 'DPM2.2':
                            i_ind, = np.where(np.array(stoptime[i_strategy]) > Simduration)
                        else:
                            i_ind = None

                        i_dosage = [i_dose.split(';')[0] for i_dose in dosage[i_strategy]]
                        i_dosage_mis_ = [i_dose.split(';')[0] for i_dose in i_dosage_mis[i_strategy]]
                        i_ind_dose = [x == y for x, y in zip(i_dosage, i_dosage_mis_)]
                        i_ind_dose = np.array([i for i, x in enumerate(i_ind_dose) if x])
                        i_ind = np.intersect1d(i_ind, i_ind_dose)
                        if i_strategy == 'CPM':
                            # Select the virtual patient who received drug 1 as the initial treatment.#
                            i_ind_dose = np.array([i for i, i_dose in enumerate(i_dosage) if i_dose == '(1.0,0.0)'])
                            i_ind = np.intersect1d(i_ind, i_ind_dose)

                        i_paramID = np.array(stoptime['paramID'])[i_ind]
                        i_stoptime_mis_np = np.array(i_stoptime_mis[i_strategy])[i_ind]
                        # 'i_ind' is used to select a virtual patient as an exmaple for plotting the treatment results.#
                        if i_strategy == 'CPM':
                            # For the CPM treatment example, the aim is to illustrate a case without misspecification shorter than the total #
                            # simulation duration. When drug 2 efficacy is misspecified by a factor of 30x, the resulting survival time becomes #
                            # longer than the total simulation duration.#
                            i_ind, = np.where(i_stoptime_mis_np > Simduration)
                        elif i_strategy == 'DPM2.2':
                            # For the DPM treatment example, the aim is to illustrate a case without misspecification longer than the total#
                            # simulation duration. When drug 2 efficacy is misspecified by a factor of 1/30, the resulting survival time becomes#
                            # shorter than the total simulation duration.#
                            i_ind, = np.where(i_stoptime_mis_np < Simduration)
                        # i_ind = random.sample(list(np.arange(len(i_paramID))), 1)
                        if len(i_ind) > 0:
                            i_paramID = i_paramID[i_ind]
                            i_stoptime_mis_np = i_stoptime_mis_np[i_ind]
                            if i_strategy == 'CPM':
                                i_ind = np.argmin(np.array(stoptime[i_strategy])[i_ind])
                            elif i_strategy == 'DPM2.2':
                                i_ind = np.argmin(i_stoptime_mis_np)

                            i_paramID = i_paramID[i_ind]
                            i_ind, = np.where(np.array(stoptime['paramID']) == i_paramID)[0]
                            i_para_sel = {key: value[i_ind] for key, value in para.items()}
                            i_para_sel = DPM_miscellaneous_fillful(i_para_sel)
                        else:
                            i_para_sel = None
                        para_sel[i_mis][i_strategy] = i_para_sel
                else:
                    i_setind = setind[i_set]
                    i_stoptime_mis_sel = deepcopy(i_stoptime_mis)
                    for i_key in i_stoptime_mis_sel:
                        i_stoptime_mis_sel[i_key] = [x for x, y in zip(i_stoptime_mis_sel[i_key], i_setind) if y]

                    hz, p = DPM_analysis_hazard_ratio(i_stoptime_mis_sel, Strategy_name, Simduration)
                    km_mis[i_set][i_mis]['hz_ratio'] = hz
                    km_mis[i_set][i_mis]['p'] = p
                    km_mis[i_set][i_mis]['num'] = len(i_stoptime_mis_sel['paramID'])
                    for i_strategy in Strategy_name:
                        km_mis[i_set][i_mis][i_strategy] = DPM_analysis_KM(i_stoptime_mis_sel[i_strategy], Simduration)

        result_total = {i_strategy: [] for i_strategy in Strategy_name}

        hz_ratio_Ben = []
        hz_ratio_NoBen = []
        p_Ben = []
        p_NoBen = []
        hz_Ben = deepcopy(result_total)
        hz_NoBen = deepcopy(result_total)
        firstuse_Ben = deepcopy(result_total)
        firstuse_NoBen = deepcopy(result_total)
        average_dose_intensity_Ben = deepcopy(result_total)
        average_dose_intensity_NoBen = deepcopy(result_total)
        average_move_num_Ben = deepcopy(result_total)
        average_move_num_NoBen = deepcopy(result_total)

        hz_ratio_Ben_mis = {i_mis: [] for i_mis in mis}
        hz_ratio_NoBen_mis = {i_mis: [] for i_mis in mis}
        p_Ben_mis = {i_mis: [] for i_mis in mis}
        p_NoBen_mis = {i_mis: [] for i_mis in mis}

        hz_Ben_mis = {i_mis: deepcopy(result_total) for i_mis in mis}
        hz_NoBen_mis = {i_mis: deepcopy(result_total) for i_mis in mis}
        firstuse_Ben_mis = {i_mis: deepcopy(result_total) for i_mis in mis}
        firstuse_NoBen_mis = {i_mis: deepcopy(result_total) for i_mis in mis}
        average_dose_intensity_Ben_mis = {i_mis: deepcopy(result_total) for i_mis in mis}
        average_dose_intensity_NoBen_mis = {i_mis: deepcopy(result_total) for i_mis in mis}
        average_move_num_Ben_mis = {i_mis: deepcopy(result_total) for i_mis in mis}
        average_move_num_NoBen_mis = {i_mis: deepcopy(result_total) for i_mis in mis}

        result = {**{i_strategy: [] for i_strategy in Strategy_name}}

        with tqdm(total=nsample, ncols=150, desc='Runing...') as pbar:
            for i in range(nsample):
                i_stoptime_Ben = deepcopy(result)
                i_stoptime_NoBen = deepcopy(result)
                i_dosage_Ben = deepcopy(result)
                i_dosage_NoBen = deepcopy(result)

                i_stoptime_Ben_mis = {i_mis: deepcopy(result) for i_mis in mis}
                i_stoptime_NoBen_mis = {i_mis: deepcopy(result) for i_mis in mis}
                i_dosage_Ben_mis = {i_mis: deepcopy(result) for i_mis in mis}
                i_dosage_NoBen_mis = {i_mis: deepcopy(result) for i_mis in mis}

                i_Benidx = random.sample(Benidx, 2*samplesize)
                i_NoBenidx = random.sample(NoBenidx, 2*samplesize)
                assert set(i_Benidx).intersection(set(i_NoBenidx)) == set()

                i_Benidx_CPM = random.sample(i_Benidx, samplesize)
                i_Benidx_DPM = set(i_Benidx).difference(set(i_Benidx_CPM))
                i_NoBenidx_CPM = random.sample(i_NoBenidx, samplesize)
                i_NoBenidx_DPM = set(i_NoBenidx).difference(set(i_NoBenidx_CPM))

                assert set(i_Benidx_CPM).union(set(i_Benidx_DPM)) == set(i_Benidx)
                assert set(i_NoBenidx_CPM).union(set(i_NoBenidx_DPM)) == set(i_NoBenidx)
                assert len(set(i_Benidx_CPM)) == len(set(i_Benidx_DPM)) == len(set(i_NoBenidx_CPM)) == len(set(i_NoBenidx_DPM)) == samplesize
                assert set(i_Benidx_CPM).intersection(set(i_Benidx_DPM)) == set()
                assert set(i_NoBenidx_CPM).intersection(set(i_NoBenidx_DPM)) == set()

                i_Benidx = {Strategy_name[0]: i_Benidx_CPM, Strategy_name[1]: i_Benidx_DPM}
                i_NoBenidx = {Strategy_name[0]: i_NoBenidx_CPM, Strategy_name[1]: i_NoBenidx_DPM}

                for i_strategy in Strategy_name:
                    i_stoptime_Ben[i_strategy] = [stoptime[i_strategy][ii] for ii in i_Benidx[i_strategy]]
                    i_stoptime_NoBen[i_strategy] = [stoptime[i_strategy][ii] for ii in i_NoBenidx[i_strategy]]

                    i_dosage_Ben[i_strategy] = [dosage[i_strategy][ii] for ii in i_Benidx[i_strategy]]
                    i_dosage_NoBen[i_strategy] = [dosage[i_strategy][ii] for ii in i_NoBenidx[i_strategy]]

                    for i_mis in mis:
                        i_stoptime_Ben_mis[i_mis][i_strategy] = [stoptime_mis[i_mis][i_strategy][ii] for ii in i_Benidx[i_strategy]]
                        i_stoptime_NoBen_mis[i_mis][i_strategy] = [stoptime_mis[i_mis][i_strategy][ii] for ii in i_NoBenidx[i_strategy]]

                        i_dosage_Ben_mis[i_mis][i_strategy] = [dosage_mis[i_mis][i_strategy][ii] for ii in i_Benidx[i_strategy]]
                        i_dosage_NoBen_mis[i_mis][i_strategy] = [dosage_mis[i_mis][i_strategy][ii] for ii in i_NoBenidx[i_strategy]]
                # Case without misspecification.#
                # Benefit group.#
                hz_ratio_Ben, p_Ben, hz_Ben, firstuse_Ben, average_dose_intensity_Ben, average_move_num_Ben = \
                    DPM_run_drugmis_analysis_2(i_stoptime_Ben, i_dosage_Ben, hz_ratio_Ben, p_Ben, hz_Ben, firstuse_Ben,
                                               average_dose_intensity_Ben, average_move_num_Ben)
                # NoBenefit group.#
                hz_ratio_NoBen, p_NoBen, hz_NoBen, firstuse_NoBen, average_dose_intensity_NoBen, average_move_num_NoBen = \
                    DPM_run_drugmis_analysis_2(i_stoptime_NoBen, i_dosage_NoBen, hz_ratio_NoBen, p_NoBen, hz_NoBen, firstuse_NoBen,
                                               average_dose_intensity_NoBen, average_move_num_NoBen)

                # Case with misspecification.#
                for i_mis in mis:
                    hz_ratio_Ben_mis[i_mis], p_Ben_mis[i_mis], hz_Ben_mis[i_mis], firstuse_Ben_mis[i_mis], \
                        average_dose_intensity_Ben_mis[i_mis], average_move_num_Ben_mis[i_mis] = \
                        DPM_run_drugmis_analysis_2(i_stoptime_Ben_mis[i_mis], i_dosage_Ben_mis[i_mis], hz_ratio_Ben_mis[i_mis], p_Ben_mis[i_mis],
                                                   hz_Ben_mis[i_mis], firstuse_Ben_mis[i_mis], average_dose_intensity_Ben_mis[i_mis],
                                                   average_move_num_Ben_mis[i_mis])

                    hz_ratio_NoBen_mis[i_mis], p_NoBen_mis[i_mis], hz_NoBen_mis[i_mis], firstuse_NoBen_mis[i_mis], \
                        average_dose_intensity_NoBen_mis[i_mis], average_move_num_NoBen_mis[i_mis] = \
                        DPM_run_drugmis_analysis_2(i_stoptime_NoBen_mis[i_mis], i_dosage_NoBen_mis[i_mis], hz_ratio_NoBen_mis[i_mis],
                                                   p_NoBen_mis[i_mis], hz_NoBen_mis[i_mis], firstuse_NoBen_mis[i_mis],
                                                   average_dose_intensity_NoBen_mis[i_mis], average_move_num_NoBen_mis[i_mis])
                pbar.update(1)
        result_Ben = hz_ratio_Ben, p_Ben, hz_Ben, firstuse_Ben, average_dose_intensity_Ben, average_move_num_Ben
        result_NoBen = hz_ratio_NoBen, p_NoBen, hz_NoBen, firstuse_NoBen, average_dose_intensity_NoBen, average_move_num_NoBen
        result_Ben_mis = hz_ratio_Ben_mis, p_Ben_mis, hz_Ben_mis, firstuse_Ben_mis, average_dose_intensity_Ben_mis, average_move_num_Ben_mis
        result_NoBen_mis = hz_ratio_NoBen_mis, p_NoBen_mis, hz_NoBen_mis, firstuse_NoBen_mis, average_dose_intensity_NoBen_mis, \
            average_move_num_NoBen_mis
        with bz2.BZ2File(filename_save, 'wb') as f:
            pickle.dump((survival_num, survival_diff, para_sel, km, km_mis, result_Ben, result_NoBen, result_Ben_mis, result_NoBen_mis), f)
    else:
        with bz2.BZ2File(filename_save, 'rb') as f:
            survival_num, survival_diff, para_sel, km, km_mis, result_Ben, result_NoBen, result_Ben_mis, result_NoBen_mis = pickle.load(f)

    pathsave_fig = os.path.join(pathsave, 'figure')
    if not os.path.exists(pathsave_fig):
        os.makedirs(pathsave_fig)

    # Values #
    hz_ratio_Ben, hz_ratio_ci_Ben, falseneg_Ben, hz_Ben, hz_ci_Ben, firstuse_Ben, average_dose_intensity_Ben, average_move_num_Ben = \
        DPM_run_drugmis_analysis_3(result_Ben, 'fn')
    hz_ratio_NoBen, hz_ratio_ci_NoBen, falsepos_NoBen, hz_NoBen, hz_ci_NoBen, firstuse_NoBen, average_dose_intensity_NoBen, \
        average_move_num_NoBen = DPM_run_drugmis_analysis_3(result_NoBen, 'fp')

    result = {i_mis: [] for i_mis in mis}

    hz_ratio_Ben_mis = deepcopy(result)
    hz_ratio_ci_Ben_mis = deepcopy(result)
    falseneg_Ben_mis = deepcopy(result)
    hz_Ben_mis = deepcopy(result)
    hz_ci_Ben_mis = deepcopy(result)
    firstuse_Ben_mis = deepcopy(result)
    average_dose_intensity_Ben_mis = deepcopy(result)
    average_move_num_Ben_mis = deepcopy(result)

    hz_ratio_NoBen_mis = deepcopy(result)
    hz_ratio_ci_NoBen_mis = deepcopy(result)
    falsepos_NoBen_mis = deepcopy(result)
    hz_NoBen_mis = deepcopy(result)
    hz_ci_NoBen_mis = deepcopy(result)
    firstuse_NoBen_mis = deepcopy(result)
    average_dose_intensity_NoBen_mis = deepcopy(result)
    average_move_num_NoBen_mis = deepcopy(result)

    for i_mis in mis:
        hz_ratio_Ben_mis[i_mis], hz_ratio_ci_Ben_mis[i_mis], falseneg_Ben_mis[i_mis], hz_Ben_mis[i_mis], hz_ci_Ben_mis[i_mis], \
            firstuse_Ben_mis[i_mis], average_dose_intensity_Ben_mis[i_mis], average_move_num_Ben_mis[i_mis] \
            = DPM_run_drugmis_analysis_3(result_Ben_mis, 'fn', mis_=i_mis)

        hz_ratio_NoBen_mis[i_mis], hz_ratio_ci_NoBen_mis[i_mis], falsepos_NoBen_mis[i_mis], hz_NoBen_mis[i_mis], hz_ci_NoBen_mis[i_mis], \
            firstuse_NoBen_mis[i_mis], average_dose_intensity_NoBen_mis[i_mis], average_move_num_NoBen_mis[i_mis] \
            = DPM_run_drugmis_analysis_3(result_NoBen_mis, 'fp', mis_=i_mis)

    # Plot #
    # km #
    pathsave_km = os.path.join(pathsave_fig, 'km')
    if not os.path.exists(pathsave_km):
        os.makedirs(pathsave_km)
    setname = ['total', 'first 2 diff', 'first 2 same']
    km_cpm = {i_set: {i: None for i in mis+['control']} for i_set in setname}
    km_dpm = {i_set: {i: None for i in mis + ['control']} for i_set in setname}
    color = ['#D95F02', '#D95F02', '#D95F02', '#1B9E77', '#1B9E77', '#1B9E77', 'k']
    linestyle = ['-', '--', ':', '-', '--', ':', '-']
    xtickstep = 300
    par_multi = {'color': color, 'linestyle': linestyle, 'duration': Simduration, 'xtick step': xtickstep}
    for i_setname in setname:
        par_multi['totalnum'] = km[i_setname]['num']
        km_ref = km[i_setname][Strategy_name[0]]
        km_treat = km[i_setname][Strategy_name[1]]
        hzr = list(km[i_setname]['hz_ratio'].values())[0]
        par = {'color': ['#1B9E77', '#D95F02'], '0': Strategy_name[0], '1': Strategy_name[1],
               'hzr': hzr, 'p': 0, 'duration': Simduration, 'xtick step': xtickstep, 'totalnum': km[i_setname]['num']}
        if i_setname == 'total':
            i_title = 'Fig3 (left) Total'
        elif i_setname == 'first 2 diff':
            i_title = 'Fig3 (middle) Benefit group'
        elif i_setname == 'first 2 same':
            i_title = 'Fig3 (right) Non-Benefit group'
        else:
            raise ValueError('Error.')
        DPM_plot_KM(km_ref, km_treat, par, i_title)
        i_pathsave = os.path.join(pathsave_km, i_title + '.pdf')
        # i_pathsave = os.path.join(pathsave_km, i_setname+'.pdf')
        plt.savefig(i_pathsave, format='pdf')
        plt.close()

        km_cpm[i_setname]['control'] = km_ref
        km_dpm[i_setname]['control'] = km_treat
        for i_mis in mis:
            if km_mis[i_setname][i_mis]['CPM']['median_survival'] == np.inf:
                km_mis[i_setname][i_mis]['CPM']['median_survival'] = Simduration
            if km_mis[i_setname][i_mis]['DPM2.2']['median_survival'] == np.inf:
                km_mis[i_setname][i_mis]['DPM2.2']['median_survival'] = Simduration
            km_cpm[i_setname][i_mis] = km_mis[i_setname][i_mis]['CPM']
            km_dpm[i_setname][i_mis] = km_mis[i_setname][i_mis]['DPM2.2']

        if i_setname == 'total':
            i_title_CPM = 'FigS2A (left) Total,CPM'
            i_title_DPM = 'FigS2A (right) Total,DPM'
        elif i_setname == 'first 2 diff':
            i_title_CPM = 'FigS2B (left) Benefit group,CPM'
            i_title_DPM = 'FigS2B (right) Benefit group,DPM'
        elif i_setname == 'first 2 same':
            i_title_CPM = 'FigS2C (left) Non-Benefit group,CPM'
            i_title_DPM = 'FigS2C (right) Non-Benefit group,DPM'
        else:
            raise ValueError('Error.')
        # FigS2 (left), CPM#
        DPM_plot_KM_multi(km_cpm[i_setname], par_multi, i_title_CPM)
        i_pathsave = os.path.join(pathsave_km, i_title_CPM+'.pdf')
        plt.savefig(i_pathsave, format='pdf')
        plt.close()
        # FigS2 (right), DPM#
        DPM_plot_KM_multi(km_dpm[i_setname], par_multi, i_title_DPM)
        i_pathsave = os.path.join(pathsave_km, i_title_DPM+'.pdf')
        plt.savefig(i_pathsave, format='pdf')
        plt.close()

    # Using one parameter to illustrate the situation misspeciciation make CPM better with drug 2 efficacy increase to 30x. #
    # When drug 2 efficacy is increased by 30×, CPM shows improved performance.#
    ratio = 30
    para = para_sel['./30x_atsim/']['CPM']
    para['Num_drug'] = NUM_DRUG_DEFAULT_VAL
    mis_para = deepcopy(para)
    mis_para['Sa.S.D2.'] = mis_para['Sa.S.D2.'] * ratio
    mis_para['Sa.R1.D2.'] = mis_para['Sa.R1.D2.'] * ratio
    mis_para['Sa.R2.D2.'] = mis_para['Sa.R2.D2.'] * ratio
    mis_para['Sa.R12.D2.'] = mis_para['Sa.R12.D2.'] * ratio

    pathsave_example = os.path.join(pathsave_fig, 'example 30x')
    if not os.path.exists(pathsave_example):
        os.makedirs(pathsave_example)
    DPM_run_plot_1PAR_drugmis(par=para, mis_par=mis_para, Strategy_name=Strategy_name, misspecification_atsim=True, pathsave=pathsave_example)

    # When the efficacy of drug 2 is reduced to 1/30 of its original value, DPM shows diminished performance.#
    ratio = 1/30
    para = para_sel['./div30_atsim/']['DPM2.2']
    para['Num_drug'] = NUM_DRUG_DEFAULT_VAL
    mis_para = deepcopy(para)
    mis_para['Sa.S.D2.'] = mis_para['Sa.S.D2.'] * ratio
    mis_para['Sa.R1.D2.'] = mis_para['Sa.R1.D2.'] * ratio
    mis_para['Sa.R2.D2.'] = mis_para['Sa.R2.D2.'] * ratio
    mis_para['Sa.R12.D2.'] = mis_para['Sa.R12.D2.'] * ratio

    pathsave_example = os.path.join(pathsave_fig, 'example div30')
    if not os.path.exists(pathsave_example):
        os.makedirs(pathsave_example)
    DPM_run_plot_1PAR_drugmis(par=para, mis_par=mis_para, Strategy_name=Strategy_name, misspecification_atsim=True, pathsave=pathsave_example)

    color = ['#D04F3A', '#D04F3A', '#227D88', '#227D88']
    linestyle = ['-', '--', '-', '--']
    # Average drug 2 introduction timestep.#
    name = 'Average drug 2 introduction timestep'
    filename_drug2_intromove = os.path.join(pathsave, f'{name}.pdf')
    col_name = 'Average drug 2 introduction timestep'
    fig_name = 'Fig5A '
    DPM_run_drugmis_analysis_4(firstuse_Ben, firstuse_NoBen, firstuse_Ben_mis, firstuse_NoBen_mis, mis, col_name, filename_drug2_intromove,
                               color, linestyle, name, fig_name)

    # Average drug 2 intensity.#
    name = 'Average drug 2 intensity'
    filename_drug2_intensity = os.path.join(pathsave, f'{name}.pdf')
    col_name = 'Average drug2 intensity'
    fig_name = 'Fig5B '
    DPM_run_drugmis_analysis_4(average_dose_intensity_Ben, average_dose_intensity_NoBen, average_dose_intensity_Ben_mis,
                               average_dose_intensity_NoBen_mis, mis, col_name, filename_drug2_intensity, color, linestyle, name, fig_name)

    # Drug 2 time weighted average dose intensity.#
    name = 'Drug 2 time weighted average dose intensity'
    filename_drug2_intensity_weight_ave = os.path.join(pathsave, f'{name}.pdf')
    col_name = 'Drug 2 time weighted average dose intensity.'
    fig_name = 'Fig5C '
    DPM_run_drugmis_analysis_4(average_move_num_Ben, average_move_num_NoBen, average_move_num_Ben_mis, average_move_num_NoBen_mis, mis, col_name,
                               filename_drug2_intensity_weight_ave, color, linestyle, name, fig_name)

    # Hazard value.#
    name = 'Hazard'
    filename_hz = os.path.join(pathsave, f'{name}.pdf')
    col_name = 'Hazard'
    fig_name = 'Fig4C '
    DPM_run_drugmis_analysis_4(hz_Ben, hz_NoBen, hz_Ben_mis, hz_NoBen_mis, mis, col_name, filename_hz, color, linestyle, name, fig_name,
                               format_sci=True)

    # Confidence interval for the hazard.#
    name = 'CI of Hazard'
    filename_hz_ci = os.path.join(pathsave, f'{name}.pdf')
    col_name = 'CI of Hazard'
    fig_name = 'not used'
    DPM_run_drugmis_analysis_4(hz_ci_Ben, hz_ci_NoBen, hz_ci_Ben_mis, hz_ci_NoBen_mis, mis, col_name, filename_hz_ci, color, linestyle, name,
                               fig_name, format_sci=True, plot=False)

    color = ['#D04F3A', '#227D88']
    linestyle = ['-', '-', '--']
    # Hazard ratio.#
    name = 'Hazard ratio'
    filename_hzratio = os.path.join(pathsave, f'{name}.pdf')
    col_name = 'Hazard ratio'
    fig_name = 'Fig4B '
    DPM_run_drugmis_analysis_4(hz_ratio_Ben, hz_ratio_NoBen, hz_ratio_Ben_mis, hz_ratio_NoBen_mis, mis, col_name, filename_hzratio, color,
                               linestyle, name, fig_name,)

    # Confidence interval for the hazard ratio.#
    name = 'CI of Hazard ratio'
    filename_hzratio = os.path.join(pathsave, f'TableS1 {name}.pdf')
    col_name = 'CI of Hazard ratio'
    fig_name = 'not used'
    DPM_run_drugmis_analysis_4(hz_ratio_ci_Ben, hz_ratio_ci_NoBen, hz_ratio_ci_Ben_mis, hz_ratio_ci_NoBen_mis, mis, col_name, filename_hzratio,
                               color, linestyle, name, fig_name, plot=False)

    # False positive and negative rate.#
    name = 'False positive and negative rate'
    filename_falserate = os.path.join(pathsave, f'{name}.pdf')
    col_name = 'False positive and negative rate'
    fig_name = 'Fig4A '
    DPM_run_drugmis_analysis_4(falseneg_Ben, falsepos_NoBen, falseneg_Ben_mis, falsepos_NoBen_mis, mis, col_name, filename_falserate,
                               color, linestyle, name, fig_name, format_sci=False)
    return


def DPM_run_drugmis_info(stoptime, dosage, Strategy_name, Simduration):
    hz_ratio, p = DPM_analysis_hazard_ratio(stoptime, Strategy_name, Simduration)
    hz_ratio, p = list(hz_ratio.values())[0], list(p.values())[0]
    firstuse = {i_strategy: None for i_strategy in Strategy_name}
    average_dose_intensity = {i_strategy: None for i_strategy in Strategy_name}
    average_move_num = {i_strategy: None for i_strategy in Strategy_name}
    hz = {i_strategy: None for i_strategy in Strategy_name}
    km = {i_strategy: None for i_strategy in Strategy_name}
    for i_strategy in Strategy_name:
        firstuse[i_strategy], average_dose_intensity[i_strategy], average_move_num[i_strategy] = \
            DPM_run_drugmis_analysisdose(dosage[i_strategy], i_strategy)

        km[i_strategy] = DPM_analysis_KM(stoptime[i_strategy], Simduration)
        hz[i_strategy] = km[i_strategy]['hazard']['mean'].tolist()[0]
    return hz_ratio, p, km, hz, firstuse, average_dose_intensity, average_move_num


def DPM_run_drugmis_analysisdose(dose, strategyname, inddrug=1):
    firstuse, max_num_change, num_change, average_move_num, average_dose_intensity = [], [], [], [], []
    use_atbegin = 0
    for i, i_dose in enumerate(dose):
        i_dose = i_dose.split(';')
        if '-1' in i_dose:
            i_dose = i_dose[:i_dose.index('-1')]
        i_firstuse, i_average_move_num, drugovermove, drugtotal, i_num_change,  i_current = None, None, 0, 0, 0, i_dose[0]
        if i_current == '(0.0,1.0)':
            use_atbegin += 1

        for j, j_step in enumerate(i_dose):
            i_val = [float(i_val) for i_val in j_step[1:-2].split(',')]
            drugtotal = drugtotal + i_val[inddrug]
            drugovermove = drugovermove + (j+1) * i_val[inddrug]
            if strategyname == 'CPM' and i_current not in ['(0.0,1.0)', '(0.0,0.0)', '(1.0,0.0)']:
                assert strategyname == 'CPM' and i_current not in ['(0.0,1.0)', '(0.0,0.0)', '(1.0,0.0)']
            if i_firstuse is None and i_val[inddrug] != 0:
                i_firstuse = j+1
            if j_step != i_current and j_step != '(0.0,0.0)':
                i_num_change += 1
                i_current = j_step

        i_average_move_num = drugovermove/drugtotal if drugtotal != 0 else None
        i_firstuse = len(i_dose) + 1 if i_firstuse is None else i_firstuse
        i_average_move_num = len(i_dose) + 1 if i_average_move_num is None else i_average_move_num
        i_average_dose_intensity = drugtotal/len(i_dose)

        firstuse.append(i_firstuse)
        average_move_num.append(i_average_move_num)
        num_change.append(i_num_change/len(i_dose))
        max_num_change.append(i_num_change)
        average_dose_intensity.append(i_average_dose_intensity)

    return np.mean(firstuse), np.mean(average_dose_intensity), np.mean(average_move_num)


def DPM_run_drugmis_ci(data, n_iterations=N_ITERATIONS, alpha=ALPHA_CI):
    n_size, means = int(len(data)), []
    with tqdm(total=n_iterations, ncols=150, desc='Runing...') as pbar:
        for i in range(n_iterations):
            s = resample(data, n_samples=n_size)
            means.append(np.mean(s))
            pbar.update(1)

    p = ((1.0 - alpha)/2.0) * 100
    lower = np.percentile(means, p)

    p = (alpha + ((1.0 - alpha)/2.0)) * 100
    upper = np.percentile(means, p)
    assert lower < upper
    return lower, upper
