from DPM_lib import deepcopy, plt, pd, KaplanMeierFitter, CoxPHFitter, statistics, ExponentialFitter, tabulate, itertools, os, np, \
    bz2, pickle, tqdm, random, sns, supervenn, rv_continuous, loguniform, pearsonr
from DPM_plot import DPM_plot_KM, DPM_plot_KM_multi, DPM_plot_KM_multi2, DPM_plot_contour, DPM_plot_LOD_multi, DPM_plot_hz_ratio, DPM_plot_hz, \
    DPM_plot_misspec_pop
from DPM_constant import FIG_DPI, PLT_INTERACTIVE, LOD_LIST, MISSPECIFICATION_LOD_STR, LOGBASE, MUTATION_RATE, ALL_POSSIBLE_CELLTYPE_2DRUG
from DPM_generate import DPM_generate_pdf


def DPM_analysis_stat(stoptime, Strategy_name, Simduration, bool_select=None, titlestr='All'):
    num_betterall, num_sigbetterall, num_sigbetter, num_sigworse, paramID_sigbetter, paramID_sigworse = \
        dict(), dict(), dict(), dict(), dict(), dict()
    if bool_select is not None:
        stoptime = {item: list(itertools.compress(value, bool_select)) for (item, value) in stoptime.items()}
    df_stoptime = pd.DataFrame.from_dict(stoptime, dtype=int)
    for i_strategy in Strategy_name:
        df_stoptime[i_strategy] = df_stoptime[i_strategy].apply(lambda x: x if x < Simduration else Simduration)
    survival_median = df_stoptime.median()[Strategy_name].astype(int)
    survival_mean = df_stoptime.mean()[Strategy_name].apply(np.ceil).astype(int)
    survival_5y = df_stoptime[Strategy_name].apply(lambda x: x >= Simduration).apply(sum) / df_stoptime.shape[0] * 100
    survival_5y = survival_5y.round(2)
    for i, i_strategy in enumerate(df_stoptime[Strategy_name]):
        sub_Strategy_name = [ii for ii in Strategy_name if ii is not i_strategy]
        # No. of cases strategy numerically better than all others
        num_betterall[i_strategy] = \
            np.array(df_stoptime[i_strategy] > df_stoptime[sub_Strategy_name].max(axis=1)).sum()
        # No. of cases strategy significantly better than all others
        num_sigbetterall[i_strategy] = \
            np.array((df_stoptime[i_strategy] > df_stoptime[sub_Strategy_name].max(axis=1) + 8 * 7) &
                     (df_stoptime[i_strategy] > 1.25 * df_stoptime[sub_Strategy_name].max(axis=1))).sum()
        i_sigbetter, i_sigworse, i_paramID_sigbetter, i_paramID_sigworse = dict(), dict(), dict(), dict()
        i_stoptime = df_stoptime[i_strategy]
        i_sigbetter[i_strategy], i_sigworse[i_strategy], i_paramID_sigbetter[i_strategy], i_paramID_sigworse[i_strategy] = \
            'N.A.', 'N.A.', 'N.A.', 'N.A.'
        for i_substrategy in sub_Strategy_name:
            bool_sigbetter = (df_stoptime[i_substrategy] > i_stoptime + 8 * 7) & (df_stoptime[i_substrategy] > 1.25 * i_stoptime)
            i_sigbetter[i_substrategy] = np.array(bool_sigbetter).sum()
            i_paramID_sigbetter[i_substrategy] = bool_sigbetter.index[bool_sigbetter]

            bool_sigworse = (i_stoptime > df_stoptime[i_substrategy] + 8 * 7) & (i_stoptime > 1.25 * df_stoptime[i_substrategy])
            i_sigworse[i_substrategy] = np.array(bool_sigworse).sum()
            i_paramID_sigworse[i_substrategy] = bool_sigworse.index[bool_sigworse]

        i_sigbetter = {k: i_sigbetter[k] for k in Strategy_name}
        i_sigworse = {k: i_sigworse[k] for k in Strategy_name}
        i_paramID_sigbetter = {k: i_paramID_sigbetter[k] for k in Strategy_name}
        i_paramID_sigworse = {k: i_paramID_sigworse[k] for k in Strategy_name}

        num_sigbetter[i_strategy] = i_sigbetter
        num_sigworse[i_strategy] = i_sigworse
        paramID_sigbetter[i_strategy] = i_paramID_sigbetter
        paramID_sigworse[i_strategy] = i_paramID_sigworse

    data = list()
    data.append(['Median survival, day'] + list(survival_median.values))
    data.append(['Mean survival, day'] + list(survival_mean.values))
    data.append(['Survival at 5y, %'] + list(survival_5y.values))
    data.append(['No. of cases strategy numerically better than all others'] + list(num_betterall.values()))
    data.append(['No. of cases strategy significantly better than all others'] + list(num_sigbetterall.values()))
    # for i_strategy in Strategy_name:
    #     i_num = num_sigbetter[i_strategy]
    #     data.append([f'No. of cases significantly better than {i_strategy}'] + list(i_num.values()))
    # data.append(['Median survival, day'] + list(survival_median.values))
    # data.append(['Mean survival, day'] + list(survival_mean.values))
    # data.append(['Survival at 5y, %'] + list(survival_5y.values))
    # data.append(['No. of cases strategy numerically better than all others'] + list(num_betterall.values()))
    # data.append(['No. of cases strategy significantly better than all others'] + list(num_sigbetterall.values()))
    for i_strategy in Strategy_name:
        i_num_sigbetter, i_num_sigworse = num_sigbetter[i_strategy], num_sigworse[i_strategy]
        data.append([f'No. of cases significantly better than {i_strategy}'] + list(i_num_sigbetter.values()))
        data.append([f'No. of cases significantly worse than {i_strategy}'] + list(i_num_sigworse.values()))

    col_names = [f'Patients ({titlestr})']
    col_names.extend(Strategy_name)
    print(tabulate(data, headers=col_names, tablefmt='rst', numalign='left'))

    # df_data = pd.DataFrame(data, columns=[f'Patients ({titlestr})']+Strategy_name)

    # Significantly better means at least 8 week of absolute improvement
    par = {'duration': Simduration, 'binsize': 14, 'xtick step': 300}
    for i in list(itertools.combinations(Strategy_name, 2)):
        if 'strategy0' in i:
            ref = stoptime['strategy0']
            treat = stoptime[list(filter(lambda x: x != 'strategy0', i))[0]]
            name = ('strategy0', list(filter(lambda x: x != 'strategy0', i))[0])
        else:
            ref, treat = stoptime[i[0]], stoptime[i[1]]
            name = i
        # i_ID_sigbetter = paramID_sigbetter[name[0]][name[1]]
        # i_ID_sigworese = paramID_sigworse[name[0]][name[1]]
        par['name'] = name
        # par['id_sigbetter'] = i_ID_sigbetter
        # par['id_sigworse'] = i_ID_sigworese
        DPM_plot_contour(ref, treat, par, titlestr)
        # DPM_plot_density2d_surface(ref, treat, par)

    km = dict()
    for i_strategy in Strategy_name:
        i_stoptime = stoptime[i_strategy]
        i_stoptime = [i_val if i_val <= Simduration else Simduration + 1 for i_val in i_stoptime]
        km[i_strategy] = DPM_analysis_KM(i_stoptime, Simduration)

    # plot all KM
    par = {'duration': Simduration, 'xtick step': 300, 'totalnum': df_stoptime.shape[0]}
    DPM_plot_KM_multi(km, par, titlestr)

    p = DPM_analysis_pairwise_logrank_test(stoptime, Strategy_name, Simduration)
    hz = dict()
    for i in list(itertools.combinations(Strategy_name, 2)):
        if 'strategy0' in i:
            ref = stoptime['strategy0']
            treat = stoptime[list(filter(lambda x: x != 'strategy0', i))[0]]
            name = ('strategy0', list(filter(lambda x: x != 'strategy0', i))[0])
        else:
            ref, treat = stoptime[i[0]], stoptime[i[1]]
            name = i
        hz[name] = DPM_analysis_HZ(ref, treat, Simduration)

    for i in hz.keys():
        idx_p = [idx for idx, i_val in enumerate(p.keys()) if set(i) == set(i_val)][0]
        km_ref = km[i[0]]
        km_treat = km[i[1]]
        par = {'color': ['k', 'b'], '0': i[0], '1': i[1], 'hzr': hz[i], 'p': p[list(p.keys())[idx_p]],
               'duration': Simduration, 'xtick step': 300, 'totalnum': df_stoptime.shape[0]}
        DPM_plot_KM(km_ref, km_treat, par, titlestr)
        plt.close('all')
    return


def DPM_analysis_misspec_pop(para, pop, info, LOD, pathsave, pathload):
    class MutationFraction(rv_continuous):
        def _pdf(self, x, k, LOD, const):
            return (1.0/const)*DPM_generate_pdf(x, k, LOD)

    def DPM_analysis_pop_1(pop_, mutationFraction_distribution_, LOD_, pathsave_):
        bins = np.insert(np.logspace(-9, np.log10(i_LOD), num=int(1e3)), 0, 0)
        for i, i_cell in enumerate(ALL_POSSIBLE_CELLTYPE_2DRUG):
            plt.rcParams['font.size'] = 9
            fig, ax = plt.subplots(6, 7)
            fig.set_size_inches(15, 10)
            plt.tight_layout(pad=0.9, h_pad=0.9, w_pad=0.9)
            pearson_pdf = list()
            pearson_loguni = list()
            with tqdm(total=len(pop_), ncols=150) as pbar:
                for j, j_pop_ in enumerate(pop_):
                    j_ax = ax[j // 7, j % 7]
                    j_pop_ = np.asarray(j_pop_)
                    j_total = np.sum(j_pop_, axis=1)
                    j_pop_ = j_pop_[:, i]/j_total
                    ind = np.flatnonzero(j_pop_ < LOD_)
                    j_pop_ = j_pop_[ind]

                    pdf_samples = mutationFraction_distribution_.rvs(k=MUTATION_RATE, LOD=i_LOD, const=1, size=len(j_pop_))
                    loguni_samples = loguniform.rvs(MUTATION_RATE, i_LOD, size=len(j_pop_))

                    hist, _ = np.histogram(j_pop_, bins=bins, weights=[1/len(j_pop_)]*len(j_pop_))
                    hist_pdf, _ = np.histogram(pdf_samples, bins, weights=[1/len(j_pop_)]*len(j_pop_))
                    hist_loguni, _ = np.histogram(loguni_samples, bins, weights=[1/len(j_pop_)]*len(j_pop_))

                    # plt.stairs(hist, bins, color='r', lw=1, label='sim')
                    # plt.stairs(hist_pdf, bins, color='r', lw=1, label='sim')
                    # plt.stairs(hist_loguni, bins, color='r', lw=1, label='sim')

                    pearson_pdf.append(pearsonr(hist, hist_pdf)[0])
                    pearson_loguni.append(pearsonr(hist, hist_loguni)[0])

                    DPM_plot_misspec_pop(j_ax, j, hist, hist_pdf, hist_loguni, bins, i_cell, legend_fontsize=8)
                    pbar.update(1)

                ax[-1, -1].plot(np.arange(len(pearson_pdf)), pearson_pdf, color='k', lw=1, label='pdf')
                ax[-1, -1].plot(np.arange(len(pearson_loguni)), pearson_loguni, color='b', lw=1, label='loguni')
                ax[-1, -1].set_ylabel('Pearson corr')
                ax[-1, -1].set_xlabel('Step')
                ax[-1, -1].legend(loc='best', prop={'size': 8}, frameon=False, ncol=1)
            i_pathsave_ = os.path.join(pathsave_, f"{i_cell} cell.pdf")
            plt.savefig(i_pathsave_, format='pdf', bbox_inches='tight')
            plt.close('all')
        return

    paramID_sigbetter = info['total']['paramID sigbetter']
    pathsave = os.path.join(pathsave, 'pop')
    for i_LOD in LOD:
        i_LOD = float(i_LOD)
        mutationFraction_distribution = MutationFraction(name='mutationFraction_distribution', a=MUTATION_RATE, b=i_LOD)
        # i_pdf = mutationFraction_distribution.pdf(x=bins, k=MUTATION_RATE, LOD=i_LOD, const=1)
        # i_loguni = loguniform.pdf(bins, MUTATION_RATE, i_LOD)
        for _, i_strategy in enumerate(pop):
            i_pathsave = os.path.join(pathsave, str(i_LOD), i_strategy)
            if not os.path.exists(i_pathsave):
                os.makedirs(i_pathsave)
            i_pop = pop[i_strategy]
            DPM_analysis_pop_1(i_pop, mutationFraction_distribution, i_LOD, i_pathsave)

    para_df = pd.DataFrame.from_dict(para)
    para_df['X'] = para_df.loc[:, ['Spop', 'R1pop', 'R2pop', 'R12pop']].sum(axis=1)
    print('{:.5e}'.format(para_df['X'].min()))
    print('{:.5e}'.format(para_df['X'].max()))
    para_df['Spop per'] = para_df['Spop'] / para_df['X']
    para_df['R1pop per'] = para_df['R1pop'] / para_df['X']
    para_df['R2pop per'] = para_df['R2pop'] / para_df['X']
    para_df['R12pop per'] = para_df['R12pop'] / para_df['X']
    assert set(para_df['R12pop per']) == {0}
    para_df['R1pop per'].min()
    para_df['R1pop per'].max()
    para_df['R2pop per'].min()
    para_df['R2pop per'].max()
    return


def DPM_analysis_misspec(km, info, km_mis, info_mis, LOD, setname, Strategy_name, Simduration, pathsave, plot):
    km_mis = dict(zip(MISSPECIFICATION_LOD_STR, km_mis))
    info_mis = dict(zip(MISSPECIFICATION_LOD_STR, info_mis))

    LOD_plot = LOD_LIST
    LOD_all = LOD + ['nomis']

    keys = []
    for i in MISSPECIFICATION_LOD_STR:
        keys.extend([i_LOD + ' ' + i for i_LOD in LOD if i_LOD in LOD_plot])
    keys.extend([LOD_all[-1]])

    for i, i_LOD in enumerate(LOD):
        for i_setname in setname:
            if plot:
                path_folder = os.path.join(pathsave, i_setname)
                if not os.path.exists(path_folder):
                    os.makedirs(path_folder)
                # True
                filename = os.path.join(path_folder, 'nomis.pdf')
                if not os.path.exists(filename):
                    titlestr = i_setname + ', ' + LOD_all[-1]
                    DPM_analysis_misspec_plot(Strategy_name, Simduration, km[i_setname], titlestr)
                    plt.savefig(filename, dpi=FIG_DPI)
                    plt.close('all')
                for i_mis in MISSPECIFICATION_LOD_STR:
                    if i_LOD in LOD_plot:
                        filename = os.path.join(path_folder, 'LOD ' + i_LOD + ' ' + i_mis + '.pdf')
                        titlestr = i_setname + ', set by ' + i_mis + ', LOD: ' + i_LOD
                        DPM_analysis_misspec_plot(Strategy_name, Simduration, km_mis[i_mis][i][i_setname], titlestr)
                        plt.savefig(filename, dpi=FIG_DPI)
                        plt.close('all')

    for i_setname in setname:
        result = {i_strategy: [] for i_strategy in Strategy_name}
        median_survial, hazard, firstd2, movenum_d2, numdchange = [], [], [], [], []
        for _ in range(len(MISSPECIFICATION_LOD_STR)):
            median_survial.append(deepcopy(result))
            hazard.append(deepcopy(result))
            firstd2.append(deepcopy(result))
            movenum_d2.append(deepcopy(result))
            numdchange.append(deepcopy(result))

        km_set_i = km[i_setname]
        info_i_set = info[i_setname]

        km_mis_set_i, hz_ratio_set_i, km_mis_set_i_num, info_mis_set_i = [], [], [], []
        for i, (key_i, km_mis_i) in enumerate(km_mis.items()):
            km_mis_i_set_i = [i_val[i_setname] for i_val in km_mis_i]
            km_mis_set_i.append(km_mis_i_set_i)

            info_mis_i = info_mis[key_i]
            info_mis_i_set_i = [i_val[i_setname] for i_val in info_mis_i]
            info_mis_set_i.append(info_mis_i_set_i)

            hz_ratio_mis_i_set_i = [i_km_mis_i_set_i['hz_ratio'] for i_km_mis_i_set_i in km_mis_i_set_i]
            hz_ratio_mis_i_set_i.append(km_set_i['hz_ratio'])
            hz_ratio_set_i.append(hz_ratio_mis_i_set_i)

            km_mis_i_set_i_num = [i_km_mis_i_set_i['num'] for i_km_mis_i_set_i in km_mis_i_set_i]
            km_mis_set_i_num.append(km_mis_i_set_i_num)

        path_folder = os.path.join(pathsave, i_setname)
        if not os.path.exists(path_folder):
            os.makedirs(path_folder)

        for i_strategy in Strategy_name:
            km_i_set_strategy, km_i_set_strategy_plot, info_i_set_strategy = [], [], []
            for i in range(len(MISSPECIFICATION_LOD_STR)):
                i_km_i_set_strategy = [i_val[i_strategy] for i_val in km_mis_set_i[i]]
                km_i_set_strategy.append(i_km_i_set_strategy)
                i_km_i_set_strategy_plot = [i_km_i_set_strategy[ii] for ii in range(len(i_km_i_set_strategy)) if LOD[i] in LOD_plot]
                km_i_set_strategy_plot.append(i_km_i_set_strategy_plot)
                info_i_set_strategy.append([i_info_0_i_set[i_strategy] for i_info_0_i_set in info_mis_set_i[i]])

            ms_i_set_strategy, hazard_i_set_strategy, firstd2_i_set_strategy, movenum_d2_i_set_strategy, \
                ave_numdrugchange_i_set_strategy = [], [], [], [], []

            for i in range(len(MISSPECIFICATION_LOD_STR)):
                i_ms_i_set_strategy = [i_km_i_set_strategy['median_survival'] for i_km_i_set_strategy in km_i_set_strategy[i]]
                i_ms_i_set_strategy.append(km_set_i[i_strategy]['median_survival'])
                ms_i_set_strategy.append(i_ms_i_set_strategy)
                median_survial[i][i_strategy] = i_ms_i_set_strategy

                i_hazard_i_set_strategy = [i_km_i_set_strategy['hazard'] for i_km_i_set_strategy in km_i_set_strategy[i]]
                i_hazard_i_set_strategy.append(km_set_i[i_strategy]['hazard'])
                hazard_i_set_strategy.append(i_hazard_i_set_strategy)
                hazard[i][i_strategy] = i_hazard_i_set_strategy

                i_firstd2_i_set_strategy = [i_info_i_set_strategy['first drug 2'] for i_info_i_set_strategy in info_i_set_strategy[i]]
                i_firstd2_i_set_strategy.append(info_i_set[i_strategy]['first drug 2'])
                firstd2_i_set_strategy.append(i_firstd2_i_set_strategy)
                firstd2[i][i_strategy] = i_firstd2_i_set_strategy

                i_movenum_d2_i_set_strategy = [i_info_i_set_strategy['average move number'] for i_info_i_set_strategy in info_i_set_strategy[i]]
                i_movenum_d2_i_set_strategy.append(info_i_set[i_strategy]['average move number'])
                movenum_d2_i_set_strategy.append(i_movenum_d2_i_set_strategy)
                movenum_d2[i][i_strategy] = i_movenum_d2_i_set_strategy

                i_ave_numdrugchange_i_set_strategy = [i_info_i_set_strategy['drug changes'] for i_info_i_set_strategy in info_i_set_strategy[i]]
                i_ave_numdrugchange_i_set_strategy.append(info_i_set[i_strategy]['drug changes'])
                ave_numdrugchange_i_set_strategy.append(i_ave_numdrugchange_i_set_strategy)
                numdchange[i][i_strategy] = i_ave_numdrugchange_i_set_strategy

            i_km = {i_key: None for i_key in keys}
            for i, i_LOD in enumerate(LOD):
                if i_LOD in LOD_plot:
                    for ii, i_mis in enumerate(MISSPECIFICATION_LOD_STR):
                        i_key = i_LOD + ' ' + i_mis
                        i_km[i_key] = {**km_i_set_strategy_plot[ii][i], **{'num': km_mis_set_i_num[ii][i]}}

            i_km[keys[-1]] = {**km_set_i[i_strategy], 'num': km_set_i['num']}

            titlestr = i_setname + ', ' + i_strategy
            par = {'duration': Simduration, 'xtick step': 300}
            DPM_plot_KM_multi2(i_km, par, titlestr)
            plt.savefig(os.path.join(path_folder, i_strategy + '.pdf'), dpi=FIG_DPI)
            plt.close('all')

        titlestr = i_setname
        DPM_plot_LOD_multi(median_survial, LOD_all, titlestr)
        plt.savefig(os.path.join(path_folder, 'median survival.pdf'), dpi=FIG_DPI)
        plt.close('all')

        DPM_plot_hz_ratio(hz_ratio_set_i, LOD_all, titlestr)
        plt.savefig(os.path.join(path_folder, 'hz ratio.pdf'), dpi=FIG_DPI)
        plt.close('all')

        DPM_plot_hz(hazard, LOD_all, titlestr)
        plt.savefig(os.path.join(path_folder, 'hz.pdf'), dpi=FIG_DPI)
        plt.close('all')

        DPM_plot_LOD_multi(firstd2, LOD_all, titlestr, ylablestr='first drug2')
        plt.savefig(os.path.join(path_folder, 'first drug2.pdf'), dpi=FIG_DPI)
        plt.close('all')

        DPM_plot_LOD_multi(movenum_d2, LOD_all, titlestr, ylablestr='Average move number normalized to dose of drug 2')
        plt.savefig(os.path.join(path_folder, 'ave movenum drug2.pdf'), dpi=FIG_DPI)
        plt.close('all')

        DPM_plot_LOD_multi(numdchange, LOD_all, titlestr, ylablestr='Average number of drug changes divided by the length of the clinical course')
        plt.savefig(os.path.join(path_folder, 'ave num drug change.pdf'), dpi=FIG_DPI)
        plt.close('all')
    return


def DPM_analysis_misspec_par(para, info, info_mis, LOD, pathsave, pathload):
    def DPM_analysis_misspec_par_1(i_value, i_paraval):
        num = [i_value.count(ii) for ii in i_paraval]
        assert sum(num) == len(i_value)
        return [i_num / len(i_value) * 100 for i_num in num]

    def DPM_analysis_misspec_par_2(name):
        bins = [0, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.9, 1]
        ind = para_df[name] == 0
        ind = np.flatnonzero(ind)
        para_df.loc[ind, name] = 0
        for ii in range(len(bins)-1):
            ind = (bins[ii] < para_df[name]) & (para_df[name] <= bins[ii+1])
            ind = np.flatnonzero(ind)
            para_df.loc[ind, name] = bins[ii+1]

    paramID = para['paramID']
    # paramID_first2diff = [i for i, _ in enumerate(para['paramID']) if setind['first 2 diff'][i]]

    if not os.path.exists(os.path.join(pathsave, 'par')):
        os.makedirs(os.path.join(pathsave, 'par'))

    try:
        with bz2.BZ2File(os.path.join(pathload, 'para_index.pckl'), 'rb') as f:
            paramID_sigbetter, paramID_nosigbetter, paramID_mis_losesig, paramID_mis_gainsig = pickle.load(f)
    except FileNotFoundError:
        paramID_sigbetter, paramID_nosigbetter = [], []
        paramID_mis_losesig = {i_mis: {i_LOD: [] for i_LOD in LOD} for i_mis in MISSPECIFICATION_LOD_STR}
        paramID_mis_gainsig = deepcopy(paramID_mis_losesig)
        with tqdm(total=len(paramID), ncols=150) as pbar:
            for i, i_paramID in enumerate(paramID):
                if i_paramID in info['total']['paramID sigbetter']:
                    paramID_sigbetter.append(i)
                    for j, j_mis in enumerate(MISSPECIFICATION_LOD_STR):
                        for k, k_LOD in enumerate(LOD):
                            jk_info_mis = info_mis[j][k]['total']['paramID sigbetter']
                            if i_paramID not in jk_info_mis:
                                paramID_mis_losesig[j_mis][k_LOD].append(i)
                else:
                    paramID_nosigbetter.append(i)
                    for j, j_mis in enumerate(MISSPECIFICATION_LOD_STR):
                        for k, k_LOD in enumerate(LOD):
                            jk_info_mis = info_mis[j][k]['total']['paramID sigbetter']
                            if i_paramID in jk_info_mis:
                                paramID_mis_gainsig[j_mis][k_LOD].append(i)
                pbar.update(1)
        assert set(paramID_sigbetter + paramID_nosigbetter) == set(range(len(paramID)))
        with bz2.BZ2File(os.path.join(pathload, 'para_index.pckl'), 'wb') as f:
            pickle.dump((paramID_sigbetter, paramID_nosigbetter, paramID_mis_losesig, paramID_mis_gainsig), f)

    para['X'] = [sum([para['Spop'][i], para['R1pop'][i], para['R2pop'][i], para['R12pop'][i]]) for i in range(len(para['Spop']))]
    para_df = pd.DataFrame.from_dict(para)
    para_df['Spop'], para_df['R1pop'], para_df['R2pop'],  para_df['R12pop'] = \
        para_df['Spop']/para_df['X'], para_df['R1pop']/para_df['X'], para_df['R2pop']/para_df['X'], \
        para_df['R12pop']/para_df['X']

    path_folder = os.path.join(pathsave, 'set')
    if not os.path.exists(path_folder):
        os.makedirs(path_folder)
    set_total = set(range(len(para['paramID'])))
    set_sigbetter = set(paramID_sigbetter)
    for i_mis in MISSPECIFICATION_LOD_STR:
        sets = [set_total, set_sigbetter]
        labels = ['total', 'sigbetter']
        for i_LOD in LOD:
            sets.append(set(paramID_mis_losesig[i_mis][i_LOD]))
            labels.append('lose sigbetter ' + i_LOD)
        plt.figure(figsize=(18, 8))
        supervenn(sets, labels, widths_minmax_ratio=0.1, col_annotations_area_height=1.2, side_plots='right', chunks_ordering='minimize gaps')
        plt.title('mis ' + i_mis, fontsize=18)
        plt.ylabel('Sets', fontsize=18)
        plt.xlabel('Items', fontsize=18)
        plt.tight_layout()
        plt.savefig(os.path.join(path_folder, 'set ' + i_mis + '.pdf'), dpi=FIG_DPI)
        plt.close('all')

    ind_par = 1
    DPM_analysis_misspec_par_2('Spop'), DPM_analysis_misspec_par_2('R1pop'), DPM_analysis_misspec_par_2('R2pop')

    result = {'total': None, 'sigbetter': None, **dict(zip(LOD, [None]*len(LOD)))}
    val = {i_mis: {i_para: deepcopy(result) for i_para in list(para.keys())[ind_par:]} for i_mis in MISSPECIFICATION_LOD_STR}
    for i_mis in MISSPECIFICATION_LOD_STR:
        for i, i_para in enumerate(list(para.keys())[ind_par:]):
            if i_para in ['X', 'Spop', 'R1pop', 'R2pop', 'R12pop', 'g0_S', 'T.R1..S.', 'T.R2..S.']:

                i_paraval_sorted = sorted(set(para_df[i_para]))
                val[i_mis][i_para]['total'] = DPM_analysis_misspec_par_1(para_df[i_para].values.tolist(), i_paraval_sorted)

                i_paraval_sigbetter = para_df.iloc[paramID_sigbetter, para_df.columns.get_loc(i_para)].values.tolist()
                val[i_mis][i_para]['sigbetter'] = DPM_analysis_misspec_par_1(i_paraval_sigbetter, i_paraval_sorted)

                for i_LOD in LOD:
                    i_parval_mis = para_df.iloc[paramID_mis_losesig[i_mis][i_LOD], para_df.columns.get_loc(i_para)].values.tolist()
                    val[i_mis][i_para][i_LOD] = DPM_analysis_misspec_par_1(i_parval_mis, i_paraval_sorted)
            else:
                val[i_mis][i_para]['total'] = para[i_para]
                val[i_mis][i_para]['sigbetter'] = para_df.iloc[paramID_sigbetter, para_df.columns.get_loc(i_para)].values.tolist()
                for i_LOD in LOD:
                    val[i_mis][i_para][i_LOD] = para_df.iloc[paramID_mis_losesig[i_mis][i_LOD], para_df.columns.get_loc(i_para)].values.tolist()

    color = ['r', 'b', 'g', 'c', 'y', 'k', 'm', 'blueviolet']
    title = ['Spop', 'R1pop', 'R2pop', 'R12pop', 'g0', 'Sa.S.D1', 'Sa.S.D2', 'Sa.R1.D1', 'Sa.R2.D2', 'T.StoR1', 'T.StoR2']
    x = list(result.keys())
    for i_mis in MISSPECIFICATION_LOD_STR:
        for i, i_para in enumerate(list(para.keys())[ind_par:]):
            print(i_para)
            if i_para == 'X':
                continue
            plt.rcParams['font.size'] = 21
            plt.figure()
            fig = plt.gcf()
            fig.set_size_inches(24, 11)
            i_val = val[i_mis][i_para]
            if i_para in ['Spop', 'R1pop', 'R2pop', 'R12pop', 'g0_S', 'T.R1..S.', 'T.R2..S.']:
                i_paraval_sorted = sorted(set(para_df[i_para]))
                if len(i_paraval_sorted) > len(color):
                    color = ["#" + ''.join([random.choice('0123456789ABCDEF') for _ in range(6)]) for i in range(len(i_paraval_sorted))]
                i_val = pd.DataFrame.from_dict(i_val)
                baseline = np.zeros(i_val.shape[1])

                for j in range(i_val.shape[0]):
                    i_row = i_val.iloc[j, :]
                    plt.plot(x, i_row, color=color[j], linewidth=2, marker='d', markersize=7)
                    # plt.bar(x, i_row.values, bottom=baseline, color=color[j])
                    # baseline += i_row.values

                plt.ylim([-10, 80])
                plt.yticks(list(range(0, 80, 10)))
                plt.ylabel('percentage')
                plt.legend(['{:0.2e}'.format(i) for i in i_paraval_sorted], loc='upper center', frameon=False, ncol=int(len(i_paraval_sorted)))
            else:
                df = pd.DataFrame(columns=['group', 'value'])
                for key, value in i_val.items():
                    i_df = pd.DataFrame({'group': np.repeat(key, len(value)), 'value': value})
                    df = pd.concat([df, i_df])

                df = pd.DataFrame(df.to_dict('records'))
                sns.violinplot(x='group', y='value', data=df, order=i_val.keys())
                # plt.xticks(labels=i_val.keys())
                # plt.yscale('log', base=10)
                # plt.ylim([1e-15, 1e3])

            plt.title(title[i] + '  ' + i_mis)
            plt.xlabel('Group')
            plt.show()
            plt.savefig(os.path.join(pathsave, 'par', title[i] + '  ' + i_mis + '.pdf'), dpi=FIG_DPI)
            plt.close('all')
    return


def DPM_analysis_misspec_plot(Strategy_name, Simduration, data, titlestr):
    plt.ioff()
    km = {key: data[key] for key in Strategy_name}
    hz, p = data['hz_ratio'], data['p']
    for i in hz.keys():
        idx_p = [idx for idx, i_val in enumerate(p.keys()) if set(i) == set(i_val)][0]
        km_ref = km[i[0]]
        km_treat = km[i[1]]
        par = {'color': ['k', 'b'], '0': i[0], '1': i[1], 'hzr': hz[i], 'p': p[list(p.keys())[idx_p]], 'duration': Simduration,
               'xtick step': 300, 'totalnum': data['num']}
        DPM_plot_KM(km_ref, km_treat, par, titlestr)
        if PLT_INTERACTIVE:
            plt.show()
        # plt.close('all')
    return

def DPM_analysis_hazard_ratio(stoptime, Strategy_name, Simduration):
    p = DPM_analysis_pairwise_logrank_test(stoptime, Strategy_name, Simduration)
    hz = dict()
    for i in list(itertools.combinations(Strategy_name, 2)):
        if 'strategy0' in i:
            ref = stoptime['strategy0']
            treat = stoptime[list(filter(lambda x: x != 'strategy0', i))[0]]
            name = ('strategy0', list(filter(lambda x: x != 'strategy0', i))[0])
        else:
            ref, treat = stoptime[i[0]], stoptime[i[1]]
            name = i
        hz[name] = DPM_analysis_HZ(ref, treat, Simduration)
    return hz, p


def DPM_analysis_KM(data, duration):
    E = [1 if i_val <= duration else 0 for i_val in data]
    epf = ExponentialFitter().fit(data, E)
    hazard = {'mean': epf.hazard_.mean(), 'ci': epf.confidence_interval_hazard_.mean()}

    kmf = KaplanMeierFitter()
    kmf.fit(data, E)
    km_interval = kmf.confidence_interval_survival_function_
    km = kmf.survival_function_

    t = km.index.values
    val = km.iloc[:].values.flatten()
    interval_lower = km_interval.iloc[:, 0].values.flatten()
    interval_upper = km_interval.iloc[:, 1].values.flatten()
    median_survival = kmf.median_survival_time_

    idx = np.where(t <= duration)[0]
    t, val, interval_lower, interval_upper = t[idx], val[idx], interval_lower[idx], interval_upper[idx]
    t, val, interval_lower, interval_upper = np.append(t, duration), \
        np.append(val, val[-1]), \
        np.append(interval_lower, interval_lower[-1]), \
        np.append(interval_upper, interval_upper[-1])

    return {'t': t, 'val': val, 'median_survival': median_survival, 'interval_lower': interval_lower, 'interval_upper': interval_upper,
            'hazard': hazard}


def DPM_analysis_HZ(data_ref, data, Simduration):
    if (len(data_ref) != 0) & (len(data) != 0):
        treat = np.concatenate((np.zeros(len(data_ref)), np.ones(len(data))))
        val = data_ref + data
        E = [1 if i_val <= Simduration else 0 for i_val in data_ref]
        E.extend([1 if i_val <= Simduration else 0 for i_val in data])
        E = np.array(E)
        d = {'val': val, 'E': E, 'treat': treat}
        df = pd.DataFrame(data=d)
        cph = CoxPHFitter()
        cph.fit(df, duration_col='val', event_col='E')
        hz_ratio = cph.hazard_ratios_.values[0]
    else:
        hz_ratio = None
    return hz_ratio


def DPM_analysis_pairwise_logrank_test(stoptime, Strategy_name, Simduration):
    G, T, E = tuple(), tuple(), tuple()
    flag_empty = False
    for i_strategy in Strategy_name:
        i_stop = stoptime[i_strategy]
        if not i_stop:
            flag_empty = True
            break
        G = G + tuple(len(i_stop) * [str(i_strategy)])
        E = E + tuple([1 if i_val <= Simduration else 0 for i_val in i_stop])
        T = T + tuple(i_stop)
    p = statistics.pairwise_logrank_test(T, G, E) if not flag_empty else None
    p_out = dict(zip(p.name, p.p_value)) if not flag_empty else None
    return p_out


def DPM_analysis_dose(dose, strategyname, inddrug=1):
    firstuse, max_num_change, num_change, average_move_num = [], [], [], []
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
            if strategyname == 'strategy0' and i_current not in ['(0.0,1.0)', '(0.0,0.0)', '(1.0,0.0)']:
                assert strategyname == 'strategy0' and i_current not in ['(0.0,1.0)', '(0.0,0.0)', '(1.0,0.0)']
            if i_firstuse is None and i_val[inddrug] != 0:
                i_firstuse = j+1
            if j_step != i_current and j_step != '(0.0,0.0)':
                i_num_change += 1
                i_current = j_step

        i_average_move_num = drugovermove/drugtotal if drugtotal != 0 else None
        i_firstuse = len(i_dose) + 1 if i_firstuse is None else i_firstuse
        i_average_move_num = len(i_dose) + 1 if i_average_move_num is None else i_average_move_num
        firstuse.append(i_firstuse)
        average_move_num.append(i_average_move_num)
        num_change.append(i_num_change/len(i_dose))
        max_num_change.append(i_num_change)

        # if strategyname == 'strategy0' and i_num_change > 1:
        #     print(i_num_change)
        #     print(i_dose)

        # firstuse = np.array(firstuse)
        # num_change = np.array(num_change)
        # a = firstuse[num_change==0]

    return dict(zip(['first drug 2', 'average move number', 'drug changes', 'max drug changes'],
                    [np.mean(firstuse), np.mean(average_move_num), np.mean(num_change), np.max(max_num_change)]))


def DPM_analysis_sigbetter(stoptime, Strategy):
    paramID, stoptime_ref, stoptime_test = stoptime['paramID'], stoptime[Strategy[0]], stoptime[Strategy[1]]
    ind = np.logical_and(np.array(stoptime_test) > np.array(stoptime_ref) + 30*2, np.array(stoptime_test) > 1.25 * np.array(stoptime_ref))
    return list(itertools.compress(paramID, ind))
