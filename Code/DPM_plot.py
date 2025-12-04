from DPM_lib import os, plt, mpl, math, PdfPages
from DPM_miscellaneous import DPM_miscellaneous_treatment_change_time
from DPM_constant import *
''' This script plots the simulation results of the DPM (Dynamic Precision Medicine) model. '''
'''
CHECKED
'''


if not PLT_INTERACTIVE:
    mpl.use('pdf')  # 'pdf'
else:
    # mpl.use('TkAgg')
    mpl.use("macosx")
print('Plot backend is ', mpl.get_backend())


def DPM_plot_all(Num_drug, strategy, Simduration, Limit_molecular_detection, Limit_radiologicdetection, Limit_mortality, pathsave):
    X_all = dict()
    for i, (i_key, i_strategy) in enumerate(strategy.items()):
        if i_strategy is not None:
            title = i_key
            DPM_plot_1strategy(Num_drug, i_strategy, Simduration, Limit_molecular_detection, Limit_radiologicdetection, Limit_mortality, title,
                               pathsave)
            t, X, _ = i_strategy
            i_strategy = (t, X.sum(axis=0))
            X_all[i_key] = i_strategy
        else:
            X_all[i_key] = None
    DPM_plot_allstrategy(X_all, Simduration, Limit_radiologicdetection, Limit_mortality, pathsave)
    return


def DPM_plot_1strategy(Num_drug, strategy, Simduration, Limit_moleculardetection, Limit_radiologicdetection, Limit_mortality, title, pathsave):
    # Plot simulation result for a single strategy.#
    color_X = ('#45936D', '#edc165', '#009ee9', '#b95559') if Num_drug == 2 else COLOR_X_3DRUG if Num_drug == 3 else None
    label_X = LABEL_X_2DRUG if Num_drug == 2 else LABEL_X_3DRUG if Num_drug == 3 else None
    color_drug = ('#D04F3A', '#227D88') if Num_drug == 2 else COLOR_DRUG_3DRUG if Num_drug == 3 else None
    label_drug = LABEL_2DRUG if Num_drug == 2 else LABEL_3DRUG if Num_drug == 3 else None
    legend_order = LEGEND_ORDER_2DRUG if Num_drug == 2 else LEGEND_ORDER_3DRUG if Num_drug == 3 else None
    legend_colnum = LEGEND_COLNUM_2DRUG if Num_drug == 2 else LEGEND_COLNUM_3DRUG if Num_drug == 3 else None

    t, X, d = strategy
    Xtotal = np.sum(X, axis=0)
    # X = np.vstack((Xtotal, X))
    diffpts = DPM_miscellaneous_treatment_change_time(d)

    # plt.rcParams['figure.dpi'] = FIG_DPI
    plt.rcParams['font.size'] = 14
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(10, 6)
    ax = plt.axes()
    for i_type in range(X.shape[0]):
        X1 = X[i_type, :]
        # X1[X1 < MINX] = MINX
        plt.plot(t/T_NORM, X1, color=color_X[i_type], label=label_X[i_type])

    ymaxval = 16  # math.ceil(math.log(np.max(X), LOGBASE)) + YINCREASE
    # xmaxval = t[-1]/T_NORM

    X_Limit_moleculardetection = Xtotal * Limit_moleculardetection
    plt.plot(t/T_NORM, X_Limit_moleculardetection, color='k', linestyle='dotted')
    if ymaxval > math.log(Limit_radiologicdetection, LOGBASE):
        plt.axhline(y=Limit_radiologicdetection, color='k', linestyle='dashdot')
        plt.text(np.max(t)/T_NORM/2, 2 * Limit_radiologicdetection, 'Limit of radiologic detection', fontsize=TEXTFONT)
    if ymaxval > math.log(Limit_mortality, LOGBASE):
        plt.axhline(y=Limit_mortality, color='k', linestyle='dashed')
        plt.text(np.max(t)/T_NORM/2, 2 * Limit_mortality, 'Limit of mortality', fontsize=TEXTFONT)

    diffpts = np.append(diffpts, (t[-1]))
    diffpts = diffpts.astype(int)
    for i_treat in range(len(diffpts)-1):
        d_i = d[:, diffpts[i_treat]:diffpts[i_treat+1]]
        d_i = np.unique(d_i, axis=1)
        d_i_sum = np.sum(d_i)
        y_basal = ymaxval - YINCREASE
        for i_drug in range(Num_drug):
            if d_i[i_drug] != 0:
                plt.fill_between(t[diffpts[i_treat]:diffpts[i_treat + 1]]/T_NORM, LOGBASE ** y_basal,
                                 LOGBASE ** (y_basal + YINCREASE * d_i[i_drug]/d_i_sum), color=color_drug[i_drug])
                plt.hlines(y=LOGBASE ** (y_basal + YINCREASE * d_i[i_drug] / d_i_sum), xmin=t[diffpts[i_treat]] / T_NORM,
                           xmax=t[diffpts[i_treat + 1] - 1]/T_NORM, color='black')
                # if i_treat != 0:
                plt.vlines(x=t[diffpts[i_treat]]/T_NORM, ymin=LOGBASE ** (ymaxval-YINCREASE), ymax=LOGBASE ** ymaxval, color='black')
                y_basal = y_basal + YINCREASE * d_i[i_drug]/d_i_sum

    # Plot the artificial label.#
    y_basal = ymaxval - YINCREASE
    d_lab = np.full(Num_drug, 1/Num_drug, dtype=float)
    t_lab = np.array(range(-int(t[-1]), -int(t[-1]/2), 1))
    for i in range(Num_drug):
        plt.fill_between(t_lab.ravel(), LOGBASE ** y_basal, LOGBASE ** (y_basal + YINCREASE * d_lab[i]), color=color_drug[i],
                         label=label_drug[i])
        y_basal = y_basal + YINCREASE * d_lab[i]
    # plt.plot(t_lab, t_lab * Limit_moleculardetection, color='k', linestyle='dotted', label='Limit of molecular detection')

    plt.xlabel('Time[days]')
    plt.ylabel('Number of cells')
    plt.title(title)
    handles, labels = plt.gca().get_legend_handles_labels()
    leg = plt.legend([handles[idx] for idx in legend_order], [labels[idx] for idx in legend_order],
                     bbox_to_anchor=LEGEND_BOXANCHOR, loc='upper center', ncol=legend_colnum)
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_linewidth(0)
    ax.set_yscale('log', base=LOGBASE)
    # plt.xlim([XMINVAL, xmaxval + XINCREASE])
    plt.xlim([XMINVAL, Simduration + XINCREASE])
    plt.ylim([YMINVAL, LOGBASE ** ymaxval])
    ylocmaj = mpl.ticker.LogLocator(base=LOGBASE, numticks=len(range(-2, int(ymaxval), 2)))
    ax.yaxis.set_major_locator(ylocmaj)
    ylocmin = mpl.ticker.LogLocator(base=LOGBASE, subs=np.arange(0, LOGBASE) * 1/LOGBASE, numticks=100)
    ax.yaxis.set_minor_locator(ylocmin)
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_yticks([10**-2, 10**0, 10**3, 10**6, 10**9, 10 ** 13])
    ax.set_xticks(np.arange(0, math.ceil(Simduration/STEPSIZE_DEFAULT_VAL)+1)*STEPSIZE_DEFAULT_VAL)
    ax.grid()
    plt.tight_layout()
    plt.show()
    return


def DPM_plot_allstrategy(X_all, Simduration, Limit_radiologicdetection, Limit_mortality, pathsave):
    # plt.rcParams['figure.dpi'] = FIG_DPI
    color = ['#D04F3A', '#D04F3A', '#227D88', '#227D88']
    linestyle = ['solid', 'dashed', 'solid', 'dashed']
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(9, 5)
    ax = plt.axes()
    Xmax = 0
    tmax = 0
    for i, (i_key, i_strategy) in enumerate(X_all.items()):
        if i_strategy is not None:
            t, i_X_all = i_strategy
            i_X_all[i_X_all < MINX] = MINX
            Xmax = max(Xmax, max(i_X_all))
            tmax = max(tmax, max(t))
            plt.plot(t/T_NORM, i_X_all, color=color[i], linestyle=linestyle[i], label=i_key)

    ymaxval = 16  # math.ceil(math.log(Xmax, LOGBASE)) + YINCREASE/2
    # xmaxval = tmax/T_NORM
    if ymaxval > math.log(Limit_radiologicdetection, LOGBASE):
        plt.axhline(y=Limit_radiologicdetection, color='k', linestyle='dashdot')
        plt.text(tmax/T_NORM/2, 2 * Limit_radiologicdetection, 'Limit of radiologic detection', fontsize=TEXTFONT+1)
    if ymaxval > math.log(Limit_mortality, LOGBASE):
        plt.axhline(y=Limit_mortality, color='k', linestyle='dashed')
        plt.text(tmax/T_NORM/2, 2 * Limit_mortality, 'Limit of mortality', fontsize=TEXTFONT+1)

    plt.xlabel('Time[Months]')
    plt.ylabel('Number of cells')
    plt.title('Total cell number')
    handles, labels = plt.gca().get_legend_handles_labels()
    leg = plt.legend(handles, labels, bbox_to_anchor=LEGEND_BOXANCHOR, loc='upper center', ncol=LEGEND_COLNUM_XALL-1)
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_linewidth(0)
    ax.set_yscale('log', base=LOGBASE)
    plt.xlim([XMINVAL, Simduration + XINCREASE])
    plt.ylim([LOGBASE ** YMINVAL_XALL, LOGBASE ** ymaxval])
    ylocmaj = mpl.ticker.LogLocator(base=LOGBASE, numticks=len(range(-2, int(ymaxval), 2)))
    ax.yaxis.set_major_locator(ylocmaj)
    ylocmin = mpl.ticker.LogLocator(base=LOGBASE, subs=np.arange(0, LOGBASE) * 1 / LOGBASE, numticks=100)
    ax.yaxis.set_minor_locator(ylocmin)
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xticks(np.arange(0, math.ceil(Simduration/STEPSIZE_DEFAULT_VAL)+1)*STEPSIZE_DEFAULT_VAL)
    ax.set_yticks([10**-2, 10**0, 10**3, 10**6, 10**9, 10 ** 13])
    plt.show()
    plt.tight_layout()
    ax.grid()
    # plt.savefig(os.path.join(pathsave, 'totalcellnum' + '.pdf'), dpi=FIG_DPI)
    # plt.close('all')
    return


def DPM_plot_KM(data0, data, par, titlestr):
    if 'p' in par.keys():
        p_str = DPM_plot_pstr(par['p'])
    else:
        p_str = ''

    plt.rcParams['font.size'] = 20
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(12.25, 9)
    ax = plt.axes()
    l1, = plt.plot(data0['t'], data0['val'], color=par['color'][0])
    ax.fill_between(data0['t'], data0['interval_lower'], data0['interval_upper'],
                    color=par['color'][0], alpha=.3)
    l2, = plt.plot(data['t'], data['val'], color=par['color'][1])
    ax.fill_between(data['t'], data['interval_lower'], data['interval_upper'],
                    color=par['color'][1], alpha=.3)
    l3, = plt.plot([], [], ' ')
    l4, = plt.plot([], [], ' ')
    l5, = plt.plot([], [], ' ')
    if p_str:
        l6, = plt.plot([], [], ' ')
    else:
        l6 = None
    # plt.title(f"{int(par['totalnum'])} of patients ({titlestr})")
    plt.title(titlestr)
    plt.xticks(list(range(0, par['duration']+par['xtick step'], par['xtick step'])))
    plt.yticks(np.round(np.arange(0, 1.2, 0.2), 1).tolist())
    plt.xlabel('Survival time[days]')
    plt.ylabel('Fraction of surviving patients')
    if p_str:
        plt.legend([l1, l2, l3, l4, l5, l6],
                   [par['0'], par['1'],
                    f"median survival time of {par['0']} (days): {int(data0['median_survival'])}",
                    f"median survival time of {par['1']} (days): {int(data['median_survival'])}",
                    f"hazard ratio {par['1']} over {par['0']}: {round(par['hzr'], 3)}",
                    p_str],
                   loc='upper right', frameon=False)
    else:
        plt.legend([l1, l2, l3, l4, l5],
                   [par['0'], par['1'],
                    f"median survival time of {par['0']} (days): {int(data0['median_survival'])}",
                    f"median survival time of {par['1']} (days): {int(data['median_survival'])}",
                    f"hazard ratio {par['1']} over {par['0']}: {round(par['hzr'], 3)}"],
                   loc='upper right', frameon=False)
    return


def DPM_plot_KM_multi(km, par, titlestr):
    plt.rcParams['font.size'] = 14
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(12.25, 9)
    ax = plt.axes()
    keys = km.keys()
    legh = []
    legstr = []
    color = par['color']
    linestyle = par['linestyle']
    for i, i_key in enumerate(keys):
        i_km = km[i_key]
        l_1i, = plt.plot(i_km['t'], i_km['val'], color=color[i], linestyle=linestyle[i])
        ax.fill_between(i_km['t'], i_km['interval_lower'], i_km['interval_upper'], color=color[i], alpha=.3)
        l_2i, = plt.plot([], [], ' ')
        legh.extend([l_1i, l_2i])
        legstr.extend([i_key, f"median survival time of {i_key} (days): {int(i_km['median_survival'])}"])
        # l_i, = plt.plot([], [], ' ')
        # legh.append(l_i)
        # legstr.extend([f"# of patients {par['totalnum']:,}, ({titlestr})"])

    # plt.title(f"{int(par['totalnum'])} of patients ({titlestr})")
    plt.title(titlestr)
    plt.xticks(list(range(0, par['duration']+par['xtick step'], par['xtick step'])))
    plt.yticks(np.round(np.arange(0, 1.2, 0.2), 1).tolist())
    plt.xlabel('Survival time[days]')
    plt.ylabel('Fraction of surviving patients')
    plt.legend(legh, legstr, loc='upper right', frameon=False)
    return


def DPM_plot_multi(x, data, leg, name, color, linestyle):
    plt.rcParams['font.size'] = 17
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(13, 7)
    legh, legstr = [], []
    xx = np.arange(len(x))
    for i in range(len(leg)):

        j_l, = plt.plot(xx.tolist(), data[i, :], color=color[i], linestyle=linestyle[i], linewidth=2, marker='d', markersize=7)
        legh.append(j_l)

    plt.ylabel(name)
    plt.xticks(xx.tolist(), x)
    plt.legend(legh, leg, bbox_to_anchor=(0.5, 1.17), loc='upper center', frameon=False, ncol=4)
    return


def DPM_plot_pstr(p):
    if p > 0.05:
        p_str = 'n.s.'
    elif 0.01 < p <= 0.05:
        p_str = '*p ≤ 0.05'
    elif 0.001 < p <= 0.01:
        p_str = '**p ≤ 0.01'
    elif 0.0001 < p <= 0.01:
        p_str = '***p ≤ 0.001'
    else:
        p_str = '****p ≤ 0.0001'
    return p_str


def DPM_plot_df2pdf(df, filename, color_, numpages=(1, 1), pagesize=(11, 8.5)):
    def DPM_plot_df2pdf_1(df_):
        # draw as table
        alternating_colors = [color_[0] * len(df_.columns)]
        alternating_colors = alternating_colors * len(df_)
        # plt.rc('font', size=30)
        fig_, ax = plt.subplots(figsize=pagesize)
        ax.axis('tight')
        ax.axis('off')
        the_table = ax.table(cellText=df_.values,
                             rowLabels=None,
                             colLabels=df_.columns,
                             rowColours=color_[1]*len(df_),
                             colColours=color_[1]*len(df_.columns),
                             cellColours=alternating_colors,
                             loc='center')
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(8)
        the_table.scale(1, 4)
        return fig_

    with PdfPages(filename) as pdf:
        nh, nv = numpages
        rows_per_page = len(df)//nh
        cols_per_page = len(df.columns)//nv
        for i in range(0, nh):
            for j in range(0, nv):
                page = df.iloc[(i*rows_per_page):min((i+1)*rows_per_page, len(df)),
                               (j*cols_per_page):min((j+1)*cols_per_page, len(df.columns))]
                fig = DPM_plot_df2pdf_1(page)
                if nh > 1 or nv > 1:
                    # Add part/page number at bottom-center of page
                    fig.text(0.5, 0.5/pagesize[0],
                             'Part-{}x{}: Page-{}'.format(i+1, j+1, i*nv + j + 1), ha='center', fontsize=9)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close()
    return
