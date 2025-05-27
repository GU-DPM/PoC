import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import matplotlib as mpl


def DPM_miscellaneous_treatment_change_time(d):
    d_diff = np.diff(d)
    d_diff_boolen = d_diff.astype(bool)
    if d_diff_boolen.ndim != 1:
        pos_changedose = np.sum(d_diff_boolen, axis=0).astype(bool)
    else:
        pos_changedose = d_diff_boolen
    pos_changedose, = np.where(pos_changedose)
    # Position where drug dose changes.
    pos_changedose += 1
    # At the beginning, index=0, drug dose always changes. It can be seen as the start of applying treatment, previous drug doses are all 0.
    pos_changedose = np.insert(pos_changedose, 0, 0)
    return pos_changedose


def DPM_miscellaneous_cutbydose(drug, Spop, R1pop, R2pop, R12pop):
    index, = np.where(drug < 0)
    if index.any():
        index = index[0]
        drug = drug[:index]
        Spop = Spop[:index+1]
        R1pop = R1pop[:index+1]
        R2pop = R2pop[:index+1]
        R12pop = R12pop[:index+1]
    return drug, Spop, R1pop, R2pop, R12pop


def DPM_miscellaneous_cutbycell(drug, Spop, R1pop, R2pop, R12pop, Limit_mortality):
    X = np.vstack((Spop, R1pop, R2pop, R12pop))
    total = X.sum(axis=0)
    maxcell = np.amax(X, axis=0)
    index_total, = np.where(total > Limit_mortality)
    index_maxcell, = np.where(maxcell < 1)
    index_total = index_total[0] if index_total.any() else index_total
    index_maxcell = index_maxcell[0] if index_maxcell.any() else index_maxcell
    index = []
    if index_total.any() and index_maxcell.any():
        index = min([index_total, index_maxcell])
    elif index_total.any():
        index = index_total
    elif index_maxcell.any():
        index = index_maxcell
    if index.any():
        drug = drug[:index-1]
        Spop = Spop[:index]
        R1pop = R1pop[:index]
        R2pop = R2pop[:index]
        R12pop = R12pop[:index]
    return drug, Spop, R1pop, R2pop, R12pop


def DPM_miscellaneous_dataplot(drug1, drug2, Spop, R1pop, R2pop, R12pop, timestep):
    t = np.arange(0, Spop.shape[0]*timestep, timestep)
    X = np.vstack((Spop, R1pop, R2pop, R12pop))
    d = np.vstack((drug1, drug2))
    data = (t, X, d)
    return data


def DPM_plot_1strategy(data, title='', savename='ex', Num_drug=2, Simduration=1800, Limit_moleculardetection=1/1e4,
                       Limit_radiologicdetection=1e9, Limit_mortality=1e13):
    mpl.use('pdf')
    # Plot simulation result of one strategy.
    LOGBASE = 10
    MINX = .1
    XMINVAL = -1
    XINCREASE = 45
    YINCREASE = 2
    YMINVAL = 1e-2
    TEXTFONT = 10
    LEGEND_BOXANCHOR = (0.5, 1.15)
    color_X = ('b', 'g', 'c', 'm', 'r')
    label_X = ('Total cells', 'S cells', 'R1 cells', 'R2 cells', 'R12 cells')
    color_drug = ('g', 'b')
    label_drug = ('Drug 1', 'Drug 2')
    legend_colnum = 4
    legend_order = np.reshape(np.array(range(0, 8)), (-1, legend_colnum)).flatten('F')
    legend_order = legend_order[:-1]

    t, X, d = data
    Xtotal = 1.5 * np.sum(X, axis=0)
    X = np.vstack((Xtotal, X))
    diffpts = DPM_miscellaneous_treatment_change_time(d)

    plt.rcParams['font.size'] = 14
    # plt.figure()
    # ax = plt.axes()
    fig, ax = plt.subplots()
    # fig.set_size_inches(fig_width-1, fig_height+1)
    fig.set_size_inches(14, 9)

    ax.set_yscale('log', base=LOGBASE)
    for i_type in range(X.shape[0]):
        X1 = X[i_type, :]
        X1[X1 < MINX] = MINX
        plt.plot(t, X1, color=color_X[i_type], label=label_X[i_type])

    ymaxval = math.ceil(math.log(np.max(X), LOGBASE)) + YINCREASE
    # xmaxval = t[-1]/T_NORM

    X_Limit_moleculardetection = X[0, ] * Limit_moleculardetection
    plt.plot(t, X_Limit_moleculardetection, color='k', linestyle='dotted')
    if ymaxval > math.log(Limit_radiologicdetection, LOGBASE):
        plt.axhline(y=Limit_radiologicdetection, color='k', linestyle='dashdot')
        plt.text(np.max(t)/2, 2 * Limit_radiologicdetection, 'Limit of radiologic detection', fontsize=TEXTFONT)
    if ymaxval > math.log(Limit_mortality, LOGBASE):
        plt.axhline(y=Limit_mortality, color='k', linestyle='dashed')
        plt.text(np.max(t)/2, 2 * Limit_mortality, 'Limit of mortality', fontsize=TEXTFONT)

    diffpts = np.append(diffpts, (len(t)-1))
    diffpts = diffpts.astype(int)
    for i_treat in range(len(diffpts)-1):
        d_i = d[:, diffpts[i_treat]:diffpts[i_treat+1]]
        d_i = np.unique(d_i, axis=1)
        d_i_sum = d_i.sum()
        y_basal = ymaxval - YINCREASE
        for i_drug in range(Num_drug):
            if d_i.shape[1] != 0 and d_i[i_drug] != 0:
                i_begin = diffpts[i_treat]  # if i_treat == 0 else diffpts[i_treat]-1
                plt.fill_between(t[i_begin:diffpts[i_treat + 1]+1], LOGBASE ** y_basal,
                                 LOGBASE ** (y_basal + YINCREASE * d_i[i_drug]/d_i_sum), color=color_drug[i_drug])
                plt.hlines(y=LOGBASE ** (y_basal + YINCREASE * d_i[i_drug] / d_i_sum), xmin=t[i_begin],
                           xmax=t[diffpts[i_treat + 1]], color='black')
                if i_treat != 0:
                    plt.vlines(x=t[i_begin], ymin=LOGBASE ** (ymaxval-YINCREASE), ymax=LOGBASE ** ymaxval, color='black')
                y_basal = y_basal + YINCREASE * d_i[i_drug]/d_i_sum

    # Plot artificial label.
    y_basal = ymaxval - YINCREASE
    d_lab = np.full(Num_drug, 1/Num_drug, dtype=float)
    t_lab = np.array(range(-int(t[-1]), -int(t[-1]/2), 1))
    for i in range(Num_drug):
        plt.fill_between(t_lab, LOGBASE ** y_basal, LOGBASE ** (y_basal + YINCREASE * d_lab[i]), color=color_drug[i], label=label_drug[i])
        y_basal = y_basal + YINCREASE * d_lab[i]
    plt.plot(t_lab, t_lab * Limit_moleculardetection, color='k', linestyle='dotted')
    plt.xlabel('Time[Months]')
    plt.ylabel('Number of cells')
    plt.title(title)
    handles, labels = plt.gca().get_legend_handles_labels()
    leg = plt.legend([handles[idx] for idx in legend_order], [labels[idx] for idx in legend_order],
                     bbox_to_anchor=LEGEND_BOXANCHOR, loc='upper center', ncol=legend_colnum)
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_linewidth(0)
    # plt.xlim([XMINVAL, xmaxval + XINCREASE])
    plt.xlim([XMINVAL, Simduration + XINCREASE])
    plt.ylim([YMINVAL, LOGBASE ** ymaxval])
    ylocmaj = mpl.ticker.LogLocator(base=LOGBASE, numticks=len(range(-2, int(ymaxval), 2)))
    ax.yaxis.set_major_locator(ylocmaj)
    ylocmin = mpl.ticker.LogLocator(base=LOGBASE, subs=np.arange(0, LOGBASE) * 1/LOGBASE, numticks=100)
    ax.yaxis.set_minor_locator(ylocmin)
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    # ax.set_xticks(t)
    plt.show()
    plt.savefig(savename+'.pdf', dpi=300)
    return


def DPM_plot_1csv(filename='./singleResults.csv'):
    savename = 'ex'
    df = pd.read_csv(filename)
    colnum_drug1 = -1
    colnum_Spop = 1
    colnum_R1pop = 2
    colnum_R2pop = 3
    colnum_R12pop = 4
    timestep = 45
    drug1 = df.iloc[:, colnum_drug1].to_numpy()
    Spop = df.iloc[:, colnum_Spop].to_numpy()
    R1pop = df.iloc[:, colnum_R1pop].to_numpy()
    R2pop = df.iloc[:, colnum_R2pop].to_numpy()
    R12pop = df.iloc[:, colnum_R12pop].to_numpy()

    drug1, Spop, R1pop, R2pop, R12pop = DPM_miscellaneous_cutbydose(drug1, Spop, R1pop, R2pop, R12pop)

    drug2 = 1 - drug1
    t = np.arange(0, Spop.shape[0]*timestep, timestep)
    X = np.vstack((Spop, R1pop, R2pop, R12pop))
    d = np.vstack((drug1, drug2))
    data = (t, X, d)
    title = ''
    DPM_plot_1strategy(data)


def DPM_plot_csv(filename='./individualResults_31731026.csv'):
    savename = 'ex'
    df = pd.read_csv(filename)
    timestep = 45
    colnum_drug1_cpm = 13
    colnum_drug1_dpm = 14
    colnum_Spop_cpm, colnum_R1pop_cpm, colnum_R2pop_cpm, colnum_R12pop_cpm = 1, 2, 3, 4
    colnum_Spop_dpm, colnum_R1pop_dpm, colnum_R2pop_dpm, colnum_R12pop_dpm = 5, 6, 7, 8
    colnum_Spop_dpmtrial, colnum_R1pop_dpmtrial, colnum_R2pop_dpmtrial, colnum_R12pop_dpmtrial = 9, 10, 11, 12

    drug1_cpm = df.iloc[:, colnum_drug1_cpm].to_numpy()
    drug1_dpm = df.iloc[:, colnum_drug1_dpm].to_numpy()

    Spop_cpm = df.iloc[:, colnum_Spop_cpm].to_numpy()
    R1pop_cpm = df.iloc[:, colnum_R1pop_cpm].to_numpy()
    R2pop_cpm = df.iloc[:, colnum_R2pop_cpm].to_numpy()
    R12pop_cpm = df.iloc[:, colnum_R12pop_cpm].to_numpy()

    Spop_dpm = df.iloc[:, colnum_Spop_dpm].to_numpy()
    R1pop_dpm = df.iloc[:, colnum_R1pop_dpm].to_numpy()
    R2pop_dpm = df.iloc[:, colnum_R2pop_dpm].to_numpy()
    R12pop_dpm = df.iloc[:, colnum_R12pop_dpm].to_numpy()

    Spop_dpmtrial, R1pop_dpmtrial, R2pop_dpmtrial, R12pop_dpmtrial = \
        (df.iloc[:, colnum_Spop_dpmtrial].to_numpy(),
         df.iloc[:, colnum_R1pop_dpmtrial].to_numpy(),
         df.iloc[:, colnum_R2pop_dpmtrial].to_numpy(),
         df.iloc[:, colnum_R12pop_dpmtrial].to_numpy())

    drug1_dpmtrial = np.concatenate((drug1_dpm[:2], drug1_cpm[2:]))
    drug1_cpm, Spop_cpm, R1pop_cpm, R2pop_cpm, R12pop_cpm = \
        DPM_miscellaneous_cutbydose(drug1_cpm, Spop_cpm, R1pop_cpm, R2pop_cpm, R12pop_cpm)
    drug1_dpm, Spop_dpm, R1pop_dpm, R2pop_dpm, R12pop_dpm = \
        DPM_miscellaneous_cutbydose(drug1_dpm, Spop_dpm, R1pop_dpm, R2pop_dpm, R12pop_dpm)

    Limit_mortality = 1e13
    drug1_dpmtrial, Spop_dpmtrial, R1pop_dpmtrial, R2pop_dpmtrial, R12pop_dpmtrial = \
        DPM_miscellaneous_cutbycell(drug1_dpmtrial, Spop_dpmtrial, R1pop_dpmtrial, R2pop_dpmtrial, R12pop_dpmtrial, Limit_mortality)

    drug2_cpm, drug2_dpm, drug2_dpmtrial = 1 - drug1_cpm, 1 - drug1_dpm, 1- drug1_dpmtrial
    data_cpm = DPM_miscellaneous_dataplot(drug1_cpm, drug2_cpm, Spop_cpm, R1pop_cpm, R2pop_cpm, R12pop_cpm, timestep)
    data_dpm = DPM_miscellaneous_dataplot(drug1_dpm, drug2_dpm, Spop_dpm, R1pop_dpm, R2pop_dpm, R12pop_dpm, timestep)
    data_dpmtrial = DPM_miscellaneous_dataplot(drug1_dpmtrial, drug2_dpmtrial, Spop_dpmtrial, R1pop_dpmtrial, R2pop_dpmtrial, R12pop_dpmtrial,
                                               timestep)

    DPM_plot_1strategy(data_cpm, title='cpm', savename='cpm')
    DPM_plot_1strategy(data_dpm, title='dpm', savename='dpm')
    DPM_plot_1strategy(data_dpmtrial, title='dpmtrail', savename='dpmtrail')
    return


DPM_plot_csv(filename='./individualResults_31731026.csv')
DPM_plot_1csv(filename='./singleResults.csv')

