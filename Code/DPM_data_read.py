from DPM_lib import plt, sm, pd, np, sns, lognorm, tabulate, stat, math
from DPM_read_save import DPM_save_fcs, DPM_read_fcspckl
from DPM_plot import DPM_plot_linearreg
from DPM_data_analysis import DPM_data_analysis_cellcount

flag_read_1 = True
filename = ['../Data/Cell count/Effect of Combination Rx on MM361 par transduced w GFPmut_Exp1.xlsx',
            '../Data/Cell count/Effect of Combination Rx on MM361 par transduced w GFPmut_Exp2.xlsx',
            ]
sheetname = ['EXP1', 'EXP2']
numrep = len(filename)
par = [(-1e3, 1e3, 0.93, 0.80),  # (0)
       (-1e3, 1e3, 0.94, 0.70),  # (1)
       ]
index = [0, 3, 6, 9]
condition = ['Veh',
             'TP(1ug/ml)+ICI(100nM)',
             'Palbo(100nM)',
             'TP(1ug)+ICI(100nM)+Palbo(100nM)',
             'TP(1ug)+ICI(100nM)+Palbo(50nM)']
colname = [numrep * [i_col] for i_col in condition]
colnamepd = [item for sublist in colname for item in sublist]
MM361_GFPmut = pd.DataFrame(data=None, index=index, columns=colnamepd, dtype=float)
##
if flag_read_1:
    ylim = [10, 15]
    data_all_rep1 = pd.read_excel(filename[0], sheet_name=sheetname[0], engine='openpyxl')

    num_t0_rep1 = data_all_rep1.iloc[2, 3:6].to_numpy(dtype=int)
    num_t0_rep1 = np.append(num_t0_rep1, data_all_rep1.iloc[10, 3:6].to_numpy(dtype=int))
    num_t0_rep1 = np.append(num_t0_rep1, data_all_rep1.iloc[18, 3:7].to_numpy(dtype=int))
    num_t0_rep1 = np.append(num_t0_rep1, data_all_rep1.iloc[26, 3:6].to_numpy(dtype=int))
    num_t0_rep1 = np.append(num_t0_rep1, data_all_rep1.iloc[34, 3:6].to_numpy(dtype=int))

    Veh_rep1 = data_all_rep1.iloc[3:6, 3:7].to_numpy(dtype=int)
    Veh_rep1_3d, Veh_rep1_6d, Veh_rep1_9d = Veh_rep1[0, :],  Veh_rep1[1, :],  Veh_rep1[2, :]
    TP_ICI_rep1 = data_all_rep1.iloc[11:14, 3:7].to_numpy(dtype=int)
    TP_ICI_rep1_3d, TP_ICI_rep1_6d, TP_ICI_rep1_9d = TP_ICI_rep1[0, :], TP_ICI_rep1[1, :], TP_ICI_rep1[2, :]
    Palbo_rep1 = data_all_rep1.iloc[19:22, 3:7].to_numpy(dtype=int)
    Palbo_rep1_3d, Palbo_rep1_6d, Palbo_rep1_9d = Palbo_rep1[0, :], Palbo_rep1[1, :], Palbo_rep1[2, :]
    TP_ICI_Palbofull_rep1 = data_all_rep1.iloc[27:30, 3:7].to_numpy(dtype=int)
    TP_ICI_Palbofull_rep1_3d, TP_ICI_Palbofull_rep1_6d, TP_ICI_Palbofull_rep1_9d = \
        TP_ICI_Palbofull_rep1[0, :], TP_ICI_Palbofull_rep1[1, :], TP_ICI_Palbofull_rep1[2, :]
    TP_ICI_Palbohalf_rep1 = data_all_rep1.iloc[35:, 3:7].to_numpy(dtype=int)
    TP_ICI_Palbohalf_rep1_3d, TP_ICI_Palbohalf_rep1_6d, TP_ICI_Palbohalf_rep1_9d = \
        TP_ICI_Palbohalf_rep1[0, :], TP_ICI_Palbohalf_rep1[1, :], TP_ICI_Palbohalf_rep1[2, :]

    Veh_rep1 = [num_t0_rep1, Veh_rep1_3d, Veh_rep1_6d, Veh_rep1_9d]
    DPM_data_analysis_cellcount(index, Veh_rep1, condition[0] + f' rep{1}', ylim)
    plt.close('all')

    TP_ICI_rep1 = [num_t0_rep1, TP_ICI_rep1_3d, TP_ICI_rep1_6d, TP_ICI_rep1_9d]
    DPM_data_analysis_cellcount(index, TP_ICI_rep1, condition[1] + f' rep{1}', ylim)
    plt.close('all')

    Palbo_rep1 = [num_t0_rep1, Palbo_rep1_3d, Palbo_rep1_6d, Palbo_rep1_9d]
    DPM_data_analysis_cellcount(index, Palbo_rep1, condition[2] + f' rep{1}', ylim)
    plt.close('all')

    TP_ICI_Palbofull_rep1 = [num_t0_rep1, TP_ICI_Palbofull_rep1_3d, TP_ICI_Palbofull_rep1_6d, TP_ICI_Palbofull_rep1_9d]
    DPM_data_analysis_cellcount(index, TP_ICI_Palbofull_rep1, condition[3] + f' rep{1}', ylim)
    plt.close('all')

    TP_ICI_Palbohalf_rep1 = [num_t0_rep1, TP_ICI_Palbohalf_rep1_3d, TP_ICI_Palbohalf_rep1_6d, TP_ICI_Palbohalf_rep1_9d]
    DPM_data_analysis_cellcount(index, TP_ICI_Palbohalf_rep1, condition[4] + f' rep{1}', ylim)
    plt.close('all')

    del Veh_rep1, TP_ICI_rep1, Palbo_rep1, TP_ICI_Palbofull_rep1, TP_ICI_Palbohalf_rep1
    ## rep2
    data_all_rep2 = pd.read_excel(filename[1], sheet_name=sheetname[1], engine='openpyxl')

    num_t0_rep2 = data_all_rep2.iloc[2, 3:6].to_numpy(dtype=int)
    num_t0_rep2 = np.append(num_t0_rep2, data_all_rep2.iloc[10, 3:6].to_numpy(dtype=int))
    num_t0_rep2 = np.append(num_t0_rep2, data_all_rep2.iloc[18, 3:6].to_numpy(dtype=int))
    num_t0_rep2 = np.append(num_t0_rep2, data_all_rep2.iloc[26, 3:6].to_numpy(dtype=int))
    num_t0_rep2 = np.append(num_t0_rep2, data_all_rep2.iloc[34, 3:6].to_numpy(dtype=int))

    Veh_rep2 = data_all_rep2.iloc[3:6, 3:7].to_numpy(dtype=int)
    Veh_rep2_3d, Veh_rep2_6d, Veh_rep2_9d = Veh_rep2[0, :],  Veh_rep2[1, :],  Veh_rep2[2, :]
    TP_ICI_rep2 = data_all_rep2.iloc[11:14, 3:7].to_numpy(dtype=int)
    TP_ICI_rep2_3d, TP_ICI_rep2_6d, TP_ICI_rep2_9d = TP_ICI_rep2[0, :], TP_ICI_rep2[1, :], TP_ICI_rep2[2, :]
    Palbo_rep2 = data_all_rep2.iloc[19:22, 3:7].to_numpy(dtype=int)
    Palbo_rep2_3d, Palbo_rep2_6d, Palbo_rep2_9d = Palbo_rep2[0, :], Palbo_rep2[1, :], Palbo_rep2[2, :]
    TP_ICI_Palbofull_rep2 = data_all_rep2.iloc[27:30, 3:7].to_numpy(dtype=int)
    TP_ICI_Palbofull_rep2_3d, TP_ICI_Palbofull_rep2_6d, TP_ICI_Palbofull_rep2_9d = \
        TP_ICI_Palbofull_rep2[0, :], TP_ICI_Palbofull_rep2[1, :], TP_ICI_Palbofull_rep2[2, :]
    TP_ICI_Palbohalf_rep2 = data_all_rep2.iloc[35:, 3:7].to_numpy(dtype=int)
    TP_ICI_Palbohalf_rep2_3d, TP_ICI_Palbohalf_rep2_6d, TP_ICI_Palbohalf_rep2_9d = \
        TP_ICI_Palbohalf_rep2[0, :], TP_ICI_Palbohalf_rep2[1, :], TP_ICI_Palbohalf_rep2[2, :]

    Veh_rep2 = [num_t0_rep2, Veh_rep2_3d, Veh_rep2_6d, Veh_rep2_9d]
    DPM_data_analysis_cellcount(index, Veh_rep2, condition[0] + f' rep{2}', ylim)
    plt.close('all')

    TP_ICI_rep2 = [num_t0_rep2, TP_ICI_rep2_3d, TP_ICI_rep2_6d, TP_ICI_rep2_9d]
    DPM_data_analysis_cellcount(index, TP_ICI_rep2, condition[1] + f' rep{2}', ylim)
    plt.close('all')

    Palbo_rep2 = [num_t0_rep2, Palbo_rep2_3d, Palbo_rep2_6d, Palbo_rep2_9d]
    DPM_data_analysis_cellcount(index, Palbo_rep2, condition[2] + f' rep{2}', ylim)
    plt.close('all')

    TP_ICI_Palbofull_rep2 = [num_t0_rep2, TP_ICI_Palbofull_rep2_3d, TP_ICI_Palbofull_rep2_6d, TP_ICI_Palbofull_rep2_9d]
    DPM_data_analysis_cellcount(index, TP_ICI_Palbofull_rep2, condition[3] + f' rep{2}', ylim)
    plt.close('all')

    TP_ICI_Palbohalf_rep2 = [num_t0_rep2, TP_ICI_Palbohalf_rep2_3d, TP_ICI_Palbohalf_rep2_6d, TP_ICI_Palbohalf_rep2_9d]
    DPM_data_analysis_cellcount(index, TP_ICI_Palbohalf_rep2, condition[4] + f' rep{2}', ylim)
    plt.close('all')

    # print(tabulate(MM361_GFPmut, headers='keys', tablefmt='psql'))

    Veh = [np.append(num_t0_rep1, num_t0_rep2), np.append(Veh_rep1_3d, Veh_rep2_3d),
           np.append(Veh_rep1_6d, Veh_rep2_6d), np.append(Veh_rep1_9d, Veh_rep2_9d)]
    DPM_data_analysis_cellcount(index, Veh, condition[0], ylim)
    plt.close('all')

    TP_ICI = [np.append(num_t0_rep1, num_t0_rep2), np.append(TP_ICI_rep1_3d, TP_ICI_rep2_3d),
              np.append(TP_ICI_rep1_6d, TP_ICI_rep2_6d), np.append(TP_ICI_rep1_9d, TP_ICI_rep2_9d)]
    DPM_data_analysis_cellcount(index, TP_ICI, condition[1], ylim)
    plt.close('all')

    Palbo = [np.append(num_t0_rep1, num_t0_rep2), np.append(Palbo_rep1_3d, Palbo_rep2_3d), np.append(Palbo_rep1_6d, Palbo_rep2_6d),
             np.append(Palbo_rep1_9d, Palbo_rep2_9d)]
    DPM_data_analysis_cellcount(index, Palbo, condition[2], ylim)
    plt.close('all')

    TP_ICI_Palbofull = [np.append(num_t0_rep1, num_t0_rep2), np.append(TP_ICI_Palbofull_rep1_3d, TP_ICI_Palbofull_rep2_3d),
                        np.append(TP_ICI_Palbofull_rep1_6d, TP_ICI_Palbofull_rep2_6d), np.append(TP_ICI_Palbofull_rep1_9d, TP_ICI_Palbofull_rep2_9d)]
    DPM_data_analysis_cellcount(index, TP_ICI_Palbofull, condition[3], ylim)
    plt.close('all')

    TP_ICI_Palbohalf = [np.append(num_t0_rep1, num_t0_rep2), np.append(TP_ICI_Palbohalf_rep1_3d, TP_ICI_Palbohalf_rep2_3d),
                        np.append(TP_ICI_Palbohalf_rep1_6d, TP_ICI_Palbohalf_rep2_6d),
                        np.append(TP_ICI_Palbohalf_rep1_9d, TP_ICI_Palbohalf_rep2_9d)]
    DPM_data_analysis_cellcount(index, TP_ICI_Palbohalf, condition[4], ylim)
    plt.close('all')

flag_read_2 = True
filename = ['../Data/Cell count/Combination Rx on BT474 par.xlsx']
sheetname = ['P102']
numrep = len(filename)
time = [0, 3, 6, 9]
condition = ['BT474 Veh',
             'BT474 TP(1ug/ml)',
             'BT474 ICI(100nM)',
             'BT474 Palbo(100nM)',
             'BT474 TP(1ug/ml)+ICI(100nM)',
             'BT474 TP(1ug)+ICI(100nM)+Palbo(100nM)',
             'BT474 TP(1ug)+ICI(100nM)+Palbo(50nM)']
colname = [numrep * [i_col] for i_col in condition]
colnamepd = [item for sublist in colname for item in sublist]
BT474 = pd.DataFrame(data=None, index=index, columns=colnamepd, dtype=float)
##
if flag_read_2:
    ylim = [9, 15]
    for i, i_filename in enumerate(filename):
        data_all_i = pd.read_excel(i_filename, sheet_name=sheetname[i], engine='openpyxl')
        num_t0 = data_all_i.iloc[1, 3:7].to_numpy(dtype=int)
        num_t0 = np.append(num_t0, data_all_i.iloc[9, 3:8].to_numpy(dtype=int))
        num_t0 = np.append(num_t0, data_all_i.iloc[40, 3:7].to_numpy(dtype=int))
        num_t0 = np.append(num_t0, data_all_i.iloc[48, 3:7].to_numpy(dtype=int))
        Veh_3d = data_all_i.iloc[2, 3:9].to_numpy(dtype=int)
        Veh_6d = data_all_i.iloc[3, 3:9].to_numpy(dtype=int)
        Veh_9d = data_all_i.iloc[4, 3:9].to_numpy(dtype=int)

        TP_3d = data_all_i.iloc[10, 3:9].to_numpy(dtype=int)
        TP_6d = data_all_i.iloc[11, 3:7].to_numpy(dtype=int)
        TP_9d = data_all_i.iloc[12, 3:7].to_numpy(dtype=int)

        ICI_3d = data_all_i.iloc[17, 3:7].to_numpy(dtype=int)
        ICI_6d = data_all_i.iloc[18, 3:7].to_numpy(dtype=int)
        ICI_9d = data_all_i.iloc[19, 3:7].to_numpy(dtype=int)

        Palbo_3d = data_all_i.iloc[25, 3:7].to_numpy(dtype=int)
        Palbo_6d = data_all_i.iloc[26, 3:7].to_numpy(dtype=int)
        Palbo_9d = data_all_i.iloc[27, 3:7].to_numpy(dtype=int)

        TP_ICI_3d = data_all_i.iloc[33, 3:7].to_numpy(dtype=int)
        TP_ICI_6d = data_all_i.iloc[34, 3:7].to_numpy(dtype=int)
        TP_ICI_9d = data_all_i.iloc[35, 3:7].to_numpy(dtype=int)

        TP_ICI_Palbofull_3d = data_all_i.iloc[41, 3:7].to_numpy(dtype=int)
        TP_ICI_Palbofull_6d = data_all_i.iloc[42, 3:7].to_numpy(dtype=int)
        TP_ICI_Palbofull_9d = data_all_i.iloc[43, 3:7].to_numpy(dtype=int)

        TP_ICI_Palbohalf_3d = data_all_i.iloc[49, 3:9].to_numpy(dtype=int)
        TP_ICI_Palbohalf_6d = data_all_i.iloc[50, 3:7].to_numpy(dtype=int)
        TP_ICI_Palbohalf_9d = data_all_i.iloc[51, 3:7].to_numpy(dtype=int)

        Veh_i = [num_t0, Veh_3d, Veh_6d, Veh_9d]
        TP_i = [num_t0, TP_3d, TP_6d, TP_9d]
        ICI_i = [num_t0, ICI_3d, ICI_6d, ICI_9d]
        Palbo_i = [num_t0, Palbo_3d, Palbo_6d, Palbo_9d]
        TP_ICI_i = [num_t0, TP_ICI_3d, TP_ICI_6d, TP_ICI_9d]
        TP_ICI_Palbofull_i = [num_t0, TP_ICI_Palbofull_3d, TP_ICI_Palbofull_6d, TP_ICI_Palbofull_9d]
        TP_ICI_Palbohalf_i = [num_t0, TP_ICI_Palbohalf_3d, TP_ICI_Palbohalf_6d, TP_ICI_Palbohalf_9d]
        DPM_data_analysis_cellcount(time, Veh_i, condition[0], ylim)
        plt.close('all')
        DPM_data_analysis_cellcount(time, TP_i, condition[1], ylim)
        plt.close('all')
        DPM_data_analysis_cellcount(time, ICI_i, condition[2], ylim)
        plt.close('all')
        DPM_data_analysis_cellcount(time, Palbo_i, condition[3], ylim)
        plt.close('all')
        DPM_data_analysis_cellcount(time, TP_ICI_i, condition[4], ylim)
        plt.close('all')
        DPM_data_analysis_cellcount(time, TP_ICI_Palbofull_i, condition[5], ylim)
        plt.close('all')
        DPM_data_analysis_cellcount(time, TP_ICI_Palbohalf_i, condition[6], ylim)
        plt.close('all')

##
# p100
filename_p100_ctrl = ('../Data/Exp6/Riggins Bahnassy 052021 test/b052021_301.pckl',  # day 0
                      '../Data/Exp6/Riggins Bahnassy 052321/b052322_101.pckl',   # day 3
                      '../Data/Exp6/Riggins Bahnassy 052622/a052622_201.pckl',   # day 6
                      '../Data/Exp6/Riggins Bahnassy 053122/a053122_101.pckl',   # day 11
                      '../Data/Exp6/Riggins Bahnassy 053122/D11 CTL P100.pckl',  # day 11
                      '../Data/Exp6/Riggins Bahnassy 060122/D12 CTL P100.pckl'   # day 12
                      )
filename_p100_CTV = ('../Data/Exp6/Riggins Bahnassy 052021 test/b052021_302.pckl',  # day 0
                     '../Data/Exp6/Riggins Bahnassy 052321/b052322_102.pckl',  # day 3
                     '../Data/Exp6/Riggins Bahnassy 052622/a052622_202.pckl',  # day 6
                     '../Data/Exp6/Riggins Bahnassy 053122/a053122_102.pckl',  # day 11
                     '../Data/Exp6/Riggins Bahnassy 053122/D11 CTV P100.pckl',  # day 11
                     '../Data/Exp6/Riggins Bahnassy 060122/D12 CTV P100.pckl'  # day 12
                     )
# p117
filename_p117_ctrl = ('../Data/Exp6/Riggins Bahnassy 052021 test/b052021_303.pckl',  # day 0
                      '../Data/Exp6/Riggins Bahnassy 052321/b052322_103.pckl',   # day 3
                      '../Data/Exp6/Riggins Bahnassy 052622/a052622_203.pckl',   # day 6
                      '../Data/Exp6/Riggins Bahnassy 053122/a053122_103.pckl',   # day 11
                      '../Data/Exp6/Riggins Bahnassy 053122/D11 CTL P117.pckl',  # day 11
                      '../Data/Exp6/Riggins Bahnassy 060122/D12 CTL P117.pckl'   # day 12
                      )
filename_p117_CTV = ('../Data/Exp6/Riggins Bahnassy 052021 test/b052021_304.pckl',  # day 0
                     '../Data/Exp6/Riggins Bahnassy 052321/b052322_104.pckl',  # day 3
                     '../Data/Exp6/Riggins Bahnassy 052622/a052622_204.pckl',  # day 6
                     '../Data/Exp6/Riggins Bahnassy 053122/a053122_104.pckl',  # day 11
                     '../Data/Exp6/Riggins Bahnassy 053122/D11 CTV P117.pckl',  # day 11
                     '../Data/Exp6/Riggins Bahnassy 060122/D12 CTV P117.pckl'  # day 12
                     )
# GFPwt
filename_GFPwt_ctrl = ('../Data/Exp6/Riggins Bahnassy 052021 test/b052021_305.pckl',  # day 0
                       '../Data/Exp6/Riggins Bahnassy 052321/b052322_105.pckl',  # day 3
                       '../Data/Exp6/Riggins Bahnassy 052622/a052622_205.pckl',  # day 6
                       '../Data/Exp6/Riggins Bahnassy 053122/a053122_105.pckl',  # day 11
                       '../Data/Exp6/Riggins Bahnassy 053122/D11 CTL GFPwt P3.pckl',  # day 11
                       '../Data/Exp6/Riggins Bahnassy 060122/D12 CTL GFPwt P3.pckl'  # day 12
                       )
filename_GFPwt_CTV = ('../Data/Exp6/Riggins Bahnassy 052021 test/b052021_306.pckl',  # day 0
                      '../Data/Exp6/Riggins Bahnassy 052321/b052322_106.pckl',  # day 3
                      '../Data/Exp6/Riggins Bahnassy 052622/a052622_206.pckl',  # day 6
                      '../Data/Exp6/Riggins Bahnassy 053122/a053122_106.pckl',  # day 11
                      '../Data/Exp6/Riggins Bahnassy 053122/D11 CTV GFPwt P3.pckl',  # day 11
                      '../Data/Exp6/Riggins Bahnassy 060122/D12 CTV GFPwt P3.pckl'  # day 12
                      )
# test data
filename_test = ('./fc_data/Day 0 CTV_2uM_(b031821_304).pckl',  # day 0
                 './fc_data/Day 5 CTV_2uM_(b032421_102).pckl',  # day 5
                 './fc_data/Day 5 unstained ctrl (b032421_101).pckl',  # day 5
                 )

# p100 ctrl
p100_ctrl_live, p100_ctrl_dead = DPM_read_fcspckl(filename_p100_ctrl)
del filename_p100_ctrl
# p100 CTV
p100_CTV_live, p100_CTV_dead = DPM_read_fcspckl(filename_p100_CTV)
del filename_p100_CTV
# p117 ctrl
p117_ctrl_live, p117_ctrl_dead = DPM_read_fcspckl(filename_p117_ctrl)
del filename_p117_ctrl
# p117 CTV
p117_CTV_live, p117_CTV_dead = DPM_read_fcspckl(filename_p117_CTV)
del filename_p117_CTV
# GFPwt ctrl
GFPwt_ctrl_live, GFPwt_ctrl_dead = DPM_read_fcspckl(filename_GFPwt_ctrl)
del filename_GFPwt_ctrl
# GFPwt CTV
GFPwt_CTV_live, GFPwt_CTV_dead = DPM_read_fcspckl(filename_GFPwt_CTV)
del filename_GFPwt_CTV
# test data
test_live, test_dead = DPM_read_fcspckl(filename_test)
del filename_test

all_ctrl = []
all_ctrl.extend(p100_ctrl_live)
all_ctrl.extend(p117_ctrl_live)
all_ctrl.extend(GFPwt_ctrl_live)
palette = sns.color_palette(None, len(all_ctrl))
leg_str = ('p100_ctrl_d0', 'p100_ctrl_d3', 'p100_ctrl_d6', 'p100_ctrl_d11', 'p100_ctrl_d11', 'p100_ctrl_d12',
           'p117_ctrl_d0', 'p117_ctrl_d3', 'p117_ctrl_d6', 'p117_ctrl_d11', 'p117_ctrl_d11', 'p117_ctrl_d12',
           'GFPwt_ctrl_d0', 'GFPwt_ctrl_d3', 'GFPwt_ctrl_d6', 'GFPwt_ctrl_d11', 'GFPwt_ctrl_d11', 'GFPwt_ctrl_d12',)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.gca().set_xscale('log')
for i, i_ctrl in enumerate(all_ctrl):
    # i_ctrl = [x for x in i_ctrl if x > 0]
    i_leg = leg_str[i]
    n, bins, _ = ax.hist(i_ctrl, bins=np.logspace(np.log10(1e0), np.log10(1e3), int(2e2)), histtype='step', density=True,
                         color=palette[i], label=i_leg)
    ax.legend(loc='best', prop={'size': 'small'}, frameon=False)
plt.close()
all_ctrl_toge = []
[all_ctrl_toge.extend(x) for x in all_ctrl]
# n, bins, _ = plt.hist(all_ctrl_toge, bins=np.logspace(np.log10(1e0), np.log10(1e3), int(2e2)), histtype='step', density=True, color=palette[0])
# plt.gca().set_xscale('log')
ctrl_mean = stat.mean(all_ctrl_toge)
ctrl_median = stat.median(all_ctrl_toge)

# p100_CTV_live
# p100_CTV_range = [(1e3, 260000),  # day0
#                   (300, np.inf),  # day3
#                   (260, np.inf),  # day6
#                   (10, np.inf),   # day11
#                   (10, np.inf),   # day11
#                   (-np.inf, np.inf)  # day12
#                   ]
# p100_CTV_leg = ['day00', 'day03', 'day06', 'day11', 'day11', 'day12', 'ctrl']
# stat.mean(p100_CTV_live[0])
# stat.median(p100_CTV_live[0])
# p100_CTV_live_sel = []
# for i, i_p100_CTV_live in enumerate(p100_CTV_live):
#     p100_CTV_live_sel.append([x for x in i_p100_CTV_live if p100_CTV_range[i][0] < x < p100_CTV_range[i][1]])
#
# n, bins, _ = plt.hist(p100_CTV_live_sel[0], bins=np.logspace(np.log10(1e0), np.log10(1e7), int(5e2)), histtype='step',
#                       density=True, color='k', label='day 0')
# p100_CTV_live_shape, p100_CTV_live_loc, p100_CTV_live_scale = lognorm.fit(p100_CTV_live_sel[0], loc=0)
# pdf_fit = lognorm.pdf(bins, p100_CTV_live_shape, p100_CTV_live_loc, p100_CTV_live_scale)
# plt.plot(bins, pdf_fit, label='pdf fit day0')
# plt.legend(loc='best', prop={'size': 'small'}, frameon=False)
# plt.gca().set_xscale('log')
#
# p100_ctrl_liveall = []
# [p100_ctrl_liveall.extend(x) for x in p100_ctrl_live]
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# plt.gca().set_xscale('log')
# for i, i_p100_CTV_live_sel in enumerate(p100_CTV_live_sel):
#     print(len(i_p100_CTV_live_sel))
#     n, bins, _ = plt.hist(i_p100_CTV_live_sel, bins=np.logspace(np.log10(1e0), np.log10(1e6), int(6e2)), histtype='step',
#                           density=False, color=palette[i], label=p100_CTV_leg[i]+' '+str(len(i_p100_CTV_live_sel)))
# # plot ctrl
# for i, i_p100_ctrl_live in enumerate(p100_ctrl_live):
#     i_lab = 'ctrl' if i == 0 else None
#     n, bins, _ = plt.hist(i_p100_ctrl_live, bins=np.logspace(np.log10(1e0), np.log10(1e6), int(5e2)), histtype='step',
#                           density=False, color='k', label=i_lab)
# ax.legend(loc='best', prop={'size': 'small'}, frameon=False)
# mean_val = stat.median(p100_CTV_live[0])  # stat.mean(p100_CTV_live[0]) stat.median(p100_CTV_live[0])
# for i in range(math.floor(math.log2(mean_val))):
#     plt.axvline(mean_val/(2**i), 0, 1e3, color='b', linestyle='-')
#     plt.text(mean_val/(2**i), 1e3, str(i), fontdict=None)
# plt.ylim([0, 1e3])
# plt.title('p100', pad=20)
# fig.set_size_inches(15, 8)
# fig.tight_layout()
# plt.close('all')
# del p100_CTV_live, p100_ctrl_live
#
# # p117_CTV_live
# p117_CTV_range = [(1e3, 260000),  # day0
#                   (-np.inf, np.inf),  # day3
#                   (-np.inf, np.inf),  # day6
#                   (-np.inf, np.inf),   # day11
#                   (-np.inf, np.inf),   # day11
#                   (-np.inf, np.inf)  # day12
#                   ]
# p117_CTV_leg = ['day00', 'day03', 'day06', 'day11', 'day11', 'day12', 'ctrl']
# stat.mean(p117_CTV_live[0])
# stat.median(p117_CTV_live[0])
# p117_CTV_live_sel = []
# for i, i_p117_CTV_live in enumerate(p117_CTV_live):
#     p117_CTV_live_sel.append([x for x in i_p117_CTV_live if p117_CTV_range[i][0] < x < p117_CTV_range[i][1]])
#
# n, bins, _ = plt.hist(p117_CTV_live_sel[0], bins=np.logspace(np.log10(1e0), np.log10(1e7), int(5e2)), histtype='step',
#                       density=True, color='k', label='day 0')
# p117_CTV_live_shape, p117_CTV_live_loc, p117_CTV_live_scale = lognorm.fit(p117_CTV_live_sel[0], loc=0)
# pdf_fit = lognorm.pdf(bins, p117_CTV_live_shape, p117_CTV_live_loc, p117_CTV_live_scale)
# plt.plot(bins, pdf_fit, label='pdf fit day0')
# plt.legend(loc='best', prop={'size': 'small'}, frameon=False)
# plt.gca().set_xscale('log')
#
# p117_ctrl_liveall = []
# [p117_ctrl_liveall.extend(x) for x in p117_ctrl_live]
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# plt.gca().set_xscale('log')
# for i, i_p117_CTV_live_sel in enumerate(p117_CTV_live_sel):
#     print(len(i_p117_CTV_live_sel))
#     n, bins, _ = plt.hist(i_p117_CTV_live_sel, bins=np.logspace(np.log10(1e0), np.log10(1e6), int(6e2)), histtype='step',
#                           density=False, color=palette[i], label=p117_CTV_leg[i]+' '+str(len(i_p117_CTV_live_sel)))
# # plot ctrl
# for i, i_p117_ctrl_live in enumerate(p117_ctrl_live):
#     i_lab = 'ctrl' if i == 0 else None
#     n, bins, _ = plt.hist(i_p117_ctrl_live, bins=np.logspace(np.log10(1e0), np.log10(1e6), int(5e2)), histtype='step',
#                           density=False, color='k', label=i_lab)
# ax.legend(loc='best', prop={'size': 'small'}, frameon=False)
# mean_val = stat.median(p117_CTV_live[0])  # stat.mean(p117_CTV_live[0]) stat.median(p117_CTV_live[0])
# for i in range(math.floor(math.log2(mean_val))):
#     plt.axvline(mean_val/(2**i), 0, 1e3, color='b', linestyle='-')
#     plt.text(mean_val/(2**i), 1e3, str(i), fontdict=None)
# plt.ylim([0, 1e3])
# plt.title('p117', pad=20)
# fig.set_size_inches(15, 8)
# fig.tight_layout()
# plt.close('all')
# del p117_CTV_live, p117_ctrl_live
#
# # GFPwt_CTV_live
# GFPwt_CTV_range = [(1e3, 260000),  # day0
#                    (-np.inf, np.inf),  # day3
#                    (-np.inf, np.inf),  # day6
#                    (-np.inf, np.inf),   # day11
#                    (-np.inf, np.inf),   # day11
#                    (-np.inf, np.inf)  # day12
#                   ]
# GFPwt_CTV_leg = ['day00', 'day03', 'day06', 'day11', 'day11', 'day12', 'ctrl']
# stat.mean(GFPwt_CTV_live[0])
# stat.median(GFPwt_CTV_live[0])
# GFPwt_CTV_live_sel = []
# for i, i_GFPwt_CTV_live in enumerate(GFPwt_CTV_live):
#     GFPwt_CTV_live_sel.append([x for x in i_GFPwt_CTV_live if GFPwt_CTV_range[i][0] < x < GFPwt_CTV_range[i][1]])
#
# n, bins, _ = plt.hist(GFPwt_CTV_live_sel[0], bins=np.logspace(np.log10(1e0), np.log10(1e7), int(5e2)), histtype='step',
#                       density=True, color='k', label='day 0')
# GFPwt_CTV_live_shape, GFPwt_CTV_live_loc, GFPwt_CTV_live_scale = lognorm.fit(GFPwt_CTV_live_sel[0], loc=0)
# pdf_fit = lognorm.pdf(bins, GFPwt_CTV_live_shape, GFPwt_CTV_live_loc, GFPwt_CTV_live_scale)
# plt.plot(bins, pdf_fit, label='pdf fit day0')
# plt.legend(loc='best', prop={'size': 'small'}, frameon=False)
# plt.gca().set_xscale('log')
#
# GFPwt_ctrl_liveall = []
# [GFPwt_ctrl_liveall.extend(x) for x in GFPwt_ctrl_live]
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# plt.gca().set_xscale('log')
# for i, i_GFPwt_CTV_live_sel in enumerate(GFPwt_CTV_live_sel):
#     print(len(i_GFPwt_CTV_live_sel))
#     n, bins, _ = plt.hist(i_GFPwt_CTV_live_sel, bins=np.logspace(np.log10(1e0), np.log10(1e6), int(6e2)), histtype='step',
#                           density=False, color=palette[i], label=GFPwt_CTV_leg[i]+' '+str(len(i_GFPwt_CTV_live_sel)))
# # plot ctrl
# for i, i_GFPwt_ctrl_live in enumerate(GFPwt_ctrl_live):
#     i_lab = 'ctrl' if i == 0 else None
#     n, bins, _ = plt.hist(i_GFPwt_ctrl_live, bins=np.logspace(np.log10(1e0), np.log10(1e6), int(5e2)), histtype='step',
#                           density=False, color='k', label=i_lab)
# ax.legend(loc='best', prop={'size': 'small'}, frameon=False)
# mean_val = stat.median(GFPwt_CTV_live[0])  # stat.mean(GFPwt_CTV_live[0]) stat.median(GFPwt_CTV_live[0])
# for i in range(math.floor(math.log2(mean_val))):
#     plt.axvline(mean_val/(2**i), 0, 1e3, color='b', linestyle='-')
#     plt.text(mean_val/(2**i), 1e3, str(i), fontdict=None)
# plt.ylim([0, 1e3])
# plt.title('GFPwt', pad=20)
# fig.set_size_inches(15, 8)
# fig.tight_layout()
# plt.close('all')
# del GFPwt_CTV_live, GFPwt_ctrl_live

# test
test_range = [(1e3, 260000),  # day0
              (-np.inf, np.inf),  # day5
              ]
test_ctrl_live = test_live[-1]
test_live.pop(-1)
test_leg = ['day0', 'day5', 'ctrl']
stat.mean(test_live[0])
stat.median(test_live[0])
test_live_sel = []
for i, i_test_live in enumerate(test_live):
    test_live_sel.append([x for x in i_test_live if test_range[i][0] < x < test_range[i][1]])

n, bins, _ = plt.hist(test_live_sel[0], bins=np.logspace(np.log10(1e0), np.log10(1e7), int(5e2)), histtype='step',
                      density=True, color='k', label='day 0')
test_live_shape, test_live_loc, test_live_scale = lognorm.fit(test_live_sel[0], loc=0)
pdf_fit = lognorm.pdf(bins, test_live_shape, test_live_loc, test_live_scale)
plt.plot(bins, pdf_fit, label='pdf fit day0')
plt.legend(loc='best', prop={'size': 'small'}, frameon=False)
plt.gca().set_xscale('log')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.gca().set_xscale('log')
for i, i_test_live_sel in enumerate(test_live_sel):
    print(len(i_test_live_sel))
    # n, bins, _ = plt.hist(i_test_live_sel, bins=np.logspace(np.log10(1e0), np.log10(1e6), int(6e2)), histtype='step',
    #                       density=False, color=palette[i], label=test_leg[i]+' '+str(len(i_test_live_sel)))
    n, bin_edges = np.histogram(i_test_live_sel, bins=np.logspace(np.log10(1e0), np.log10(1e6), int(6e2)), density=False)
    x_bins = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(bin_edges) - 1)]
    plt.plot(x_bins, n, label=test_leg[i])
# plot ctrl
# n, bins, _ = plt.hist(test_ctrl_live, bins=np.logspace(np.log10(1e0), np.log10(1e6), int(5e2)), histtype='step',
#                       density=False, color='k', label='ctrl')
n, bin_edges = np.histogram(test_ctrl_live, bins=np.logspace(np.log10(1e0), np.log10(1e6), int(5e2)), density=False)
x_bins = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(bin_edges) - 1)]
plt.plot(x_bins, n, label='ctrl')
ax.legend(loc='best', prop={'size': 'small'}, frameon=False)
mean_val = test_live_scale  # stat.mean(GFPwt_CTV_live[0]) stat.median(GFPwt_CTV_live[0])
for i in range(math.floor(math.log2(mean_val))):
    plt.axvline(mean_val/(2**i), 0, 1e3, color='b', linestyle='-')
    plt.text(mean_val/(2**i), 1e3, str(i), fontdict=None)
plt.ylim([0, 1e3])
plt.title('test', pad=20)
fig.set_size_inches(15, 8)
fig.tight_layout()
plt.close('all')
del test_live, test_ctrl_live

test_live = dict({0:  dict({'cellnum': len(test_live_sel[0]), 'gen': 0})})
plt_flag = True
method = 'OP'  # OP/EM
generation = [4, 5, 6]  # 4 5 6
mingen, maxgen = 4, 6
time = 5  # days
fix = False
cellnum, gen, scale, shape, loc, weight, aic, pdf, pdf_each, data, bins, x_bins, comb = \
    DPM_fcs_analysis_fitlognorm(test_live_sel[1], mean_val, generation, mingen, maxgen, method, fix)
# 1, 2: (4,5), 3: (4,6), 4: (5,6), 5: (4,5,6)

if plt_flag:
    plt.rcParams['figure.figsize'] = (7, 6)
    plt.rcParams['font.size'] = 14
    plt.gca().set_xscale('log')
    n, bin_edges = np.histogram(data, bins=np.exp(bins), density=True)
    plt.plot(np.exp(x_bins), n, label='pdf data')
    plt.plot(np.exp(x_bins), pdf, label='pdf fit total case 4')
    for i, i_pdf in enumerate(pdf_each):
        plt.plot(np.exp(x_bins), i_pdf, label='pdf fit '+str(i))
    plt.ylabel('Probability density')
    plt.xlabel('Cell trace violet value')
    plt.legend(loc='upper left', prop={'size': 'large'}, frameon=False)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
test_live.update({time: dict({'cellnum': cellnum, 'gen': gen}), 'maxgen': max(gen)})
plt.close()

plt.rcParams['figure.figsize'] = (4, 3)
plt.rcParams['font.size'] = 12
plt.plot(aic[2:])
plt.title('AIC for different cases')
# plt.ylabel('AIC [Akaike information criterion]')
plt.xticks([0, 1, 2, 3, 4], labels=['case 1', 'case 2', 'case 3', 'case 4', 'case 5'])
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.tight_layout()

model = {3}  # 1,2,3
model1_par, model2_par, model3_par = DPM_fcs_model_fit(test_live, model)
if model1_par is not None:
    DPM_fcs_plot_model_12(model1_par, test_live)
    plt.text(3, 2.3e4, 'model 1', horizontalalignment='center')
plt.close()
if model2_par is not None:
    DPM_fcs_plot_model_12(model2_par, test_live)
    plt.text(3, 2.3e4, 'model 2', horizontalalignment='center')
plt.close()
if model3_par is not None:
    DPM_fcs_plot_model_3(model3_par, test_live)
    fig.axes[0].text(3, 2e4, 'model 3', horizontalalignment='center')
