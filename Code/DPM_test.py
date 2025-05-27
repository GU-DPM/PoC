from DPM_run import *

# DPM_run_generate_R1R2_pdftrue()
# DPM_run_analysis(pathload=pathload, pathsave=pathsave, filename_pattern='result', Strategy_name=['strategy0', 'strategy2.2'])
##
filename_pattern = '_para'
subcolone_LOD = [1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01]
misspecification_LOD = ['pdf', 'loguni']
pathload = './pnas'
for i_subcolone_LOD in subcolone_LOD:
    for i_misspecification_LOD in misspecification_LOD:
        pathsave = f'./misspecification/PNAS_mis_LOD/{i_subcolone_LOD:.0e}/{i_misspecification_LOD}/'
        DPM_run_par_csv_folder(pathload=pathload,
                               pathsave=pathsave,
                               filename_pattern=filename_pattern,
                               run_sim=True,
                               subcolone_LOD=i_subcolone_LOD,
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


# pathload = './misspecification/PNAS/'
# pathsave = './misspecification/PNAS/figure'
# for i_LOD in ['1e-06', '1e-05', '1e-04', '1e-03', '1e-02', '1e-01']:
#     pathloadmis = './csv_result/2drug3million/PNAS_mis_LOD/' + i_LOD
#     DPM_run_processing(Num_drug=2,
#                        pathload=pathload,
#                        pathsave=pathsave,
#                        Num_stepdiff=2,
#                        Strategy_name=['strategy0', 'strategy2.2'],
#                        pathloadmis='',
#                        use_parallel=False)

# subcolone_LOD = 1e-6
filename_pattern = '_para'

# misspecification_LOD = 0
# pathload = './csv_data/2drug3million/now'
# pathsave = f'./csv_result/2drug3million/mis_LOD/{subcolone_LOD:.0e}/{misspecification_LOD}/'

subcolone_LOD = [1e-05, 1e-04, 1e-03, 1e-02, 1e-01]
misspecification_LOD = [0, 'pdf', 'loguni', 'max']
pathload = './pnas'
# for i_subcolone_LOD in subcolone_LOD:
#     for i_misspecification_LOD in misspecification_LOD:
#         pathsave = f'./misspecification/PNAS_mis_LOD/{i_subcolone_LOD:.0e}/{i_misspecification_LOD}/'
i_subcolone_LOD = 1e-03
i_misspecification_LOD = 'loguni'
pathsave = f'./misspecification/PNAS_mis_LOD/{i_subcolone_LOD:.0e}/{i_misspecification_LOD}/'
# pathsave = f'./misspecification/test/'
DPM_run_par_csv_folder(pathload=pathload,
                       pathsave=pathsave,
                       filename_pattern=filename_pattern,
                       run_sim=True,
                       subcolone_LOD=i_subcolone_LOD,
                       misspecification_LOD=i_misspecification_LOD,
                       Strategy_name=['strategy0', 'strategy2.2'],
                       fullinput=True,
                       save_filename_param=True,
                       save_filename_stopt=True,
                       save_filename_dosage=True,
                       save_filename_pop=False,
                       save_filename_eachtimepoint=False,
                       misspecification_sim_only=True,
                       # misspecification=False,
                       use_parallel=False)

misspecification_LOD = 'pdf loguni'
pathload = './csv_data/2drug3million/now'
pathsave = f'./csv_result/2drug3million/mis_LOD/{subcolone_LOD:.0e}/{misspecification_LOD}/'
DPM_run_par_csv_folder(pathload=pathload,
                       pathsave=pathsave,
                       filename_pattern=filename_pattern,
                       Num_drug=2,
                       dose_method='d',
                       Simduration=1800,
                       Stepsize=45,
                       run_sim=False,
                       subcolone_LOD=subcolone_LOD,
                       misspecification_LOD=misspecification_LOD,
                       mutation_rate=7.1e-7,
                       Strategy_name=['strategy0', 'strategy2.2'],
                       fullinput=True,
                       save_filename_param=True,
                       save_filename_stopt=True,
                       save_filename_dosage=True,
                       save_filename_pop=False,
                       save_filename_eachtimepoint=False,
                       misspecification_sim_only=True,
                       misspecification=False,
                       use_parallel=False)


# 2 drug case
filename_2drug_1 = './csv_data/example_DPM_parameter_2drug.csv'
filename_2drug_2 = './csv_data/example_DPM_parameter_2drug2.csv'
filename_2drug_1_mis = './csv_data/example_DPM_parameter_2drug_mis.csv'
filename_2drug_2_mis = './csv_data/example_DPM_parameter_2drug_mis2.csv'
filename_2drugtest = './csv_data/example_DPM_parameter_2drugtest.csv'
filename_2drugtest_mis = './csv_data/example_DPM_parameter_2drugtest_mis.csv'
# 3 drug case
filename_3drug_1 = './csv_data/example_DPM_parameter_3drug.csv'
filename_3drug_2 = './csv_data/example_DPM_parameter_3drug2.csv'
filename_3drug_1_mis = './csv_data/example_DPM_parameter_3drug_mis.csv'
filename_3drug_2_mis = './csv_data/example_DPM_parameter_3drug_mis2.csv'
# misspecification_filename_csv=filename_2drugtest_mis,
misspecification_filename_csv = filename_2drug_1_mis,
DPM_run_par_csv(filename_csv=filename_2drugtest,
                par_ind=0,
                Num_drug=2,
                pathsave='./csv_result/',
                dose_method='d',
                dose_interval=0.1,
                use_input_dose_combination_only=False,
                Simduration=1800,
                Limit_mortality=1e13,
                Stepsize=45,
                Limit_radiologicdetection=1e9,
                run_sim=True,
                erase_preresult=True,
                PAR_criterisa=[],
                Strategy_name=['strategy0', 'strategy2.2'],
                fullinput=True,
                save_filename_param=True,
                save_filename_stopt=True,
                save_filename_pop=False,
                save_filename_dosage=True,
                save_filename_eachtimepoint=True,
                misspecification_sim_only=False,
                lookahead_step=5,
                Maxnum_subseq=500, subtreedepth=5)

# DPM_run_par_csv(filename_csv=[filename_2drug_1, filename_2drug_2], misspecification_filename_csv=[filename_2drug_1_mis, filename_2drug_2_mis],
#                 Num_drug=2, dose_method='d', dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.5, 0.5),
#                 Simduration=1800, Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True,
#                 PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True,
#                 save_filename_pop=True, save_filename_dosage=True, save_filename_eachtimepoint=True,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=500, subtreedepth=5)

# DPM_run_par_csv(filename_csv=filename_2drug_1, Num_drug=2, dose_method='d',
#                 dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.5, 0.5), Simduration=1800,
#                 Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True, PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True,
#                 save_filename_pop=True, save_filename_dosage=True, save_filename_eachtimepoint=True,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=500, subtreedepth=5)

# DPM_run_par_csv(filename_csv=[filename_2drug_1, filename_2drug_2],
#                 Num_drug=2, dose_method='d',
#                 dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.5, 0.5), Simduration=1800,
#                 Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True, PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True,
#                 save_filename_pop=True, save_filename_dosage=True, save_filename_eachtimepoint=True,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=500, subtreedepth=5)


# DPM_run_par_csv(filename_csv=filename_2drug_1, misspecification_filename_csv=filename_2drug_1_mis, Num_drug=2, dose_method='d',
#                 dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.5, 0.5), Simduration=1800,
#                 Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True, PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True, save_filename_pop=True, save_filename_dosage=True,
#                 save_filename_eachtimepoint=False, misspecification_sim_only=True,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=500, subtreedepth=5)
#

# DPM_run_par_csv(filename_csv=[filename_2drug_1, filename_2drug_2], misspecification_filename_csv=[filename_2drug_1_mis, filename_2drug_2_mis],
#                 Num_drug=2, dose_method='d', dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.5, 0.5),
#                 Simduration=1800, Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True,
#                 PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True, save_filename_pop=True, save_filename_dosage=True,
#                 save_filename_eachtimepoint=True, misspecification_sim_only=True,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=500, subtreedepth=5)

# # (1) 3 drugs.
# DPM_run_par_csv(filename_csv=filename_3drug_1, misspecification_filename_csv=filename_3drug_1_mis, Num_drug=3, dose_method='d',
#                 dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.2, 0.3, 0.5), Simduration=1800,
#                 Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True, PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True, save_filename_pop=True, save_filename_dosage=True,
#                 save_filename_eachtimepoint=False, misspecification_sim_only=False,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=50, subtreedepth=5)
# # (2) 3 drugs.
# DPM_run_par_csv(filename_csv=[filename_3drug_1, filename_3drug_2], misspecification_filename_csv=[filename_3drug_1_mis, filename_3drug_2_mis],
#                 Num_drug=3, dose_method='d',
#                 dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.2, 0.3, 0.5), Simduration=1800,
#                 Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True, PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True, save_filename_pop=True, save_filename_dosage=True,
#                 save_filename_eachtimepoint=False, misspecification_sim_only=False,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=50, subtreedepth=5)
# # (3) 3 drugs.
# DPM_run_par_csv(filename_csv=filename_3drug_1, Num_drug=3, dose_method='d',
#                 dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.2, 0.3, 0.5), Simduration=1800,
#                 Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True, PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True, save_filename_pop=True, save_filename_dosage=True,
#                 save_filename_eachtimepoint=False, misspecification_sim_only=False,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=50, subtreedepth=5)
# # (4) 3 drugs.
# DPM_run_par_csv(filename_csv=[filename_3drug_1, filename_3drug_2], Num_drug=3, dose_method='d',
#                 dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.2, 0.3, 0.5), Simduration=1800,
#                 Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True, PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True, save_filename_pop=True, save_filename_dosage=True,
#                 save_filename_eachtimepoint=False, misspecification_sim_only=False,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=50, subtreedepth=5)
# # (5) 3 drugs.
# DPM_run_par_csv(filename_csv=filename_3drug_1, misspecification_filename_csv=filename_3drug_1_mis, Num_drug=3, dose_method='d',
#                 dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.2, 0.3, 0.5), Simduration=1800,
#                 Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True, PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True, save_filename_pop=True, save_filename_dosage=True,
#                 save_filename_eachtimepoint=False, misspecification_sim_only=True,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=50, subtreedepth=5)
# # (6) 3 drugs.
# DPM_run_par_csv(filename_csv=[filename_3drug_1, filename_3drug_2], misspecification_filename_csv=[filename_3drug_1_mis, filename_3drug_2_mis],
#                 Num_drug=3, dose_method='d',
#                 dose_interval=0.01, use_input_dose_combination_only=False, dose_combination=(0.2, 0.3, 0.5), Simduration=1800,
#                 Limit_mortality=1e13, Stepsize=45, Limit_radiologicdetection=1e9, run_sim=True, erase_preresult=True, PAR_criterisa=[],
#                 Strategy_name=['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4', 'strategy5', 'strategy6',
#                                'strategy7', 'strategy8', 'strategy9'],
#                 fullinput=True, save_filename_param=True, save_filename_stopt=True, save_filename_pop=True, save_filename_dosage=True,
#                 save_filename_eachtimepoint=False, misspecification_sim_only=True,
#                 pathsave='./csv_result', lookahead_step=5, Maxnum_subseq=50, subtreedepth=5)


# par = {'Num_drug': 2, 'Spop': 4.995e9, 'R1pop': 500, 'R2pop': 5e6, 'R12pop': 0, 'g0_S': 0.1287, 'g0_R1': 0.1287, 'g0_R2': 0.1287, 'g0_R12': 0.1287,
#        'Sa.S.D1.': 5.895, 'Sa.S.D2.': 0.1179, 'Sa.R1.D1.': 0.0053917, 'Sa.R1.D2.': 0.1179, 'Sa.R2.D1.': 5.895, 'Sa.R2.D2.': 0.094321,
#        'Sa.R12.D1.': 0.0053917, 'Sa.R12.D2.': 0.094321, 'T.S..S.': 0, 'T.S..R1.': 0, 'T.S..R2.': 0, 'T.S..R12.': 0, 'T.R1..S.': 1e-11,
#        'T.R1..R1.': 0, 'T.R1..R2.': 0, 'T.R1..R12.': 0, 'T.R2..S.': 2.15e-10, 'T.R2..R1.': 0, 'T.R2..R2.': 0, 'T.R2..R12.': 0, 'T.R12..S.': 0,
#        'T.R12..R1.': 2.15e-10, 'T.R12..R2.': 1e-11, 'T.R12..R12.': 0}
# DPM_run_plot_1PAR(par=par, Strategy_name=['strategy0', 'strategy2.2'], pathsave='./1par_result/')

# DPM_save_default_PARset(Num_drug=2, par_save_block_size=1e4, pathsave='./parset_data/2drug')
# DPM_save_default_PARset(Num_drug=3, par_save_block_size=1e4, pathsave='./parset_data/3drug')

if __name__ == '__main__':
    # DPM_run_generate_R1R2_pdftrue()
    pathload = './csv_result/2drug3million/pdftrue/'
    pathsave = './csv_result/2drug3million/pdftrue/figure'
    pathload_misLOD = './csv_result/2drug3million/pdftrue_mis_LOD/'

    LOD = ['1e-06',
           '1e-05',
           '1e-04',
           '1e-03',
           '1e-02',
           '1e-01']
    filename_pattern = 'result'
    DPM_run_analysis_LOD(Num_drug=2,
                         LOD=LOD,
                         filename_pattern=filename_pattern,
                         pathload=pathload,
                         pathsave=pathsave,
                         pathloadmis=pathload_misLOD,
                         Strategy_name=['strategy0', 'strategy2.2'],
                         plot=True)
