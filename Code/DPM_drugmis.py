from DPM_run_drugmis import *
from DPM_run import DPM_run_par_csv_folder, DPM_run_generate_misDrug2, DPM_run_processing

'''generate parameter sets for different misspecified drug ratios.'''
# DPM_run_generate_misDrug2(pathload='./pnas/', pathsave='./pnas_div5/', ratio=1/5)
# DPM_run_generate_misDrug2(pathload='./pnas/', pathsave='./pnas_div10/', ratio=1/10)
# DPM_run_generate_misDrug2(pathload='./pnas/', pathsave='./pnas_div30/', ratio=1/30)
# DPM_run_generate_misDrug2(pathload='./pnas/', pathsave='./pnas_5x/', ratio=5)
# DPM_run_generate_misDrug2(pathload='./pnas/', pathsave='./pnas_10x/', ratio=10)
# DPM_run_generate_misDrug2(pathload='./pnas/', pathsave='./pnas_30x/', ratio=30)

# pathload = './pnas/'
# misspecification_fileload = ['./pnas_div5/', './pnas_div10/', './pnas_div30/', './pnas_5x/', './pnas_10x/', './pnas_30x/']
# pathsave = ['./div5_atsim/', './div10_atsim/', './div30_atsim/', './5x_atsim/', './10x_atsim/', './30x_atsim/']
# filename_pattern = '_para.csv'
# for i in [1]:   #range(len(misspecification_fileload)):
#     DPM_run_par_csv_folder(pathload=pathload,
#                            pathsave=pathsave[i],
#                            misspecification_fileload=misspecification_fileload[i],
#                            filename_pattern=filename_pattern,
#                            run_sim=True,
#                            Strategy_name=['strategy0', 'strategy2.2'],
#                            fullinput=True,
#                            save_filename_param=True,
#                            save_filename_stopt=True,
#                            save_filename_dosage=True,
#                            save_filename_pop=True,
#                            save_filename_eachtimepoint=False,
#                            misspecification_sim_only=True,
#                            misspecification_ofdecision=False,
#                            misspecification_atsim=True,
#                            use_parallel=False)

# pathload = ['./div5_atsim/', './div10_atsim/', './div30_atsim/', './5x_atsim/', './10x_atsim/', './30x_atsim/']
# for i_path in pathload:
#     DPM_run_processing(Num_drug=2,
#                        pathload=i_path,
#                        Strategy_name=['strategy0', 'strategy2.2'],
#                        use_parallel=False)

# mis = ['pnas_5x_oftrue', 'pnas_10x_oftrue', 'pnas_30x_oftrue', 'pnas_div5_oftrue', 'pnas_div10_oftrue', 'pnas_div30_oftrue']
# orders matters, first three, 5x, 10x and 30x, then div5, div10, div30
mis = ['./5x_atsim/', './10x_atsim/', './30x_atsim/', './div5_atsim/', './div10_atsim/', './div30_atsim/']
pathload = './pnas/'
pathsave = './mis_atsim/'
pathload_mis = './'
filename_pattern = 'result'
strategy_name = ['strategy0', 'strategy2.2']
DPM_run_drugmis_analysis(mis=mis,
                         filename_pattern=filename_pattern,
                         pathload=pathload,
                         pathsave=pathsave,
                         pathloadmis=pathload_mis,
                         Strategy_name=strategy_name)
