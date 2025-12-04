from DPM_run_drugmis import *
from DPM_run import DPM_run_par_csv_folder, DPM_run_generate_misDrug2, DPM_run_processing

'''
This is the main function used to run the analyses for Deepak's paper,#
"A Proof-of-Concept Clinical Trial Design for Evolutionary Guided Precision Medicine for Cancer"#
(doi: 10.1101/2025.05.23.25328210)#
Run this function to execute the study's workflow.#
'''

# (1)
'''
The first step is to generate the misspecified drug 2 sensitivies. There are six scenarios, in which the true potency of drug 2#
differs from the assumed (values in the pnas paper) potency by the following factors:#
1. Actual potency is 5x lower than assumed (pathsave='./pnas_div5', ratio=1/5),#
2. Actual potency is 10x lower than the assumed (pathsave='./pnas_div10', ratio=1/10),#
3. Actual potency is 30x lower than the assumed (pathsave='./pnas_div30', ratio=1/30),#
4. Actual potency is 5x higher than the thought (pathsave='./pnas_5x', ratio=5),#
5. Actual potency is 10x higher than the thought (pathsave='./pnas_10x', ratio=10),#
6. Actual potency is 30x higher than the thought (pathsave='./pnas_30x', ratio=30),#
'''
flag_1thstep = False
if flag_1thstep:
    pathsave = ['./div5/', './div10/', './div30/', './5x/', './10x/', './30x/']
    ratio = [1/5, 1/10, 1/30, 5, 10, 30]
    for i_pathsave, i_ratio in zip(pathsave, ratio):
        DPM_run_generate_misDrug2(pathload_='./pnas/', pathsave_=i_pathsave, ratio=i_ratio)

# (2)
'''
The second step is to run the simulations under the different misspecification scenarios.#
misspecification_atsim=True means that during simulation, the parameter values generated in step 1, which are saved in the folders:#
('./div5', './div10/', './div30/', './5x/', './10x/', './30x/')#
will be used. At the drug selection phase, the parameters from the './pnas/' folder are used.#
These scenarios correspond to those decribed in Deepak's paper.#

'misspecification_atdecision = False' (which is the default setting and should remain so) means that during decision making,# 
the parameter values always come from the './pnas/' folder. These parameters represent the estimated potency of drug 2,#
and this estimation should not be divided into different parameter sets such as those in the folders:#  
('./div5', './div10/', './div30/', './5x/', './10x/', './30x/')#
Thus, the parameter estimation happens once, and the misspecification occurs during simulation (true situations),#
where different scenarios or parameter sets from above folders are applied.#
'''
flag_2thstep = False
if flag_2thstep:
    pathload = './pnas/'
    misspecification_fileload = ['./div5/', './div10/', './div30/', './5x/', './10x/', './pnas_30x/']
    pathsave = ['./div5_atsim/', './div10_atsim/', './div30_atsim/', './5x_atsim/', './10x_atsim/', './30x_atsim/']
    filename_pattern = '_para.csv'
    for i in range(len(misspecification_fileload)):
        DPM_run_par_csv_folder(pathload=pathload,
                               pathsave=pathsave[i],
                               misspecification_fileload=misspecification_fileload[i],
                               filename_pattern=filename_pattern,
                               run_sim=True,
                               Strategy_name=['CPM', 'DPM2.2'],
                               save_filename_param=True,
                               save_filename_stopt=True,
                               save_filename_dosage=True,
                               save_filename_pop=True,
                               save_filename_eachtimepoint=False,
                               misspecification_sim_only=True,
                               misspecification_atsim=True,
                               use_parallel=False)

# (3)
'''
The third step is to run the simulations without misspecification.#
'''
flag_3thstep = False
if flag_3thstep:
    pathload = './pnas/'
    pathsave = './pnas_sim/'
    filename_pattern = '_para.csv'
    DPM_run_par_csv_folder(pathload=pathload,
                           pathsave=pathsave,
                           filename_pattern=filename_pattern,
                           run_sim=True,
                           Strategy_name=['CPM', 'DPM2.2'],
                           save_filename_param=True,
                           save_filename_stopt=True,
                           save_filename_dosage=True,
                           save_filename_pop=True,
                           save_filename_eachtimepoint=False,
                           misspecification_sim_only=False,
                           misspecification_atsim=False,
                           use_parallel=False)
# (4)
'''
The fourth step involves processing the simulation outputs, covering both misspecification and non-misspecification scenarios.
'''
flag_4thstep = True
if flag_4thstep:
    pathload = ['pnas_sim', './div5_atsim/', './div10_atsim/', './div30_atsim/', './5x_atsim/', './10x_atsim/', './30x_atsim/']
    for i_path in pathload:
        DPM_run_processing(Num_drug=2,
                           pathload=i_path,
                           Strategy_name=['CPM', 'DPM2.2'],
                           use_parallel=False)

# (5)
'''
The fifth step is to analyze the results from the fourth step and genrate the raw figures (Fig.3,4,5,S2, Tab.S1) used in Deepak's paper.#
The order must be preserved: first the three values 5x, 10x and 30x, followed by div5, div10, div30.#
'''
flag_5thstep = True
if flag_5thstep:
    mis = ['./5x_atsim/', './10x_atsim/', './30x_atsim/', './div5_atsim/', './div10_atsim/', './div30_atsim/']
    pathload = './pnas_sim/'
    pathsave = './mis_atsim/'
    pathload_mis = './'
    filename_pattern = 'result'
    strategy_name = ['CPM', 'DPM2.2']
    DPM_run_drugmis_analysis(mis=mis,
                             filename_pattern=filename_pattern,
                             pathload=pathload,
                             pathsave=pathsave,
                             pathloadmis=pathload_mis,
                             Strategy_name=strategy_name)
