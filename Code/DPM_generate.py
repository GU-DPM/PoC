from DPM_lib import itertools, sys, rv_continuous, loguniform, simps, plt, deepcopy, math
from DPM_constant import *
""" This script checks arguments used in the DPM (Dynamic Precision Medicine) model. """


# Generate the default discrete dose combination based on drug number.
def DPM_generate_discrete_dose_combination(num_drug):
    dose_list = [float(0)]
    for i in range(int(num_drug)):
        dose_list.append(float(1 / (i + 1)))
    dose_list.sort()
    dose_combination = itertools.product(dose_list, repeat=int(num_drug))
    dose_combination = [x for x in dose_combination if sum(x) == 1]
    return dose_combination


# Generate the continuous dose combination based on drug number and dose interval.
def DPM_generate_continuous_dose_combination(num_drug, dose_interval):
    dose_list = np.arange(0, 1+dose_interval, dose_interval, dtype=float)
    dose_list.sort()
    dose_combination = itertools.product(dose_list, repeat=int(num_drug))
    dose_combination = [x for x in dose_combination if sum(x) == 1]
    return dose_combination


# Generate default parameters in 2 drug case according to PNAS(2012).
def DPM_generate_default_par_2drug():
    # g0: basel proliferation rate for 4 kinds for cells: S, R1, R2, R12.
    g0 = [0.001, 0.0026, 0.007, 0.0184, 0.0487, 0.1287, 0.34]
    # R1 to X0 ratio.
    ratioR1toX0 = [0, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 9e-1]
    # R2 to X0 ratio.
    ratioR2toX0 = [0, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 9e-1]
    # Sa:4×2 matrix of drug sensitivities.
    # Sa(S,D1)/g0.
    Sa_ratio_S_D1tog0 = [5.6e-4, 0.0054, 0.0517, 0.4964, 4.7683, 45.8045, 440]
    # Sa(S,D2)/Sa(S,D1).
    Sa_ratio_S_D2toS_D1 = [4e-4, 0.0015, 0.0054, 0.02, 0.0737, 0.2714, 1e0]
    # Sa(R1,D1)/Sa(S,D1).
    Sa_ratio_R1_D1toS_D1 = [0, 1e-5, 9.5635e-5, 9.1461e-4, 0.0087, 0.0837, 0.8]
    # Sa(R2,D2)/Sa(S,D2).
    Sa_ratio_R2_D2toS_D2 = [0, 1e-5, 9.5635e-5, 9.1461e-4, 0.0087, 0.0837, 0.8]
    # T:4×4 transition rate matrix, R1->R12=S->R2; R2->R12=S->R1.
    T_StoR1 = [1e-11, 2.154e-10, 4.642e-9, 1e-7, 2.154e-6, 4.642e-5, 1e-3]
    T_StoR2 = [1e-11, 2.154e-10, 4.642e-9, 1e-7, 2.154e-6, 4.642e-5, 1e-3]
    return g0, ratioR1toX0, ratioR2toX0, Sa_ratio_S_D1tog0, Sa_ratio_S_D2toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, T_StoR1, T_StoR2


# Generate default parameters in 3 drug case according to Biology Direct(2016).
def DPM_generate_default_par_3drug():
    # g0: basel proliferation rate for 4 kinds for cells: S, R1, R2, R3, R12, R13, R23, R123.
    g0 = [0.001, 0.0026, 0.007, 0.0184, 0.0487, 0.1287, 0.34]
    # S to X0 ratio.
    ratioStoX0 = [1e-1, 9e-1]
    # R1 to X0 ratio.
    ratioR1toX0 = [0e0, 1e-5, 1e-3, 1e-1, 9e-1]
    # R2 to X0 ratio.
    ratioR2toX0 = [0e0, 1e-5, 1e-3, 1e-1, 9e-1]
    # R3 to X0 ratio.
    ratioR3toX0 = [0e0, 1e-5, 1e-3, 1e-1, 9e-1]
    # R12 to X0 ratio.
    ratioR12toX0 = [0e0, 1e-5, 1e-3]
    # R23 to X0 ratio.
    ratioR23toX0 = [0e0, 1e-5, 1e-3]
    # R13 to X0 ratio.
    ratioR13toX0 = [0e0, 1e-5, 1e-3]
    # Sa(S,D1)/g0.
    Sa_ratio_S_D1tog0 = [5e-1, 5e0]
    # Sa(S,D2)/Sa(S,D1).
    Sa_ratio_S_D2toS_D1 = [1e-1, 0.3, 1e0]
    # Sa(S,D3)/Sa(S,D1).
    Sa_ratio_S_D3toS_D1 = [1e-1, 0.3, 1e0]
    # Sa(R1,D1)/Sa(S,D1).
    Sa_ratio_R1_D1toS_D1 = [0e0, 3e-1, 1e0]
    # Sa(R2,D2)/Sa(S,D2).
    Sa_ratio_R2_D2toS_D2 = [0e0, 3e-1, 1e0]
    # Sa(R3,D3)/Sa(S,D3).
    Sa_ratio_R3_D3toS_D3 = [0e0, 3e-1, 1e0]
    # T:8×8 transition rate matrix.
    T_StoR1 = [1e-11, 1e-9, 1e-7, 1e-5, 1e-3]
    T_StoR2 = [1e-11, 1e-9, 1e-7, 1e-5, 1e-3]
    T_StoR3 = [1e-11, 1e-9, 1e-7, 1e-5, 1e-3]
    return g0, ratioStoX0, ratioR1toX0, ratioR2toX0, ratioR3toX0, ratioR12toX0, ratioR23toX0, ratioR13toX0, Sa_ratio_S_D1tog0, \
        Sa_ratio_S_D2toS_D1, Sa_ratio_S_D3toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, \
        Sa_ratio_R3_D3toS_D3, T_StoR1, T_StoR2, T_StoR3


# Generate default parameter criterisa index based on drug number.
def DPM_generate_PAR_criterisa_list(Num_drug):
    if Num_drug == 2:
        PAR_criterisa_list = [0, 1, 2, 3, 4, 5, 6]
    elif Num_drug == 3:
        PAR_criterisa_list = [0, 7, 8]
    else:
        PAR_criterisa_list = None
    return PAR_criterisa_list


# Generate heading for .csv saving files.
def DPM_generate_heading_csv(Num_drug, Simduration, Stepsize):
    # Heading of param csv file.
    Heading_param_csv: list[str] = HEADING_2DRUG_PARAM_CSV if Num_drug == 2 else HEADING_3DRUG_PARAM_CSV if Num_drug == 3 else None
    # Heading of stopt csv file.
    Heading_stopt_csv = HEADING_STOPT_CSV
    # Heading of dosage csv file.
    Heading_dosage_csv = HEADING_DOSAGE_CSV
    Heading_dosage_csv_str = 'Drug1 dosage,Drug2 dosage' if Num_drug == 2 else 'Drug1 dosage,Drug2 dosage,Drug3 dosage' if Num_drug == 3 else None
    Heading_dosage_csv.extend([f'({Heading_dosage_csv_str}) at t={t_i}' for t_i in np.arange(0, Simduration, Stepsize)])
    # Heading of pop csv file.
    Heading_pop_csv = HEADING_POP_CSV
    Heading_pop_csv_str = ','.join(ALL_POSSIBLE_CELLTYPE_2DRUG) if Num_drug == 2 else ','.join(ALL_POSSIBLE_CELLTYPE_3DRUG) if Num_drug == 3 \
        else None
    Heading_pop_csv.extend(f'({Heading_pop_csv_str}) at t={t_i}' for t_i in np.arange(Stepsize, Simduration + Stepsize, Stepsize))
    # Heading of eachtimepoint csv file.
    Heading_eachtimepoint_csv = HEADING_EACHTIMEPOINT_CSV
    Heading_eachtimepoint_csv_str = [f'({Heading_pop_csv_str}) at t={t_i}, ({Heading_dosage_csv_str}) at t={t_i}'.split(', ')
                                     for t_i in np.arange(0, Simduration + SIMTIMESTEP_DEFAULT_VAL, SIMTIMESTEP_DEFAULT_VAL)]
    Heading_eachtimepoint_csv.extend(list(itertools.chain.from_iterable(Heading_eachtimepoint_csv_str)))
    return Heading_param_csv, Heading_stopt_csv, Heading_dosage_csv, Heading_pop_csv, Heading_eachtimepoint_csv


# Generate truncated decision tree. Used in strategy 6-9.
def DPM_generate_truncate_decisiontree(decision_tree, tree_terminate, index):
    # Delete the trees in the variable "decision_tree" which have the same value from index 0 to the same depth in the variable
    # "tree_terminate".
    truncate_boolen = [True if x[:len(tree_terminate)] == tree_terminate else False for x in decision_tree]
    # Find the index of True
    truncate_index = [i for i, x in enumerate(truncate_boolen) if x]
    # Find the index of False
    keep_index = [i for i, x in enumerate(truncate_boolen) if not x]
    decision_tree = [x for i, x in enumerate(decision_tree) if i in keep_index]
    if (np.array(truncate_index, dtype=int) < index).any():
        '''This should not happen becasue the tuple in decision_tree is sorted element by element from index 0 to the end. If it happens, it will 
        delete the trees whose index are smaller than the current index number. The shorten of the decision_tree will make the length of the 
        decision_tree smaller than the current index. The trees whose index are smaller than the current index will move forward and their new index 
        are smaller than the current index. The simulation of these trees are missed and the result will be wrong. But because the tuple in the 
        decision_tree is orderd, e.g., 5 steps look forward, 3 drug combination: (0,0,0,0,0)->(0,0,0,0,1)->(0,0,0,0,2)->(0,0,0,1,0)->(0,0,0,1,1)->
        (0,0,0,1,2)->(0,0,0,2,0)->(0,0,0,2,1)->(0,0,0,2,2)->(0,0,1,0,0)->(0,0,1,0,0)->(0,0,1,0,1)->(0,0,1,0,2)->...->...->(2,2,2,2,2). If the 
        simulation result shows we should exclude the nodes begin with (0,0,1,_,_), it would happens at the first time simulating the node begin with 
        (0,0,1,_,_). Because when the first time the nodes begin with (0,0,1,_,_) show up, we either keep it or exclude it. If we exclude it, then the
        nodes begin with (0,0,1,_,_) will not appear again. If we keep it, because the nodes are ordered, we begin to explore all the nodes begin with
        (0,0,1,_,_). We will decide to whether keep the nodes begin with (0,0,1,_,_) based on the simulation results of elements x and y in the tuple 
        (0,0,1,x,y) but they begin with (0,0,1,_,_). After we explored all the nodes begin with (0,0,1,_,_), there are no nodes begin with (0,0,1,_,_) 
        left. I think that it is safe because the nodes are ordered, we will explore all the nodes begin with (0,0,1,_,_) after it first happened. 
        If we didn't explore all the nodes begin with (0,0,1,_,_), instead we explore (0,0,2,_,_) and come back to (0,0,1,_,_). If (0,0,2,_,_) has a 
        smaller values than all the nodes begin with (0,0,1,_,_), we will exclude the nodes begin with (0,0,1,_,_) which didn't excluded before. 
        That's the situation we don't want. I add a check here just make sure this situation never happen.'''
        print('Truncate tree error. delete the trees whose index are smaller than the current index number.')
        DPM_print_errorandclose()
        sys.exit()
    return decision_tree


# Generate initial cell number X0 based on input parameter.
def DPM_generate_X0(PAR):
    return DPM_generate_X0_2drug(PAR) if PAR['Num_drug'] == 2 else DPM_generate_X0_3drug(PAR) if PAR['Num_drug'] == 3 else None


# Generate initial cell number X0 based on input parameter, 2 drug case.
def DPM_generate_X0_2drug(PAR):
    return np.array([PAR['Spop'], PAR['R1pop'], PAR['R2pop'], PAR['R12pop']], dtype=float)


# Generate initial cell number X0 based on input parameter, 3 drug case.
def DPM_generate_X0_3drug(PAR):
    return np.array([PAR['Spop'], PAR['R1pop'], PAR['R2pop'], PAR['R3pop'], PAR['R12pop'], PAR['R13pop'], PAR['R23pop'], PAR['R123pop']], dtype=float)


# Generate g0 based on input parameter.
def DPM_generate_g0(PAR):
    return DPM_generate_g0_2drug(PAR) if PAR['Num_drug'] == 2 else DPM_generate_g0_3drug(PAR) if PAR['Num_drug'] == 3 else None


# Generate g0 based on input parameter, 2 drug case.
def DPM_generate_g0_2drug(PAR):
    return np.array([PAR['g0_S'], PAR['g0_R1'], PAR['g0_R2'], PAR['g0_R12']], dtype=float)


# Generate g0 based on input parameter, 3 drug case.
def DPM_generate_g0_3drug(PAR):
    return np.array([PAR['g0_S'], PAR['g0_R1'], PAR['g0_R2'], PAR['g0_R3'], PAR['g0_R12'], PAR['g0_R13'], PAR['g0_R23'], PAR['g0_R123']], dtype=float)


# Generate Sa based on input parameter.
def DPM_generate_Sa(PAR):
    return DPM_generate_Sa_2drug(PAR) if PAR['Num_drug'] == 2 else DPM_generate_Sa_3drug(PAR) if PAR['Num_drug'] == 3 else None


# Generate Sa based on input parameter, 2 drug case.
def DPM_generate_Sa_2drug(PAR):
    Sa = np.zeros((PAR['Num_cell_type'], PAR['Num_drug']), dtype=float)

    # Sensitivity of S cell on drug 1 and drug 2.
    Sa[0, :] = np.array([PAR['Sa.S.D1.'], PAR['Sa.S.D2.']])
    # Sensitivity of R1 cell on drug 1 and drug 2.
    Sa[1, :] = np.array([PAR['Sa.R1.D1.'], PAR['Sa.R1.D2.']])
    # Sensitivity of R2 cell on drug 1 and drug 2.
    Sa[2, :] = np.array([PAR['Sa.R2.D1.'], PAR['Sa.R2.D2.']])
    # Sensitivity of R12 cell on drug 1 and drug 2.
    Sa[3, :] = np.array([PAR['Sa.R12.D1.'], PAR['Sa.R12.D2.']])

    return Sa


# Generate Sa based on input parameter, 3 drug case.
def DPM_generate_Sa_3drug(PAR):
    Sa = np.zeros((PAR['Num_cell_type'], PAR['Num_drug']), dtype=float)

    # Sensitivity of S cell on drug 1, drug 2 and drug 3.
    Sa[0, :] = np.array([PAR['Sa.S.D1.'], PAR['Sa.S.D2.'], PAR['Sa.S.D3.']])
    # Sensitivity of R1 cell on drug 1, drug 2 and drug 3.
    Sa[1, :] = np.array([PAR['Sa.R1.D1.'], PAR['Sa.R1.D2.'], PAR['Sa.R1.D3.']])
    # Sensitivity of R2 cell on drug 1, drug 2 and drug 3.
    Sa[2, :] = np.array([PAR['Sa.R2.D1.'], PAR['Sa.R2.D2.'], PAR['Sa.R2.D3.']])
    # Sensitivity of R3 cell on drug 1, drug 2 and drug 3.
    Sa[3, :] = np.array([PAR['Sa.R3.D1.'], PAR['Sa.R3.D2.'], PAR['Sa.R3.D3.']])
    # Sensitivity of R12 cell on drug 1, drug 2 and drug 3.
    Sa[4, :] = np.array([PAR['Sa.R12.D1.'], PAR['Sa.R12.D2.'], PAR['Sa.R12.D3.']])
    # Sensitivity of R13 cell on drug 1, drug 2 and drug 3.
    Sa[5, :] = np.array([PAR['Sa.R13.D1.'], PAR['Sa.R13.D2.'], PAR['Sa.R13.D3.']])
    # Sensitivity of R23 cell on drug 1, drug 2 and drug 3.
    Sa[6, :] = np.array([PAR['Sa.R23.D1.'], PAR['Sa.R23.D2.'], PAR['Sa.R23.D3.']])
    # Sensitivity of R123 cell on drug 1, drug 2 and drug 3.
    Sa[7, :] = np.array([PAR['Sa.R123.D1.'], PAR['Sa.R123.D2.'], PAR['Sa.R123.D3.']])

    return Sa


# Generate T based on input parameter.
def DPM_generate_T(PAR):
    T = DPM_generate_T_2drug(PAR) if PAR['Num_drug'] == 2 else DPM_generate_T_3drug(PAR) if PAR['Num_drug'] == 3 else None
    return T


# Generate T based on input parameter, 2 drug case.
def DPM_generate_T_2drug(PAR):
    T = np.zeros((PAR['Num_cell_type'], PAR['Num_cell_type']), dtype=float)
    # Transition to S.
    T[0, :] = np.array([PAR['T.S..S.'], PAR['T.S..R1.'], PAR['T.S..R2.'], PAR['T.S..R12.']])
    # Transition to R1.
    T[1, :] = np.array([PAR['T.R1..S.'], PAR['T.R1..R1.'], PAR['T.R1..R2.'], PAR['T.R1..R12.']])
    # Transition to R2.
    T[2, :] = np.array([PAR['T.R2..S.'], PAR['T.R2..R1.'], PAR['T.R2..R2.'], PAR['T.R2..R12.']])
    # Transition to R12.
    T[3, :] = np.array([PAR['T.R12..S.'], PAR['T.R12..R1.'], PAR['T.R12..R2.'], PAR['T.R12..R12.']])
    return T


# Generate T based on input parameter, 3 drug case.
def DPM_generate_T_3drug(PAR):
    T = np.zeros((PAR['Num_cell_type'], PAR['Num_cell_type']), dtype=float)
    keys = [i for i in PAR.keys() if 'T.' in i]
    # Transition to S.
    T[0, :] = np.array([PAR['T.S..S.'], PAR['T.S..R1.'], PAR['T.S..R2.'], PAR['T.S..R3.'], PAR['T.S..R12.'], PAR['T.S..R13.'],
                        PAR['T.S..R23.'], PAR['T.S..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.S..' in i_key]
    # Transition to R1.
    T[1, :] = np.array([PAR['T.R1..S.'], PAR['T.R1..R1.'], PAR['T.R1..R2.'], PAR['T.R1..R3.'], PAR['T.R1..R12.'], PAR['T.R1..R13.'],
                        PAR['T.R1..R23.'], PAR['T.R1..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R1..' in i_key]
    # Transition to R2.
    T[2, :] = np.array([PAR['T.R2..S.'], PAR['T.R2..R1.'], PAR['T.R2..R2.'], PAR['T.R2..R3.'], PAR['T.R2..R12.'], PAR['T.R2..R13.'],
                        PAR['T.R2..R23.'], PAR['T.R2..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R2..' in i_key]
    # Transition to R3.
    T[3, :] = np.array([PAR['T.R3..S.'], PAR['T.R3..R1.'], PAR['T.R3..R2.'], PAR['T.R3..R3.'], PAR['T.R3..R12.'], PAR['T.R3..R13.'],
                        PAR['T.R3..R23.'], PAR['T.R3..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R3..' in i_key]
    # Transition to R12.
    T[4, :] = np.array([PAR['T.R12..S.'], PAR['T.R12..R1.'], PAR['T.R12..R2.'], PAR['T.R12..R3.'], PAR['T.R12..R12.'], PAR['T.R12..R13.'],
                        PAR['T.R12..R23.'], PAR['T.R12..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R12..' in i_key]
    # Transition to R13.
    T[5, :] = np.array([PAR['T.R13..S.'], PAR['T.R13..R1.'], PAR['T.R13..R2.'], PAR['T.R13..R3.'], PAR['T.R13..R12.'], PAR['T.R13..R13.'],
                        PAR['T.R13..R23.'], PAR['T.R13..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R13..' in i_key]
    # Transition to R23.
    T[6, :] = np.array([PAR['T.R23..S.'], PAR['T.R23..R1.'], PAR['T.R23..R2.'], PAR['T.R23..R3.'], PAR['T.R23..R12.'], PAR['T.R23..R13.'],
                        PAR['T.R23..R23.'], PAR['T.R23..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R23..' in i_key]
    # Transition to R123.
    T[7, :] = np.array([PAR['T.R123..S.'], PAR['T.R123..R1.'], PAR['T.R123..R2.'], PAR['T.R123..R3.'], PAR['T.R123..R12.'], PAR['T.R123..R13.'],
                        PAR['T.R123..R23.'], PAR['T.R123..R123.']])

    return T


def DPM_generate_pdf(x, k, LOD):
    return (k/x**2)/(1-k/LOD+k)


def DPM_generate_misspecification_subcolone(par, subcolone_LOD, misspecification_LOD, mutation_rate, celltype, mis_specfiy_pop, Strategy_name,
                                            norm_pdf=1, Xtotal=None):
    # def p(x, k, LOD):
    #     return (k/x**2)/(1-k/LOD+k)

    class MutationFraction(rv_continuous):
        def _pdf(self, x, k, LOD, const):
            return (1.0/const)*DPM_generate_pdf(x, k, LOD)

    # x = np.linspace(mutation_rate, subcolone_LOD, int(1e7))
    # norm_pdf = simps(DPM_generate_pdf(x, mutation_rate, subcolone_LOD), x)
    mutationFraction_distribution = MutationFraction(name='mutationFraction_distribution', a=mutation_rate, b=subcolone_LOD)
    # pdf = mutationFraction_distribution.pdf(k=mutation_rate, LOD=subcolone_LOD, x=x, const=norm_pdf)
    # cdf = mutationFraction_distribution.cdf(k=mutation_rate, LOD=subcolone_LOD, x=x, const=norm_pdf)
    # plt.plot(x, cdf)
    # plt.xscale("log")
    # plt.yscale("log")

    for i, i_par in enumerate(par):
        # print(i)
        i_par_ori = deepcopy(i_par)
        flag_mis = False
        if Xtotal is None:
            i_X = i_par['Spop'] + i_par['R1pop'] + i_par['R2pop'] + i_par['R12pop']
        else:
            i_X = Xtotal
        for i_sub in celltype:
            if i_par[i_sub]/i_X < 2*subcolone_LOD:
                flag_mis = True
                mis_specfiy_pop[Strategy_name][i_sub]['total'].append(i_X)
                mis_specfiy_pop[Strategy_name][i_sub]['percent'].append(i_par[i_sub]/i_X)
                if misspecification_LOD == 0:
                    i_par[i_sub] = 0.0
                    mis_specfiy_pop[Strategy_name][i_sub]['est'].append(i_par[i_sub]/i_X)
                elif misspecification_LOD == 'pdf':
                    if i_sub in ['R1pop', 'R2pop']:
                        val = i_X * 2 * mutationFraction_distribution.rvs(k=mutation_rate, LOD=subcolone_LOD, const=norm_pdf, size=1)
                        i_par[i_sub] = np.ceil(val).item() if val >= 1 else 0
                    else:
                        i_par[i_sub] = 0.0
                    mis_specfiy_pop[Strategy_name][i_sub]['est'].append(i_par[i_sub]/i_X)
                elif misspecification_LOD == 'max':
                    if i_sub in ['R1pop', 'R2pop']:
                        val = i_X * (2 * subcolone_LOD - 1e-9)
                        i_par[i_sub] = np.ceil(val).item() if val >= 1 else 0
                    else:
                        i_par[i_sub] = 0.0
                    mis_specfiy_pop[Strategy_name][i_sub]['est'].append(i_par[i_sub] / i_X)
                elif misspecification_LOD == 'loguni':
                    if i_sub in ['R1pop', 'R2pop']:
                        val = i_X * 2 * loguniform.rvs(mutation_rate, subcolone_LOD, size=1)
                        i_par[i_sub] = np.ceil(val).item() if val >= 1 else 0
                    else:
                        i_par[i_sub] = 0.0
                    mis_specfiy_pop[Strategy_name][i_sub]['est'].append(i_par[i_sub] / i_X)
        if misspecification_LOD == 'pdf':
            if i_par['R1pop'] > 0:
                assert i_par['R1pop'] >= 2 * mutation_rate * i_X
            if i_par['R2pop'] > 0:
                assert i_par['R2pop'] >= 2 * mutation_rate * i_X
        if flag_mis:
            i_par['Spop'] = i_X - i_par['R1pop'] - i_par['R2pop'] - i_par['R12pop']
            assert math.isclose(i_par['Spop'] + i_par['R1pop'] + i_par['R2pop'] + i_par['R12pop'],
                                i_par_ori['Spop'] + i_par_ori['R1pop'] + i_par_ori['R2pop'] + i_par_ori['R12pop'], rel_tol=1e-6)
        else:
            assert i_par == i_par_ori
            # for i_key, _ in i_par.items():
            #     assert math.isclose(i_par[i_key], i_par_ori[i_key], rel_tol=1e-5)
        par[i] = i_par
    return par, mis_specfiy_pop
