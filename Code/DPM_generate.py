from DPM_lib import itertools
from DPM_constant import *
'''This script checks arguments used in the DPM (Dynamic Precision Medicine) model.'''
'''
CHECKED
'''


# Generate the default discrete dose combination based on the number of drugs.#
def DPM_generate_discrete_dose_combination(num_drug):
    dose_list = [float(0)]
    for i in range(int(num_drug)):
        dose_list.append(float(1 / (i + 1)))
    dose_list.sort()
    dose_combination = itertools.product(dose_list, repeat=int(num_drug))
    dose_combination = [x for x in dose_combination if sum(x) == 1]
    return dose_combination


# Generate the continuous dose combination based on the number of drugs number and the dose interval.#
def DPM_generate_continuous_dose_combination(num_drug, dose_interval):
    dose_list = np.arange(0, 1+dose_interval, dose_interval, dtype=float)
    dose_list = np.sort(dose_list)
    dose_combination = itertools.product(dose_list.tolist(), repeat=int(num_drug))
    dose_combination = [x for x in dose_combination if sum(x) == 1]
    return dose_combination


# Generate default parameters for the two-drug case following PNAS (2012).#
def DPM_generate_default_par_2drug():
    # g0: basel proliferation rate for the 4 cell types: S, R1, R2, R12.#
    g0 = [0.001, 0.0026, 0.007, 0.0184, 0.0487, 0.1287, 0.34]
    # Ratio of R1 to X0.#
    ratioR1toX0 = [0, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 9e-1]
    # Ratio of R2 to X0.#
    ratioR2toX0 = [0, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 9e-1]
    # Sa:4×2 matrix of drug sensitivities.#
    # Sa(S,D1)/g0.#
    Sa_ratio_S_D1tog0 = [5.6e-4, 0.0054, 0.0517, 0.4964, 4.7683, 45.8045, 440]
    # Sa(S,D2)/Sa(S,D1).#
    Sa_ratio_S_D2toS_D1 = [4e-4, 0.0015, 0.0054, 0.02, 0.0737, 0.2714, 1e0]
    # Sa(R1,D1)/Sa(S,D1).#
    Sa_ratio_R1_D1toS_D1 = [0, 1e-5, 9.5635e-5, 9.1461e-4, 0.0087, 0.0837, 0.8]
    # Sa(R2,D2)/Sa(S,D2).#
    Sa_ratio_R2_D2toS_D2 = [0, 1e-5, 9.5635e-5, 9.1461e-4, 0.0087, 0.0837, 0.8]
    # T:4×4 transition rate matrix, R1->R12=S->R2, and R2->R12=S->R1.#
    T_StoR1 = [1e-11, 2.154e-10, 4.642e-9, 1e-7, 2.154e-6, 4.642e-5, 1e-3]
    T_StoR2 = [1e-11, 2.154e-10, 4.642e-9, 1e-7, 2.154e-6, 4.642e-5, 1e-3]
    return g0, ratioR1toX0, ratioR2toX0, Sa_ratio_S_D1tog0, Sa_ratio_S_D2toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, T_StoR1, T_StoR2


# Generate default parameters for the 3 drug case based on Biology Direct(2016).#
def DPM_generate_default_par_3drug():
    # g0: baseline proliferation rate for the eight cell types: S, R1, R2, R3, R12, R13, R23, R123.#
    g0 = [0.001, 0.0026, 0.007, 0.0184, 0.0487, 0.1287, 0.34]
    # Ratio of S to X0.#
    ratioStoX0 = [1e-1, 9e-1]
    # Ratio of R1 to X0.#
    ratioR1toX0 = [0e0, 1e-5, 1e-3, 1e-1, 9e-1]
    # Ratio of R2 to X0.#
    ratioR2toX0 = [0e0, 1e-5, 1e-3, 1e-1, 9e-1]
    # Ratio of R3 to X0.#
    ratioR3toX0 = [0e0, 1e-5, 1e-3, 1e-1, 9e-1]
    # Ratio of R12 to X0.#
    ratioR12toX0 = [0e0, 1e-5, 1e-3]
    # Ratio of R23 to X0 ratio.#
    ratioR23toX0 = [0e0, 1e-5, 1e-3]
    # Ratio of R13 to X0 ratio.#
    ratioR13toX0 = [0e0, 1e-5, 1e-3]
    # Sa(S,D1)/g0.#
    Sa_ratio_S_D1tog0 = [5e-1, 5e0]
    # Sa(S,D2)/Sa(S,D1).#
    Sa_ratio_S_D2toS_D1 = [1e-1, 0.3, 1e0]
    # Sa(S,D3)/Sa(S,D1).#
    Sa_ratio_S_D3toS_D1 = [1e-1, 0.3, 1e0]
    # Sa(R1,D1)/Sa(S,D1).#
    Sa_ratio_R1_D1toS_D1 = [0e0, 3e-1, 1e0]
    # Sa(R2,D2)/Sa(S,D2).#
    Sa_ratio_R2_D2toS_D2 = [0e0, 3e-1, 1e0]
    # Sa(R3,D3)/Sa(S,D3).#
    Sa_ratio_R3_D3toS_D3 = [0e0, 3e-1, 1e0]
    # T:8×8 transition rate matrix.#
    T_StoR1 = [1e-11, 1e-9, 1e-7, 1e-5, 1e-3]
    T_StoR2 = [1e-11, 1e-9, 1e-7, 1e-5, 1e-3]
    T_StoR3 = [1e-11, 1e-9, 1e-7, 1e-5, 1e-3]
    return g0, ratioStoX0, ratioR1toX0, ratioR2toX0, ratioR3toX0, ratioR12toX0, ratioR23toX0, ratioR13toX0, Sa_ratio_S_D1tog0, \
        Sa_ratio_S_D2toS_D1, Sa_ratio_S_D3toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, \
        Sa_ratio_R3_D3toS_D3, T_StoR1, T_StoR2, T_StoR3


# Generate the default parameter criteria index based on the number of drugs.#
def DPM_generate_PAR_criterisa_list(Num_drug):
    if Num_drug == 2:
        PAR_criterisa_list = [0, 1, 2, 3, 4, 5, 6]
    elif Num_drug == 3:
        PAR_criterisa_list = [0, 7, 8]
    else:
        PAR_criterisa_list = None
    return PAR_criterisa_list


# Generate the header for .csv output files.#
def DPM_generate_header_csv(Num_drug, Simduration, Stepsize):
    # Header of the parameter csv file.#
    Header_param_csv = HEADER_2DRUG_PARAM_CSV if Num_drug == 2 else HEADER_3DRUG_PARAM_CSV if Num_drug == 3 else None
    # Header of the stopt csv file.#
    Header_stopt_csv = HEADER_STOPT_CSV
    # Header of the dosage csv file.#
    Header_dosage_csv = HEADER_DOSAGE_CSV
    Header_dosage_csv_str = 'Drug1 dosage,Drug2 dosage' if Num_drug == 2 else 'Drug1 dosage,Drug2 dosage,Drug3 dosage' if Num_drug == 3 else None
    Header_dosage_csv.extend([f'({Header_dosage_csv_str}) at t={t_i}' for t_i in np.arange(0, Simduration, Stepsize).tolist()])
    # Header of the pop csv file.#
    Header_pop_csv = HEADER_POP_CSV
    Header_pop_csv_str = ','.join(ALL_POSSIBLE_CELLTYPE_2DRUG) if Num_drug == 2 else ','.join(ALL_POSSIBLE_CELLTYPE_3DRUG) if Num_drug == 3 \
        else None
    Header_pop_csv.extend(f'({Header_pop_csv_str}) at t={t_i}' for t_i in np.arange(Stepsize, Simduration + Stepsize, Stepsize).tolist())
    # Header of the eachtimepoint csv file.#
    Header_eachtimepoint_csv = HEADER_EACHTIMEPOINT_CSV
    Header_eachtimepoint_csv_str = [f'({Header_pop_csv_str}) at t={t_i}, ({Header_dosage_csv_str}) at t={t_i}'.split(', ')
                                    for t_i in np.arange(0, Simduration + SIMTIMESTEP_DEFAULT_VAL, SIMTIMESTEP_DEFAULT_VAL).tolist()]
    Header_eachtimepoint_csv.extend(list(itertools.chain.from_iterable(Header_eachtimepoint_csv_str)))
    return Header_param_csv, Header_stopt_csv, Header_dosage_csv, Header_pop_csv, Header_eachtimepoint_csv


# Generate initial cell number X0 based on the input parameter.#
def DPM_generate_X0(PAR):
    return DPM_generate_X0_2drug(PAR) if PAR['Num_drug'] == 2 else DPM_generate_X0_3drug(PAR) if PAR['Num_drug'] == 3 else None


# Generate initial cell number X0 based on the input parameter, 2 drug case.#
def DPM_generate_X0_2drug(PAR):
    return np.array([PAR['Spop'], PAR['R1pop'], PAR['R2pop'], PAR['R12pop']], dtype=float)


# Generate initial cell number X0 based on the input parameter, 3 drug case.#
def DPM_generate_X0_3drug(PAR):
    return np.array([PAR['Spop'], PAR['R1pop'], PAR['R2pop'], PAR['R3pop'], PAR['R12pop'], PAR['R13pop'], PAR['R23pop'], PAR['R123pop']],
                    dtype=float)


# Generate g0 based on the input parameter.#
def DPM_generate_g0(PAR):
    return DPM_generate_g0_2drug(PAR) if PAR['Num_drug'] == 2 else DPM_generate_g0_3drug(PAR) if PAR['Num_drug'] == 3 else None


# Generate g0 based on the input parameter, 2 drug case.#
def DPM_generate_g0_2drug(PAR):
    return np.array([PAR['g0_S'], PAR['g0_R1'], PAR['g0_R2'], PAR['g0_R12']], dtype=float)


# Generate g0 based on the input parameter, 3 drug case.#
def DPM_generate_g0_3drug(PAR):
    return np.array([PAR['g0_S'], PAR['g0_R1'], PAR['g0_R2'], PAR['g0_R3'], PAR['g0_R12'], PAR['g0_R13'], PAR['g0_R23'], PAR['g0_R123']],
                    dtype=float)


# Generate Sa based on the input parameter.#
def DPM_generate_Sa(PAR):
    return DPM_generate_Sa_2drug(PAR) if PAR['Num_drug'] == 2 else DPM_generate_Sa_3drug(PAR) if PAR['Num_drug'] == 3 else None


# Generate Sa based on the input parameter, 2 drug case.#
def DPM_generate_Sa_2drug(PAR):
    Sa = np.zeros((PAR['Num_cell_type'], PAR['Num_drug']), dtype=float)

    # Sensitivity of S cell on drug 1 and drug 2.#
    Sa[0, :] = np.array([PAR['Sa.S.D1.'], PAR['Sa.S.D2.']])
    # Sensitivity of R1 cell on drug 1 and drug 2.#
    Sa[1, :] = np.array([PAR['Sa.R1.D1.'], PAR['Sa.R1.D2.']])
    # Sensitivity of R2 cell on drug 1 and drug 2.#
    Sa[2, :] = np.array([PAR['Sa.R2.D1.'], PAR['Sa.R2.D2.']])
    # Sensitivity of R12 cell on drug 1 and drug 2.#
    Sa[3, :] = np.array([PAR['Sa.R12.D1.'], PAR['Sa.R12.D2.']])
    return Sa


# Generate Sa based on the input parameter, 3 drug case.#
def DPM_generate_Sa_3drug(PAR):
    Sa = np.zeros((PAR['Num_cell_type'], PAR['Num_drug']), dtype=float)

    # Sensitivity of S cell on drug 1, drug 2 and drug 3.#
    Sa[0, :] = np.array([PAR['Sa.S.D1.'], PAR['Sa.S.D2.'], PAR['Sa.S.D3.']])
    # Sensitivity of R1 cell on drug 1, drug 2 and drug 3.#
    Sa[1, :] = np.array([PAR['Sa.R1.D1.'], PAR['Sa.R1.D2.'], PAR['Sa.R1.D3.']])
    # Sensitivity of R2 cell on drug 1, drug 2 and drug 3.#
    Sa[2, :] = np.array([PAR['Sa.R2.D1.'], PAR['Sa.R2.D2.'], PAR['Sa.R2.D3.']])
    # Sensitivity of R3 cell on drug 1, drug 2 and drug 3.#
    Sa[3, :] = np.array([PAR['Sa.R3.D1.'], PAR['Sa.R3.D2.'], PAR['Sa.R3.D3.']])
    # Sensitivity of R12 cell on drug 1, drug 2 and drug 3.#
    Sa[4, :] = np.array([PAR['Sa.R12.D1.'], PAR['Sa.R12.D2.'], PAR['Sa.R12.D3.']])
    # Sensitivity of R13 cell on drug 1, drug 2 and drug 3.#
    Sa[5, :] = np.array([PAR['Sa.R13.D1.'], PAR['Sa.R13.D2.'], PAR['Sa.R13.D3.']])
    # Sensitivity of R23 cell on drug 1, drug 2 and drug 3.#
    Sa[6, :] = np.array([PAR['Sa.R23.D1.'], PAR['Sa.R23.D2.'], PAR['Sa.R23.D3.']])
    # Sensitivity of R123 cell on drug 1, drug 2 and drug 3.#
    Sa[7, :] = np.array([PAR['Sa.R123.D1.'], PAR['Sa.R123.D2.'], PAR['Sa.R123.D3.']])
    return Sa


# Generate T based on the input parameter.#
def DPM_generate_T(PAR):
    T = DPM_generate_T_2drug(PAR) if PAR['Num_drug'] == 2 else DPM_generate_T_3drug(PAR) if PAR['Num_drug'] == 3 else None
    return T


# Generate T based on the input parameter, 2 drug case.#
def DPM_generate_T_2drug(PAR):
    T = np.zeros((PAR['Num_cell_type'], PAR['Num_cell_type']), dtype=float)
    # Transition to S.#
    T[0, :] = np.array([PAR['T.S..S.'], PAR['T.S..R1.'], PAR['T.S..R2.'], PAR['T.S..R12.']])
    # Transition to R1.#
    T[1, :] = np.array([PAR['T.R1..S.'], PAR['T.R1..R1.'], PAR['T.R1..R2.'], PAR['T.R1..R12.']])
    # Transition to R2.#
    T[2, :] = np.array([PAR['T.R2..S.'], PAR['T.R2..R1.'], PAR['T.R2..R2.'], PAR['T.R2..R12.']])
    # Transition to R12.#
    T[3, :] = np.array([PAR['T.R12..S.'], PAR['T.R12..R1.'], PAR['T.R12..R2.'], PAR['T.R12..R12.']])
    return T


# Generate T based on the input parameter, 3 drug case.#
def DPM_generate_T_3drug(PAR):
    T = np.zeros((PAR['Num_cell_type'], PAR['Num_cell_type']), dtype=float)
    keys = [i for i in PAR.keys() if 'T.' in i]
    # Transition to S.#
    T[0, :] = np.array([PAR['T.S..S.'], PAR['T.S..R1.'], PAR['T.S..R2.'], PAR['T.S..R3.'], PAR['T.S..R12.'], PAR['T.S..R13.'],
                        PAR['T.S..R23.'], PAR['T.S..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.S..' in i_key]
    # Transition to R1.#
    T[1, :] = np.array([PAR['T.R1..S.'], PAR['T.R1..R1.'], PAR['T.R1..R2.'], PAR['T.R1..R3.'], PAR['T.R1..R12.'], PAR['T.R1..R13.'],
                        PAR['T.R1..R23.'], PAR['T.R1..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R1..' in i_key]
    # Transition to R2.#
    T[2, :] = np.array([PAR['T.R2..S.'], PAR['T.R2..R1.'], PAR['T.R2..R2.'], PAR['T.R2..R3.'], PAR['T.R2..R12.'], PAR['T.R2..R13.'],
                        PAR['T.R2..R23.'], PAR['T.R2..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R2..' in i_key]
    # Transition to R3.#
    T[3, :] = np.array([PAR['T.R3..S.'], PAR['T.R3..R1.'], PAR['T.R3..R2.'], PAR['T.R3..R3.'], PAR['T.R3..R12.'], PAR['T.R3..R13.'],
                        PAR['T.R3..R23.'], PAR['T.R3..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R3..' in i_key]
    # Transition to R12.#
    T[4, :] = np.array([PAR['T.R12..S.'], PAR['T.R12..R1.'], PAR['T.R12..R2.'], PAR['T.R12..R3.'], PAR['T.R12..R12.'], PAR['T.R12..R13.'],
                        PAR['T.R12..R23.'], PAR['T.R12..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R12..' in i_key]
    # Transition to R13.#
    T[5, :] = np.array([PAR['T.R13..S.'], PAR['T.R13..R1.'], PAR['T.R13..R2.'], PAR['T.R13..R3.'], PAR['T.R13..R12.'], PAR['T.R13..R13.'],
                        PAR['T.R13..R23.'], PAR['T.R13..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R13..' in i_key]
    # Transition to R23.#
    T[6, :] = np.array([PAR['T.R23..S.'], PAR['T.R23..R1.'], PAR['T.R23..R2.'], PAR['T.R23..R3.'], PAR['T.R23..R12.'], PAR['T.R23..R13.'],
                        PAR['T.R23..R23.'], PAR['T.R23..R123.']])
    _ = [PAR.pop(i_key) for i_key in keys if 'T.R23..' in i_key]
    # Transition to R123.#
    T[7, :] = np.array([PAR['T.R123..S.'], PAR['T.R123..R1.'], PAR['T.R123..R2.'], PAR['T.R123..R3.'], PAR['T.R123..R12.'], PAR['T.R123..R13.'],
                        PAR['T.R123..R23.'], PAR['T.R123..R123.']])

    return T
