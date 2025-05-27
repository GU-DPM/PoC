from random import randrange
from DPM_assign_check import *
from DPM_generate import DPM_generate_truncate_decisiontree
import matplotlib.pyplot as plt


def DPM_strategy_pro(X0, g0, Sa, T, Num_drug, tnow, X_Strategy_sim, X_Strategy_true, t_Strategy_sim, d_Strategy_sim, Stepsize, maxthreshold,
                     Simtimestep, Limit_mortality, Limit_radiologicdetection,  subcolone_LOD, misspecification_LOD, mutation_rate, celltype,
                     mis_specfiy_pop, Strategy_name):
    flag_treat, flag_break = False, False
    # no treat, let proliferate
    treat_i = np.zeros([Num_drug, 1], dtype=int)
    t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
    d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))

    X_true_i, t_true_i, _ = DPM_sim_model(X0, T, g0, Sa, d_i, t_i, maxthreshold, True)
    X_true_total_i = np.sum(X_true_i, axis=0)

    # True dead, stop sim.
    if np.any(X_true_total_i >= Limit_mortality):
        X_sim, flag_break = X_true_i, True
    else:
        X_sim = np.tile(X_Strategy_sim[:, [-1]], (1, X_true_i.shape[1]))

    # True cell number of all types are smaller than 1 (cured), should not happen, error.
    assert not np.all(X_true_i[:, -1] < 1)

    X_Strategy = np.append(X_Strategy_sim, X_sim, 1) if tnow == 0 else np.append(X_Strategy_sim, X_sim[:, 1:], 1)
    t_Strategy = np.append(t_Strategy_sim, t_i + tnow) if tnow == 0 else np.append(t_Strategy_sim, t_true_i[1:] + tnow)
    d_Strategy = np.append(d_Strategy_sim, d_i, 1)

    X_Strategy_true = np.append(X_Strategy_true, X_true_i, 1) if tnow == 0 else np.append(X_Strategy_true, X_true_i[:, 1:], 1)

    # true total cell population reemerges from radiologic detection, do subcolone detection and treat again
    if X_true_total_i[-1] >= Limit_radiologicdetection:
        flag_treat = True
        X0, mis_specfiy_pop = \
            DPM_strategy_subcolone_detection(X_true_i[:, -1], subcolone_LOD, Limit_radiologicdetection, misspecification_LOD, mutation_rate,
                                             celltype, mis_specfiy_pop, Strategy_name)
        X_Strategy[:, -1] = X0
    else:
        X0 = X_true_i[:, -1]
    return X0, X_Strategy, X_Strategy_true, t_Strategy, d_Strategy, flag_treat,  flag_break, mis_specfiy_pop


def DPM_strategy_LOD(X_sim, t_sim, d_sim, X_true, t_true, tnow, X_Strategy_sim, X_Strategy_true, t_Strategy_sim, d_Strategy_sim,
                     Stepsize, Limit_mortality, Limit_radiologicdetection, subcolone_LOD, misspecification_LOD, mutation_rate, celltype,
                     mis_specfiy_pop, Strategy_name):
    flag_break, subcolone_detection, flag_treat = False, False, True
    X_sim_total = np.sum(X_sim, axis=0)
    X_true_total = np.sum(X_true, axis=0)

    # If total true cell population is >= Limit_mortality (death) or
    # the true cell number of all types are smaller than 1 (cured), stop.
    if np.any(X_true_total >= Limit_mortality) or np.all(X_true[:, -1] < 1):
        X0_sim = X_sim[:, 0]
        X_sim, t_sim, flag_break = X_true, t_true, True
        X_sim[:, 0] = X0_sim

    X_Strategy = np.append(X_Strategy_sim, X_sim, 1) if tnow == 0 else np.append(X_Strategy_sim, X_sim[:, 1:], 1)
    t_Strategy = np.append(t_Strategy_sim, t_sim + tnow) if tnow == 0 else np.append(t_Strategy_sim, t_sim[1:] + tnow)
    d_Strategy = np.append(d_Strategy_sim, d_sim, 1)

    X_Strategy_true = np.append(X_Strategy_true, X_true, 1) if tnow == 0 else np.append(X_Strategy_true, X_true[:, 1:], 1)

    # If the simulation cell number of all types are smaller than 1
    if np.all(X_sim[:, -1] < 1):
        # If true cell population is >= radiologicdetection
        if X_true_total[-1] >= Limit_radiologicdetection:
            X_true_total_atStepsize = np.sum(X_Strategy_true[:, range(0, X_Strategy_true.shape[1], Stepsize)][:, :-1], axis=0)
            replapse = np.any(X_true_total_atStepsize < Limit_radiologicdetection)
            # If it is a relapse (true total cell number reemerge from radiologic LOD), do subcolone detection
            if replapse:
                subcolone_detection = True
                X_LOD, mis_specfiy_pop = \
                    DPM_strategy_subcolone_detection(X_true[:, -1], subcolone_LOD, Limit_radiologicdetection, misspecification_LOD,
                                                     mutation_rate, celltype, mis_specfiy_pop, Strategy_name)
                X_sim[:, -1] = X_LOD
                # # Because did subcolone detection, use the drug treat the majority cells.
                # treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_LOD_i, PAR['Num_drug'], Sa_i)
                # drug_used[ind_drug] = 1
        # If true cell population is < radiologicdetection stop treating because doctors cannot not detect the tumor.
        elif X_true_total[-1] < Limit_radiologicdetection:
            flag_treat = False

    # If simulation cell population is bigger or equal than Limit_mortality and true cell population is not, will do subcolone detection.
    if np.any(X_sim_total >= Limit_mortality) and np.all(X_true_total < Limit_mortality):
        subcolone_detection = True
        X_LOD, mis_specfiy_pop = \
            DPM_strategy_subcolone_detection(X_true[:, -1], subcolone_LOD, Limit_radiologicdetection, misspecification_LOD,
                                             mutation_rate, celltype, mis_specfiy_pop, Strategy_name)
        X_sim[:, -1] = X_LOD
        # # Because did subcolone detection, use the drug treat the majority cells.
        # treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_LOD_i, PAR['Num_drug'], Sa_i)
        # drug_used[ind_drug] = 1
    return X_Strategy, X_Strategy_true, t_Strategy, d_Strategy, X_sim, flag_break, subcolone_detection, flag_treat, mis_specfiy_pop


def DPM_strategy_subcolone_detection(x, subcolone_LOD, Limit_radiologicdetection, misspecification_LOD, mutation_rate, celltype,
                                     mis_specfiy_pop, Strategy_name):
    # If true cell population is bigger or radiologic detection, can do a tissue biopsy otherwise a liquid biopsy, the LOD will be
    # higher in liquid biopsy.
    liquid_biopsy_ratio = LIQUID_BIOPSY_LOD_RATIO_DEFAULT_VAl  # if subcolone_LOD < max(np.array(LOD_LIST).astype(float)) else 1
    subcolone_LOD_i = subcolone_LOD * liquid_biopsy_ratio if np.sum(x) < Limit_radiologicdetection else subcolone_LOD
    par = dict(zip(celltype, x))
    X_LOD, mis_specfiy_pop = \
        DPM_generate_misspecification_subcolone([par], subcolone_LOD_i, misspecification_LOD, mutation_rate, celltype[1:], mis_specfiy_pop,
                                                Strategy_name)
    X_LOD = np.fromiter(X_LOD[0].values(), dtype=float)
    return X_LOD, mis_specfiy_pop


def DPM_strategy_select_drug_majority(x, Num_drug, Sa):
    treat = np.zeros([Num_drug, 1], dtype=int)
    # Index of the majority cell type.
    ind_drug = np.argmax(Sa[np.argmax(x), :])
    treat[ind_drug] = 1
    return treat, ind_drug


# Define Strategy 0.
def DPM_strategy_0(PAR, Simduration, Stepsize, Simtimestep, Limit_mortality, Limit_radiologicdetection, LSsim, misspecification_ofdecision,
                   misspecification_atsim, mis_PAR):
    # Strategy 0: Current personalized medicine:
    maxthreshold = Limit_mortality
    X0_i = DPM_generate_X0(PAR)
    g0_i = DPM_generate_g0(PAR)
    Sa_i = DPM_generate_Sa(PAR)
    T_i = DPM_generate_T(PAR)
    mis_Sa_i = DPM_generate_Sa(mis_PAR) if misspecification_ofdecision or misspecification_atsim else None

    nadir = X0_i.sum()
    tnow = 0.0
    t_Strategy0_i = np.zeros([0], dtype=float)
    X_Strategy0_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy0_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    drug_used = np.zeros([PAR['Num_drug']], dtype=int)

    '''Determine the drug usage.'''
    if misspecification_ofdecision:
        treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_i, PAR['Num_drug'], mis_Sa_i)
    elif misspecification_atsim:
        treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_i, PAR['Num_drug'], Sa_i)
    else:
        treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_i, PAR['Num_drug'], Sa_i)

    drug_used[ind_drug] = 1

    while tnow < Simduration:
        t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
        d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))

        if misspecification_atsim:
            X_i, t_i, d_i = DPM_sim_model(X0_i, T_i, g0_i, mis_Sa_i, d_i, t_i, maxthreshold, LSsim)
        else:
            X_i, t_i, d_i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)

        X_Strategy0_i = np.append(X_Strategy0_i, X_i, 1) if tnow == 0 else np.append(X_Strategy0_i, X_i[:, 1:], 1)
        t_Strategy0_i = np.append(t_Strategy0_i, t_i + tnow) if tnow == 0 else np.append(t_Strategy0_i, t_i[1:] + tnow)
        d_Strategy0_i = np.append(d_Strategy0_i, d_i, 1)

        # If total cell population is bigger or equal than Limit_mortality, stop. Mortality happens.
        # or if the cell number of all types are smaller than 1, cured. Stop.
        if (X_i[:, -1].sum() >= Limit_mortality) or np.all(X_i[:, -1] < 1):
            break

        # (1) If total cell population is bigger than 2 times of nadir, total cell population is bigger than Limit_radiologicdetection
        # (tumor can be detected), there are drugs have not been used (each drug is used only once) or
        # (2) If the total cell population reemerges from a level below the detection,
        # and there are drugs have not been used (each drug is used only once), switch to the other unused drugs.
        if ((X_i[:, -1].sum() >= 2 * nadir and X_i[:, -1].sum() >= Limit_radiologicdetection) or
                (X_i[:, 0].sum() < Limit_radiologicdetection <= X_i[:, -1].sum())) and any(drug_used == 0):
            treat_i = np.zeros([PAR['Num_drug'], 1], dtype=int)
            index_drug_not_used = [i for i, x in enumerate(list(drug_used)) if x == 0]
            treat_i[index_drug_not_used] = 1
            drug_used[index_drug_not_used] = 1
            nadir = X_i[:, -1].sum()

        # If the total cell population is smaller than the current nadir and bigger than the Limit_radiologicdetection (tumor can be detected),
        # update the nadir.
        if Limit_radiologicdetection <= X_i[:, -1].sum() < nadir:
            nadir = X_i[:, -1].sum()
        # Update X0_i and tnow.
        X0_i = X_i[:, -1]
        tnow += Stepsize
    return t_Strategy0_i, X_Strategy0_i, d_Strategy0_i


# Define Strategy 0, LOD
def DPM_strategy_LOD_0(PAR,Simduration, Stepsize, Simtimestep, Limit_mortality, Limit_radiologicdetection, PAR_LOD, subcolone_LOD,
                       misspecification_LOD, mutation_rate, mis_specfiy_pop, Strategy_name):
    # Strategy 0: Current personalized medicine:
    maxthreshold = Limit_mortality
    X0_true_i, X0_LOD_i, g0_i, Sa_i, T_i = \
        DPM_generate_X0(PAR), DPM_generate_X0(PAR_LOD), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)

    celltype = ['Spop', 'R1pop', 'R2pop', 'R12pop']
    nadir = X0_true_i.sum()
    if nadir < Limit_radiologicdetection:
        raise Exception('Initial tumor cell number is smaller than limit of radiologic detction.')
    tnow = 0.0
    t_Strategy0_i = np.zeros([0], dtype=float)
    X_Strategy0_i, X_Strategy0_true_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float), np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy0_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    drug_used = np.zeros([PAR['Num_drug']], dtype=int)

    treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_LOD_i, PAR['Num_drug'], Sa_i)
    drug_used[ind_drug] = 1

    flag_treat, subcolone_detection = True, False
    while tnow < Simduration:
        if flag_treat:
            t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))

            X_i, t_i, d_i = DPM_sim_model(X0_LOD_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, False)
            X_true_i, t_true_i, _ = DPM_sim_model(X0_true_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, True)
            X_true_total_i = np.sum(X_true_i, axis=0)

            X_Strategy0_i, X_Strategy0_true_i, t_Strategy0_i, d_Strategy0_i, X_i, flag_break, subcolone_detection, flag_treat, mis_specfiy_pop =\
                DPM_strategy_LOD(X_i, t_i, d_i, X_true_i, t_true_i, tnow, X_Strategy0_i, X_Strategy0_true_i, t_Strategy0_i, d_Strategy0_i,
                                 Stepsize, Limit_mortality, Limit_radiologicdetection, subcolone_LOD, misspecification_LOD, mutation_rate,
                                 celltype, mis_specfiy_pop, Strategy_name)
            if flag_break:
                break

            if subcolone_detection:
                treat_i, ind_drug = DPM_strategy_select_drug_majority(X_i[:, -1], PAR['Num_drug'], Sa_i)
                drug_used[ind_drug] = 1

            # (1) If the true total cell population is bigger than 2 times of nadir and total cell population is bigger than
            # Limit_radiologicdetection (tumor can be detected), there are drugs have not been used (each drug is used only once) or
            # (2) If the true total cell population reemerges from radiologic detection,
            # and there are drugs have not been used (each drug is used only once), switch to the other unused drugs.
            if ((X_true_total_i[-1] >= 2 * nadir and X_true_total_i[-1] >= Limit_radiologicdetection) or
                    (X_true_total_i[0] < Limit_radiologicdetection <= X_true_total_i[-1])) \
                    and any(drug_used == 0) and (not subcolone_detection):
                treat_i = np.zeros([PAR['Num_drug'], 1], dtype=int)
                index_drug_not_used = [i for i, x in enumerate(list(drug_used)) if x == 0]
                treat_i[index_drug_not_used], drug_used[index_drug_not_used] = 1, 1
                nadir = X_true_total_i[-1]
            # If the total cell population is smaller than the current nadir and bigger than the Limit_radiologicdetection
            # (tumor can be detected), update the nadir.
            if X_true_total_i[-1] < nadir:
                nadir = X_true_i[:, -1].sum() if X_true_i[:, -1].sum() >= Limit_radiologicdetection else 0
            # Update X0_LOD_i, X0_true_i and tnow.
            X0_LOD_i = X_i[:, -1]
            X0_true_i = X_true_i[:, -1]
        else:
            # no treat, let proliferate
            X0_true_i, X_Strategy0_i, X_Strategy0_true_i, t_Strategy0_i, d_Strategy0_i, flag_treat, flag_break, mis_specfiy_pop = \
                DPM_strategy_pro(X0_true_i, g0_i, Sa_i, T_i, PAR['Num_drug'], tnow, X_Strategy0_i, X_Strategy0_true_i, t_Strategy0_i, d_Strategy0_i,
                                 Stepsize, maxthreshold, Simtimestep, Limit_mortality, Limit_radiologicdetection, subcolone_LOD,
                                 misspecification_LOD, mutation_rate, celltype, mis_specfiy_pop, Strategy_name)
            if flag_break:
                break

            if flag_treat:
                X0_LOD_i = X0_true_i
                treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_true_i, PAR['Num_drug'], Sa_i)
                drug_used[ind_drug] = 1

        tnow += Stepsize
    return t_Strategy0_i, X_Strategy0_i, d_Strategy0_i, mis_specfiy_pop


# Define Strategy 1.
def DPM_strategy_1(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, mis_specification, mis_PAR):
    # Strategy 1: Minimize the total cell population.
    # In each Stepsize, select the d_i that minimizes the total cell population.
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_Strategy1_i = np.zeros([0], dtype=float)
    X_Strategy1_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy1_i = np.zeros([PAR['Num_drug'], 0], dtype=float)

    while tnow < Simduration:
        X_i_total, t_i_total, d_i_total, X_i_end_total, t_sim = [], [], [], [], []
        t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
        for i_dose_combination in dose_combination:
            treat_i = np.array([i_dose_combination], dtype=float).T
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
            d_i_total.append(d_i)
            if mis_specification:
                mis_X_i, mis_t_i, _ = DPM_sim_model(X0_i, mis_T_i, mis_g0_i, mis_Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(mis_X_i)
                X_i_end_total.append(mis_X_i[:, -1].sum())
                t_sim.append(mis_t_i[-1])
            else:
                X_i, t__i, d__i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(X_i)
                t_i_total.append(t__i)
                d_i_total[-1] = d__i
                X_i_end_total.append(X_i[:, -1].sum())
                t_sim.append(t__i[-1])

        # Find the d_i that minimizes total cell population.
        if np.argmin(X_i_end_total) < Limit_mortality:
            index_select = np.argmin(X_i_end_total)
        else:
            # This means simulation of all of the dose combination reach mortality. Then select the longest survival time.
            index_select = np.argmax(t_sim)
        d_i_minimum = d_i_total[index_select]

        if mis_specification:
            X_i_minimum, t_i_minimum, d_i_minimum = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, LSsim)
        else:
            X_i_minimum, t_i_minimum = X_i_total[index_select], t_i_total[index_select]

        X_Strategy1_i = np.append(X_Strategy1_i, X_i_minimum, 1) if tnow == 0 else np.append(X_Strategy1_i, X_i_minimum[:, 1:], 1)
        t_Strategy1_i = np.append(t_Strategy1_i, t_i_minimum + tnow) if tnow == 0 else np.append(t_Strategy1_i, t_i_minimum[1:] + tnow)
        d_Strategy1_i = np.append(d_Strategy1_i, d_i_minimum, 1)

        # If total cell population is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i_minimum[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i_minimum[:, -1] < 1):
            break
        # Update X0_i and tnow.
        X0_i = X_i_minimum[:, -1]
        tnow += Stepsize
    return t_Strategy1_i, X_Strategy1_i, d_Strategy1_i


# Define Strategy 2.
def DPM_strategy_2(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, Strategy2threshold,
                   misspecification_ofdecision, misspecification_atsim, mis_PAR):
    # Strategy 2: Minimize the risk of incurable cells developing unless there is an immediate threat of mortality.
    # Strategy 2.1: threshold is 1e9
    # Strategy 2.2: threshold is 1e11
    maxthreshold = Limit_mortality
    if misspecification_atsim:
        X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(mis_PAR), DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)
        g0_true_i, Sa_true_i, T_true_i = DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    else:
        X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
        g0_true_i, Sa_true_i, T_true_i = None, None, None

    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        misspecification_ofdecision else (None, None, None)

    tnow = 0.0
    t_Strategy2_i = np.zeros([0], dtype=float)
    X_Strategy2_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy2_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    while tnow < Simduration:
        X_i_total, t_i_total, d_i_total, X_i_end_total, X_i_end_multi_resis, t_sim = [], [], [], [], [], []
        t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
        for i_dose_combination in dose_combination:
            treat_i = np.array([i_dose_combination], dtype=float).T
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
            d_i_total.append(d_i)
            if misspecification_ofdecision:
                mis_X_i, mis_t_i, _ = DPM_sim_model(X0_i, mis_T_i, mis_g0_i, mis_Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(mis_X_i)
                X_i_end_total.append(mis_X_i[:, -1].sum())
                X_i_end_multi_resis.append(mis_X_i[-1, -1])
                t_sim.append(mis_t_i[-1])
            elif misspecification_atsim:
                mis_X_i, mis_t_i, _ = DPM_sim_model(X0_i, T_true_i, g0_true_i, Sa_true_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(mis_X_i)
                X_i_end_total.append(mis_X_i[:, -1].sum())
                X_i_end_multi_resis.append(mis_X_i[-1, -1])
                t_sim.append(mis_t_i[-1])
            else:
                X_i, t__i, d__i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, True)
                X_i_total.append(X_i)
                t_i_total.append(t__i)
                d_i_total[-1] = d__i
                X_i_end_total.append(X_i[:, -1].sum())
                X_i_end_multi_resis.append(X_i[-1, -1])
                t_sim.append(t__i[-1])
        # If one or more doses cause mortality while others do not, exclude the dose(s) that result in mortality.

        assert max(t_sim) <= Stepsize
        if min(X_i_end_total) < Limit_mortality <= max(X_i_end_total):
            ind = [i for i, x in enumerate(X_i_end_total) if x < Limit_mortality]
            X_i_total = [X_i_total[i] for i in ind]
            if not (misspecification_ofdecision or misspecification_atsim):
                t_i_total = [t_i_total[i] for i in ind]
            d_i_total = [d_i_total[i] for i in ind]
            X_i_end_total = [X_i_end_total[i] for i in ind]
            X_i_end_multi_resis = [X_i_end_multi_resis[i] for i in ind]
            t_sim = [t_sim[i] for i in ind]
            assert min(X_i_end_total) < Limit_mortality

        if min(X_i_end_total) < Limit_mortality:
            # If the current total cell population does not exceed the threshold, minimize the multiply-resistant population.
            if X0_i.sum() <= Strategy2threshold:
                index_minimum_end_multi_resis = np.argmin(X_i_end_multi_resis)
                index_same_minimum_end_multi_resis = [i for i, x in enumerate(X_i_end_multi_resis) if x ==
                                                      X_i_end_multi_resis[index_minimum_end_multi_resis]]
                # If there is only one d_i that gives the minimum multiply-resistant population.
                if len(index_same_minimum_end_multi_resis) == 1:
                    index_select = index_minimum_end_multi_resis
                # If there are multiple d_i that give the same minimum multiply-resistant population.
                else:
                    X_i_end_total_same_minimum_end_multi_resis = [X_i_end_total[i] for i in index_same_minimum_end_multi_resis]
                    # Find the d_i gives the minimum total cell population in all the d_i giving the same minimum multiply-resistant population.
                    index_minimum_end_total_same_minimum_end_multi_resis = np.argmin(X_i_end_total_same_minimum_end_multi_resis)
                    index_select = index_same_minimum_end_multi_resis[index_minimum_end_total_same_minimum_end_multi_resis]
            # If the current total cell population exceeds the threshold, minimize the total population.
            else:
                index_select = np.argmin(X_i_end_total)
        # If all dose options result in mortality, select the dose that extends survival the longest.
        else:
            index_select = np.argmax(t_sim)

        d_i_minimum = d_i_total[index_select]
        if misspecification_ofdecision or misspecification_atsim:
            X_i_minimum, t_i_minimum, d_i_minimum = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, LSsim)
        else:
            X_i_minimum, t_i_minimum = X_i_total[index_select], t_i_total[index_select]

        X_Strategy2_i = np.append(X_Strategy2_i, X_i_minimum, 1) if tnow == 0 else np.append(X_Strategy2_i, X_i_minimum[:, 1:], 1)
        t_Strategy2_i = np.append(t_Strategy2_i, t_i_minimum + tnow) if tnow == 0 else np.append(t_Strategy2_i, t_i_minimum[1:] + tnow)
        d_Strategy2_i = np.append(d_Strategy2_i, d_i_minimum, 1)

        # If total cell population is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i_minimum[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i_minimum[:, -1] < 1):
            break
        # Update X0_i and tnow
        X0_i = X_i_minimum[:, -1]
        tnow += Stepsize
    return t_Strategy2_i, X_Strategy2_i, d_Strategy2_i


# Define Strategy 2_LOD.
def DPM_strategy_LOD_2(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, Limit_radiologicdetection, Strategy2threshold,
                       PAR_LOD, subcolone_LOD, misspecification_LOD, mutation_rate, mis_specfiy_pop, Strategy_name):
    # Strategy 2: Minimize the risk of incurable cells developing unless there is an immediate threat of mortality.
    # Strategy 2.1: threshold is 1e9
    # Strategy 2.2: threshold is 1e11
    maxthreshold = Limit_mortality
    X0_true_i, X0_LOD_i, g0_i, Sa_i, T_i = \
        DPM_generate_X0(PAR), DPM_generate_X0(PAR_LOD), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    if X0_true_i.sum() < Limit_radiologicdetection:
        raise Exception('Initial tumor cell number is smaller than limit of radiologic detction.')

    celltype = ['Spop', 'R1pop', 'R2pop', 'R12pop']
    tnow = 0.0
    t_Strategy2_i = np.zeros([0], dtype=float)
    X_Strategy2_i, X_Strategy2_true_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float), np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy2_i = np.zeros([PAR['Num_drug'], 0], dtype=float)

    flag_treat, subcolone_detection = True, False
    while tnow < Simduration:
        if flag_treat:
            X_i_total, t_i_total, d_i_total, X_i_end_total, X_i_end_multi_resis, t_sim = [], [], [], [], [], []
            for i_dose_combination in dose_combination:
                t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
                treat_i = np.array([i_dose_combination], dtype=float).T
                d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
                d_i_total.append(d_i)
                X_i, t_i, d_i = DPM_sim_model(X0_LOD_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, True)
                X_i_total.append(X_i)
                t_i_total.append(t_i)
                if not d_i.any():
                    d_i = np.array(treat_i)
                d_i_total[-1] = d_i
                X_i_end_total.append(X_i[:, -1].sum())
                X_i_end_multi_resis.append(X_i[-1, -1])
                t_sim.append(t_i[-1])
            # If the current total cell population does not exceed the threshold, minimize the multiply-resistant population.
            if X0_true_i.sum() <= Strategy2threshold:
                index_minimum_end_multi_resis = np.argmin(X_i_end_multi_resis)
                index_same_minimum_end_multi_resis = [i for i, x in enumerate(X_i_end_multi_resis) if x ==
                                                      X_i_end_multi_resis[index_minimum_end_multi_resis]]
                # If there is only one d_i that gives the minimum multiply-resistant population.
                if len(index_same_minimum_end_multi_resis) == 1:
                    index_select = index_minimum_end_multi_resis
                # If there are multiple d_i that give the same minimum multiply-resistant population.
                else:
                    X_i_end_total_same_minimum_end_multi_resis = [X_i_end_total[i] for i in index_same_minimum_end_multi_resis]
                    # Find the d_i gives the minimum total cell population in all the d_i giving the same minimum multiply-resistant population.
                    index_minimum_end_total_same_minimum_end_multi_resis = np.argmin(X_i_end_total_same_minimum_end_multi_resis)
                    index_select = index_same_minimum_end_multi_resis[index_minimum_end_total_same_minimum_end_multi_resis]
            # If the current total cell population exceeds the threshold, minimize the total population.
            else:
                index_select = np.argmin(X_i_end_total)

            t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
            X_i_minimum, t_i_minimum,  d_i_minimum = X_i_total[index_select], t_i_total[index_select], d_i_total[index_select]
            X_true_i, t_true_i, _ = DPM_sim_model(X0_true_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, True)

            X_Strategy2_i, X_Strategy2_true_i, t_Strategy2_i, d_Strategy2_i, X_i_minimum, flag_break, subcolone_detection, flag_treat, \
                mis_specfiy_pop = \
                DPM_strategy_LOD(X_i_minimum, t_i_minimum,  d_i_minimum, X_true_i, t_true_i, tnow, X_Strategy2_i, X_Strategy2_true_i,
                                 t_Strategy2_i, d_Strategy2_i, Stepsize, Limit_mortality, Limit_radiologicdetection, subcolone_LOD,
                                 misspecification_LOD, mutation_rate, celltype, mis_specfiy_pop, Strategy_name)
            if flag_break:
                break

            # Update X0_LOD_i, X0_true_i and tnow
            X0_LOD_i = X_i_minimum[:, -1]
            X0_true_i = X_true_i[:, -1]
        else:
            # no treat, let proliferate
            X0_true_i, X_Strategy2_i, X_Strategy2_true_i, t_Strategy2_i, d_Strategy2_i, flag_treat, flag_break, mis_specfiy_pop = \
                DPM_strategy_pro(X0_true_i, g0_i, Sa_i, T_i, PAR['Num_drug'], tnow, X_Strategy2_i, X_Strategy2_true_i, t_Strategy2_i,
                                 d_Strategy2_i, Stepsize, maxthreshold, Simtimestep, Limit_mortality, Limit_radiologicdetection, subcolone_LOD,
                                 misspecification_LOD, mutation_rate, celltype, mis_specfiy_pop, Strategy_name)
            if flag_break:
                break

            if flag_treat:
                X0_LOD_i = X0_true_i

        tnow += Stepsize
    return t_Strategy2_i, X_Strategy2_i, d_Strategy2_i, mis_specfiy_pop


# Define Strategy 3.
def DPM_strategy_3(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, mis_specification, mis_PAR):
    # Strategy 3: Minimize the predicted total cell population unless the first multiply-resistant cell will arise by the selection of the d_i
    # which gives the minimum total cell population.

    # At each Stepsize:
    # If the predicted multiply-resistant population < 1 or it is curable, select d_i to minimize the total cell population.
    # If the selected d_i rises the first multiply-resistant cell, re-select d_i to minimize the multiply-resistant population.
    # else if the current multiply-resistant >= 1 and multiply-resistant is not curable, minimize the total cell population.
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_Strategy3_i = np.zeros([0], dtype=float)
    X_Strategy3_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy3_i = np.zeros([PAR['Num_drug'], 0], dtype=float)

    while tnow < Simduration:
        X_i_total, t_i_total, d_i_total, X_i_end_total, X_i_end_multi_resis, t_sim = [], [], [], [], [], []
        for i_dose_combination in dose_combination:
            t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
            treat_i = np.array([i_dose_combination], dtype=float).T
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
            d_i_total.append(d_i)
            if mis_specification:
                mis_X_i, mis_t_i, _ = DPM_sim_model(X0_i, mis_T_i, mis_g0_i, mis_Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(mis_X_i)
                t_i_total.append(mis_t_i)
                X_i_end_total.append(mis_X_i[:, -1].sum())
                X_i_end_multi_resis.append(mis_X_i[-1, -1])
                t_sim.append(mis_t_i[-1])
            else:
                X_i, t__i, d__i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(X_i)
                t_i_total.append(t__i)
                d_i_total[-1] = d__i
                X_i_end_total.append(X_i[:, -1].sum())
                X_i_end_multi_resis.append(X_i[-1, -1])
                t_sim.append(t__i[-1])

        # If the current multiply-resistant < 1 or the multiply-resistant is curable (curable means any Sa[-1, :] >= g0),
        # minimize the total cell population.
        index_minimum_end_total = np.argmin(X_i_end_total)
        index_select = index_minimum_end_total
        if (X0_i[-1] < 1 or (np.any(np.greater_equal(mis_Sa_i[-1, :], mis_g0_i[-1]))) if mis_specification
                else np.any(np.greater_equal(Sa_i[-1, :], g0_i[-1]))):
            X_i_minimum = X_i_total[index_minimum_end_total]
            # If the first multiply-resistant cell will arise under the selected d_i, minimize the multiply-resistant population.
            if X0_i[-1] < 1 <= X_i_minimum[-1, -1]:
                index_minimum_end_multi_resis = np.argmin(X_i_end_multi_resis)
                index_same_minimum_end_multi_resis = [i for i, x in enumerate(X_i_end_multi_resis) if x ==
                                                      X_i_end_multi_resis[index_minimum_end_multi_resis]]
                # If there is only 1 d_i that gives the minimum multiply-resistant population.
                if len(index_same_minimum_end_multi_resis) == 1:
                    index_select = index_minimum_end_multi_resis
                # If there are multiple d_i that give the same minimum multiply-resistant populaiton, minimize the total cell population in all
                # the d_i giving the same minimum multiply-resistant populaiton.
                else:
                    X_i_end_total_same_minimum_end_multi_resis = [X_i_end_total[i] for i in index_same_minimum_end_multi_resis]
                    # Find the d_i gives the minimum total cell population in all the d_i giving the same minimum multiply-resistant population.
                    index_minimum_end_total_same_minimum_end_multi_resis = np.argmin(X_i_end_total_same_minimum_end_multi_resis)
                    index_select = index_same_minimum_end_multi_resis[index_minimum_end_total_same_minimum_end_multi_resis]
                # If the minimize of multiply-resistant populaiton reach mortality, still minimize total cell population.
                if X_i_end_total[index_select] >= Limit_mortality:
                    index_select = index_minimum_end_total
        # If current multiply-resistant >= 1 and multiply-resistant is not curable, minimize the total cell population.
        else:
            pass
        # This means simulation of all of the dose combination reach mortality. Then select the longest survival time.
        if X_i_end_total[index_select] >= Limit_mortality:
            index_select = np.argmax(t_sim)
        d_i_minimum = d_i_total[index_select]
        if mis_specification:
            X_i_minimum, t_i_minimum, d_i_minimum = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, LSsim)
        else:
            X_i_minimum, t_i_minimum = X_i_total[index_select], t_i_total[index_select]

        X_Strategy3_i = np.append(X_Strategy3_i, X_i_minimum, 1) if tnow == 0 else np.append(X_Strategy3_i, X_i_minimum[:, 1:], 1)
        t_Strategy3_i = np.append(t_Strategy3_i, t_i_minimum + tnow) if tnow == 0 else np.append(t_Strategy3_i, t_i_minimum[1:] + tnow)
        d_Strategy3_i = np.append(d_Strategy3_i, d_i_minimum, 1)

        # If total cell populatin is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i_minimum[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i_minimum[:, -1] < 1):
            break
        # Update X0_i and tnow.
        X0_i = X_i_minimum[:, -1]
        tnow += Stepsize
    return t_Strategy3_i, X_Strategy3_i, d_Strategy3_i


# Define Strategy 4.
def DPM_strategy_4(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, mis_specification, mis_PAR):
    # Strategy 4: Estimate the time to either incurability or death, and react to the most proximal threat as long as there is a chance of cure.
    # In each Stepsize, evaluate the predicted durations toward incurability (multiply-resistant population >= 1) and
    # mortality (population >= 1e13) dictated by the growth of S, R1, R2, ... and multiply-resistant populations.
    # For each dosage combination d_i, define τ_inc(d) as the predicted time to incurability (multiply-resistant >= 1), given the currently
    # observed population and d_i fixed.
    # Define τ_S(d) as the predicted time to S causing mortality (S > 1e13), given the currently observed population and d_i fixed.
    # τ_R1(d), R1 causing mortality (R1 > 1e13).
    # τ_R2(d), R2 causing mortality (R2 > 1e13).
    # ...
    # τ_X(d), X type causing mortality (X > 1e13).
    # τ_multi_resis(d), multiply-resistant population causing mortality (multiply_resistant > 1e13).
    # If the current multiply-resistant population < 1 or multiply-resistant population is curable, i.e., there exists some d_i such that each each
    # component of diag(Sa*d) > g0, vary d to maximize min(τ_inc,τ_s,τ_R1,τ_R2,...,τ_multi_resis) with the constraint that
    # min(τ_S,τ_R1,τ_R2,...,τ_multi_resis) > Stepsize. If such a dosage combination does not exist, maximize min(τ_S,τ_R1,τ_R2,..,τ_multi_resis).
    # If the current R_multi_resis >= 1, and multiply-resistant population is not curable, maximize min(τ_S,τ_R1,τ_R2,...,τ_multi_resis).
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_Strategy4_i = np.zeros([0], dtype=float)
    X_Strategy4_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy4_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    curable_multi_resis = False

    # If the current multiply-resistant populations < 1 or multiply-resistant populations is curable (curable means any Sa[-1, :] >= g0).
    if np.any(np.greater_equal(mis_Sa_i[-1, :], mis_g0_i[-1])) if mis_specification else np.any(np.greater_equal(Sa_i[-1, :], g0_i[-1])):
        curable_multi_resis = True

    while tnow < Simduration:
        τ_total = np.zeros([PAR['Num_cell_type'] + 1, 0], dtype=float)
        d_τ_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
        timeleft = Simduration - tnow
        t_i = np.arange(0, timeleft + Simtimestep, Simtimestep)
        for i_dose_combination in dose_combination:
            treat_i = np.array([i_dose_combination], dtype=float).T
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
            τ_i_drug = np.zeros(PAR['Num_cell_type'] + 1, dtype=float)
            if mis_specification:
                mis_X_τ_i, _, _ = DPM_sim_model(X0_i, mis_T_i, mis_g0_i, mis_Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_τ_i_use = mis_X_τ_i
            else:
                X_τ_i, _, _ = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_τ_i_use = X_τ_i

            # Calculate the τ_inc,τ_s,τ_R1,τ_R2,...,τ_multi_resis.
            index_celltype = [-1]
            index_celltype.extend(list(range(PAR['Num_cell_type'])))
            limit_celltype = [1]
            limit_celltype.extend([Limit_mortality] * PAR['Num_cell_type'])
            for i in range(len(index_celltype)):
                occurrences_τ, = np.where(X_τ_i_use[index_celltype[i], :] >= limit_celltype[i])
                if not occurrences_τ.size == 0:
                    τ_i_drug[i] = occurrences_τ[0]
                else:
                    τ_i_drug[i] = timeleft

            τ_total = np.append(τ_total, np.reshape(τ_i_drug, (len(τ_i_drug), 1)), 1)
            d_τ_i = np.append(d_τ_i, treat_i, 1)

        τ_S_to_multi_resis = τ_total[1:, :]
        min_τ_S_to_multi_resis = np.amin(τ_S_to_multi_resis, 0)
        index_min_τ_S_to_multi_resis_biggerthanStepsize, = np.where(min_τ_S_to_multi_resis > Stepsize)
        # If the current multiply-resistant population < 1 or it is curable.
        if X0_i[-1] < 1 or curable_multi_resis:
            # If there exists min(τ_S, τ_R1, τ_R2,...,τ_multi_resis) > Stepsize, then maximize min(τ_inc, τ_S, τ_R1, τ_R2, τ_R12) among the d_i that
            # meet the criteria min(τ_S, τ_R1, τ_R2,...,τ_multi_resis) > Stepsize.
            if index_min_τ_S_to_multi_resis_biggerthanStepsize.size:
                τ_total = τ_total[:, index_min_τ_S_to_multi_resis_biggerthanStepsize]
                d_τ_i = d_τ_i[:, index_min_τ_S_to_multi_resis_biggerthanStepsize]
                min_τ_inc_S_to_multi_resis = np.amin(τ_total, 0)
                index_maxmin_τ_inc_S_to_multi_resis = np.argmax(min_τ_inc_S_to_multi_resis)
                d_τ = np.array([d_τ_i[:, index_maxmin_τ_inc_S_to_multi_resis]]).T
            # If there doesn't exist exists min(τ_S, τ_R1, τ_R2,...,τ_multi_resis) > Stepsize, then maximize min(τ_S, τ_R1, τ_R2,...,τ_multi_resis).
            else:
                index_maxmin_τ_S_to_multi_resis = np.argmax(min_τ_S_to_multi_resis)
                d_τ = np.array([d_τ_i[:, index_maxmin_τ_S_to_multi_resis]]).T
        # If the current multiply-resistant population >= 1 and it is not curable, maximize min(τ_S, τ_R1, τ_R2,...,τ_multi_resis).
        elif X0_i[-1] >= 1 and (not curable_multi_resis):
            index_maxmin_τ_S_to_multi_resis = np.argmax(min_τ_S_to_multi_resis)
            d_τ = np.array([d_τ_i[:, index_maxmin_τ_S_to_multi_resis]]).T
        else:
            raise ValueError('Wrong situation')

        t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
        d_i = np.tile(d_τ, (1, t_i.shape[0] - 1))
        X_i, t_i, d_i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)

        X_Strategy4_i = np.append(X_Strategy4_i, X_i, 1) if tnow == 0 else np.append(X_Strategy4_i, X_i[:, 1:], 1)
        t_Strategy4_i = np.append(t_Strategy4_i, t_i + tnow) if tnow == 0 else np.append(t_Strategy4_i, t_i[1:] + tnow)
        d_Strategy4_i = np.append(d_Strategy4_i, d_i, 1)

        # If total cell population is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i[:, -1] < 1):
            break
        # Update X0_i and tnow.
        X0_i = X_i[:, -1]
        tnow += Stepsize
    return t_Strategy4_i, X_Strategy4_i, d_Strategy4_i


# Define Strategy 5.
def DPM_strategy_5(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, lookahead_step, LSsim, mis_specification, mis_PAR):
    # Multistep extension of strategy1.
    # In each Stepsize, select the d_i that minimizes the total cell population.
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_Strategy5_i = np.zeros([0], dtype=float)
    X_Strategy5_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy5_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    # Generate the decision tree using the variable value of 'lookahead_step' as the value of depth.
    decision_tree_all = list(itertools.product(list(range(len(dose_combination))), repeat=lookahead_step))
    decision_tree_all.sort()
    # Terminate one node if the total populations exceed some other nodes of the same depth,
    # terminate one node if the multiply-resistant populations exceed some other nodes of the same depth or
    # terminate one node if the total and multiply-resistant populations both exceed some other nodes of the same depth.
    terminate = 'total_pop'  # 'multi_resis' or 'both'
    while tnow < Simduration:
        if mis_specification:
            node, total_pop, multi_resis, _, t_sim = DPM_strategy_recurse_multistep(
                X0_i, mis_T_i, mis_g0_i, mis_Sa_i, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, lookahead_step,
                maxthreshold, LSsim, tnow, terminate, decision_tree_all)
        else:
            node, total_pop, multi_resis, _, t_sim = DPM_strategy_recurse_multistep(
                X0_i, T_i, g0_i, Sa_i, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, lookahead_step, maxthreshold, LSsim,
                tnow, terminate, decision_tree_all)

        node_minimum = None
        for index in reversed(range(total_pop.shape[1])):
            X_i_end_total = total_pop[:, index]
            index_minimum_end_total = np.argmin(X_i_end_total)
            if X_i_end_total.min() < np.inf:
                if X_i_end_total.min() < Limit_mortality:
                    node_minimum = node[index_minimum_end_total, :]
                    break
                else:
                    # This means simulation at this depth all nodes reach mortality. Then select the longest survival time.
                    index_maximum_survival = np.argmax(t_sim[:, index])
                    node_minimum = node[index_maximum_survival, :]
                    break

        d_i = np.empty([PAR['Num_drug'], 0], dtype=float)
        for i_deep in node_minimum:
            d_i = np.append(d_i, np.tile(np.array([dose_combination[i_deep]], dtype=float).T, (1, Stepsize)), 1)
        t_i = np.arange(0, node_minimum.shape[0] * Stepsize + Simtimestep, Simtimestep)
        X_i, t_i, d_i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)

        X_Strategy5_i = np.append(X_Strategy5_i, X_i, 1) if tnow == 0 else np.append(X_Strategy5_i, X_i[:, 1:], 1)
        t_Strategy5_i = np.append(t_Strategy5_i, t_i + tnow) if tnow == 0 else np.append(t_Strategy5_i, t_i[1:] + tnow)
        d_Strategy5_i = np.append(d_Strategy5_i, d_i, 1)

        # If total cell population is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i[:, -1] < 1):
            break
        # Update X0_i and tnow.
        X0_i = X_i[:, -1]
        tnow += node_minimum.shape[0] * Stepsize
    return t_Strategy5_i, X_Strategy5_i, d_Strategy5_i


# Define Strategy 6.
def DPM_strategy_6(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, lookahead_step, LSsim, Strategy6threshold,
                   mis_specification, mis_PAR):
    # Multistep extension of strategy2.
    # Minimize the risk of incurable cells developing unless there is an immediate threat of mortality.
    # Strategy 6: threshold is 1e9
    # Strategy 7: threshold is 1e11
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_Strategy6_i = np.zeros([0], dtype=float)
    X_Strategy6_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy6_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    # Generate the decision tree using the variable value of 'lookahead_step' as the value of depth.
    decision_tree_all = list(itertools.product(list(range(len(dose_combination))), repeat=lookahead_step))
    decision_tree_all.sort()
    # Terminate one node if both the total populations and multiply-resistant populations exceed some other nodes of the same depth.
    terminate = 'both'
    while tnow < Simduration:
        if mis_specification:
            node, total_pop, multi_resis, _, t_sim = DPM_strategy_recurse_multistep(
                X0_i, mis_T_i, mis_g0_i, mis_Sa_i, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, lookahead_step,
                maxthreshold, LSsim, tnow, terminate, decision_tree_all)
        else:
            node, total_pop, multi_resis, _, t_sim = DPM_strategy_recurse_multistep(
                X0_i, T_i, g0_i, Sa_i, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, lookahead_step, maxthreshold, LSsim,
                tnow, terminate, decision_tree_all)

        index_select = None
        for index in reversed(range(total_pop.shape[1])):
            X_i_end_total = total_pop[:, index]
            multi_resis_end = multi_resis[:, index]
            if multi_resis_end.min() < np.inf and X_i_end_total.min() < np.inf:
                # If the current total cell population does not exceed the threshold, minimize the multiply-resistant population.
                if X0_i.sum() <= Strategy6threshold and np.argmin(multi_resis_end) < Limit_mortality:
                    index_minimum_end_multi_resis = np.argmin(multi_resis_end)
                    index_same_minimum_end_multi_resis = [i for i, x in enumerate(multi_resis_end) if x ==
                                                          multi_resis_end[index_minimum_end_multi_resis]]
                    # If there is only one d_i that gives the minimum multiply-resistant population.
                    if len(index_same_minimum_end_multi_resis) == 1:
                        index_select = index_minimum_end_multi_resis
                        break
                    else:
                        # If there multiple d_i that give the same minimum multiply-resistant population.
                        X_i_end_total_same_minimum_end_multi_resis = [X_i_end_total[i] for i in index_same_minimum_end_multi_resis]
                        # Find the d_i gives the minimum total cell population in all the d_i giving the same minimum multiply-resistant population.
                        index_minimum_end_total_same_minimum_end_multi_resis = np.argmin(X_i_end_total_same_minimum_end_multi_resis)
                        index_select = index_same_minimum_end_multi_resis[index_minimum_end_total_same_minimum_end_multi_resis]
                        break
                # If the current total cell population exceeds the threshold, minimize the total population.
                elif np.argmin(X_i_end_total) < Limit_mortality:
                    index_select = np.argmin(X_i_end_total)
                    break
                else:
                    # This means simulation at this depth all nodes reach mortality. Then select the longest survival time.
                    index_select = np.argmax(t_sim[:, index])
                    break

        node_minimum = node[index_select, :]
        d_i = np.empty([PAR['Num_drug'], 0], dtype=float)
        for i_deep in node_minimum:
            d_i = np.append(d_i, np.tile(np.array([dose_combination[i_deep]]).T, (1, Stepsize)), 1)
        t_i = np.arange(0, node_minimum.shape[0] * Stepsize + Simtimestep, Simtimestep)
        X_i, t_i, d_i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)

        X_Strategy6_i = np.append(X_Strategy6_i, X_i, 1) if tnow == 0 else np.append(X_Strategy6_i, X_i[:, 1:], 1)
        t_Strategy6_i = np.append(t_Strategy6_i, t_i + tnow) if tnow == 0 else np.append(t_Strategy6_i, t_i[1:] + tnow)
        d_Strategy6_i = np.append(d_Strategy6_i, d_i, 1)

        # If total cell population is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i[:, -1] < 1):
            break
        # Update X0_i and tnow
        X0_i = X_i[:, -1]
        tnow += node_minimum.shape[0] * Stepsize
    return t_Strategy6_i, X_Strategy6_i, d_Strategy6_i


# Define Strategy 8.
def DPM_strategy_8(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, lookahead_step, LSsim, mis_specification, mis_PAR):
    # Multistep extension of strategy 3.
    # Strategy 3: Minimize the predicted total cell population unless the first multiply-resistant cell will arise by the selection of the d_i
    # which gives the minimum total cell population.
    # At each Stepsize:
    # If the predicted multiply-resistant < 1 or multiply-resistant is curable, select d_i to minimize the total cell population.
    # # # If the selected d_i rises the first multiply-resistant cell, re-select d_i to minimize the multiply-resistant population.
    # else if the current multiply-resistant >= 1 and multiply-resistant is not curable, minimize the total cell population.
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_Strategy8_i = np.zeros([0], dtype=float)
    X_Strategy8_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy8_i = np.zeros([PAR['Num_drug'], 0], dtype=float)

    # Generate the decision tree using the variable value of 'lookahead_step' as the value of depth.
    decision_tree_all = list(itertools.product(list(range(len(dose_combination))), repeat=lookahead_step))
    decision_tree_all.sort()
    # Terminate one node if both the total populations and multiply-resistant populations exceed some other nodes of the same depth.
    terminate = 'both'
    while tnow < Simduration:
        if mis_specification:
            node, total_pop, multi_resis, _, t_sim = DPM_strategy_recurse_multistep(
                X0_i, mis_T_i, mis_g0_i, mis_Sa_i, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, lookahead_step,
                maxthreshold, LSsim, tnow, terminate, decision_tree_all)
        else:
            node, total_pop, multi_resis, _, t_sim = DPM_strategy_recurse_multistep(
                X0_i, T_i, g0_i, Sa_i, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, lookahead_step, maxthreshold, LSsim,
                tnow, terminate, decision_tree_all)

        index_select = None
        for index in reversed(range(total_pop.shape[1])):
            X_i_end_total = total_pop[:, index]
            multi_resis_end = multi_resis[:, index]
            if multi_resis_end.min() < np.inf and X_i_end_total.min() < np.inf:
                index_minimum_end_total = np.argmin(X_i_end_total)
                index_select = index_minimum_end_total
                # If the current multiply-resistant < 1 or the multiply-resistant is curable (curable means any Sa[-1, :] >= g0),
                # minimize the total cell population.
                if X0_i[-1] < 1 or (np.any(np.greater_equal(mis_Sa_i[-1, :], mis_g0_i[-1])) if mis_specification else np.any(
                        np.greater_equal(Sa_i[-1, :], g0_i[-1]))):  # g0_i[0].item() - Sa_i[-1, :].sum() <= 0:
                    multi_resis_end_minimum = multi_resis_end[index_minimum_end_total]
                    # If the first multiply-resistant cell will arise under the selected d_i, minimize the multiply-resistant population.
                    if X0_i[-1] < 1 <= multi_resis_end_minimum:
                        index_minimum_end_multi_resis = np.argmin(multi_resis_end)
                        index_same_minimum_end_multi_resis = [i for i, x in enumerate(multi_resis_end) if x ==
                                                              multi_resis_end[index_minimum_end_multi_resis]]
                        # If there is only 1 d_i that gives the minimum multiply-resistant population.
                        if len(index_same_minimum_end_multi_resis) == 1:
                            index_select = index_minimum_end_multi_resis
                        # If there are multiple d_i that give the same minimum multiply-resistant populaiton, minimize the total cell population in
                        # all the d_i giving the same minimum multiply-resistant populaiton.
                        else:
                            X_i_end_total_same_minimum_end_multi_resis = [X_i_end_total[i] for i in index_same_minimum_end_multi_resis]
                            # Find the d_i gives the minimum total cell population in all the d_i giving the same minimum multiply-resistant
                            # population.
                            index_minimum_end_total_same_minimum_end_multi_resis = np.argmin(X_i_end_total_same_minimum_end_multi_resis)
                            index_select = index_same_minimum_end_multi_resis[index_minimum_end_total_same_minimum_end_multi_resis]
                        # If the minimize of multiply-resistant populaiton reach mortality, still minimize total cell population.
                        if X_i_end_total[index_select] >= Limit_mortality:
                            index_select = index_minimum_end_total
                # If current multiply-resistant >= 1 and multiply-resistant is not curable, minimize the total cell population.
                else:
                    pass
                # This means simulation of all of the dose combination reach mortality. Then select the longest survival time.
                if X_i_end_total[index_select] >= Limit_mortality:
                    index_select = np.argmax(t_sim[:, index])

        node_minimum = node[index_select, :]
        d_i = np.empty([PAR['Num_drug'], 0], dtype=float)
        for i_deep in node_minimum:
            d_i = np.append(d_i, np.tile(np.array([dose_combination[i_deep]]).T, (1, Stepsize)), 1)
        t_i = np.arange(0, node_minimum.shape[0] * Stepsize + Simtimestep, Simtimestep)
        X_i, t_i, d_i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)

        X_Strategy8_i = np.append(X_Strategy8_i, X_i, 1) if tnow == 0 else np.append(X_Strategy8_i, X_i[:, 1:], 1)
        t_Strategy8_i = np.append(t_Strategy8_i, t_i + tnow) if tnow == 0 else np.append(t_Strategy8_i, t_i[1:] + tnow)
        d_Strategy8_i = np.append(d_Strategy8_i, d_i, 1)

        # If total cell populatin is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i[:, -1] < 1):
            break
        # Update X0_i and tnow.
        X0_i = X_i[:, -1]
        tnow += node_minimum.shape[0] * Stepsize
    return t_Strategy8_i, X_Strategy8_i, d_Strategy8_i


# Define Strategy 9.
def DPM_strategy_9(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, Maxnum_subseq, subtreedepth, LSsim,
                   mis_specification, mis_PAR):
    # ALTO (Adaptive Long-Term Optimization) in the Biology Direct paper.
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow: int = 0

    # Gennerate the decision tree using the variable value of 'subtreedepth' as the value of depth.
    decision_tree_all = list(itertools.product(list(range(len(dose_combination))), repeat=subtreedepth))
    decision_tree_all.sort()
    X_init_total = [X0_i]
    node_pre_total = [np.empty([0], dtype=int)]
    # No terminate in strategy 9. Didn't terminate even if the total and multiply-resistant populations both exceed other nodes of the same depth.
    # Because we didn't use strategy 1-3, which only use total and multiply-resistant populations to decide which strategy will be selected.
    # Strategy 9 will select the dose combination giving the longest survival time, it could not be fully decided by the two values which are total
    # and multiply-resistant populations.
    node_whole, node_select, stop, terminate = None, None, None, 'no'
    while tnow < Simduration:
        depth = min(subtreedepth, -(-(Simduration - tnow)//Stepsize))
        tnow += depth * Stepsize
        timeleft = Simduration - tnow
        # Record the total population for each initial X0.
        total_pop_allinitial = np.empty(0, dtype=float)
        # Record the each population size at the ending time for each initial X0.
        each_pop_allinitial = np.empty([0, X0_i.shape[0]], dtype=float)
        # Record the simulation times for each used nodes for each initial X0.
        t_sim_allinitial = np.empty(0, dtype=float)
        # Record the used nodes for each initial X0.
        node_allinitial = np.empty([0, depth], dtype=int)
        # All previous nodes.
        node_pre = []
        for i in range(len(X_init_total)):
            X0_init_i = X_init_total[i]
            node_pre_i = node_pre_total[i]
            if mis_specification:
                node, total_pop, each_pop, t_sim = DPM_strategy_9_node(X0_init_i, mis_T_i, mis_g0_i, mis_Sa_i, dose_combination, Stepsize,
                                                                       Simtimestep, Limit_mortality, subtreedepth, LSsim, depth, decision_tree_all)

            else:
                node, total_pop, each_pop, t_sim = DPM_strategy_9_node(X0_init_i, T_i, g0_i, Sa_i, dose_combination, Stepsize, Simtimestep,
                                                                       Limit_mortality, subtreedepth, LSsim, depth, decision_tree_all)

            # Check the first dimension of the variable 'node', 'total_pop', 'each_pop' and 't_sim' are the same.
            if node.shape[0] != total_pop.shape[0] != each_pop.shape[0] != t_sim.shape[0]:
                print('Dimension error. The first dimension of the variable "node", "total_pop", "each_pop" and "t_sim" are not same.')
                DPM_print_errorandclose()
                sys.exit()
            node_allinitial = np.append(node_allinitial, node, 0)
            total_pop_allinitial = np.append(total_pop_allinitial, total_pop, 0)
            each_pop_allinitial = np.append(each_pop_allinitial, each_pop, 0)
            t_sim_allinitial = np.append(t_sim_allinitial, t_sim, 0)
            node_pre.extend([node_pre_i] * node.shape[0])
        if mis_specification:
            node_whole, X_init, stop = DPM_strategy_9_nodeselection(
                dose_combination, node_allinitial, total_pop_allinitial, each_pop_allinitial, t_sim_allinitial, node_pre, mis_g0_i, mis_T_i, mis_Sa_i,
                Limit_mortality, Simtimestep, Maxnum_subseq, timeleft)
        else:
            node_whole, X_init, stop = DPM_strategy_9_nodeselection(
                dose_combination, node_allinitial, total_pop_allinitial, each_pop_allinitial, t_sim_allinitial, node_pre, g0_i, T_i, Sa_i,
                Limit_mortality, Simtimestep, Maxnum_subseq, timeleft)

        if stop:
            node_select = node_whole
            break
        X_init_total = X_init
        node_pre_total = node_whole

    d_i = np.empty([PAR['Num_drug'], 0], dtype=int)
    for i_deep in node_select.astype(int):
        d_i = np.append(d_i, np.tile(np.array([dose_combination[i_deep]]).T, (1, Stepsize)), 1)
    t_i = np.arange(0, node_select.shape[0] * Stepsize + Simtimestep, Simtimestep)
    X_Strategy9_i, t_Strategy9_i, d_Strategy9_i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)

    return t_Strategy9_i, X_Strategy9_i, d_Strategy9_i


# Calculate the node results of Strategy 9.
def DPM_strategy_9_node(X0_i, T_i, g0_i, Sa_i, dose_combination, Stepsize, Simtimestep, Limit_mortality, subtreedepth, LSsim, depth,
                        decision_tree_all):
    if depth < subtreedepth:
        decision_tree = list(itertools.product(list(range(len(dose_combination))), repeat=depth))
        decision_tree.sort()
    else:
        decision_tree = decision_tree_all
    # Record the total populations for each used nodes.
    total_pop = np.empty(0, dtype=float)
    # Record the each population size at the ending time.
    each_pop = np.empty([0, X0_i.shape[0]], dtype=float)
    # Record the simulation times for each used nodes.
    t_sim = np.empty(0, dtype=float)
    # Record the used nodes.
    node = np.empty([0, depth], dtype=int)
    for index in range(len(decision_tree)):
        # Select one node.
        i_node = decision_tree[index]

        d_i = np.empty([len(dose_combination[0]), 0], dtype=float)
        for i_deep in i_node:
            d_i = np.append(d_i, np.tile(np.array([dose_combination[i_deep]]).T, (1, Stepsize)), 1)
        t_i = np.arange(0, len(i_node) * Stepsize + Simtimestep, Simtimestep)
        X_node_i, t_slow_i, _ = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, Limit_mortality, LSsim)
        t_node_i = t_slow_i[-1]

        # # Approximate solution.
        # pos_changedose = DPM_treatment_change_time(np.array(i_node))
        # tnow_node = 0
        # X0_constant = X0_i
        # for i in range(len(pos_changedose)):
        #     i_constant = i_node[pos_changedose[i]]
        #     if pos_changedose[i] != pos_changedose[-1]:
        #         i_duration = pos_changedose[i + 1] - pos_changedose[i]
        #     else:
        #         i_duration = len(i_node) - pos_changedose[i]
        #     d_i = np.tile(np.array([dose_combination[i_constant]]).T, (1, Stepsize * i_duration))
        #     t_i = np.arange(0, i_duration * Stepsize + Simtimestep, Simtimestep)
        #     X_node_constant_i, t_node_constant_i = DPM_strategy_9_approxsol(
        #     d_i, t_i, X0_constant, T_i, g0_i, Sa_i, Simtimestep, Limit_mortality, LSsim)
        #     tnow_node += t_node_constant_i[-1]
        #     X0_constant = X_node_constant_i[:, -1]
        #     if X_node_constant_i[:, -1].sum() >= Limit_mortality:
        #         break
        # X_node_i = X_node_constant_i
        # t_node_i = tnow_node

        node = np.append(node, np.reshape(np.array(i_node, dtype=int), (1, depth)), 0)
        total_pop = np.append(total_pop, X_node_i[:, -1].sum())
        each_pop = np.append(each_pop, np.reshape(X_node_i[:, -1], (1, X_node_i.shape[0])), 0)
        t_sim = np.append(t_sim, t_node_i)
    return node, total_pop, each_pop, t_sim


# Run node selection procedure of Strategy 9.
def DPM_strategy_9_nodeselection(dose_combination, node, total_pop, each_pop, t_sim, node_pre, g0_i, T_i, Sa_i, Limit_mortality, Simtimestep,
                                 Maxnum_subseq, timeleft):
    stop, node_whole, X_init = False, None, None
    # If all dose the combination reach mortality before or at the last step. Select the longest survival time.
    if total_pop.min() >= Limit_mortality:
        index_select = np.argmax(t_sim)
        node_pre_select = node_pre[index_select]
        node_select = node[index_select]
        node_whole = np.append(node_pre_select, node_select)
        stop = True
        return node_whole, None, stop
    # If cured happens, stop.
    elif any([all(i < 1) for i in each_pop]):
        cured_truefalse = [all(i < 1) for i in each_pop]
        cured_index = [i for i, x in enumerate(cured_truefalse) if x]
        # Select any cured node because they all cured.
        index_select = cured_index[0]
        node_pre_select = node_pre[index_select]
        node_select = node[index_select]
        node_whole = np.append(node_pre_select, node_select)
        stop = True
        return node_whole, None, stop
    # If time left is 0, survive at the end of the simulation, any node is fine.
    elif timeleft <= 0:
        node_whole = [np.concatenate((node_pre[i], node[i])) for i in range(len(node_pre))]
        node_whole = node_whole[randrange(len(node_whole))]
        stop = True
        return node_whole, None, stop
    # Not all nodes reach mortality or cured.
    else:
        # Because not all nodes reach mortality, first exclude the nodes that reach mortality.
        if np.greater_equal(total_pop, Limit_mortality).any():
            mortality_truefalse = [i < Limit_mortality for i in total_pop]
            node, total_pop, each_pop = node[mortality_truefalse, :], total_pop[mortality_truefalse], each_pop[mortality_truefalse, :]
            t_sim, node_pre = t_sim[mortality_truefalse], [node_pre[i] for i, x in enumerate(mortality_truefalse) if x]
        # Check if there are nodes that have the same cell numbers of each cell type population, delete the duplicate nodes and keep the nodes unique.
        _, unique_index, _ = np.unique(each_pop, return_index=True, return_counts=True, axis=0)
        if len(unique_index) != each_pop.shape[0]:
            unique_index.sort()
            node, total_pop, each_pop = node[unique_index, :], total_pop[unique_index], each_pop[unique_index, :]
            t_sim, node_pre = t_sim[unique_index], [node_pre[i] for i in unique_index]

        # Variable 'X_l': the total populations at Tmax (Simduration, e.g., 1800 days) by administering full dosage of each drug simultaneously
        # (although this might not be allowed due to toxicity) for each nodes.
        # Variable 'X_u': the minimum total populations at the maximum monitoring time Tmax (Simduration, e.g., 1800 days) over all valid static
        # treatments for each nodes.
        # Variable 'geometric_mean_XlXu': geometric means of variable 'X_l' and 'X_u' for each nodes.
        X_l, X_u = np.full(node.shape[0], -1, dtype=float), np.full(node.shape[0], -1, dtype=float)
        geometric_mean_XlXu = np.full(node.shape[0], -1, dtype=float)
        t_i = np.arange(0, timeleft + Simtimestep, Simtimestep)
        d_fulldose = np.full([len(dose_combination[0]), len(t_i) - 1], 1, dtype=float)
        # Make variable 'maxthreshold' large and set variable 'LSsim' to 'False'.
        # The simulation won't stop if it reaches the 'maxthreshold' value or cure happens.
        maxthreshold = Limit_mortality * 1e100
        LSsim = False
        for i in range(each_pop.shape[0]):
            X0_i = each_pop[i, :]
            # X_fulldose_i, t_fulldose_i = DPM_run_approx_sol(d_fulldose, t_i, X0_i, T_i, g0_i, Sa_i, Simtimestep, maxthreshold, LSsim)

            X_fulldose_i, t_fulldose_i, _ = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_fulldose, t_i, maxthreshold, LSsim)

            X_l[i] = X_fulldose_i[:, -1].sum()
            X_staticdose_i = np.inf
            for i_dose_combination in dose_combination:
                treat_i = np.array([i_dose_combination], dtype=float).T
                d_static_i = np.tile(treat_i, (1, t_i.shape[0] - 1))

                # X_static_i, t_static_i = DPM_run_approx_sol(d_static_i, t_i, X0_i, T_i, g0_i, Sa_i, Simtimestep, maxthreshold, LSsim)
                X_static_i, t_static_i, _ = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_static_i, t_i, maxthreshold, LSsim)
                X_static_i, t_static_i, _ = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_static_i, t_i, maxthreshold, LSsim)

                X_staticdose_i = min(X_staticdose_i, X_static_i[:, -1].sum())
                # print(X_static_i[:, -1].sum())

            X_u[i] = X_staticdose_i
            geometric_mean_XlXu[i] = (np.log10(X_l[i])+np.log10(X_u[i])) * (1 / 2)
            pass
        # Select nodes based on two categorys.
        # (1) Based on the proof from Biology Direct paper, d_i yields a longer survival time than d_j if X_l_j > X_u_i and X_l_j > Limit_mortality.
        # Find the X_l that exceed the upper bound of any other nodes and the value of itself is larger than'Limit_mortality'.
        index_remove = (X_l > min(X_u)) * (X_l > Limit_mortality)
        node = node[np.logical_not(index_remove), :]
        each_pop = each_pop[np.logical_not(index_remove), :]
        node_pre = [node_pre[i] for i, x in enumerate(np.logical_not(index_remove)) if x]
        if node.shape[0] > Maxnum_subseq:
            geometric_mean_XlXu = geometric_mean_XlXu[np.logical_not(index_remove)]
            index_sort_geometric_mean_XlXu = np.argsort(geometric_mean_XlXu)
            node = node[index_sort_geometric_mean_XlXu[:Maxnum_subseq], :]
            each_pop = each_pop[index_sort_geometric_mean_XlXu[:Maxnum_subseq], :]
            node_pre = [node_pre[i] for i in index_sort_geometric_mean_XlXu[:Maxnum_subseq]]
        # connect nodes.
        node_whole = [np.concatenate((node_pre[i], node[i])) for i in range(len(node_pre))]
        X_init = list(each_pop)
    return node_whole, X_init, stop


# # Run approximation solution used in Strategy 9.
# def DPM_strategy_9_approxsol(d_input, t_input, X0_i, T_i, g0_i, Sa_i, Simtimestep, Limit_mortality, LSsim):
#     d_cut = d_input[:, 0:1]
#     t_cut = t_input[[0, -1]]
#     X_cut_i, t_cut_i, _ = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_cut, t_cut, Limit_mortality, LSsim, Simtimestep_control=False)
#     # Find any cell type change from smaller than 1 to bigger than 1.
#     index_celltype_alter = [True if X0_i[i] < 1 <= X_cut_i[i, -1] else False for i in range(len(X0_i))]
#     X0_alter = X0_i
#     t_alter = t_input
#     d_alter = d_input
#     stop = False
#     tnow = 0
#     X_out = t_out = X_alter_i = t_alter_i = None
#     if any(index_celltype_alter):
#         while any(index_celltype_alter):
#             X_alter_i, t_alter_i, _ = DPM_sim_model(X0_alter, T_i, g0_i, Sa_i, d_alter, t_alter, Limit_mortality, LSsim, True, None, None,
#                                                     index_celltype_alter=index_celltype_alter, index_celltype_alter_type=1)
#
#             # If reach mortality, the stop the simulation and return output.
#             if X_alter_i[:, -1].sum() >= Limit_mortality and LSsim:
#                 X_out = X_alter_i
#                 t_out = t_alter_i + tnow
#                 stop = True
#                 break
#             index_celltype_alter_i = [False if X0_alter[i] < 1 <= X_alter_i[i, -1] else True for i in range(len(X0_alter))]
#             index_celltype_alter = [bool(index_celltype_alter * index_celltype_alter_i) for index_celltype_alter, index_celltype_alter_i
#                                     in zip(index_celltype_alter, index_celltype_alter_i)]
#             t_alter = np.arange(0, t_input[-1] - t_alter_i[-1] + Simtimestep, Simtimestep)
#             d_alter = np.tile(d_cut, (1, t_alter.shape[0] - 1))
#             X0_alter = X_alter_i[:, -1]
#             tnow += t_alter_i[-1]
#         if not stop:
#             if t_alter.size != 0:
#                 X_out, t_out, _ = DPM_sim_model(
#                 X0_alter, T_i, g0_i, Sa_i, d_cut, t_alter[[0, -1]], Limit_mortality, LSsim, Simtimestep_control=False)
#                 t_out += tnow
#                 # If reach mortality, then re-simulate using the default Simtimestep in order to get the exact time when mortality happens.
#                 if X_out[:, -1].sum() >= Limit_mortality and LSsim:
#                     d_alter = np.tile(d_cut, (1, t_alter.shape[0] - 1))
#                     X_out, t_out, _ = DPM_sim_model(X0_alter, T_i, g0_i, Sa_i, d_alter, t_alter, Limit_mortality, LSsim)
#                     t_out += tnow
#             else:
#                 X_out = X_alter_i
#                 t_out = t_alter_i + tnow
#     else:
#         if X_cut_i[:, -1].sum() >= Limit_mortality and LSsim:
#             X_out, t_out, _ = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_input, t_input, Limit_mortality, LSsim)
#         else:
#             X_out = X_cut_i
#             t_out = t_cut_i
#
#     return X_out, t_out


# Define recursive simulation used in multistep Strategies.
def DPM_strategy_recurse_multistep(X0_i, T_i, g0_i, Sa_i, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, lookahead_step,
                                   maxthreshold, LSsim, tnow, terminate, decision_tree_all):
    depth = min(lookahead_step, (Simduration - tnow) / Stepsize)
    if depth < lookahead_step:
        decision_tree = list(itertools.product(list(range(len(dose_combination))), repeat=depth))
        decision_tree.sort()
    else:
        decision_tree = decision_tree_all
    # Record the total populations for each used nodes.
    total_pop = np.empty([0, depth], dtype=float)
    # Record the multiply-resistant populations for each used nodes.
    multi_resis = np.empty([0, depth], dtype=float)
    # Record the each population size at the ending time.
    each_pop = np.empty([0, X0_i.shape[0], depth], dtype=float)
    # Record the simulation times for each used nodes.
    t_sim = np.empty([0, depth], dtype=float)
    # Record the used nodes.
    node = np.empty([0, depth], dtype=int)
    index = 0
    while index < len(decision_tree):
        # Select one node.
        i_node = decision_tree[index]
        # Initialize the total population for this node. All equal to -1 for convenient.
        i_total_pop_depth = np.full(depth, -1, dtype=float)
        # Initialize the multiply-resistant populations for this node. All equal to -1 for convenient.
        i_multi_resis_depth = np.full(depth, -1, dtype=float)
        # Initialize each population size at the ending time.
        i_eachpop = np.full([X0_i.shape[0], depth], -1, dtype=float)
        # Initialize the simulation times for this node. All equal to -1 for convenient.
        i_t_sim = np.full(depth, -1, dtype=float)
        node, total_pop, multi_resis, each_pop, t_sim, i_tree_depthnow = DPM_strategy_recurse_multistep_run(
            X0_i, T_i, g0_i, Sa_i, i_node, i_total_pop_depth, i_multi_resis_depth, i_eachpop, i_t_sim, dose_combination, Stepsize, Simtimestep,
            Limit_mortality, maxthreshold, LSsim, total_pop, multi_resis, node, each_pop, t_sim, terminate, depth, depthnow=0)
        if i_tree_depthnow is not None:
            decision_tree = DPM_generate_truncate_decisiontree(decision_tree, i_tree_depthnow, index)
        else:
            index += 1
    return node, total_pop, multi_resis, each_pop, t_sim


# Run recursive simulation used in multistep Strategies.
def DPM_strategy_recurse_multistep_run(X0_i, T_i, g0_i, Sa_i, i_node, i_total_pop_depth, i_multi_resis_depth, i_eachpop, i_t_sim, dose_combination,
                                       Stepsize, Simtimestep, Limit_mortality, maxthreshold, LSsim, total_pop, multi_resis, node, each_pop, t_sim,
                                       terminate, depth, depthnow):
    if depthnow > depth - 1:
        total_pop = np.append(total_pop, np.reshape(i_total_pop_depth, (1, depth)), 0)
        multi_resis = np.append(multi_resis, np.reshape(i_multi_resis_depth, (1, depth)), 0)
        node = np.append(node, np.reshape(np.array(i_node, dtype=int), (1, depth)), 0)
        each_pop = np.append(each_pop, np.reshape(i_eachpop, (1, i_eachpop.shape[0], i_eachpop.shape[1])), 0)
        t_sim = np.append(t_sim, np.reshape(np.array(i_t_sim, dtype=float), (1, depth)), 0)
        return node, total_pop, multi_resis, each_pop, t_sim, None

    i_dose_combination = dose_combination[i_node[depthnow]]
    treat_i = np.array([i_dose_combination], dtype=float).T
    t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
    d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
    X_i, t_i, d_i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)

    # Total populations.
    total_pop_i = X_i[:, -1].sum()
    # Multiply-resistant populations.
    multi_resis_i = X_i[-1, -1]

    if total_pop.any() and multi_resis.any():
        total_pop_samedepth = total_pop[:, depthnow]
        multi_resis_samedepth = multi_resis[:, depthnow]
        if terminate == 'total_pop':
            if np.less(total_pop_samedepth, total_pop_i).any() and not np.greater(total_pop_samedepth, Limit_mortality).all():
                i_tree_depthnow = i_node[:depthnow + 1]
                return node, total_pop, multi_resis, each_pop, t_sim, i_tree_depthnow
        elif terminate == 'multi_resis':
            if np.less(multi_resis_samedepth, multi_resis_i).any() and not np.greater(total_pop_samedepth, Limit_mortality).all():
                i_tree_depthnow = i_node[:depthnow + 1]
                return node, total_pop, multi_resis, each_pop, t_sim, i_tree_depthnow
        elif terminate == 'both':
            if (np.less(total_pop_samedepth, total_pop_i) * np.less(multi_resis_samedepth, multi_resis_i)).any() and \
                    not np.greater(total_pop_samedepth, Limit_mortality).all():
                i_tree_depthnow = i_node[:depthnow + 1]
                return node, total_pop, multi_resis, each_pop, t_sim, i_tree_depthnow
        else:
            pass

    # Stop if reach motality.
    stop = False
    if total_pop_i >= Limit_mortality:
        i_total_pop_depth[depthnow:] = total_pop_i
        i_total_pop_depth[depthnow + 1:] = np.inf
        i_multi_resis_depth[depthnow:] = multi_resis_i
        i_multi_resis_depth[depthnow + 1:] = np.inf
        i_eachpop[:, depthnow] = X_i[:, -1]
        i_eachpop[:, depthnow + 1:] = np.inf
        i_t_sim[depthnow] = t_i[-1]
        stop = True
    # Stop if cured.
    if all(X_i[:, -1] < 1):
        i_total_pop_depth[depthnow:] = total_pop_i
        i_total_pop_depth[depthnow + 1:] = 0
        i_multi_resis_depth[depthnow:] = multi_resis_i
        i_multi_resis_depth[depthnow + 1:] = 0
        i_eachpop[:, depthnow] = X_i[:, -1]
        i_eachpop[:, depthnow + 1:] = 0
        i_t_sim[depthnow] = t_i[-1]
        stop = True
    if stop:
        total_pop = np.append(total_pop, np.reshape(i_total_pop_depth, (1, depth)), 0)
        multi_resis = np.append(multi_resis, np.reshape(i_multi_resis_depth, (1, depth)), 0)
        node = np.append(node, np.reshape(np.array(i_node, dtype=int), (1, depth)), 0)
        each_pop = np.append(each_pop, np.reshape(i_eachpop, (1, i_eachpop.shape[0], i_eachpop.shape[1])), 0)
        t_sim = np.append(t_sim, np.reshape(np.array(i_t_sim, dtype=float), (1, depth)), 0)
        return node, total_pop, multi_resis, each_pop, t_sim, None

    i_total_pop_depth[depthnow] = total_pop_i
    i_multi_resis_depth[depthnow] = multi_resis_i
    i_eachpop[:, depthnow] = X_i[:, -1]
    i_t_sim[depthnow] = t_i[-1]

    depthnow += 1
    X0_i = X_i[:, -1]
    node, total_pop, multi_resis, each_pop, t_sim, i_tree_depthnow = DPM_strategy_recurse_multistep_run(
        X0_i, T_i, g0_i, Sa_i, i_node, i_total_pop_depth, i_multi_resis_depth, i_eachpop, i_t_sim, dose_combination, Stepsize, Simtimestep,
        Limit_mortality, maxthreshold, LSsim, total_pop, multi_resis, node, each_pop, t_sim, terminate, depth, depthnow)
    return node, total_pop, multi_resis, each_pop, t_sim, i_tree_depthnow


# Define optimal strategy
def DPM_strategy_optimal(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, mis_specification, mis_PAR):
    # Optimal strategy by Dymos package
    PARopt = PAR if not mis_PAR else mis_PAR
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PARopt), DPM_generate_g0(PARopt), DPM_generate_Sa(PARopt), DPM_generate_T(PARopt)
    import openmdao.api as om
    import dymos as dm

    class DPMstrategyoptimalmodel(om.ExplicitComponent):
        """
        The is the same ODE model used in 2drug case, redefine here to use the Dymos package
        """
        def initialize(self):
            self.options.declare('num_nodes', types=int)
            self.options.declare('g0')
            self.options.declare('Sa')
            self.options.declare('T')

        def setup(self):
            nn = self.options['num_nodes']
            g0 = self.options['g0']
            Sa = self.options['Sa']
            T = self.options['T']
            Num_state = Sa.shape[0]

            # Input state variables
            # self.add_input('time', shape=(nn,), units='d')
            self.add_input('S', shape=(nn,), units='N')
            self.add_input('R1', shape=(nn,), units='N')
            self.add_input('R2', shape=(nn,), units='N')
            self.add_input('R12', shape=(nn,), units='N')
            if Num_state == 3:
                self.add_input('R3', shape=(nn,), units='N')
                self.add_input('R13', shape=(nn,), units='N')
                self.add_input('R23', shape=(nn,), units='N')
                self.add_input('R123', shape=(nn,), units='N')

            # Input parameters
            # self.add_input('d', shape=(nn,))   #
            self.add_input('d', shape=1, tags=['dymos.static_target'])

            # Outputs
            self.add_output('Sdot', shape=(nn,), units='N/d', tags=['dymos.state_rate_source:S'])
            self.add_output('R1dot', shape=(nn,), units='N/d', tags=['dymos.state_rate_source:R1'])
            self.add_output('R2dot', shape=(nn,), units='N/d', tags=['dymos.state_rate_source:R2'])
            self.add_output('R12dot', shape=(nn,), units='N/d', tags=['dymos.state_rate_source:R12'])
            if Num_state == 3:
                self.add_output('R3dot', val=np.zeros(nn), tags=['dymos.state_rate_source:R3'])
                self.add_output('R13dot', val=np.zeros(nn), tags=['dymos.state_rate_source:R13'])
                self.add_output('R23dot', val=np.zeros(nn), tags=['dymos.state_rate_source:R23'])
                self.add_output('R123dot', val=np.zeros(nn), tags=['dymos.state_rate_source:R123'])

            # self.add_output('total', shape=(nn,))

        def compute(self, inputs, outputs):
            nn = self.options['num_nodes']
            # time = inputs['time']
            # pprint(time)
            g0, Sa, T = self.options['g0'], self.options['Sa'], self.options['T']
            Num_state = Sa.shape[0]

            d1 = inputs['d']

            S_op, R1_op, R2_op, R12_op = inputs['S'], inputs['R1'], inputs['R2'], inputs['R12']
            X = np.array((S_op, R1_op, R2_op, R12_op))

            # d = np.append(d1.reshape((1, len(d1))), 1.0 - d1.reshape((1, len(d1))), 0)

            I_identity = np.identity(Num_state)
            Sdot, R1dot, R2dot, R12dot = np.zeros(nn), np.zeros(nn), np.zeros(nn), np.zeros(nn)

            # Mask the terms involved in the initial populations below 1.
            # for i_node in range(self.options['num_nodes']):
            #     x = np.array((S[i_node], R1[i_node], R2[i_node], R12[i_node]))
            #     mask = np.ones((Num_state, Num_state), dtype=float)
            #     d_i_node = d[:, i_node]
            #     for i in range(Num_state):
            #         mask[:, i] = 0 if x[i] < 1 else 1
            #     f = (I_identity + T) @ np.diag(g0) - np.diag(Sa @ d_i_node)
            #     # f = f * mask
            #     dx = f @ x
            #     Sdot[i_node], R1dot[i_node], R2dot[i_node], R12dot[i_node] = dx[0], dx[1], dx[2], dx[3]

            d_op = np.array([inputs['d'], 1-inputs['d']]).reshape((Sa.shape[1],))
            f = (I_identity + T) @ np.diag(g0) - np.diag(Sa @ d_op)
            dx = f @ X

            # pprint(S)
            outputs['Sdot'] = dx[0, :]  # Sdot
            outputs['R1dot'] = dx[1, :]  # R1dot
            outputs['R2dot'] = dx[2, :]  # R2dot
            outputs['R12dot'] = dx[3, :]  # R12dot
            # outputs['total'] = total

    p = om.Problem(model=om.Group())

    optimizer = 'IPOPT'
    p.driver = om.pyOptSparseDriver(print_results=False)
    p.driver.options['optimizer'] = optimizer
    p.driver.opt_settings['max_iter'] = int(1e4)
    if optimizer == 'IPOPT':
        p.driver.opt_settings['print_level'] = 4
    p.driver.declare_coloring()

    # p.model.add_subsystem('sum', om.ExecComp('X = S + R1 + R2 + R12'))

    traj = p.model.add_subsystem('traj', dm.Trajectory())
    # transcription = dm.Radau(num_segments=17, order=3, compressed=True)
    transcription = dm.GaussLobatto(num_segments=17, order=3, compressed=True)
    # transcription = dm.ExplicitShooting(num_segments=20, num_steps_per_segment=20, method='rk4')

    """ phase0 """
    phase0 = dm.Phase(ode_class=DPMstrategyoptimalmodel, transcription=transcription, ode_init_kwargs={'g0': g0_i, 'Sa': Sa_i, 'T': T_i})
    phase0 = traj.add_phase('phase0', phase0)
    # add time
    phase0.set_time_options(fix_initial=True, fix_duration=True, units='d')
    # add state
    phase0.add_state('S', val=X0_i[0], fix_initial=True, lower=0, units='N', rate_source='Sdot')
    phase0.add_state('R1', val=X0_i[1], fix_initial=True, lower=0, units='N', rate_source='R1dot')
    phase0.add_state('R2', val=X0_i[2], fix_initial=True, lower=0, units='N', rate_source='R2dot')
    phase0.add_state('R12', val=X0_i[3], fix_initial=True, lower=0, units='N', rate_source='R12dot')
    # add control variable
    phase0.add_parameter('d', val=0, opt=True, lower=0, upper=1, targets='d', static_target=True)
    # phase0.add_control('d', continuity=False, rate_continuity=False, lower=0, upper=1)
    # add constraint
    # phase0.add_path_constraint('total', lower=0, upper=Limit_mortality)
    # add objective
    # phase0.add_objective('time', loc='final', scaler=-1)
    phase0.add_objective('S', loc='final', scaler=-1)
    # p.model.add_objective('sum.X')

    # Connections
    # p.model.connect('traj.phase0.states:S', 'sum.S')
    # p.model.connect('traj.phase0.states:R1', 'sum.R1')
    # p.model.connect('traj.phase0.states:R2', 'sum.R2')
    # p.model.connect('traj.phase0.states:R12', 'sum.R12')

    # Setup the Problem
    p.setup()

    p['traj.phase0.t_initial'] = 0.0
    p['traj.phase0.t_duration'] = Stepsize

    p.set_val('traj.phase0.states:S', val=X0_i[0])
    p.set_val('traj.phase0.states:R1', val=X0_i[1])
    p.set_val('traj.phase0.states:R2', val=X0_i[2])
    p.set_val('traj.phase0.states:R12', val=X0_i[3])
    # p.set_val('traj.phase0.controls:d', val=0)
    p.set_val('traj.phase0.parameters:d', val=0)

    # phase0.set_refine_options(tol=1.0E-7)  # refine_iteration_limit=3
    dm.run_problem(p, simulate=True)
    t = p.get_val('traj.phase0.timeseries.time')
    # d = p.get_val('traj.phase0.timeseries.controls:d')
    d = p.get_val('traj.phase0.timeseries.parameters:d')
    total = p.get_val('traj.phase0.timeseries.total')
    plt.plot(t, d)
    plt.ylim([0, 1])

    plt.plot(t, total)
    plt.ylim([0, 1e10])

    S = p.get_val('traj.phase0.timeseries.states:S')
    plt.plot(t, S)
    R1 = p.get_val('traj.phase0.timeseries.states:R1')
    plt.plot(t, R1)
    R2 = p.get_val('traj.phase0.timeseries.states:R2')
    plt.plot(t, R2)
    R12 = p.get_val('traj.phase0.timeseries.states:R12')
    plt.plot(t, R12)

    from dymos.examples.plotting import plot_results
    def vanderpol_plots():
        sol = om.CaseReader('dymos_solution.db').get_case('final')
        sim = om.CaseReader('dymos_simulation.db').get_case('final')

        plot_results([('traj.phase0.timeseries.time',
                       'traj.phase0.timeseries.states:S',
                       'time (d)',
                       'S (#)'),
                      ('traj.phase0.timeseries.time',
                       'traj.phase0.timeseries.states:R1',
                       'time (d)',
                       'R1 (#)'),
                      ('traj.phase0.timeseries.time',
                       'traj.phase0.timeseries.states:R2',
                       'time (d)',
                       'R2 (#)'),
                      ('traj.phase0.timeseries.time',
                       'traj.phase0.timeseries.states:R12',
                       'time (d)',
                       'R12'),
                      ('traj.phase0.timeseries.time',
                       'traj.phase0.timeseries.parameters:d',
                       'time (d)',
                       'd'),
                      # ('traj.phase0.timeseries.time',
                      #  'traj.phase0.timeseries.total',
                      #  'time (d)',
                      #  'total (#)')
                       ],
                     title='Van Der Pol Simulation',
                     p_sol=sol, p_sim=sim)
    vanderpol_plots()

