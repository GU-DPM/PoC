from DPM_assign_check import *
'''
CHECKED
'''


def DPM_strategy_select_drug_majority(x, Num_drug, Sa):
    treat = np.zeros([Num_drug, 1], dtype=int)
    # Index of the predominant cell type.#
    ind_drug = np.argmax(Sa[np.argmax(x), :])
    treat[ind_drug] = 1
    return treat, ind_drug


# Define CPM.#
def strategy_CPM(PAR, Simduration, Stepsize, Simtimestep, Limit_mortality, Limit_radiologicdetection, LSsim, misspecification_atdecision,
                 misspecification_atsim, mis_PAR):
    # CPM: Current personalized medicine.#
    maxthreshold = Limit_mortality
    X0_i = DPM_generate_X0(PAR)
    g0_i = DPM_generate_g0(PAR)
    Sa_i = DPM_generate_Sa(PAR)
    T_i = DPM_generate_T(PAR)
    mis_Sa_i = DPM_generate_Sa(mis_PAR) if misspecification_atdecision or misspecification_atsim else None

    nadir = X0_i.sum()
    tnow = 0.0
    t_CPM_i = np.zeros([0], dtype=float)
    X_CPM_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_CPM_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    drug_used = np.zeros([PAR['Num_drug']], dtype=int)

    # Identify the drug to be used.#
    if misspecification_atdecision:
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

        X_CPM_i = np.append(X_CPM_i, X_i, 1) if tnow == 0 else np.append(X_CPM_i, X_i[:, 1:], 1)
        t_CPM_i = np.append(t_CPM_i, t_i + tnow) if tnow == 0 else np.append(t_CPM_i, t_i[1:] + tnow)
        d_CPM_i = np.append(d_CPM_i, d_i, 1)

        # If total cell population reaches or exceeds 'Limit_mortality', stop the simulation. Mortality occurs.#
        # Or if the cell count for all types falls below 1, the patient is cured. Stop simulation.#
        if (X_i[:, -1].sum() >= Limit_mortality) or np.all(X_i[:, -1] < 1):
            break

        # Switch to an unused drug if either:#
        # (1) Total cell population exceeds 2× nadir and exceeds the 'Limit_radiologicdetection' with unused drug remaining, or #
        # (2) Total cell population reemerges above the detection limit after falling below it, with unused drugs remaining.#
        # Note: Each drug is used only once.#
        if ((X_i[:, -1].sum() >= 2 * nadir and X_i[:, -1].sum() >= Limit_radiologicdetection) or
                (X_i[:, 0].sum() < Limit_radiologicdetection <= X_i[:, -1].sum())) and np.any(drug_used == 0):
            treat_i = np.zeros([PAR['Num_drug'], 1], dtype=int)
            index_drug_not_used = [i for i, x in enumerate(drug_used.tolist()) if x == 0]
            treat_i[index_drug_not_used] = 1
            drug_used[index_drug_not_used] = 1
            nadir = X_i[:, -1].sum()

        # Update nadir if total cells fall below current nadir but remain above detection limit.#
        if Limit_radiologicdetection <= X_i[:, -1].sum() < nadir:
            nadir = X_i[:, -1].sum()
        # Update 'X0_i' and 'tnow'.#
        X0_i = X_i[:, -1]
        tnow += Stepsize
    return t_CPM_i, X_CPM_i, d_CPM_i


# Define DPM1.#
def strategy_DPM1(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, mis_specification, mis_PAR):
    # DPM1: Minimize the total cell population by selecting the drug ('d_i') that achieves the greatest reduction at each time step.#
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_DPM1_i = np.zeros([0], dtype=float)
    X_DPM1_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_DPM1_i = np.zeros([PAR['Num_drug'], 0], dtype=float)

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

        # Identify the drug 'd_i' that yields the lowest total cell population.#
        if np.argmin(X_i_end_total) < Limit_mortality:
            index_select = np.argmin(X_i_end_total)
        else:
            # All dose combinations result in mortality. Select the combination with the longest survival time.#
            index_select = np.argmax(t_sim)
        d_i_minimum = d_i_total[int(index_select)]

        if mis_specification:
            X_i_minimum, t_i_minimum, d_i_minimum = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, LSsim)
        else:
            X_i_minimum, t_i_minimum = X_i_total[int(index_select)], t_i_total[int(index_select)]

        X_DPM1_i = np.append(X_DPM1_i, X_i_minimum, 1) if tnow == 0 else np.append(X_DPM1_i, X_i_minimum[:, 1:], 1)
        t_DPM1_i = np.append(t_DPM1_i, t_i_minimum + tnow) if tnow == 0 else np.append(t_DPM1_i, t_i_minimum[1:] + tnow)
        d_DPM1_i = np.append(d_DPM1_i, d_i_minimum, 1)

        # If total cell population >= 'Limit_mortality', stop simulation. Patient mortality occurs.#
        if X_i_minimum[:, -1].sum() >= Limit_mortality:
            break
        # If the cell count for all types falls below 1, the patient is cured. Stop simulation.#
        elif all(X_i_minimum[:, -1] < 1):
            break
        # Update 'X0_i' and 'tnow'.#
        X0_i = X_i_minimum[:, -1]
        tnow += Stepsize
    return t_DPM1_i, X_DPM1_i, d_DPM1_i


# Define DPM2.#
def strategy_DPM2(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, Strategy2threshold,
                  misspecification_ofdecision, misspecification_atsim, mis_PAR):
    # DPM2: Minimize the risk of incurable cells developing unless mortality is imminent.#
    # DPM2.1: threshold is 1e9.#
    # DPM2.2: threshold is 1e11.#
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
    t_DPM2_i = np.zeros([0], dtype=float)
    X_DPM2_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_DPM2_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
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

        # If some doses cause mortality and others don't, exclude the dose(s) that result in mortality.#
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
            # If the current total cell population is below the threshold, focus on minimizing the multi-resistant population.#
            if X0_i.sum() <= Strategy2threshold:
                index_minimum_end_multi_resis = np.argmin(X_i_end_multi_resis)
                index_same_minimum_end_multi_resis = [i for i, x in enumerate(X_i_end_multi_resis) if x ==
                                                      X_i_end_multi_resis[int(index_minimum_end_multi_resis)]]
                # If there is only one 'd_i' that yields the minimum multi-resistant population.#
                if len(index_same_minimum_end_multi_resis) == 1:
                    index_select = index_minimum_end_multi_resis
                # If multiple 'd_i' yield the same minimum multi-resistant population.#
                else:
                    X_i_end_total_same_minimum_end_multi_resis = [X_i_end_total[i] for i in index_same_minimum_end_multi_resis]
                    # Among the 'd_i' values that yield the same minimum multi-resistant population, select the one that produces the lowest#
                    # total cell population.#
                    index_minimum_end_total_same_minimum_end_multi_resis = np.argmin(X_i_end_total_same_minimum_end_multi_resis)
                    index_select = index_same_minimum_end_multi_resis[int(index_minimum_end_total_same_minimum_end_multi_resis)]
            # If the current total cell population exceeds the threshold, prioritie minimizing the total population.#
            else:
                index_select = np.argmin(X_i_end_total)
        # If every dose option leads to mortality, choose the one that maximizes survival time.#
        else:
            index_select = np.argmax(t_sim)

        d_i_minimum = d_i_total[index_select]
        if misspecification_ofdecision or misspecification_atsim:
            X_i_minimum, t_i_minimum, d_i_minimum = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, LSsim)
        else:
            X_i_minimum, t_i_minimum = X_i_total[index_select], t_i_total[index_select]

        X_DPM2_i = np.append(X_DPM2_i, X_i_minimum, 1) if tnow == 0 else np.append(X_DPM2_i, X_i_minimum[:, 1:], 1)
        t_DPM2_i = np.append(t_DPM2_i, t_i_minimum + tnow) if tnow == 0 else np.append(t_DPM2_i, t_i_minimum[1:] + tnow)
        d_DPM2_i = np.append(d_DPM2_i, d_i_minimum, 1)

        # If total cell population >= 'Limit_mortality', stop simulation. Patient mortality occurs.#
        if X_i_minimum[:, -1].sum() >= Limit_mortality:
            break
        # If the cell count for all types falls below 1, the patient is cured. Stop simulation.#
        elif all(X_i_minimum[:, -1] < 1):
            break
        # Update 'X0_i' and 'tnow'.#
        X0_i = X_i_minimum[:, -1]
        tnow += Stepsize
    return t_DPM2_i, X_DPM2_i, d_DPM2_i


# Define DPM3.#
def strategy_DPM3(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, mis_specification, mis_PAR):
    # DPM3: strategy balancing tumor burden and resistance emergence.#
    # Minimize total cell population, except when doing so would trigger the first multi-resistant cell.#
    # At each Stepsize:#
    # (1) If the predicted multi-resistant cells < 1 or is curable: select 'd_i' to minimize total cells.#
    # Exception: if selected 'd_i' causes first multi-resistant cell, instead select 'd_i' to minimize multi-resistant population.#
    # (2) If the current multi-resistant >= 1 and incurable: select 'd_i' to minimize total cells.#

    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_DPM3_i = np.zeros([0], dtype=float)
    X_DPM3_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_DPM3_i = np.zeros([PAR['Num_drug'], 0], dtype=float)

    while tnow < Simduration:
        X_i_total, t_i_total, d_i_total, X_i_end_total, X_i_end_multi_resis, t_sim = [], [], [], [], [], []
        t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
        for i_dose_combination in dose_combination:
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

        # If the current multi-resistant population is < 1 or is curable (curable meaning any 'Sa[-1, :] >= g0'),#
        # then minimize the total cell population.#
        index_minimum_end_total = np.argmin(X_i_end_total)
        index_select = index_minimum_end_total
        if (X0_i[-1] < 1 or (np.any(np.greater_equal(mis_Sa_i[-1, :], mis_g0_i[-1]))) if mis_specification
                else np.any(np.greater_equal(Sa_i[-1, :], g0_i[-1]))):
            X_i_minimum = X_i_total[int(index_minimum_end_total)]
            # If the first multi-resistant cell would emerge under the selected 'd_i', minimize the multi-resistant population.#
            if X0_i[-1] < 1 <= X_i_minimum[-1, -1]:
                index_minimum_end_multi_resis = np.argmin(X_i_end_multi_resis)
                index_same_minimum_end_multi_resis = [i for i, x in enumerate(X_i_end_multi_resis) if x ==
                                                      X_i_end_multi_resis[int(index_minimum_end_multi_resis)]]
                # If there is only one 'd_i' that yields the minimum multi-resistant population.#
                if len(index_same_minimum_end_multi_resis) == 1:
                    index_select = index_minimum_end_multi_resis
                # If there are multiple 'd_i' yield the same minimum multi-resistant populaiton, then among those 'd_i', choose the one that#
                # minimizes the total cell population in all.#
                else:
                    X_i_end_total_same_minimum_end_multi_resis = [X_i_end_total[i] for i in index_same_minimum_end_multi_resis]
                    # Among the 'd_i' that yield the same minimum multi-resistant population, find the one that minimizes the total cell#
                    # population.#
                    index_minimum_end_total_same_minimum_end_multi_resis = np.argmin(X_i_end_total_same_minimum_end_multi_resis)
                    index_select = index_same_minimum_end_multi_resis[int(index_minimum_end_total_same_minimum_end_multi_resis)]
                # If minimizing the multi-resistant population leads to mortality, choose the option that minimizes the total cell population.#
                if X_i_end_total[index_select] >= Limit_mortality:
                    index_select = index_minimum_end_total
        # If current multi-resistant population >= 1 and incurable, select 'd_i' to minimize the total cell population.#
        else:
            # The 'index_select' already set to 'index_minimum_end_total' above, no action needed.#
            pass
        # All dose combinations result in mortality. Select the combination with the longest survival time.#
        if X_i_end_total[index_select] >= Limit_mortality:
            index_select = np.argmax(t_sim)
        d_i_minimum = d_i_total[index_select]
        if mis_specification:
            X_i_minimum, t_i_minimum, d_i_minimum = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, LSsim)
        else:
            X_i_minimum, t_i_minimum = X_i_total[index_select], t_i_total[index_select]

        X_DPM3_i = np.append(X_DPM3_i, X_i_minimum, 1) if tnow == 0 else np.append(X_DPM3_i, X_i_minimum[:, 1:], 1)
        t_DPM3_i = np.append(t_DPM3_i, t_i_minimum + tnow) if tnow == 0 else np.append(t_DPM3_i, t_i_minimum[1:] + tnow)
        d_DPM3_i = np.append(d_DPM3_i, d_i_minimum, 1)

        # If total cell population >= 'Limit_mortality', stop simulation. Patient mortality occurs.#
        if X_i_minimum[:, -1].sum() >= Limit_mortality:
            break
        # If the cell count for all types falls below 1, the patient is cured. Stop simulation.#
        elif all(X_i_minimum[:, -1] < 1):
            break
        # Update 'X0_i' and 'tnow'.#
        X0_i = X_i_minimum[:, -1]
        tnow += Stepsize
    return t_DPM3_i, X_DPM3_i, d_DPM3_i


# Define DPM4.#
def strategy_DPM4(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, mis_specification, mis_PAR):
    # DPM4: Estimate the time to either incurability or mortality, and repond to the more imminent threat whenever cure remains possible.#
    # At each 'Stepsize', evaluate the predicted time until incurability (multi-resistant population >= 1) and predicted time until mortality#
    # (total population >= 1e13), based on the projected growth of S, R1, R2, R12.#
    # For each dosage combination 'd_i', define:#
    # τ_inc(d): predicted time to incurability (multi-resistant >= 1), given the current state and fixed 'd_i'.#
    # τ_S(d): predicted time until S reaches the mortality threshold (S > 1e13).#
    # τ_R1(d), predicted time until R1 reaches the mortality (R1 > 1e13).#
    # τ_R2(d), predicted time until R2 reaches the mortality (R2 > 1e13).#
    # τ_multi_resis(d) (τ_R12), predicted time until multi-resistant population reaches mortality threshold (multiply_resistant > 1e13).#
    # If the current multi-resistant population < 1 or multi-resistant population is curable#
    # (i.e., there exists some 'd_i' such that every component of diag(Sa*d) > g0),#
    # then adjust 'd_i' to maximize min(τ_inc,τ_s,τ_R1,τ_R2,τ_multi_resis),#
    # subject to the constraint min(τ_S,τ_R1,τ_R2,...,τ_multi_resis) > Stepsize.#
    # If no such dosage combination exists, maximize min(τ_S,τ_R1,τ_R2,..,τ_multi_resis) instead.#
    # If the current multi_resistant population >= 1 and the multi-resistant population is not curable,#
    # maximize min(τ_S,τ_R1,τ_R2,...,τ_multi_resis).#

    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_DPM4_i = np.zeros([0], dtype=float)
    X_DPM4_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_DPM4_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    curable_multi_resis = False

    # If the current multi-resistant population < 1, or if it is curable (meaning any 'Sa[-1, :] >= g0').#
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

            # Calculate the τ_inc,τ_s,τ_R1,τ_R2,τ_multi_resis.#
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

            τ_total = np.append(τ_total, np.reshape(τ_i_drug, (τ_i_drug.shape[0], 1)), 1)
            d_τ_i = np.append(d_τ_i, treat_i, 1)

        τ_S_to_multi_resis = τ_total[1:, :]
        min_τ_S_to_multi_resis = np.amin(τ_S_to_multi_resis, 0)
        index_min_τ_S_to_multi_resis_biggerthanStepsize, = np.where(min_τ_S_to_multi_resis > Stepsize)
        # If the current multi-resistant population < 1 or is curable.#
        if X0_i[-1] < 1 or curable_multi_resis:
            # If there exists any 'd_i' such that min(τ_S,τ_R1,τ_R2,τ_multi_resis) > Stepsize,#
            # then among those 'd_i', choose the one that maximizes min(τ_inc,τ_S,τ_R1,τ_R2,τ_R12).#
            if index_min_τ_S_to_multi_resis_biggerthanStepsize.size:
                τ_total = τ_total[:, index_min_τ_S_to_multi_resis_biggerthanStepsize]
                d_τ_i = d_τ_i[:, index_min_τ_S_to_multi_resis_biggerthanStepsize]
                min_τ_inc_S_to_multi_resis = np.amin(τ_total, 0)
                index_maxmin_τ_inc_S_to_multi_resis = np.argmax(min_τ_inc_S_to_multi_resis)
                d_τ = np.array([d_τ_i[:, index_maxmin_τ_inc_S_to_multi_resis]]).T
            # If no dosage combination satisfies min(τ_S,τ_R1,τ_R2,τ_multi_resis) > Stepsize,#
            # then maximize min(τ_S,τ_R1,τ_R2,τ_multi_resis).#
            else:
                index_maxmin_τ_S_to_multi_resis = np.argmax(min_τ_S_to_multi_resis)
                d_τ = np.array([d_τ_i[:, index_maxmin_τ_S_to_multi_resis]]).T
        # If the current multi-resistant population >= 1 and it is not curable, maximize min(τ_S,τ_R1,τ_R2,τ_multi_resis).#
        elif X0_i[-1] >= 1 and (not curable_multi_resis):
            index_maxmin_τ_S_to_multi_resis = np.argmax(min_τ_S_to_multi_resis)
            d_τ = np.array([d_τ_i[:, index_maxmin_τ_S_to_multi_resis]]).T
        else:
            raise ValueError('Wrong situation')

        t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
        d_i = np.tile(d_τ, (1, t_i.shape[0] - 1))
        X_i, t_i, d_i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)

        X_DPM4_i = np.append(X_DPM4_i, X_i, 1) if tnow == 0 else np.append(X_DPM4_i, X_i[:, 1:], 1)
        t_DPM4_i = np.append(t_DPM4_i, t_i + tnow) if tnow == 0 else np.append(t_DPM4_i, t_i[1:] + tnow)
        d_DPM4_i = np.append(d_DPM4_i, d_i, 1)

        # If total cell population >= 'Limit_mortality', stop simulation. Patient mortality occurs.#
        if X_i[:, -1].sum() >= Limit_mortality:
            break
        # If the cell count for all types falls below 1, the patient is cured. Stop simulation.#
        elif all(X_i[:, -1] < 1):
            break
        # Update 'X0_i' and 'tnow'.#
        X0_i = X_i[:, -1]
        tnow += Stepsize
    return t_DPM4_i, X_DPM4_i, d_DPM4_i
