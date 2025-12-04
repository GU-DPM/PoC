from DPM_lib import csv, pd, re
from DPM_assign_check import *
'''This script contains functions for reading data from files and saving data into files.'''
'''
CHECKED
'''


# Read parameters from a .cvs file.#
def DPM_read_par_csv(filename_csv, name):
    filename_text = Colored.BOLD + Colored.PURPLE + str(filename_csv) + Colored.END
    name_text = Colored.BOLD + Colored.PURPLE + str(name) + Colored.END
    # par_csv = []
    if type(filename_csv) is str:
        try:
            df = pd.read_csv(filename_csv)
            par_csv = df.to_dict(orient='records')
        except FileNotFoundError:
            print('Cannot open the file ' + filename_text + ' from the keyword input ' + name_text + '.')
            DPM_print_errorandclose()
            exit()
    else:
        DPM_print_val_type(filename_csv, name)
        print('Keyword input ' + name_text + ' should be a string containing a filename or a list of strings contaning multiple filenames.')
        DPM_print_errorandclose()
        exit()
    if not par_csv:
        print('No parameters have been input from the file ' + filename_text + '.')
        DPM_print_errorandclose()
        exit()
    else:
        return par_csv


# Read file line by line.#
def DPM_read_linebyline_from_first(filename):
    with open(filename) as f:
        for line in f:
            yield line


# Read the file length.#
def DPM_read_file_len(filename):
    with open(filename, 'r', newline='') as f:
        for i, l in enumerate(f):
            pass
    return i+1


# Open and close the file to erase existing data.#
def DPM_read_erase(filename_zip):
    for i_file in filename_zip:
        if i_file[0]:
            with open(i_file[1], 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(i_file[2])
    return


def DPM_save_default_PARset_0(i, i_par, Num_drug, Num_cell_type, X0total, PAR_criterisa, Simduration, Limit_mortality, totalcount):
    if Num_drug == 2:
        PAR_i = DPM_assign_par_default_2drug(i_par, X0total)
        PAR_i = {**{'Num_drug': Num_drug, 'Num_cell_type': Num_cell_type}, **PAR_i}
        criterisa_i = DPM_check_criterisa_2drug(PAR_criterisa, PAR_i, Simduration, Limit_mortality)
        print(round(i/totalcount, 3))
        return all(criterisa_i)
    elif Num_drug == 3:
        PAR_i = DPM_assign_par_default_3drug(i_par, X0total)
        PAR_i = {**{'Num_drug': Num_drug, 'Num_cell_type': Num_cell_type}, **PAR_i}
        criterisa_i = DPM_check_criterisa_3drug(PAR_criterisa, PAR_i, Simduration, Limit_mortality)
        print(round(i / totalcount, 3))
        return all(criterisa_i)


# Save the simulation results to .csv files.#
def DPM_save_result_csv(PAR, Stepsize, Simduration, CPM, DPM1, DPM2_1, DPM2_2, DPM3, DPM4, par_index, filename_param, filename_stopt,
                        filename_dosage, filename_pop, filename_eachtimepoint, Limit_mortality, save=True):
    g0, X0, Sa, T = DPM_generate_g0(PAR), list(DPM_generate_X0(PAR)), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    # FILENAME_result_param_yymmdd.csv format.#
    # Each line represents a parameter configuration with the following entries:#
    # (1) Parameter index.#
    # (2) Initial population compositions X0(S), X0(R1), X0(R2),...,X0(Rxyz...).#
    # (3) Growth rate g0 for each cell type g0(S), g0(R1), g0(R2),...,g0(Rxyz...).#
    # (4) Drug sensitivity matrix Sa(S,drug1), Sa(S,drug2),...,Sa(S,drugN),...,Sa(Rxyz...,drug1), Sa(Rxyz...,drug2),...,Sa(Rxyz...,drugN).#
    # (5) Transition rate matrix T(S<-S), T(S<-R1), T(S<-R2),...,T(S<-Rxyz.),...,T(Rxyz.<-S), T(Rxyz.<-R1), T(Rxyz.<->R2),...,T(Rxyz.<-Rxyz.).#
    para_PAR_index = str(par_index)
    if filename_param:
        para_g0 = ','.join([f'{i_g0:g}' for i_g0 in g0])
        para_X0composition = ','.join([f'{i_X0:.0f}' for i_X0 in X0])
        para_Sa = ','.join([f'{Sa[i,j]:g}' for i in range(Sa.shape[0]) for j in range(Sa.shape[1])])
        para_T = ','.join([f'{T[i,j]:g}' for i in range(T.shape[0]) for j in range(T.shape[1])])

        para = str(par_index) + ',' + para_X0composition + ',' + para_g0 + ',' + para_Sa + ',' + para_T
        if save:
            with open(filename_param, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(para.split(','))

    cured_survival = Simduration + Stepsize
    non_report = 'nan'
    Strategy_text = [', CPM, ', ', DPM1, ', ', DPM2.1, ', ', DPM2.2, ', ', DPM3, ', ', DPM4, ']
    Strategy = (CPM, DPM1, DPM2_1, DPM2_2, DPM3, DPM4)
    Strategy_result = {}
    for i in range(len(Strategy_text)):
        i_Strategy_result = DPM_assign_saveresult_csv(Strategy[i], Strategy_text[i], Simduration, Stepsize, cured_survival, non_report,
                                                      para_PAR_index, filename_stopt, filename_dosage, filename_pop, filename_eachtimepoint,
                                                      Limit_mortality)
        Strategy_result[Strategy_text[i].replace(', ', '')] = i_Strategy_result

    # FILENAME_result_stopt_yymmdd.csv format.#
    # Each line represents the survival times (stopping time of days) of the 7 strategies for a single parameter configuration.#
    # Line format:#
    # (1) Parameter index.#
    # (2) CPM: CPM in the PNAS paper.#
    # (3) DPM1: DPM1 in the PNAS paper.#
    # (4) DPM2.1: DPM2.1 in the PNAS paper.#
    # (5) DPM2.2: DPM2.2 in the PNAS paper.#
    # (6) DPM3: DPM3 in the PNAS paper.#
    # (7) DPM4: DPM4 in the PNAS paper.#

    # Each line has 12 entries: the parameter index and followed by the survival times for the 6 treatment strategies.#
    # Strategies that are not used are recorded as NaN.#
    # The therapy time horizon equals the value of the parameter 'Simduration' (e.g., 1800 days). Each treatment period has length 'Stepsize'#
    # (e.g., 45 days).#
    # If the patient is cured (all tumor subtype < 1), the survival time is reported as 'Simduration + Stepsize' (e.g., 1800+45=1845 days).#
    # If the survival time equals 'Simduration' (e.g., 1800 days), then at least one tumor subtype remains above 1, but the total cell#
    # population has not yet reached the mortality threshold.#

    if filename_stopt:
        survival = [Strategy_result[i_strategy.replace(', ', '')]['survival'] for i_strategy in Strategy_text]
        returnval = [i_val for i_val in survival if i_val != 'nan']
        survival = str(para_PAR_index) + ',' + ','.join(str(i) for i in survival)
        if save:
            with open(filename_stopt, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(survival.split(','))
    else:
        returnval = -1

    # FILENAME_result_dosage_yymmdd.csv format:#
    # Each line represents the dosage combination sequence of one strategy at each Stepsize interval.#
    # Line format:#
    # (1) Parameter index.#
    # (2) Strategy name (CPM, DPM1, DPM2.1, DPM2.2, DPM3, DPM4).#
    # (3) (Drug1 dosage, Drug2 dosage,...,DrugN dosage) at t = 0.#
    # (4) (Drug1 dosage, Drug2 dosage,...,DrugN dosage) at t = 1 * Stepsize.#
    # .
    # .
    # .
    # (i) (Drug1 dosage, Drug2 dosage,...,DrugN dosage) at t = (i-3) * Stepsize.#
    # If t exceeds the survival time of the strategy, then the drug dosage is set to - 1.#
    # If the strategy is not used, it will be set to NaN.#
    if filename_dosage:
        if save:
            with open(filename_dosage, 'a', newline='') as f:
                writer = csv.writer(f)
                for i_dose in [Strategy_result[i_strategy.replace(', ', '')]['dose'] for i_strategy in Strategy_text]:
                    writer.writerow(i_dose.split(', '))

    # FILENAME_result_pop_yymmdd.csv format:#
    # Each line represents the population composition dynamics of one strategy at each Stepsize interval.#
    # Line format:#
    # (1) Parameter index.#
    # (2) Strategy name (CPM, DPM1, DPM2.1, DPM2.2, DPM3, DPM4).#
    # (3) (S,R1,R2,...,Rxyz.) population size at t = 1 * Stepsize.#
    # (4) (S,R1,R2,...,Rxyz.) population size at t = 2 * Stepsize.#
    # .
    # .
    # .
    # (j) (S,R1,R2,...,Rxyz.) population size at t = Simulation stops if mortality happens.#
    # .
    # .
    # .
    # (i) (S,R1,R2,...,Rxyz.) at t = Simduration.#
    # If t exceeds the survival time of the strategy, the population size is set to -1.#
    # If the strategy is not used, its value is set to NaN.#
    # Notice that the dosage combination is recorded at the beginning of each period, whereas the population composition is recorded at the end.#

    if filename_pop:
        with open(filename_pop, 'a', newline='') as f:
            writer = csv.writer(f)
            for i_pop in [Strategy_result[i_strategy.replace(', ', '')]['pop'] for i_strategy in Strategy_text]:
                writer.writerow(i_pop.split(', '))

    # FILENAME_result_eachtimepoint_yymmdd.csv format:#
    # Each line represents the dosage and population composition of one strategy at each timepoint.#
    # Line Format:#
    # (1) Parameter index.#
    # (2) Strategy name(CPM, DPM1, DPM2.1, DPM2.2, DPM3, DPM4).#
    # (3) (Drug1 dosage, Drug2 dosage,...,DrugN dosage) at t = 0.#
    # (4) (S,R1,R2,...,Rxyz.) population size at t = 0.#
    # (5) (Drug1 dosage, Drug2 dosage,...,DrugN dosage) at t = 1.#
    # (6) (S,R1,R2,...,Rxyz.) population size at timepoint t = 1.#
    # .
    # .
    # .
    # (i-1) (Drug1 dosage, Drug2 dosage, DrugN dosage) at t = Simduration.#
    # (i) (S,R1,R2,...,Rxyz.) population size at timepoint t = Simduration.#
    # If the strategy is not used, it will be set to NaN.#

    if filename_eachtimepoint:
        if save:
            with open(filename_eachtimepoint, 'a', newline='') as f:
                writer = csv.writer(f)
                for i_eachtimepoint in [Strategy_result[i_strategy.replace(', ', '')]['eachtimepoint'] for i_strategy in Strategy_text]:
                    writer.writerow(i_eachtimepoint.split(', '))
    return returnval


def DPM_read_stoptime_csv(filename, Strategy_name):
    stoptime_key = ['paramID']
    stoptime_key.extend(Strategy_name)
    df_stoptime = pd.read_csv(filename, usecols=stoptime_key, dtype='int32')
    if df_stoptime.isnull().values.any():
        print('There is NaN values')
        DPM_print_errorandclose()
        sys.exit()
    data = df_stoptime.to_dict('list')

    return data


def DPM_read_dosage_pop_csv(filename, Strategy_name):
    ind = [i+1 for i in range(len(STRATEGY_LIST)) if STRATEGY_LIST[i] in Strategy_name]
    df_dosage = pd.read_csv(filename, header=None, dtype=str,
                            skiprows=lambda x: x % len(STRATEGY_LIST) not in ind)  # usecols=range(0, 2+Num_stepdiff),
    if df_dosage.isnull().values.any():
        print('There is NaN values')
        DPM_print_errorandclose()
        sys.exit()
    header = pd.read_csv(filename, index_col=False, nrows=0, usecols=range(0, df_dosage.shape[1])).columns.tolist()
    df_dosage = df_dosage.set_axis(header, axis=1)
    paramID = []
    dosage_out = dict()
    for i_strategy in Strategy_name:
        i_dosage = df_dosage.loc[df_dosage['Strategy name'].str.lower() == i_strategy.lower()]
        paramID.append(tuple(i_dosage['paramID'].astype(int).values))
        i_dosage = i_dosage.drop(columns=['paramID', 'Strategy name'])
        i_dosage = list(i_dosage.apply(';'.join, axis=1).astype(str).values)
        # Doses may be stored as '(np.float64(1.0),np.float64(0.0))', remove the np.float64 wrappers.#
        for j in range(len(i_dosage)):
            i_dosage[j] = re.sub(r'np\.float64\(([^)]+)\)', r'\1', i_dosage[j])
        dosage_out[i_strategy] = i_dosage
    # All paramID values are the same.#
    if not all(ele == paramID[0] for ele in paramID):
        raise Exception('Not all paramID are the same in the strategies.')

    dosage_out['paramID'] = paramID[0]
    return dosage_out
