from DPM_lib import bz2, os, csv, functools, operator, pd, pickle, tee, Parallel, delayed, multiprocessing
from DPM_assign_check import *
from DPM_generate import DPM_generate_default_par_2drug, DPM_generate_default_par_3drug
'''This script contains functions that read data from files and save data into files.'''


# Read parameters from .cvs file.
def DPM_read_par_csv(filename_csv, name):
    filename_text = Colored.BOLD + Colored.PURPLE + str(filename_csv) + Colored.END
    name_text = Colored.BOLD + Colored.PURPLE + str(name) + Colored.END
    par_csv = []
    if type(filename_csv) is str:
        try:
            df = pd.read_csv(filename_csv)
            par_csv = df.to_dict(orient='records')
        except FileNotFoundError:
            print('Can not open the file ' + filename_text + ' from the keyword input ' + name_text + '.')
            DPM_print_errorandclose()
            exit()
    else:
        DPM_print_val_type(filename_csv, name)
        print('Keyword input ' + name_text + ' should be a string of the filename or a list of strings of multiple filenames.')
        DPM_print_errorandclose()
        exit()
    if not par_csv:
        print('No parameters have been inputted from the file ' + filename_text + '.')
        DPM_print_errorandclose()
        exit()
    else:
        return par_csv


# Read file line by line.
def DPM_read_linebyline_from_first(filename):
    with open(filename) as f:
        for line in f:
            yield line


# Read file line by line in reverse.
def DPM_read_linebyline_reverse(filename):
    with open(filename, 'rb') as f:
        # move to the end of the file
        f.seek(0, os.SEEK_END)
        # get the current position
        pointer_location = f.tell()
        # create a buffer to keep the last read line
        buffer = bytearray()
        # loop till pointer reaches the top of the file
        while pointer_location >= 0:
            # move the file pointer to the location pointed by pointer_location
            f.seek(pointer_location)
            # shift pointer location -1
            pointer_location = pointer_location-1
            # read that byte/character
            new_byte = f.read(1)
            # if the read byte is new line character then it means one line is read
            if new_byte == b'\n':
                # Fetch the line from buffer and yield it
                yield buffer.decode()[::-1]
                # reinitialize the byte array to save next line
                buffer = bytearray()
            else:
                # if last read characeter is not eol then add it in buffer
                buffer.extend(new_byte)
        # as file is read completely, if there is still data in buffer, then its the first line
        if len(buffer) > 0:
            # yield the first line too
            yield buffer.decode()[::-1]


# Read file line by line in reverse backup.
def DPM_read_linebyline_reverse2(filename, buf_size=9000):
    # a generator that returns the lines of a file in reverse order
    with open(filename) as f:
        segment = None
        offset = 0
        f.seek(0, os.SEEK_END)
        file_size = remaining_size = f.tell()
        while remaining_size > 0:
            offset = min(file_size, offset+buf_size)
            f.seek(file_size-offset)
            buffer = f.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # the first line of the buffer is probably not a complete line so save it and append it to the last line of the next buffer
            if segment is not None:
                # if the previous chunk starts right from the beginning of line do not concat the segment to the last line of new chunk
                # Instead, yield the segment first
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines)-1, 0, -1):
                if lines[index]:
                    yield lines[index]
            # Don't yield None if the file was empty
            if segment is not None:
                yield segment


# Read file length.
def DPM_read_file_len(filename):
    with open(filename, 'r', newline='') as f:
        for i, l in enumerate(f):
            pass
    return i+1


# Read file last line index.
def DPM_read_lastline_index(filename):
    if os.path.exists(filename):
        if os.path.getsize(filename) != 0:
            linenum = DPM_read_file_len(filename)
            if linenum < 2:
                with open(filename, 'r', newline='') as f:
                    # csv_reader = csv.reader(f, delimiter=',')
                    # for row in csv_reader:
                    #     print(row)
                    lastline = f.readlines()[-1]
                    try:
                        firstitem = lastline.split(',')[0]
                        lastline = int(firstitem)
                    except ValueError:
                        lastline = 0
            else:
                with open(filename, 'r',) as f:
                    try:
                        end_of_file = f.seek(0, 2)
                        f.seek(end_of_file)
                        for num in range(end_of_file+1):
                            f.seek(end_of_file - num)
                            lastline = f.read()
                            if lastline.count('\n') == 2:
                                break

                        firstitem = lastline.split(',')[0]
                        lastline = int(firstitem)
                    except (ValueError, OSError, IndexError):
                        lastline = 0
        else:
            lastline = 0
    else:
        lastline = 0
    return lastline


# Open file and close in order to erase existed data.
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


# Save defualt parameter set PNAS(2012) and Biology Direct(2016) into pckl files.
def DPM_save_default_PARset(*argv, **kwargs):
    # Generation default parameter set for 2 or 3 drug case.
    # Parameter including criterias:
    # 1: PNAS(1), two drug situation
    # 2: PNAS(2), two drug situation
    # 3: PNAS(3), two drug situation
    # 4: Biology Direct(4), two drug situation
    # 5: Biology Direct(5), two drug situation
    # 6: Biology Direct(6), two drug situation
    # 7: Biology Direct(7), three drug situation
    # 8: Biology Direct(8), three drug situation

    X0total = X0TOTAL_DEFAULT_VAL
    parname_list = PARNAME_LIST_SAVE_DEFAULT_PARSET
    # If input non-keywrod arguments, error.
    if argv:
        print_quotation = True
        DPM_print_keywordargument_only(argv, parname_list, print_quotation)
        return

    # If the inputted keywords not in the allowed keywords list, error.
    print_quotation = True
    beginning_char = '\n'
    DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char)
    kwargs = {k.lower(): v for k, v in kwargs.items()}

    par = DPM_check_inputpar(kwargs, parname_list)
    Num_drug, par_save_block_size, Limit_mortality, Simduration, pathsave = \
        operator.itemgetter('Num_drug', 'par_save_block_size', 'Limit_mortality', 'Simduration', 'pathsave')(par)

    PARtotal = []
    PAR_criterisa = Num_cell_type = totalcombination = totalcount = totalcombination_copy = None
    if Num_drug == 2:
        Num_cell_type = 4
        # Default parameter criterisa used in PNAS(2012) 2 drug case.
        PAR_criterisa = [0, 1, 7, 8]
        g0, ratioR1toX0, ratioR2toX0, Sa_ratio_S_D1tog0, Sa_ratio_S_D2toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, T_StoR1, T_StoR2 \
            = DPM_generate_default_par_2drug()
        totalcombination = itertools.product(g0, ratioR1toX0, ratioR2toX0, Sa_ratio_S_D1tog0, Sa_ratio_S_D2toS_D1, Sa_ratio_R1_D1toS_D1,
                                             Sa_ratio_R2_D2toS_D2, T_StoR1, T_StoR2)
        totalcombination, totalcombination_copy = tee(totalcombination)
        iters = [g0, ratioR1toX0, ratioR2toX0, Sa_ratio_S_D1tog0, Sa_ratio_S_D2toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, T_StoR1, T_StoR2]
        totalcount = functools.reduce(operator.mul, map(len, iters), 1)
    elif Num_drug == 3:
        Num_cell_type = 8
        # Default parameter criterisa used in Biology Direct(2016) 3 drug cases.
        PAR_criterisa = [0, 1, 2, 3]
        g0, ratioStoX0, ratioR1toX0, ratioR2toX0, ratioR3toX0, ratioR12toX0, ratioR23toX0, ratioR13toX0, Sa_ratio_S_D1tog0, Sa_ratio_S_D2toS_D1, \
            Sa_ratio_S_D3toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, Sa_ratio_R3_D3toS_D3, T_StoR1, T_StoR2, T_StoR3 \
            = DPM_generate_default_par_3drug()
        totalcombination = itertools.product(g0, ratioStoX0, ratioR1toX0, ratioR2toX0, ratioR3toX0, ratioR12toX0, ratioR23toX0, ratioR13toX0,
                                             Sa_ratio_S_D1tog0, Sa_ratio_S_D2toS_D1, Sa_ratio_S_D3toS_D1, Sa_ratio_R1_D1toS_D1,
                                             Sa_ratio_R2_D2toS_D2, Sa_ratio_R3_D3toS_D3, T_StoR1, T_StoR2, T_StoR3)
        totalcombination, totalcombination_copy = tee(totalcombination)
        iters = [g0, ratioStoX0, ratioR1toX0, ratioR2toX0, ratioR3toX0, ratioR12toX0, ratioR23toX0, ratioR13toX0, Sa_ratio_S_D1tog0,
                 Sa_ratio_S_D2toS_D1, Sa_ratio_S_D3toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, Sa_ratio_R3_D3toS_D3,
                 T_StoR1, T_StoR2, T_StoR3]
        totalcount = functools.reduce(operator.mul, map(len, iters), 1)
    else:
        print('The input number of drug equals: ' + Colored.BOLD + Colored.PURPLE + str(Num_drug) + Colored.END)
        print('There is no default parameter sets for number of drug bigger than 3.')
        DPM_print_errorandclose()
        exit()

    # Filtering parameters.
    block_count = 1
    block_format = '{0:0' + str(int(-(-(np.log10(totalcount/par_save_block_size)) // 1))) + '}'
    num_cores = multiprocessing.cpu_count()
    idx_sel = Parallel(n_jobs=num_cores)(delayed(DPM_save_default_PARset_0)
                                         (i, i_par, Num_drug, Num_cell_type, X0total, PAR_criterisa, Simduration, Limit_mortality, totalcount)
                                         for i, i_par in enumerate(totalcombination))
    par_sel = [i_par for idx_sel, i_par in zip(idx_sel, totalcombination_copy) if idx_sel]
    paramID = 0
    while par_sel:
        par_block = par_sel[:int(par_save_block_size)]
        par_block = [DPM_assign_par_default_2drug(i_par, X0total) for i_par in par_block] if Num_drug == 2 else \
            [DPM_assign_par_default_3drug(i_par, X0total) for i_par in par_block] if Num_drug == 3 else None
        par_block = [{**{'paramID': i+paramID}, **i_par} for i, i_par in enumerate(par_block)]
        par_block = [{k: i_par[k] for k in HEADING_2DRUG_PARAM_CSV} for i_par in par_block] if Num_drug == 2 else \
            [{k: i_par[k] for k in HEADING_3DRUG_PARAM_CSV} for i_par in par_block] if Num_drug == 3 else None
        par_block = pd.DataFrame.from_records(par_block)

        filename = os.path.join(pathsave, 'PARset_default_PNAS_2drug_' + block_format.format(block_count) +
                                '_{:d}'.format(par_block.shape[0]) + '.csv') if Num_drug == 2 else \
            os.path.join(pathsave, 'PARset_default_Biology_Direct_3drug_' + block_format.format(block_count) +
                         '_{:d}'.format(par_block.shape[0]) + '.csv') if Num_drug == 3 else None
        par_block.to_csv(filename, index=False)

        par_sel = par_sel[int(par_save_block_size):]
        paramID += par_save_block_size
        block_count += 1
    return PARtotal


# Save the simulation results into .csv files.
def DPM_save_result_csv(PAR, Stepsize, Simduration, Strategy0, Strategy1, Strategy2_1, Strategy2_2, Strategy3, Strategy4, Strategy5,
                        Strategy6, Strategy7, Strategy8, Strategy9, par_index, filename_param, filename_stopt, filename_dosage,
                        filename_pop, filename_eachtimepoint, Limit_mortality, save=True):

    g0, X0, Sa, T = DPM_generate_g0(PAR), list(DPM_generate_X0(PAR)), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    # X0 = [i_Strategy[1][:, 0] for i_Strategy in [Strategy0, Strategy1, Strategy2_1, Strategy2_2, Strategy3, Strategy4, Strategy5, Strategy6,
    #                                              Strategy7, Strategy8, Strategy9] if i_Strategy is not None]
    # X0 = list(X0[0])

    # FILENAME_result_param_yymmdd.csv format.
    # Each line denotes a parameter configuration with the following entries:
    # (1) Parameter index.
    # (2) Initial population compositions X0(S), X0(R1), X0(R2),...,X0(Rxyz...).
    # (3) Growth rate g0 for each cell type g0(S), g0(R1), g0(R2),...,g0(Rxyz...).
    # (4) Drug sensitivity matrix Sa(S,drug1), Sa(S,drug2),...,Sa(S,drugN),...,Sa(Rxyz...,drug1), Sa(Rxyz...,drug2),...,Sa(Rxyz...,drugN).
    # (5) Transition rate matrix T(S<-S), T(S<-R1), T(S<-R2),...,T(S<-Rxyz.),...,T(Rxyz.<-S), T(Rxyz.<-R1), T(Rxyz.<->R2),...,T(Rxyz.<-Rxyz.).
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
    Strategy_text = [', strategy0, ', ', strategy1, ', ', strategy2.1, ', ', strategy2.2, ', ', strategy3, ', ', strategy4, ', ', strategy5, ',
                     ', strategy6, ', ', strategy7, ', ', strategy8, ', ', strategy9, ']
    Strategy = (Strategy0, Strategy1, Strategy2_1, Strategy2_2, Strategy3, Strategy4, Strategy5, Strategy6, Strategy7, Strategy8, Strategy9)
    Strategy_result = {}
    for i in range(len(Strategy_text)):
        i_Strategy_result = DPM_assign_saveresult_csv(Strategy[i], Strategy_text[i], Simduration, Stepsize, cured_survival, non_report,
                                                      para_PAR_index, filename_stopt, filename_dosage, filename_pop, filename_eachtimepoint,
                                                      Limit_mortality)
        Strategy_result[Strategy_text[i].replace(', ', '')] = i_Strategy_result

    # FILENAME_result_stopt_yymmdd.csv format.
    # Each line denotes the survival times (stopping time of days) of 11 strategies for one parameter configuration.
    # Line format:
    # (1) Parameter index.
    # (2) Strategy 0: strategy 0 in the PNAS paper.
    # (3) Strategy 1: strategy 1 in the PNAS paper.
    # (4) Strategy 2.1: strategy 2.1 in the PNAS paper.
    # (5) Strategy 2.2: strategy 2.2 in the PNAS paper.
    # (6) Strategy 3: strategy 3 in the PNAS paper.
    # (7) Strategy 4: strategy 4 in the PNAS paper.
    # (8) Strategy 5: multistep extension extension of strategy 1 in the Biology Direct paper.
    # (9) Strategy 6: multistep extension extension of strategy 2.1 in the Biology Direct paper.
    # (10) Strategy 7: multistep extension extension of strategy 2.2 in the Biology Direct paper.
    # (11) Strategy 8: multistep extension extension of strategy 3 in the Biology Direct paper.
    # (12) Strategy 9:  ALTO (Adaptive Long-Term Optimization) in the Biology Direct paper.

    # Each line has 12 entries: parameter index and the survival times of the 11 strategies.
    # If the strategy is not used, it will be set to nan.
    # The time horizon of the therapy equals value of the parameter Simduration (e.g., 1800 days).
    # Each period equals the value of the parameter Stepsize (e.g., 45 days).
    # If the patient is cured (every tumor subtype < 1), then the survival time is reported as the value of the parameter Simduration +
    # the value of the parameter Stepsize (e.g., 1800+45=1845) days.
    # If the survival time equals the value of the parameter Simduration (e.g., 1800 days) days, then at least one of the tumor subtype is greater
    # than 1 and the total cell population is smaller than the mortal level.

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

    # FILENAME_result_dosage_yymmdd.csv format:
    # Each line denotes the dosage combination sequence of one strategy at each Stepsize.
    # Line format:
    # (1) Parameter index.
    # (2) Strategy name (0, 1, 2.1, 2.2, 3, 4, 5, 6, 7, 8, 9).
    # (3) (Drug1 dosage, Drug2 dosage,...,DrugN dosage) at t = 0.
    # (4) (Drug1 dosage, Drug2 dosage,...,DrugN dosage) at t = 1 * Stepsize.
    # .
    # .
    # .
    # (i) (Drug1 dosage, Drug2 dosage,...,DrugN dosage) at t = (i-3) * Stepsize.
    # If t exceeds the survival time of the strategy, then the drug dosage is set to - 1.
    # If the strategy is not used, it will be set to nan.
    if filename_dosage:
        if save:
            with open(filename_dosage, 'a', newline='') as f:
                writer = csv.writer(f)
                for i_dose in [Strategy_result[i_strategy.replace(', ', '')]['dose'] for i_strategy in Strategy_text]:
                    writer.writerow(i_dose.split(', '))

    # FILENAME_result_pop_yymmdd.csv format:
    # Each line denotes the population composition dynamics of one strategy at each Stepsize.
    # Line format:
    # (1) Parameter index.
    # (2) Strategy name (0, 1, 2.1, 2.2, 3, 4, 5, 6, 7, 8, 9).
    # (3) (S,R1,R2,...,Rxyz.) population size at t = 1 * Stepsize.
    # (4) (S,R1,R2,...,Rxyz.) population size at t = 2 * Stepsize.
    # .
    # .
    # .
    # (j) (S,R1,R2,...,Rxyz.) population size at t = Simulation stops if mortality happens.
    # .
    # .
    # .
    # (i) (S,R1,R2,...,Rxyz.) at t = Simduration.
    # If t exceeds the survival time of the strategy, then the population size is set to -1.
    # If the strategy is not used, it will be set to nan.
    # Notice the dosage combination is reported at the beginning of each period, and the population composition is reported at the end of each period.

    if filename_pop:
        with open(filename_pop, 'a', newline='') as f:
            writer = csv.writer(f)
            for i_pop in [Strategy_result[i_strategy.replace(', ', '')]['pop'] for i_strategy in Strategy_text]:
                writer.writerow(i_pop.split(', '))

    # FILENAME_result_eachtimepoint_yymmdd.csv format:
    # Each line denotes the dosage and population composition of one strategy at each timepoint.
    # Line Format:
    # (1) Parameter index.
    # (2) Strategy name(0, 1, 2.1, 2.2, 3, 4, 5, 6, 7, 8, 9).
    # (3) (Drug1 dosage, Drug2 dosage,...,DrugN dosage) at t = 0.
    # (4) (S,R1,R2,...,Rxyz.) population size at t = 0.
    # (5) (Drug1 dosage, Drug2 dosage,...,DrugN dosage) at t = 1.
    # (6) (S,R1,R2,...,Rxyz.) population size at timepoint t = 1.
    # .
    # .
    # .
    # (i-1) (Drug1 dosage, Drug2 dosage, DrugN dosage) at t = Simduration.
    # (i) (S,R1,R2,...,Rxyz.) population size at timepoint t = Simduration.
    # If the strategy is not used, it will be set to nan.

    if filename_eachtimepoint:
        if save:
            with open(filename_eachtimepoint, 'a', newline='') as f:
                writer = csv.writer(f)
                for i_eachtimepoint in [Strategy_result[i_strategy.replace(', ', '')]['eachtimepoint'] for i_strategy in Strategy_text]:
                    writer.writerow(i_eachtimepoint.split(', '))
    return returnval


def DPM_save_fcs(fcs, filename):
    fc = dict()
    channel_all = FCSCHANNEL1 if 'Helix NP NIR-A' in set(fcs.channels) else FCSCHANNEL2
    for i_ch in channel_all:
        if i_ch != 'Time':
            fc[i_ch] = list(fcs[:, i_ch])
    # fc = SimpleNamespace(**fc)
    with bz2.BZ2File(filename, 'wb') as f:
        pickle.dump(fc, f)
    return


def DPM_read_fcspckl(filename):
    def DPM_read_fcspckl_nested(i_filename):
        with bz2.BZ2File(i_filename, 'rb') as f:
            i_fcs = pickle.load(f)
            val = getattr(i_fcs, 'Cell trace violet-A') if hasattr(i_fcs, 'Cell trace violet-A') else getattr(i_fcs, 'BV 421-A')
        return val

    live, dead = [], []
    for i_file in filename:
        i_file_live = i_file.replace('.pckl', '_live.pckl')
        i_file_dead = i_file.replace('.pckl', '_dead.pckl')
        live.append(DPM_read_fcspckl_nested(i_file_live))
        dead.append(DPM_read_fcspckl_nested(i_file_dead))
    return live, dead


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


def DPM_read_dosage_csv(filename, Strategy_name, Num_stepdiff):
    ind = [i+1 for i in range(len(STRATEGY_LIST[:-2])) if STRATEGY_LIST[:-2][i] in Strategy_name]
    df_dosage = pd.read_csv(filename, header=None, dtype=str,
                            skiprows=lambda x: x % len(STRATEGY_LIST[:-2]) not in ind)  # usecols=range(0, 2+Num_stepdiff),
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
        dosage_out[i_strategy] = list(i_dosage.apply(';'.join, axis=1).astype(str).values)
    # all paramID are the same
    if not all(ele == paramID[0] for ele in paramID):
        raise Exception('Not all paramID are the same in the strategies.')

    dosage_out['paramID'] = paramID[0]
    return dosage_out


def DPM_read_misspecify_LOD(pathloadmis, Strategy_name, Num_stepdiff, file_stoptime_mis,
                            file_dosage_mis, file_para_mis, i_filename_pattern, para_key):
    def DPM_read_misspecify_LOD_1(filename, name, filenameall):
        filename = filename + name
        file_out = []
        for i_val in filenameall.keys():
            file_out.extend([i_file for i_file in filenameall[i_val] if filename in i_file])
        if len(file_out) != len(filenameall.keys()):
            raise FileNotFoundError('Can not find file. ')
        return file_out

    file_stoptime = DPM_read_misspecify_LOD_1(i_filename_pattern, '_stopt', file_stoptime_mis)
    file_dosage = DPM_read_misspecify_LOD_1(i_filename_pattern, '_dosage', file_dosage_mis)
    file_para = DPM_read_misspecify_LOD_1(i_filename_pattern, '_para', file_para_mis)

    keys, stoptime, dosage, para = list(file_stoptime_mis.keys()), [], [], []
    for i, i_key in enumerate(keys):
        stoptime.append(DPM_read_stoptime_csv(os.path.join(pathloadmis, i_key, file_stoptime[i]), Strategy_name))
        dosage.append(DPM_read_dosage_csv(os.path.join(pathloadmis, i_key, file_dosage[i]), Strategy_name, Num_stepdiff))
        para.append(pd.read_csv(os.path.join(pathloadmis, i_key, file_para[i]), usecols=para_key, dtype=float).to_dict('list'))

    i_stoptime, i_dosage, i_para = dict(zip(keys, stoptime)), dict(zip(keys, dosage)), dict(zip(keys, para))

    return i_stoptime, i_dosage, i_para

