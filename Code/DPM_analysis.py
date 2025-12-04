from DPM_lib import pd, KaplanMeierFitter, CoxPHFitter, statistics, ExponentialFitter, itertools, np, re
'''
CHECKED
'''


def DPM_analysis_hazard_ratio(stoptime, Strategy_name, Simduration):
    p = DPM_analysis_pairwise_logrank_test(stoptime, Strategy_name, Simduration)
    hz = dict()
    for i in list(itertools.combinations(Strategy_name, 2)):
        if 'CPM' in i:
            ref = stoptime['CPM']
            treat = stoptime[list(filter(lambda x: x != 'CPM', i))[0]]
            name = ('CPM', list(filter(lambda x: x != 'CPM', i))[0])
        else:
            ref, treat = stoptime[i[0]], stoptime[i[1]]
            name = i
        hz[name] = DPM_analysis_HZ(ref, treat, Simduration)
    return hz, p


def DPM_analysis_KM(data, duration):
    E = [1 if i_val <= duration else 0 for i_val in data]
    epf = ExponentialFitter().fit(data, E)
    hazard = {'mean': epf.hazard_.mean(), 'ci': epf.confidence_interval_hazard_.mean()}

    kmf = KaplanMeierFitter()
    kmf.fit(data, E)
    km_interval = kmf.confidence_interval_survival_function_
    km = kmf.survival_function_

    t = km.index.values
    val = km.iloc[:].values.flatten()
    interval_lower = km_interval.iloc[:, 0].values.flatten()
    interval_upper = km_interval.iloc[:, 1].values.flatten()
    median_survival = kmf.median_survival_time_

    idx = np.where(t <= duration)[0]
    t, val, interval_lower, interval_upper = t[idx], val[idx], interval_lower[idx], interval_upper[idx]
    t, val, interval_lower, interval_upper = np.append(t, duration), \
        np.append(val, val[-1]), \
        np.append(interval_lower, interval_lower[-1]), \
        np.append(interval_upper, interval_upper[-1])

    return {'t': t, 'val': val, 'median_survival': median_survival, 'interval_lower': interval_lower, 'interval_upper': interval_upper,
            'hazard': hazard}


def DPM_analysis_HZ(data_ref, data, Simduration):
    if (len(data_ref) != 0) & (len(data) != 0):
        treat = np.concatenate((np.zeros(len(data_ref)), np.ones(len(data))))
        val = data_ref + data
        E = [1 if i_val <= Simduration else 0 for i_val in data_ref]
        E.extend([1 if i_val <= Simduration else 0 for i_val in data])
        E = np.array(E)
        d = {'val': val, 'E': E, 'treat': treat}
        df = pd.DataFrame(data=d)
        cph = CoxPHFitter()
        cph.fit(df, duration_col='val', event_col='E')
        hz_ratio = cph.hazard_ratios_.values[0]
    else:
        hz_ratio = None
    return hz_ratio


def DPM_analysis_pairwise_logrank_test(stoptime, Strategy_name, Simduration):
    G, T, E = tuple(), tuple(), tuple()
    flag_empty = False
    for i_strategy in Strategy_name:
        i_stop = stoptime[i_strategy]
        if not i_stop:
            flag_empty = True
            break
        G = G + tuple(len(i_stop) * [str(i_strategy)])
        E = E + tuple([1 if i_val <= Simduration else 0 for i_val in i_stop])
        T = T + tuple(i_stop)
    p = statistics.pairwise_logrank_test(T, G, E) if not flag_empty else None
    p_out = dict(zip(p.name, p.p_value)) if not flag_empty else None
    return p_out


def DPM_analysis_dose(dose, strategyname, inddrug=1):
    firstuse, max_num_change, num_change, average_move_num = [], [], [], []
    use_atbegin = 0
    for i, i_dose in enumerate(dose):
        i_dose = i_dose.split(';')
        if '-1' in i_dose:
            i_dose = i_dose[:i_dose.index('-1')]
        # Doses may be stored as '(np.float64(1.0),np.float64(0.0))', remove the np.float64 wrappers.#
        for j in range(len(i_dose)):
            i_dose[j] = re.sub(r'np\.float64\(([^)]+)\)', r'\1', i_dose[j])

        i_firstuse, i_average_move_num, drugovermove, drugtotal, i_num_change,  i_current = None, None, 0, 0, 0, i_dose[0]
        if i_current == '(0.0,1.0)':
            use_atbegin += 1

        for j, j_step in enumerate(i_dose):
            i_val = [float(i_val) for i_val in j_step[1:-2].split(',')]
            drugtotal = drugtotal + i_val[inddrug]
            drugovermove = drugovermove + (j+1) * i_val[inddrug]
            if strategyname == 'CPM' and i_current not in ['(0.0,1.0)', '(0.0,0.0)', '(1.0,0.0)']:
                assert strategyname == 'CPM' and i_current not in ['(0.0,1.0)', '(0.0,0.0)', '(1.0,0.0)']
            if i_firstuse is None and i_val[inddrug] != 0:
                i_firstuse = j+1
            if j_step != i_current and j_step != '(0.0,0.0)':
                i_num_change += 1
                i_current = j_step

        i_average_move_num = drugovermove/drugtotal if drugtotal != 0 else None
        i_firstuse = len(i_dose) + 1 if i_firstuse is None else i_firstuse
        i_average_move_num = len(i_dose) + 1 if i_average_move_num is None else i_average_move_num
        firstuse.append(i_firstuse)
        average_move_num.append(i_average_move_num)
        num_change.append(i_num_change/len(i_dose))
        max_num_change.append(i_num_change)

    return dict(zip(['first drug 2', 'average move number', 'drug changes', 'max drug changes'],
                    [np.mean(firstuse), np.mean(average_move_num), np.mean(num_change), np.max(max_num_change)]))


def DPM_analysis_sigbetter(stoptime, Strategy):
    paramID, stoptime_ref, stoptime_test = stoptime['paramID'], stoptime[Strategy[0]], stoptime[Strategy[1]]
    ind = np.logical_and(np.array(stoptime_test) > np.array(stoptime_ref) + 30*2, np.array(stoptime_test) > 1.25 * np.array(stoptime_ref))
    return list(itertools.compress(paramID, ind))
