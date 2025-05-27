from DPM_constant import *
from DPM_lib import itertools, plt, sm, LinearRegression, stats, lognorm, Bounds, minimize, add, solve_ivp
from DPM_plot import DPM_plot_linearreg

def DPM_data_analysis_linearreg(x, y, alpha=0.05):
    model = LinearRegression(fit_intercept=False)
    model.fit(x, y)

    yhat = model.predict(x)
    SS_Residual = sum((y-yhat)**2)
    SS_Total = sum((y-np.mean(y))**2)
    r_squared = 1-(float(SS_Residual))/SS_Total
    adjusted_r_squared = 1-(1-r_squared)*(len(y)-1)/(len(y)-x.shape[1]-1)

    # the coefficients of the regression model
    coefs = np.r_[[model.intercept_], model.coef_]
    # build an auxiliary dataframe with the constant term in it
    X_aux = x.copy()
    X_aux = np.hstack((np.ones((X_aux.shape[0], 1)), X_aux))
    # degrees of freedom
    dof = -np.diff(X_aux.shape)[0] + 1
    # Student's t-distribution table lookup
    t_val = stats.t.isf(alpha/2, dof)
    # mse of the residuals
    mse = np.sum((y - model.predict(x)) ** 2)/dof
    # inverse of the variance of the parameters
    var_params = np.diag(np.linalg.inv(X_aux.T.dot(X_aux)))
    # distance between lower and upper bound of CI
    gap = t_val * np.sqrt(mse*var_params)
    coef = coefs[1]
    coef_ci = [coef - gap[1], coef + gap[1]]

    return coef, coef_ci, r_squared, adjusted_r_squared


def DPM_data_analysis_cellcount(time, data, condition, ylim):
    data = [np.log(i) for i in data]
    numrep = [i.shape[0] for i in data]
    x = []
    [x.extend([time[i]] * i_num) for i, i_num in enumerate(numrep)]
    x = np.asarray(x).reshape(-1, 1)
    data_list = []
    [data_list.extend(list(i)) for i in data]
    data = np.asarray(data_list).reshape(-1, 1)
    x = sm.add_constant(x)
    model_i = sm.OLS(data, x).fit()
    DPM_plot_linearreg(x, data, model_i, condition, ylim)
    return

def DPM_fcs_analysis_fitlognorm(data, mean_val, generation, mingen, maxgen, method='EM', fix=False):
    all_possiable = [itertools.combinations(generation, L) for L in range(1, len(generation)+1)]
    logdata = np.array([np.log(x) for x in data if x > 0]).reshape(-1, 1)
    n, bins, _ = plt.hist(logdata, bins=np.linspace(np.min(np.array(logdata)), np.max(np.array(logdata)), int(len(logdata) / 100)),
                          histtype='step', density=True)
    plt.close()
    x_bins = [(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]
    cellnum, gen, scale, shape, loc, weight, aic, pdf, pdf_each, comb = [], [], [], [], [], [], [], [], [], []
    # filter the data
    data_use = [i_data for i_data in data if mean_val / 2 ** (maxgen + 1) < i_data < mean_val / 2 ** (mingen - 1)]
    if method.upper() == 'EM' and not fix:
        cellnum, gen, mean_best, var_best, loc_best, weight_best, aic = \
            DPM_fcs_analysis_fitlognorm_EM(data_use, logdata, mean_val, bins, all_possiable)
    elif method.upper() == 'OP':
        all_possiable = [itertools.combinations(generation, len(generation))] if fix else all_possiable
        cellnum, gen, scale, shape, loc, weight, aic, pdf, pdf_each, comb = \
            DPM_fcs_analysis_fitlognorm_OP(data_use, mean_val, bins, all_possiable, fix)

    return cellnum, gen, scale, shape, loc, weight, aic, pdf, pdf_each, data_use, bins, x_bins, comb


def DPM_fcs_analysis_fitlognorm_EM(data, logdata, mean_val, bins, all_possiable, fix):
    m, m_scale, m_shape, m_weight, m_aic, m_bic, m_gen = [], [], [], [], [], [], []
    for subset in all_possiable:
        for i, combination in enumerate(subset):
            m_gen.append(combination)
            means_init = np.array([mean_val / (2 ** i) for i in combination]).reshape(-1, 1)
            n_components = len(means_init)
            i_m = GaussianMixture(n_components, covariance_type='full', means_init=means_init, max_iter=int(1e3)).fit(logdata)
            m.append(i_m)
            m_aic.append(i_m.aic(logdata))
            m_bic.append(i_m.bic(logdata))
            m_scale.append(i_m.means_.flatten())
            m_shape.append(i_m.covariances_.flatten())
            m_weight.append(i_m.weights_.flatten())
    ind_best = np.argmin(m_aic)
    m_best, scale_best, shape_best, weight_best = m[ind_best], m_scale[ind_best], m_shape[ind_best], m_weight[ind_best]
    loc_best = [0 for _ in range(len(weight_best))]
    gen = [int(np.round(np.log2(mean_val / np.exp(i_scale)))) for i_scale in scale_best]
    ind = m_best.fit_predict(logdata)
    cellnum = [sum(ind == i_ind) for i_ind in range(len(gen))]
    x_bins = [(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]
    # best_lognormpdf = [weight_best[i] * lognorm.pdf(np.exp(x_bins), s=shape_best[i], loc=0, scale=np.exp(scale_best[i])) for i in
    #                    range(len(scale_best))]
    if len(gen) == 1:
        data_shape, data_loc, data_scale = lognorm.fit(data, loc=0)
        gen = [np.round(np.log2(mean_val / data_scale))]
        scale_best, shape_best, loc_best, weight_best = [data_scale], [data_shape], [data_loc], [1]
        best_lognormpdf = [lognorm.pdf(np.exp(x_bins), s=data_shape, loc=data_loc, scale=data_scale)]

        # plt.gca().set_xscale('log')
        # plt.plot(np.exp(x_bins), best_lognormpdf[0])
        # n, bins, _ = plt.hist(np.exp(logdata), bins=np.exp(x_bins), histtype='step', density=True)

    return cellnum, gen, scale_best, shape_best, loc_best, weight_best, m_aic


def DPM_fcs_analysis_fitlognorm_OP(data, mean_val, bins, all_possiable, fix):
    x_bins = [(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]
    m_scale, m_shape, m_loc, m_weight, m_aic, m_gen, m_pdf, m_pdf_each, m_comb = [], [], [], [], [], [], [], [], []
    for subset in all_possiable:
        for i, combination in enumerate(subset):
            m_comb.append(combination)
            if len(combination) == 1:
                i_shape, i_loc, i_scale = lognorm.fit(data, loc=0)
                i_pdf = lognorm.pdf(np.exp(x_bins), i_shape, loc=i_loc, scale=i_scale).tolist()
                i_pdf_each = [i_pdf]
                L = lognorm.pdf(data, i_shape, i_loc, i_scale)
                min_val = min(i for i in L if i > 0)
                L = [i if i > 0 else min_val for i in L]
                aic = -2*sum(np.log(L)) + 3
                # aic = -2*sum(lognorm.logpdf(data, i_shape, i_loc, i_scale)) + 3
                i_gen = [np.round(np.log2(mean_val / i_scale))]
                i_shape, i_loc, i_scale, i_weight = [i_shape], [i_loc], [i_scale], [1]
            else:
                # scale
                # scale0 = [mean_val/(2**i) for i in combination]
                # scale_lb = [mean_val/(2**(max(combination)+1)) for _ in combination]
                # scale_ub = [mean_val/(2**(min(combination)-1)) for _ in combination]

                scale0 = [mean_val/(2**combination[0])]
                scale_lb = [mean_val/(2**(combination[0]+1))]
                scale_ub = [mean_val/(2**(combination[0]-1))]

                # shape
                shape0 = [0.1 for _ in range(len(combination))]
                shape_lb = [1e-2 for _ in range(len(combination))]
                shape_ub = [1 for _ in range(len(combination))]

                # shape0 = [0.1]
                # shape_lb = [1e-2]
                # shape_ub = [1]

                # loc
                loc0 = [0 for _ in range(len(combination))]
                loc_lb = [-110 for _ in range(len(combination))]
                loc_ub = [0 for _ in range(len(combination))]

                # loc0 = [0]
                # loc_lb = [-100]
                # loc_ub = [0]

                # weight
                weight0 = [1/len(combination) for _ in range(len(combination)-1)]
                weight_lb = [0 for _ in range(len(weight0))] if not fix else [0.1 for _ in range(len(weight0))]
                weight_ub = [1 for _ in range(len(weight0))] if not fix else [1-0.1 * len(weight0) for _ in range(len(weight0))]

                x0 = np.asarray(scale0 + shape0 + loc0 + weight0)
                lb = scale_lb + shape_lb + loc_lb + weight_lb
                ub = scale_ub + shape_ub + loc_ub + weight_ub
                bounds = Bounds(lb, ub)
                A = np.concatenate((np.zeros((1, len(scale0))), np.zeros((1, len(shape0))), np.zeros((1, len(loc0))),
                                    np.ones((1, len(weight0)))), axis=1)

                x_op = minimize(DPM_fcs_analysis_fitlognorm_OP_fun, x0, args=(data, list(combination)), bounds=bounds,
                                method='trust-constr', constraints={'type': 'ineq', 'fun': lambda xval: 0.9-A.dot(xval)})
                x = x_op.x
                # x = np.append(x, 1 - sum(x[len(combination) * 3:]))
                weight = x[-len(combination) + 1:]
                weight = np.append(weight, 1 - sum(weight))
                x = x[:-len(combination) + 1]
                scale0 = x[0]
                x = np.delete(x, 0)
                shape = x[:-len(combination)]
                loc = x[-len(combination):]

                i_L, i_shape, i_loc, i_scale, i_weight, i_gen, i_pdf, i_pdf_each = \
                    [0 for _ in range(len(data))], [], [], [], [], [], [0 for _ in range(len(x_bins))], []
                for j, j_gen in enumerate(combination):
                    j_scale = scale0 if j == 0 else scale0 / (2 ** (combination[j]-combination[0]))
                    j_shape, j_loc, j_weight = shape[j], loc[j], weight[j]
                    if j_weight > 2e-2 or fix:
                        j_L = j_weight * lognorm.pdf(data, j_shape, loc=j_loc, scale=j_scale)
                        i_L = list(map(add, i_L, j_L))
                        j_pdf = j_weight * lognorm.pdf(np.exp(x_bins), j_shape, loc=j_loc, scale=j_scale)
                        i_pdf_each.append(j_pdf)
                        i_pdf = list(map(add, i_pdf, j_pdf))
                        i_shape.append(j_shape)
                        i_scale.append(j_scale)
                        i_loc.append(j_loc)
                        i_weight.append(j_weight)
                        i_gen.append(np.round(np.log2(mean_val / j_scale)))

                min_val = min(i for i in i_L if i > 0)
                i_L = [i if i > 0 else min_val for i in i_L]
                aic = -2*sum(np.log(i_L)) + len(i_shape) * 3 + (len(i_shape) - 1)

            m_scale.append(i_scale)
            m_shape.append(i_shape)
            m_loc.append(i_loc)
            m_weight.append(i_weight)
            m_gen.append(i_gen)
            m_aic.append(aic)
            m_pdf.append(i_pdf)
            m_pdf_each.append(i_pdf_each)
    ind_best = np.argmin(m_aic)
    scale_best, shape_best, loc_best, weight_best, gen_best, pdf_best, pdf_each_best = \
        m_scale[ind_best], m_shape[ind_best], m_loc[ind_best], m_weight[ind_best], m_gen[ind_best], m_pdf[ind_best], m_pdf_each[ind_best]
    if len(scale_best) == 1:
        cellnum = [len(data)]
    else:
        cellnum = [i_weight * len(data) for i_weight in weight_best]

    return cellnum, gen_best, scale_best, shape_best, loc_best, weight_best, m_aic, pdf_best, pdf_each_best, m_comb


def DPM_fcs_analysis_fitlognorm_OP_fun(x, data, combination):
    # x = np.append(x, 1-sum(x[len(combination)*3:]))

    weight = x[-len(combination)+1:]
    weight = np.append(weight, 1-sum(weight))
    x = x[:-len(combination)+1]
    scale0 = x[0]
    x = np.delete(x, 0)
    shape = x[:-len(combination)]
    loc = x[-len(combination):]

    L, scale = [0 for _ in range(len(data))], []
    for i, i_gen in enumerate(combination):
        i_scale = scale0 if i == 0 else scale0/(2**(combination[i]-combination[0]))
        i_shape, i_loc, i_weight = shape[i], loc[i], weight[i]
        # i_scale, i_shape, i_loc, i_weight = x[i::len(combination)]

        i_L = i_weight * lognorm.pdf(data, i_shape, loc=i_loc, scale=i_scale)
        L = list(map(add, L, i_L))
        scale.append(i_scale)
    pena = 0  # sum(weight[:-1]) > 0.9
    # for i in range(len(scale) - 1):
    #     i_ratio = combination[i+1]-combination[i]
    #     i_pena = 0 if 2**i_ratio * 0.9 < scale[i]/scale[i+1] < 2**i_ratio * 1.1 else 1
    #     pena += i_pena

    min_val = min(i for i in L if i > 0)
    L = [i if i > 0 else min_val for i in L]

    sumlogL = -sum(np.log(L))
    sumlogL += sumlogL + sumlogL*pena

    return sumlogL


def DPM_fcs_model_fit(data, model):
    stepsize = STEPSIZE_ODE_FCS  # sim stepsize 0.1 day
    model1_par, model2_par, model3_par = None, None, None
    if 1 in model:
        DPM_fcs_model = DPM_fcs_ode_model1
        g0, b0, D0 = [0.8], [0.9], [5e-2]
        g_lb, b_lb, D_lb = [0.1], [0.5], [1e-4]
        g_ub, b_ub, D_ub = [2.0], [1.0], [1e-1]
        x0 = np.asarray(g0 + b0 + D0)
        lb = g_lb + b_lb + D_lb
        ub = g_ub + b_ub + D_ub
        bounds = Bounds(lb, ub)
        x_op = minimize(DPM_fcs_model_fit_op, x0, args=(DPM_fcs_model, data, stepsize), bounds=bounds)
        x = x_op.x
        model1_par = {'g': x[0], 'b': x[1], 'D': x[2]}
    if 2 in model:
        DPM_fcs_model = DPM_fcs_ode_model2
        g0, b0, D0 = [0.8], [0.9], [5e-2]
        g_lb, b_lb, D_lb = [0.1], [0.5], [1e-4]
        g_ub, b_ub, D_ub = [2.0], [1.0], [1e-1]
        x0 = np.asarray(g0 + b0 + D0)
        lb = g_lb + b_lb + D_lb
        ub = g_ub + b_ub + D_ub
        bounds = Bounds(lb, ub)
        x_op = minimize(DPM_fcs_model_fit_op, x0, args=(DPM_fcs_model, data, stepsize), bounds=bounds)
        x = x_op.x
        model2_par = {'g': x[0], 'b': x[1], 'D': x[2]}
    if 3 in model:
        DPM_fcs_model = DPM_fcs_ode_model3
        g0, b0, s0 = [0.8], [0.9], [5e-2]
        g_lb, b_lb, s_lb = [0.1], [0.5], [1e-3]
        g_ub, b_ub, s_ub = [2.0], [1.0], [2e-1]
        x0 = np.asarray(g0 + b0 + s0)
        lb = g_lb + b_lb + s_lb
        ub = g_ub + b_ub + s_ub
        bounds = Bounds(lb, ub)
        x_op = minimize(DPM_fcs_model_fit_op2, x0, args=(DPM_fcs_model, data, stepsize), bounds=bounds)
        x = x_op.x
        model3_par = {'g': x[0], 'b': x[1], 's': x[2]}

    return model1_par, model2_par, model3_par


def DPM_fcs_model_fit_op(par, DPM_fcs_model, data, stepsize):
    g, b, D = par
    tpoint = [i for i in data.keys() if type(i) == int or type(i) == float]
    t_end = tpoint[-1] + 1  # last experiment day plus one
    t = np.arange(0, t_end + stepsize, stepsize)
    argsODE = (g, b, D, data[0]['cellnum'])
    max_step = stepsize
    x0 = np.zeros(FCS_maxgen)
    x0[0] = data[0]['cellnum']
    X = solve_ivp(DPM_fcs_model, t[[0, -1]], list(x0), method='RK45', rtol=SCIPY_INTEGRATE_SOLVE_IVP['RTOL'],
                  atol=SCIPY_INTEGRATE_SOLVE_IVP['ATOL'], t_eval=t, args=argsODE, max_step=max_step)
    t = X.t
    x = X.y
    mse = 0
    for i_t in tpoint[1:]:
        pos = next(x for x, val in enumerate(t) if val >= i_t)
        x_i_t = x[:, pos]
        cellnum = data[i_t]['cellnum']
        gen = [int(i_gen) for i_gen in data[i_t]['gen']]
        for i_gen, i_num in enumerate(x_i_t):
            if i_gen not in set(gen):
                mse += (i_num-100)**2
            elif i_gen in set(gen):
                mse += (i_num-cellnum[gen.index(i_gen)])**2

    return mse


def DPM_fcs_model_fit_op2(par, DPM_fcs_model, data, stepsize):
    g, b, s = par
    tpoint = [i for i in data.keys() if type(i) == int or type(i) == float]
    t_end = tpoint[-1] + 1  # last experiment day plus one
    t = np.arange(0, t_end + stepsize, stepsize)
    argsODE = (g, b, s)
    max_step = stepsize
    x0 = np.zeros(2*FCS_maxgen)
    x0[0] = data[0]['cellnum']
    X = solve_ivp(DPM_fcs_model, t[[0, -1]], list(x0), method='RK45', rtol=SCIPY_INTEGRATE_SOLVE_IVP['RTOL'],
                  atol=SCIPY_INTEGRATE_SOLVE_IVP['ATOL'], t_eval=t, args=argsODE, max_step=max_step)
    t = X.t
    x = X.y
    mse = 0
    for i_t in tpoint[1:]:
        pos = next(x for x, val in enumerate(t) if val >= i_t)
        x_i_t = x[:, pos]
        cellnum = data[i_t]['cellnum']
        gen = [int(i_gen) for i_gen in data[i_t]['gen']]
        for i_gen in range(FCS_maxgen):
            i_num = x_i_t[i_gen] + x_i_t[i_gen + FCS_maxgen]
            if i_gen not in set(gen):
                mse += (i_num-100)**2
            elif i_gen in set(gen):
                mse += (i_num-cellnum[gen.index(i_gen)])**2

    return mse


def DPM_fcs_ode_model1(t, x, g, b, D, N0):
    f = np.zeros((FCS_maxgen, FCS_maxgen), dtype=float)
    f[0, 0] = -g
    for i in range(1, len(x)):
        f[i, i-1], f[i, i] = 2*g*b, -g
    dx = f @ x
    dx[0] = dx[0] + g*D*N0
    dx[1] = dx[1] - 2*g*b*D*N0

    return list(dx)


def DPM_fcs_ode_model2(t, x, g, b, D, N0):
    f = np.zeros((FCS_maxgen, FCS_maxgen), dtype=float)
    f[0, 0] = -g*(1-D)
    for i in range(1, len(x)):
        f[i, i-1], f[i, i] = 2*g*b, -g*(1-D)
    dx = f @ x

    return list(dx)


def DPM_fcs_ode_model3(t, x, g, b, s):
    f = np.zeros((2*FCS_maxgen, 2*FCS_maxgen), dtype=float)
    f[0, 0] = -g
    f[FCS_maxgen, 0] = g*s
    for i in range(1, FCS_maxgen):
        f[i, i-1], f[i, i] = 2*g*b, -g
        f[FCS_maxgen+i, i] = g*s
    dx = f @ x

    return list(dx)


def DPM_fcs_plot_model_12(modelpar, data):
    g, b, D = modelpar['g'], modelpar['b'], modelpar['D']
    tpoint = [i for i in data.keys() if type(i) == int or type(i) == float]
    t_end = tpoint[-1] + 1  # last experiment day plus one
    t = np.arange(0, t_end + STEPSIZE_ODE_FCS, STEPSIZE_ODE_FCS)
    argsODE = (g, b, D, data[0]['cellnum'])
    max_step = STEPSIZE_ODE_FCS
    x0 = np.zeros(FCS_maxgen)
    x0[0] = data[0]['cellnum']
    X = solve_ivp(DPM_fcs_ode_model1, t[[0, -1]], list(x0), method='RK45', rtol=SCIPY_INTEGRATE_SOLVE_IVP['RTOL'],
                  atol=SCIPY_INTEGRATE_SOLVE_IVP['ATOL'], t_eval=t, args=argsODE, max_step=max_step)
    t = X.t
    x = X.y
    palette = sns.color_palette(None, x.shape[0])
    plt.rcParams['figure.figsize'] = (14, 8)
    plt.rcParams['font.size'] = 14
    for i_gen in range(x.shape[0]):
        plt.plot(t, x[i_gen, :], color=palette[i_gen], label='generation ' + str(i_gen))
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=6, prop={'size': 'small'}, frameon=False)
    tpoint = [i for i in data.keys() if type(i) == int or type(i) == float]
    for i_t in tpoint:
        cellnum = data[i_t]['cellnum']
        gen = [int(i_gen) for i_gen in data[i_t]['gen']] if hasattr(data[i_t]['gen'], '__iter__') else int(data[i_t]['gen'])
        if hasattr(gen, '__iter__'):
            for i, i_gen in enumerate(gen):
                plt.plot(i_t, cellnum[i], marker=11, color=palette[i_gen], markersize=12)
        else:
            plt.plot(i_t, cellnum, marker=11, color=palette[gen], markersize=12)
        plt.axvline(i_t, 0, 2e4, color='k', linestyle='-')
    plt.ylabel('Cell number')
    plt.xlabel('Time [days]')

    return


def DPM_fcs_plot_model_3(modelpar, data):
    g, b, s = modelpar['g'], modelpar['b'], modelpar['s']
    tpoint = [i for i in data.keys() if type(i) == int or type(i) == float]
    t_end = tpoint[-1] + 1  # last experiment day plus one
    t = np.arange(0, t_end + STEPSIZE_ODE_FCS, STEPSIZE_ODE_FCS)
    argsODE = (g, b, s)
    max_step = STEPSIZE_ODE_FCS
    x0 = np.zeros(2 * FCS_maxgen)
    x0[0] = data[0]['cellnum']
    X = solve_ivp(DPM_fcs_ode_model3, t[[0, -1]], list(x0), method='RK45', rtol=SCIPY_INTEGRATE_SOLVE_IVP['RTOL'],
                  atol=SCIPY_INTEGRATE_SOLVE_IVP['ATOL'], t_eval=t, args=argsODE, max_step=max_step)
    t = X.t
    x = X.y
    palette = sns.color_palette(None, x.shape[0])

    plt.rcParams['figure.figsize'] = (12, 8)
    plt.rcParams['font.size'] = 12
    fig, (ax1, ax2) = plt.subplots(2, 1)
    #plt.tight_layout()

    for i_gen in range(FCS_maxgen):
        x_i = x[i_gen, :] + x[i_gen + FCS_maxgen, :]
        ax1.plot(t, x_i, color=palette[i_gen], label='generation ' + str(i_gen))
        ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=6, prop={'size': 'small'}, frameon=False)

        ax2.plot(t, x[i_gen, :], color=palette[i_gen], label='generation C' + str(i_gen))
        ax2.plot(t, x[i_gen + FCS_maxgen, :], linestyle='--', color=palette[i_gen], label='generation D' + str(i_gen))
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=6, prop={'size': 'small'}, frameon=False)

    tpoint = [i for i in data.keys() if type(i) == int or type(i) == float]
    for i_t in tpoint:
        cellnum = data[i_t]['cellnum']
        gen = [int(i_gen) for i_gen in data[i_t]['gen']] if hasattr(data[i_t]['gen'], '__iter__') else int(data[i_t]['gen'])
        if hasattr(gen, '__iter__'):
            for i, i_gen in enumerate(gen):
                ax1.plot(i_t, cellnum[i], marker=11, color=palette[i_gen], markersize=12)
        else:
            ax1.plot(i_t, cellnum, marker=11, color=palette[gen], markersize=12)
        ax1.axvline(i_t, 0, 2e4, color='k', linestyle='-')
    ax1.set_ylabel('Cell number')
    ax2.set_ylabel('Cell number')
    ax2.set_xlabel('Time [days]')
    fig.text(3, 2e4, 'model 3', horizontalalignment='center')

    return