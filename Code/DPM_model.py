from DPM_lib import expm, solve_ivp, odeint
from DPM_print import *
from DPM_constant import *
from DPM_miscellaneous import DPM_miscellaneous_treatment_change_time
""" This script contains the model functions and performs the model simulation. """
'''
CHECKED
'''

solver = SOLVER


# Simulate the model.#
def DPM_sim_model(X0, T, g0, Sa, d, t, maxthreshold, LSsim, Simtimestep_control=True, index_celltype_limit=None, limit_celltype=None,
                  index_celltype_alter=None, index_celltype_alter_type=None):
    success = False
    X_i = t_i = d_i = None
    global solver
    solver_used = set()
    if not USE_MATRIXEXP:
        solver_used.add(solver)
        X_i, t_i, d_i, success = DPM_sim_model2012_scipy_integrate_odeint(
            X0, T, g0, Sa, d, t, maxthreshold, LSsim, Simtimestep_control, index_celltype_limit, limit_celltype, index_celltype_alter,
            index_celltype_alter_type)
        while not success and len(SOLVER_ALL.difference(solver_used)) != 0:
            solver = list(SOLVER_ALL.difference(solver_used))[0]
            X_i, t_i, d_i, success = DPM_sim_model2012_scipy_integrate_odeint(
                X0, T, g0, Sa, d, t, maxthreshold, LSsim, Simtimestep_control, index_celltype_limit, limit_celltype, index_celltype_alter,
                index_celltype_alter_type)
    solver = SOLVER
    if (not success) or USE_MATRIXEXP:
        X_M_i, t_M_i, d_M_i = DPM_sim_model2012(
            X0, T, g0, Sa, d, t, maxthreshold, LSsim, Simtimestep_control, index_celltype_limit, limit_celltype, index_celltype_alter,
            index_celltype_alter_type)
        X_i = X_M_i
        t_i = t_M_i
        d_i = d_M_i

    return X_i, t_i, d_i


# Simulate the model from PNAS(2012).#
def DPM_sim_model2012(X0, T, g0, Sa, d, t, maxthreshold, LSsim, Simtimestep_control=True, index_celltype_limit=None, limit_celltype=None,
                      index_celltype_alter=None, index_celltype_alter_type=None):
    if len(t) == 1:
        return np.reshape(X0, (len(X0), 1)), t, np.empty([d.shape[0], 0], dtype=float)
    if t[0] != 0:
        print('Simulation starting time point is not equal 0.')
        DPM_print_errorandclose()
        exit()
    elif np.max(np.diff(t)) > 1 and Simtimestep_control:
        print('The simulation timestep is greater than 1 day, but it must be less than 1 day. We need to detect when any cell type crosses a '  
              'population size of 1 and update the mask matrix accordingly. If the timestep exceeds 1 day, this crossing point may be detected '
              'inaccurately, which can cause the simulation to fail.')
        DPM_print_errorandclose()
        exit()

    # Simulation treatemnt solutions.#
    Num_cell_type = X0.shape[0]
    I_identity = np.identity(Num_cell_type)
    pos_changedose = DPM_miscellaneous_treatment_change_time(d)
    nintervals = 1
    intervals = np.zeros(len(t)-1)
    for i in range(len(pos_changedose)-1):
        intervals[pos_changedose[i]:pos_changedose[i+1]] = nintervals
        nintervals = nintervals + 1
    intervals[pos_changedose[-1]:] = nintervals

    # Apply the analytic exponential solution within each interval.#
    # For each constant-drug interval, identify the subintervals where any population falls below 1 or exceeds 1.#
    X = np.zeros((Num_cell_type, len(t)), dtype=float)
    X[:, 0] = X0
    mask = np.ones((Num_cell_type, Num_cell_type), dtype=float)

    # Mask the terms associated with initial populations below 1.#
    for i in range(Num_cell_type):
        mask[:, i] = 0 if X0[i] < 1 else 1

    # Initial value of X for the treatment.#
    xinit = np.full(Num_cell_type, -1, dtype=float)
    k = -1
    i = -1
    for i in range(intervals.shape[0]):
        current_interval = intervals[i]
        # Starting point of a new doseage combination.#
        if i == np.min(np.where(intervals == current_interval)):
            xinit = X[:, i]
            k = i
        # Apply the matrix exponential relative to the initial condition.#
        # Advance at least one simulation timestep (e.g., 1 day), because a treatment dose must be applied for at least one timestep:#
        # a dose given at day i is used from day i to day i + one timestep.#
        # Afterward, multiply the matrix by the mask matrix to zero out vanished terms. Columns in the matrix is masked.#
        tdiff = t[i+1] - t[k]
        f = (I_identity + T) @ np.diag(g0 * tdiff) - np.diag(Sa @ d[:, i] * tdiff)
        f = f * mask
        X[:, i+1] = expm(f) @ xinit
        # If the simulated cell number becomes too large or all cell numbers fall below 1 (stopping further doubling by definiation), #
        # terminate the simulation.#
        if (X[:, i+1].sum() >= maxthreshold or np.less(X[:, i+1], 1).all()) and LSsim:
            return X[:, :i+2], t[:i+2], d[:, :i+1]
        # If an 'argv' input is provided, compute the time at which τinc, τS, τR1, τR2,... and τR_multi_resistant reach a specific cell number.#
        if index_celltype_limit is not None:
            if X[index_celltype_limit, i+1] >= limit_celltype:
                return X[:, :i+2], t[:i+2], d[:, :i+1]

        # If 'index_celltype_alter' is provided, check whether the specified subpopulation crosses the boundary of 1.#
        # If it does, stop the simulation and return.#
        # 'index_celltype_alter' gives the index of the subpopulation to monitor; its list position must correspond to True.#
        # 'index_alter_type' specifies the direction of the crossing:#
        # 1 -> crossing from below 1 to above 1,#
        # -1 -> crossing from above 1 to below 1.#
        alter = np.zeros(Num_cell_type, dtype=int)
        for i_type in range(Num_cell_type):
            alter[i_type] = 1 if X[i_type, i] < 1 <= X[i_type, i+1] else -1 if X[i_type, i+1] < 1 <= X[i_type, i] else alter[i_type]
            if index_celltype_alter is not None:
                if np.equal(alter[index_celltype_alter], index_celltype_alter_type).any():
                    return X[:, :i+2], t[:i+2], d[:, :i+1]

        # If any of the subpopulation crosses the boundary of 1, set a new initial condition and update the mask matrix.#
        # Identify which subpopulations have crossed the boundary.#
        # A boundary-crossing event has occurred.#
        if np.any(alter):
            # Update the mask matrix.#
            for i_type in range(Num_cell_type):
                if np.equal(alter[i_type], 1):
                    mask[:, i_type] = 1
                elif np.equal(alter[i_type], -1):
                    mask[:, i_type] = 0
            # Set a new initial condition.#
            xinit = X[:, i+1]
            k = i+1

    Xout = X[:, :i+2]
    tout = t[:i+2]
    dout = d[:, :i+1]
    return Xout, tout, dout


# Simulate the solution of the PNAS(2012) model.#
def DPM_sim_model2012_scipy_integrate_odeint(X0, T, g0, Sa, d, t, maxthreshold, LSsim, Simtimestep_control=True, index_celltype_limit=None,
                                             limit_celltype=None, index_celltype_alter=None, index_celltype_alter_type=None):
    if any(elem is not None for elem in (index_celltype_limit, limit_celltype, index_celltype_alter, index_celltype_alter_type)):
        return np.reshape(X0, (len(X0), 1)), t, d, False
    if len(t) == 1:
        return np.reshape(X0, (len(X0), 1)), t, np.empty([d.shape[0], 0], dtype=float), True
    if t[0] != 0:
        print('Simulation starting time point is not equal 0.')
        DPM_print_errorandclose()
        exit()
    elif np.max(np.diff(t)) > 1 and Simtimestep_control:
        print("Simulation timestep is greater than 1 (day), but it must be less than 1 day. We need to detect when any cell type's population "
              "crosses 1 and update the mask matrix. If the timestep is larger than 1 day, the crossing time may be too inaccurate, which can "
              "cause the simulation fail.")
        DPM_print_errorandclose()
        exit()

    Num_cell_type = X0.shape[0]
    I_identity = np.identity(Num_cell_type)
    pos_changedose = DPM_miscellaneous_treatment_change_time(d)

    argsODE = (T, g0, Sa, d, Num_cell_type, I_identity, pos_changedose)

    hmax = max_step = len(t) if len(pos_changedose) == 1 else np.min(np.diff(t))

    X = ode_success = None
    # print(solver)
    if solver == 'ODEINT':
        #  hmax=SCIPY_INTEGRATE_ODEINT['hmax']
        result_odeint = odeint(DPM_sim_ODE_model2012, list(X0), t, argsODE,
                               rtol=SCIPY_INTEGRATE_ODEINT['RTOL'], atol=SCIPY_INTEGRATE_ODEINT['ATOL'],
                               full_output=True, tfirst=True, hmax=hmax)
        X = result_odeint[0].T
        ode_success = True if result_odeint[1]['message'] in ('Integration successful.', 'DPMpass')\
            else False
    elif solver == 'SOLVE_IVP':
        # tspan = float(np.diff(t[[0, -1]]))
        # max_step=SCIPY_INTEGRATE_SOLVE_IVP['max_step']
        result_solve_ivp = solve_ivp(DPM_sim_ODE_model2012, t[[0, -1]], list(X0), method='RK45',
                                     rtol=SCIPY_INTEGRATE_SOLVE_IVP['RTOL'],
                                     atol=SCIPY_INTEGRATE_SOLVE_IVP['ATOL'], t_eval=t, args=argsODE, max_step=max_step)
        X = result_solve_ivp[list(result_solve_ivp)[1]]
        ode_success = result_solve_ivp[list(result_solve_ivp)[10]]

    if (X[:, -1].sum() >= maxthreshold or np.less(X[:, -1], 1).all()) and LSsim:
        pos = np.argmax(np.sum(X, axis=0) >= maxthreshold) if X[:, -1].sum() >= maxthreshold else \
            np.argmax(np.less(X, 1).all(axis=0)) if np.less(X[:, -1], 1).all() else 0
        # if pos == 0:
        #     print("pos value should not equal to 0.")
        #     DPM_print_errorandclose()
        #     exit()
        return X[:, :pos+1], t[:pos+1], d[:, :pos], ode_success

    return X, t, d, ode_success


def DPM_sim_ODE_model2012(t, x, T, g0, Sa, d, Num_cell_type, I_identity, pos_changedose):
    d_use = DPM_sim_d_use(d, t, pos_changedose)
    mask = np.ones((Num_cell_type, Num_cell_type), dtype=float)

    # Mask the terms associated with initial population below 1.#
    for i in range(Num_cell_type):
        mask[:, i] = 0 if x[i] < 1 else 1

    f = (I_identity + T) @ np.diag(g0) - np.diag(Sa @ d_use)
    f = f * mask
    dx = f @ x
    return list(dx)


def DPM_sim_d_use(d, t, pos_changedose):
    d_use = d[:, -1]
    if len(pos_changedose) == 1:
        d_use = d[:, pos_changedose[0]]
    else:
        for i in range(len(pos_changedose) - 1):
            if pos_changedose[i] <= t < pos_changedose[i + 1]:
                d_use = d[:, pos_changedose[i]]
                break
    return d_use
