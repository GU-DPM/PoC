import numpy as np
from scipy.integrate import simps
from scipy.stats import rv_continuous, loguniform
import io
import sys
import time
import pickle
# from scipy.linalg import expm
from collections import namedtuple
from colorama import Fore
from tqdm import tqdm
# import math
import matplotlib.pyplot as plt
import matplotlib as mpl
# mpl.use('Qt5Agg')

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.problems import get_problem
from pymoo.core.evaluator import Evaluator
from pymoo.core.population import Population
from pymoo.optimize import minimize

problem = get_problem("zdt2")

# create initial data and set to the population object
X = np.random.random((10, problem.n_var))
pop = Population.new("X", X)
Evaluator().eval(problem, pop)

# effective mutation rate
k = 7.1 * 10**(-7)
k = 0
x = np.linspace(loguniform.ppf(0.01, k, 1e-1),
                loguniform.ppf(0.99, k, 1e-1), 1000)
plt.plot(x, loguniform.pdf(x, k, 1e-1),
         'r-', lw=2, alpha=0.6, label='loguniform pdf')
plt.xscale("log")

def p(x, k, LOD):
    return (k/x**2)/(1-k/LOD+k)

# LODt = [10**(-6), 10**(-5), 10**(-4), 10**(-3), 10**(-2), 10**(-1)]
# for LOD in LODt:
#     x = np.linspace(k, LOD, int(1e8))
#     print(simps(p(x, k, LOD), x))

class mutationFraction(rv_continuous):
    def _pdf(self, x, k, LOD, const):
        return (1.0/const) * p(x, k, LOD)

LOD = 10**(-1)
x = np.linspace(k, LOD, int(1e4))
mutationFraction_distribution = mutationFraction(name='mutationFraction_distribution', a=k, b=LOD)
pdf = mutationFraction_distribution.pdf(x=x, k=k, LOD=LOD, const=1)
cdf = mutationFraction_distribution.cdf(x=x,  k=k, LOD=LOD, const=1)
samples = mutationFraction_distribution.rvs(k=k, LOD=LOD, const=1, size=int(1e4))


plt.plot(x, pdf,
         'k-', lw=2, alpha=0.6, label='pdf')

ci_low = mutationFraction_distribution.ppf(.025, k=k, LOD=LOD, const=1)
ci_high = mutationFraction_distribution.ppf(.975, k=k, LOD=LOD, const=1)

plt.rcParams['figure.figsize'] = (14, 7)
plt.rcParams.update({'font.size': 16})
fig, (ax_pdf, ax_cdf) = plt.subplots(1, 2)
ax_pdf.plot(x, pdf,  lw=2, label='Probability density')
ax_pdf.hist(samples, bins=np.linspace(k, LOD, int(2e2)), density=True, label='Histogram of samples')
ax_pdf.set_xlabel('Percentage of mutation cells')
ax_pdf.set_xticks(np.linspace(k, LOD, 5))
ax_pdf.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1e'))
ax_pdf.legend(loc='best', ncol=1, frameon=0, fontsize=16)
ax_cdf.plot(x, cdf)
ax_cdf.set_xticks(np.linspace(k, LOD, 5))
ax_cdf.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1e'))
ax_cdf.set_ylabel('Cumulative distribution')
ax_cdf.set_xlabel('Percentage of mutation cells')
plt.tight_layout()
plt.suptitle(f'k={k}, LOD = {LOD:.1e}, 95% IC ({ci_low:.1e}, {ci_high:.1e})')
plt.close()


# Two drug cases
# parameters from ref "Long Range Personalized Cancer Treatment Strategies Incorporating Evolutionary Dynamics"

# S: sensitive cells, sensitive to both drugs
# R1: resistant to drug 1, sensitive to drug 2
# R2: resistant to drug 2, sensitive to drug 1
# R12: resistant to both drugs
# N: total cells, S+R1+R2+R12

# drug effect: For S, drug 1 >= drug 2, always by setting.

#              For R1:
#              (1) drug 1 < drug 2, e.g., Sa(R1,D1)/Sa(S,D1) = 1e-5 and Sa(S,D1)/g0 = 4.7683，
#              then Sa(R1,D1)/g0 = 4.7683e-5. Because Sa(R1,D2)/g0 = Sa(S,D2)/g0, if Sa(S,D2)/Sa(S,D1) = 0.2714,
#              then Sa(S,D2)/g0 = Sa(S,D2)/Sa(S,D1) * Sa(S,D1)/g0 = 0.2714 * 4.7683 = 1.294.
#              So Sa(R1,D1)/g0 = 4.7683e-5 and Sa(R1,D2)/g0 = 1.294, effect of drug 1 < drug 2
#              (2) drug 1 > drug 2, e.g., Sa(R1,D1)/Sa(S,D1) = 0.8 and Sa(S,D1)/g0 = 4.7683,
#              then Sa(R1,D1)/g0 = 3.81464. Because Sa(R1,D2)/g0 = Sa(S,D2)/g0, if Sa(S,D2)/Sa(S,D1) = 0.2714,
#              then Sa(S,D2)/g0 = Sa(S,D2)/Sa(S,D1) * Sa(S,D1)/g0 = 0.2714 * 4.7683 = 1.294.
#              So Sa(R1,D1)/g0 = 3.81464 and Sa(R1,D2)/g0 = 1.294, effect of drug 1 > drug 2.

#              For R2, drug 1 > drug 2, always.
#              Because Sa(R2,D2)/g0 = Sa(R2,D2)/Sa(S,D2) * Sa(S,D2)/Sa(S,D1) * Sa(S,D1)/g0
#              the maximum value of Sa(R2,D2)/Sa(S,D2) = 0.8, the smallest drug 2 resistant level
#              the maximum value of Sa(S,D2)/Sa(S,D1) = 1
#              then maximum value of Sa(R2,D2)/g0 = 0.8 * 1 * Sa(S,D1)/g0 = 0.8 * Sa(S,D1)/g0
#              And Sa(R2,D1)/g0 = Sa(S,D1)/g0
#              so Sa(R2,D1)/g0 > Sa(R2,D2)/g0.

#              For R12, if Sa(R12,D1) = Sa(R12,D2) = 0, then Sa(R12,D1)/g0 = Sa(R12,D2)/g0 = 0
#              Effect of drug 1 = drug 2 = 0
#              Otherwise, if Sa(R12,D1) = Sa(R1,D1) and Sa(R12,D2) = Sa(R2,D2)
#              (1) drug 1 < drug 2:
#              e.g., Sa(R1,D1)/Sa(S,D1) = 1e-5 and Sa(S,D1)/g0 = 4.7683,
#              then Sa(R1,D1)/g0 = 4.7683e-5 = Sa(R12,D1)/g0
#              If Sa(S,D2)/Sa(S,D1) = 0.2714 and Sa(R2,D2)/Sa(S,D2) = 0.8, then
#              Sa(S,D2)/g0 = Sa(S,D2)/Sa(S,D1) * Sa(S,D1)/g0 = 0.2714 * 4.7683 = 1.294. And
#              Sa(R12,D2)/g0 = Sa(R2,D2)/g0 = Sa(R2,D2)/Sa(S,D2) *  Sa(S,D2)/Sa(S,D1) * Sa(S,D1)/g0 =
#              0.8 * 0.2714 * Sa(S,D1)/g0 = 0.217 * Sa(S,D1)/g0 = 0.217 * 4.7683 = 1.0347211
#              So Sa(R12,D1)/g0 = 4.7683e-5 and Sa(R12,D2)/g0 = 1.0347211, drug 1 < drug 2
#              (2) drug 1 > drug 2:
#              e.g., Sa(R1,D1)/Sa(S,D1) = 0.8 and Sa(S,D1)/g0 = 4.7683, then
#              Sa(R1,D1)/g0 = Sa(R12,D1)/g0 = 3.81464
#              Because Sa(R12,D2)/g0 = Sa(R2,D2)/g0 = Sa(R2,D2)/Sa(S,D2) * Sa(S,D2)/Sa(S,D1) * Sa(S,D1)/g0
#              If Sa(R2,D2)/Sa(S,D2) = 1e-5, Sa(S,D2)/Sa(S,D1) = 0.2714, then
#              Sa(S,D2)/g0 = Sa(S,D2)/Sa(S,D1) * Sa(S,D1)/g0 = 0.2714 * 4.7683 = 1.294. And
#              Sa(R12,D2)/g0 = Sa(R2,D2)/g0 = 1e-5 * 0.2714 * 4.7683 = 1.294e-5
#              So Sa(R12,D1)/g0 = 3.81464 and Sa(R12,D2)/g0 = 1.294e-5, drug 1 > drug 2

# Excluding parameter combinations
# (1) PNAS(2012): (i) any combination in which one of the two drugs was completely ineffective against all cell types,
# because there are no strategic choices in this case.

# (2) PNAS(2012): (ii) any combination in which all treatment strategies resulted in survival greater than 4 y.
# All treatment strategies means strategy 0, strategy 1, strategy 2.1, strategy 2.2, strategy 3, strategy 4.
# Because the simulation was truncated at 5 y, and because a significant survival difference required at least a 25%
# improvement compared with a reference, it was not possible to ascertain if there was improvement if all strategies
# resulted in survival for longer than 4 y.

# PNAS(0) The potency of drugs on sensitive cells follows: drug 1 >= drug 2.

# According to Prof. Yeang, only 3 filtering criteria were implemented actually in PNAS(2012):
# The parameter combinations satisfying criteria (1)-(3) are included
# PNAS(1) There are S cells at t = 0, which means R1(0)/N(0) + R2(0)/N(0) < 1.

# PNAS(2) Sa(R1,D1) <= g0 and Sa(R2,D2) <= g0, which means drug 1 or drug 2 kills R1 or R2 is slower than its growth rate,
# thus no drug combinations can cure R12.

# PNAS(3) Sa(S,D1) > 0 and S(S,D2) > 0. Drug 1 and Drug 2 has at least some effects

# Sa(S,D1) < g0 and Sa(S,D2) < g0 are allowed. If Sa(S,D1) < g0, then one cannot cure patients, but DPM may still
# prolong the life span of patients relative to other treatment strategies.

# Including parameter combinations
# Biology Direct(2016): The parameter combinations satisfying criteria (4)-(6) are included

# Biology Direct(4) "Sensitive cell types can be eradicated by each drug (Sa(S,1) > g0, Sa(S,2) > g0).
# Thus in contrast Beckman, Schemman, and Yeang (2012), we consider only fundamentally curable cases."
# Sa(S,D1) > g0 and Sa(S,D2) > g0, both inequalities hold simultaneously.

# Biology Direct(5) Multiply-resistant cell types cannot be eradicated by any drug (Sa(R12,1) < g0, Sa(R12,2) < g0)
# (Sa(R12,1) < g0, Sa(R12,2) < g0) both inequality holds simultaneously.

# wrong: parameter combination satisfying (Sa(R12,1) < g0, Sa(R12,2) < g0) are included, the words below talk about
# the situation if parameter combination satisfying (Sa(R12,1) < g0, Sa(R12,2) < g0) are excluded
# --------------------------------------------------------------------------------------------------
# It means that either Sa(R12,1) > g0 or Sa(R12,2) > g0
# If Sa(R12,1) > g0, Sa(R1,1) = Sa(R12,1) > g0, Sa(S,1) > g0, and Sa(R2,1) = Sa(S,1) > g0,
# then treat with drug 1 tumor will shrink, need no strategy.
# If Sa(R12,2) > g0, Sa(R2,2) = Sa(R12,2) > g0, Sa(S,2) > Sa(R2,2) > g0, Sa(R1,2) = Sa(S,2) > g0,
# then treat with drug 2 tumor will shrink, need no strategy.
# --------------------------------------------------------------------------------------------------

# Biology Direct(6) The patient can be cured by simultaneous full dosages of all drugs
# (an invalid option due to toxicity) but cannot be cured by any valid static treatment.

# 3 drug case:
# Biology Direct(7) The potency of drugs on sensitive cells follows: drug 1 >= drug 2 >= drug 3
# based sensitivity on S cells

# Biology Direct(8) The initial size of a doubly-resistant subpopulation (R12,R13,R23) is no greater than
# the initial sizes of two singly-resistant subpopulations from which it is derived. This condition is not necessary
# for the model, but was done to reduce the computational burden.
# Translated to (R12 <= R1, R12 <= R2, R13 <= R1, R13 <= R3, R23 <= R2, R23 <= R3) all inequality holds simultaneously.

# Strategy 0: Current personalized medicine:
#   (1) Initially treat the patient with drug 1 alone if R1/(S+R1+R2+R12) <= 0.5, i.e., the R1 population does not
#   dominate the tumor. Because initial R12 = 0, so R1/(S+R1+R2+R12) = R1/(S+R1+R2) <= 0.5. S+R2 > R1, S and R2 are more
#   sensitive to drug 1, so treat patient drug 1. So the nadir is S+R1+R2 at T0?
#
#   If R1/(S+R1+R2+R12) > 0.5, treat the patient with drug 2. Don't understand, even R1 population dominate the tumor,
#   R1 can still be more sensitive to drug 1 compared to drug 2, why definitely treat drug 1？Hold this ture first
#
#   Maintain the current treatment until either one of the following events:
#   (1) The total population reaches twice the nadir population. So at the beginning nadir is S+R1+R2 at T0?
#   g0 = 0.34, Sa(R1,D1)/Sa(S,D1) = 1e-5, Sa(S,D1)/g0 = 4.7683 => Sa(R1,D1)/g0 = 4.7683e-5
#   Sa(S,D2)/Sa(S,D1) = 4e−4, Sa(R1,D2)/g0 = Sa(S,D2)/g0,
#   Sa(S,D2)/g0 = Sa(S,D2)/Sa(S,D1) * Sa(S,D1)/g0 = 4e-4 * 4.7683 = 0.00191
#   Sa(R2,D1)/g0 = Sa(S,D1)/g0 = 4.7683, R1(0)/N = 0.9, R2(0)/N = 0.1
#   R1/(S+R1+R2+R12) > 0.5
#   first use drug 2, Sa(R1,D1)/g0 = 4.7683e-5, Sa(S,D2)/g0 = 0.00191
#   About 2.1 day, tumor double, but R1 is still the majority, R1(2.1)/N(2.1) = 0.8999.
#   According to strategy 0, switch to another drug which is drug 1,
#   but R1 is still the majority and more sensitive to drug 2, why switch? Just switch

#   (2) The total population reemerges from a level below the detection threshold (1e9; relapse)
#   In Strategy 0, only using single drug, no combination

#   Strategy 1: Minimize the predicted total population. Every 45 days, adjust d to minimize the predicted total
#   population. Vary d1 and d2 between 0 and 1 with a 0.01 interval.

#   Strategy 2: Minimize the risk of incurable cells developing unless there is an immediate threat of mortality.
#   Every 45 days, adjust d to minimize the predicted R12 population if the total population dose not exceed a threshold
#   R12 is resistant to both drugs, therefore, it is often incurable.
#   But based on excluding parameter combinations (4), either Sa(R12,1) > g0 or Sa(R12,2) > g0, it's hard to say R12 is
#   incurable.
#   Vary d1 and d2 between 0 and 1 with a 0.01 interval.
#   Strategy 2.1: threshold is 1e9
#   Strategy 2.2: threshold is 1e11

#   Strategy 3: Minimize the predicted total population unless there is a prediction that the first incurable cell
#   will form within 45 days. Every 45 days, adjust d to minimize the predicted total population if the predicted
#   R12 < 1. Otherwise, adjust d to minimize the predicted R12 population.
#   However, if the current R12 >= 1 and R12 is not curable, minimize the predicted total population. That is, if R12
#   has already appeared, we no longer focus on preventing its appearance. Given that we allow for "relative"
#   resistance, it is possible that R12 is not incurable; however, in most of our parameter settings, it is incurable.


#   Strategy 4: Estimate the time to either incurability or death, and react to the most proximal threat as long as
#   there is a chance of cure. Every 45 days, evaluate the predicted durations toward incurability (R12 > 1) and
#   mortality (population >= 1e13) dictated by the growth of S, R1, R2 and R12 populations.
#   For each dosage combination d, define τ_inc(d) as the predicted time to incurability (R12 > 1), given the currently
#   observed population and d fixed.
#   Define τ_s(d) as the predicted time to S causing mortality (S > 1e13), given the currently observed population and d
#   fixed.
#   τ_R1(d), R1 causing mortality (R1 > 1e13)
#   τ_R2(d), R2 causing mortality (R2 > 1e13)
#   τ_R12(d), R12 causing mortality (R12 > 1e13)
#   If the current R12 < 1 or R12 is curable, i.e., there exists some d such that each each component of diag(Sa*d) > g0
#   vary d to maximize min(τ_inc,τ_s,τ_R1,τ_R2,τ_R12) with the constraint that min(τ_s,τ_R1,τ_R2,τ_R12) > 45 days
#   If such a dosage combination does not exist, maximize min(τ_S,τ_R1,τ_R2,τ_R12).
#   If the current R12 >= 1, and R12 is not curable, maximize min(τ_S,τ_R1,τ_R2,τ_R12)

# exploration parameter combinations
# e.g. of parameters
num_drug: int = 2
num_cell_type: int = 4
Limit_mortality: float = 1e13
Limit_radiologic_detection: float = 1e9
Limit_molecular_detection: float = 1e5
X0total: float = 1e9
Simduration: float = 1800
Simtimestep: float = 1
doseinterval: float = 0.1
Windowduration: float = 45

# parameter including criterias:
# 1: PNAS(1), two drug situation
# 2: PNAS(2), two drug situation
# 3: PNAS(3), two drug situation
# 4: Biology Direct(4), two drug situation
# 5: Biology Direct(5), two drug situation
# 6: Biology Direct(6), two drug situation
# 7: Biology Direct(7), three drug situation
# 8: Biology Direct(8), three drug situation

PARformat = namedtuple('PARformat', 'Num_drug Num_cell_type X0total ratioStoX0 ratioR1toX0 ratioR2toX0 ratioR3toX0 ratioR12toX0 ratioR23toX0 '
                                    'ratioR13toX0 g0 Sa_ratio_S_D1tog0 Sa_ratio_S_D2toS_D1 Sa_ratio_S_D3toS_D1 Sa_ratio_R1_D1toS_D1 '
                                    'Sa_ratio_R2_D2toS_D2 Sa_ratio_R3_D3toS_D3 T_StoR1 T_StoR2 T_StoR3')

PAR_i = PARformat(Num_drug=2, Num_cell_type=4, X0total=X0total, ratioStoX0=None, ratioR1toX0=1e-5, ratioR2toX0=0, ratioR3toX0=None, ratioR12toX0=None,
                  ratioR23toX0=None, ratioR13toX0=None, g0=0.05, Sa_ratio_S_D1tog0=3, Sa_ratio_S_D2toS_D1=0.3, Sa_ratio_S_D3toS_D1=None,
                  Sa_ratio_R1_D1toS_D1=0.15, Sa_ratio_R2_D2toS_D2=0.5, Sa_ratio_R3_D3toS_D3=None, T_StoR1=4e-9, T_StoR2=4e-7, T_StoR3=None)

Sa_i, X0_i, T_i = None, None, None
if PAR_i.Num_drug == 2:
    Sa_i = DPM_fun.gen_Sa_2drug(PAR_i.Num_cell_type, PAR_i.Num_drug, PAR_i.g0, PAR_i.Sa_ratio_S_D1tog0, PAR_i.Sa_ratio_S_D2toS_D1,
                                PAR_i.Sa_ratio_R1_D1toS_D1, PAR_i.Sa_ratio_R2_D2toS_D2)
    T_i = DPM_fun.gen_T_2drug(PAR_i.Num_cell_type, PAR_i.T_StoR1, PAR_i.T_StoR2)
    X0_i = DPM_fun.gen_X0_2drug(PAR_i.X0total, PAR_i.ratioR1toX0, PAR_i.ratioR2toX0)
elif PAR_i.Num_drug == 3:
    Sa_i = DPM_fun.gen_Sa_3drug(PAR_i.Num_cell_type, PAR_i.Num_drug, PAR_i.g0, PAR_i.Sa_ratio_S_D1tog0, PAR_i.Sa_ratio_S_D2toS_D1,
                                PAR_i.Sa_ratio_S_D3toS_D1, PAR_i.Sa_ratio_R1_D1toS_D1, PAR_i.Sa_ratio_R2_D2toS_D2, PAR_i.Sa_ratio_R3_D3toS_D3)
    T_i = DPM_fun.gen_T_3drug(PAR_i.Num_cell_type, PAR_i.T_StoR1, PAR_i.T_StoR2, PAR_i.T_StoR3)
    X0_i = DPM_fun.gen_X0_3drug(PAR_i.X0total, PAR_i.ratioStoX0, PAR_i.ratioR1toX0, PAR_i.ratioR2toX0, PAR_i.ratioR3toX0, PAR_i.ratioR12toX0,
                                PAR_i.ratioR23toX0, PAR_i.ratioR13toX0)

LSsim = True
g0_i = PAR_i.g0 * np.ones(PAR_i.Num_cell_type)

t_Strategy0_i, X_Strategy0_i, d_Strategy0_i = \
    DPM_fun.Strategy0_2drug(X0_i, PAR_i, T_i, g0_i, Sa_i, Simduration, Windowduration, Simtimestep, Limit_mortality,
                            Limit_radiologic_detection, LSsim)

t_Strategy1_i, X_Strategy1_i, d_Strategy1_i = \
    DPM_fun.Strategy1_2drug(X0_i, PAR_i, T_i, g0_i, Sa_i, Simduration, Windowduration, Simtimestep, doseinterval, Limit_mortality, LSsim)

Strategy2threshold: float = 1e9
t_Strategy2_1_i, X_Strategy2_1_i, d_Strategy2_1_i = DPM_fun.Strategy2_2drug(X0_i, PAR_i, T_i, g0_i, Sa_i, Simduration, Windowduration,
                                                                            Simtimestep, doseinterval, Limit_mortality, LSsim, Strategy2threshold)

Strategy2threshold: float = 1e11
t_Strategy2_2_i, X_Strategy2_2_i, d_Strategy2_2_i = DPM_fun.Strategy2_2drug(X0_i, PAR_i, T_i, g0_i, Sa_i, Simduration, Windowduration,
                                                                            Simtimestep, doseinterval, Limit_mortality, LSsim, Strategy2threshold)

t_Strategy3_i, X_Strategy3_i, d_Strategy3_i = \
    DPM_fun.Strategy3_2drug(X0_i, PAR_i, T_i, g0_i, Sa_i, Simduration, Windowduration, Simtimestep, doseinterval, Limit_mortality, LSsim)

t_Strategy4_i, X_Strategy4_i, d_Strategy4_i = \
    DPM_fun.Strategy4_2drug(X0_i, PAR_i, T_i, g0_i, Sa_i, Simduration, Windowduration, Simtimestep, doseinterval, Limit_mortality, LSsim)


Xtotal_Strategy0_i = 1.5 * np.sum(X_Strategy0_i, axis=0)
X_Strategy0_i = np.vstack((Xtotal_Strategy0_i, X_Strategy0_i))
title_Strategy = 'Strategy 0'
miny = .1  # cell number smaller than this
yminval = 1e-2
xminval = -1
legend_order = [0, 4, 1, 5, 2, 6, 3]
legend_bbox_to_anchor = (0.5, 1.15)
color_cell = ('b', 'g', 'c', 'm', 'r')
color_drug = ('g', 'b')
label_cell = ('Total cells', 'S cells', 'R1 cells', 'R2 cells', 'R12 cells')
label_drug = ('Drug 1', 'Drug 2')
DPM_fun.add_clipboard_to_figures()
DPM_fun.plot_2resis_2drug(d_Strategy0_i, t_Strategy0_i, X_Strategy0_i, Limit_mortality, miny, yminval, xminval, legend_order, legend_bbox_to_anchor,
                          Limit_radiologic_detection, Limit_molecular_detection, color_cell, color_drug, label_cell, label_drug, title_Strategy)

Xtotal_Strategy1_i = 1.5 * np.sum(X_Strategy1_i, axis=0)
X_Strategy1_i = np.vstack((Xtotal_Strategy1_i, X_Strategy1_i))
title_Strategy = 'Strategy 1'
miny = .1  # cell number smaller than this
yminval = 1e-2
xminval = -1
legend_order = [0, 4, 1, 5, 2, 6, 3]
legend_bbox_to_anchor = (0.5, 1.15)
color_cell = ('b', 'g', 'c', 'm', 'r')
color_drug = ('g', 'b')
label_cell = ('Total cells', 'S cells', 'R1 cells', 'R2 cells', 'R12 cells')
label_drug = ('Drug 1', 'Drug 2')
DPM_fun.add_clipboard_to_figures()
DPM_fun.plot_2resis_2drug(d_Strategy1_i, t_Strategy1_i, X_Strategy1_i, Limit_mortality, miny, yminval, xminval, legend_order, legend_bbox_to_anchor,
                          Limit_radiologic_detection, Limit_molecular_detection, color_cell, color_drug, label_cell, label_drug, title_Strategy)

Xtotal_Strategy2_1_i = 1.5 * np.sum(X_Strategy2_1_i, axis=0)
X_Strategy2_1_i = np.vstack((Xtotal_Strategy2_1_i, X_Strategy2_1_i))
title_Strategy = 'Strategy 2.1'
miny = .1  # cell number smaller than this
yminval = 1e-2
xminval = -1
legend_order = [0, 4, 1, 5, 2, 6, 3]
legend_bbox_to_anchor = (0.5, 1.15)
color_cell = ('b', 'g', 'c', 'm', 'r')
color_drug = ('g', 'b')
label_cell = ('Total cells', 'S cells', 'R1 cells', 'R2 cells', 'R12 cells')
label_drug = ('Drug 1', 'Drug 2')
DPM_fun.add_clipboard_to_figures()
DPM_fun.plot_2resis_2drug(d_Strategy2_1_i, t_Strategy2_1_i, X_Strategy2_1_i, Limit_mortality, miny, yminval, xminval, legend_order,
                          legend_bbox_to_anchor, Limit_radiologic_detection, Limit_molecular_detection, color_cell, color_drug,
                          label_cell, label_drug, title_Strategy)
#
Xtotal_Strategy2_2_i = 1.5 * np.sum(X_Strategy2_2_i, axis=0)
X_Strategy2_2_i = np.vstack((Xtotal_Strategy2_2_i, X_Strategy2_2_i))
title_Strategy = 'Strategy 2.2'
miny = .1  # cell number smaller than this
yminval = 1e-2
xminval = -1
legend_order = [0, 4, 1, 5, 2, 6, 3]
legend_bbox_to_anchor = (0.5, 1.15)
color_cell = ('b', 'g', 'c', 'm', 'r')
color_drug = ('g', 'b')
label_cell = ('Total cells', 'S cells', 'R1 cells', 'R2 cells', 'R12 cells')
label_drug = ('Drug 1', 'Drug 2')
DPM_fun.add_clipboard_to_figures()
DPM_fun.plot_2resis_2drug(d_Strategy2_2_i, t_Strategy2_2_i, X_Strategy2_2_i, Limit_mortality, miny, yminval, xminval, legend_order,
                          legend_bbox_to_anchor, Limit_radiologic_detection, Limit_molecular_detection, color_cell, color_drug,
                          label_cell, label_drug, title_Strategy)
#
Xtotal_Strategy3_i = 1.5 * np.sum(X_Strategy3_i, axis=0)
X_Strategy3_i = np.vstack((Xtotal_Strategy3_i, X_Strategy3_i))
title_Strategy = 'Strategy 3'
miny = .1  # cell number smaller than this
yminval = 1e-2
xminval = -1
legend_order = [0, 4, 1, 5, 2, 6, 3]
legend_bbox_to_anchor = (0.5, 1.15)
color_cell = ('b', 'g', 'c', 'm', 'r')
color_drug = ('g', 'b')
label_cell = ('Total cells', 'S cells', 'R1 cells', 'R2 cells', 'R12 cells')
label_drug = ('Drug 1', 'Drug 2')
DPM_fun.add_clipboard_to_figures()
DPM_fun.plot_2resis_2drug(d_Strategy3_i, t_Strategy3_i, X_Strategy3_i, Limit_mortality, miny, yminval, xminval, legend_order, legend_bbox_to_anchor,
                          Limit_radiologic_detection, Limit_molecular_detection, color_cell, color_drug, label_cell, label_drug, title_Strategy)

Xtotal_Strategy4_i = 1.5 * np.sum(X_Strategy4_i, axis=0)
X_Strategy4_i = np.vstack((Xtotal_Strategy4_i, X_Strategy4_i))
title_Strategy = 'Strategy 4'
miny = .1  # cell number smaller than this
yminval = 1e-2
xminval = -1
legend_order = [0, 4, 1, 5, 2, 6, 3]
legend_bbox_to_anchor = (0.5, 1.15)
color_cell = ('b', 'g', 'c', 'm', 'r')
color_drug = ('g', 'b')
label_cell = ('Total cells', 'S cells', 'R1 cells', 'R2 cells', 'R12 cells')
label_drug = ('Drug 1', 'Drug 2')
DPM_fun.add_clipboard_to_figures()
DPM_fun.plot_2resis_2drug(d_Strategy4_i, t_Strategy4_i, X_Strategy4_i, Limit_mortality, miny, yminval, xminval, legend_order, legend_bbox_to_anchor,
                          Limit_radiologic_detection, Limit_molecular_detection, color_cell, color_drug, label_cell, label_drug, title_Strategy)


def DPM_sim_ODE_model2012(t, x, d, argsODE):
    # d_use = DPM_sim_d_use(d, t, pos_changedose)

    d_t = np.arange(0, d.shape[1], 1)

    pos, = np.where(t <= d_t)
    # print(t)
    σ1, σ2 = 0, 0
    if len(pos) != 0:
        σ1 = d[0, pos[0]]
        σ2 = d[1, pos[0]]
    elif len(pos) == 0:
        σ1 = d[0, -1]
        σ2 = d[1, -1]

    θ = argsODE['θ']
    α = argsODE['α']
    κ = argsODE['κ']
    Eplus = argsODE['Eplus']
    Eminus = argsODE['Eminus']
    g0 = argsODE['g0']
    g = argsODE['g']
    d = argsODE['d']
    T_μ1 = argsODE['T_μ1']
    T_μ2 = argsODE['T_μ2']
    T_π1 = argsODE['T_π1']
    T_π2 = argsODE['T_π2']
    kr_π1 = argsODE['kr_π1']
    kr_π2 = argsODE['kr_π2']
    k_t = argsODE['k_t']

    '''4 states, no mutation happens'''
    n00_1_00 = x[0]
    n0000_1_ = x[1]
    n0100 = x[2]
    n0001 = x[3]

    '''transient mutated states'''
    n1000 = x[4]
    n0010 = x[5]
    n1001 = x[6]
    n0110 = x[7]
    n1011 = x[8]
    n1110 = x[9]

    '''non-transient mutated states'''
    n1100 = x[10]
    n0011 = x[11]
    n1101 = x[12]
    n0111 = x[13]
    n1111 = x[14]

    '''4 states, no mutation happens'''
    dn00_1_00 = - n00_1_00 * g0 * T_μ1 \
                - n00_1_00 * g0 * T_μ2 \
                + n0100 * k_t * np.exp(-κ * Eplus(σ1, θ)) \
                - n00_1_00 * k_t * np.exp(-κ * Eminus(σ1, θ, α)) \
                - n00_1_00 * kr_π2 \
                + g0 * n00_1_00 \
                - d * n00_1_00

    dn0000_1_ = - n0000_1_ * g0 * T_μ1 \
                - n0000_1_ * g0 * T_μ2 \
                + n0001 * k_t * np.exp(-κ * Eplus(σ2, θ)) \
                - n0000_1_ * k_t * np.exp(-κ * Eminus(σ2, θ, α)) \
                - n0000_1_ * kr_π1 \
                + g0 * n0000_1_ \
                - d * n0000_1_

    dn0100 = - n0100 * g * T_μ1 \
             - n0100 * g * T_μ2 \
             + n00_1_00 * k_t * np.exp(-κ * Eminus(σ1, θ, α)) \
             - n0100 * k_t * np.exp(-κ * Eplus(σ1, θ)) \
             + n0000_1_ * kr_π1 \
             + g * n0100 \
             - d * n0100

    dn0001 = - n0001 * g * T_μ1 \
             - n0001 * g * T_μ2 \
             + n0000_1_ * k_t * np.exp(-κ * Eminus(σ2, θ, α)) \
             - n0001 * k_t * np.exp(-κ * Eplus(σ2, θ)) \
             + n00_1_00 * kr_π2 \
             + g * n0001 \
             - d * n0001

    '''transient mutated states'''
    dn0010 = - n0010 * T_π2 \
             + n00_1_00 * g0 * T_μ2 \
             + n0000_1_ * g0 * T_μ2

    dn1000 = - n1000 * T_π1 \
             + n00_1_00 * g0 * T_μ1 \
             + n0000_1_ * g0 * T_μ1

    dn0110 = - n0110 * T_π2 \
             + n0100 * g * T_μ2

    dn1001 = - n1001 * T_π1 \
             + n0001 * g * T_μ1

    dn1110 = - n1110 * T_π2 \
             + n1100 * g * T_μ2

    dn1011 = - n1011 * T_π1 \
             + n0011 * g * T_μ1

    '''non-transient mutated states'''
    dn1100 = - n1100 * g * T_μ2 \
             + n0100 * g * T_μ1 \
             + n1000 * T_π1 \
             + n1101 * k_t * np.exp(-κ * Eplus(σ2, θ)) \
             - n1100 * k_t * np.exp(-κ * Eminus(σ2, θ, α)) \
             + g * n1100 \
             - d * n1100

    dn1101 = - n1101 * g * T_μ2 \
             + n1001 * T_π1 \
             + n1100 * k_t * np.exp(-κ * Eminus(σ2, θ, α)) \
             - n1101 * k_t * np.exp(-κ * Eplus(σ2, θ)) \
             + g * n1101 \
             - d * n1101

    dn0011 = - n0011 * g * T_μ1 \
             + n0001 * g * T_μ2 \
             + n0010 * T_π2 \
             + n0111 * k_t * np.exp(-κ * Eplus(σ1, θ)) \
             - n0011 * k_t * np.exp(-κ * Eminus(σ1, θ, α)) \
             + g * n0011 \
             - d * n0011

    dn0111 = - n0111 * g * T_μ1 \
             + n0110 * T_π2 \
             + n0011 * k_t * np.exp(-κ * Eminus(σ1, θ, α)) \
             - n0111 * k_t * np.exp(-κ * Eplus(σ1, θ)) \
             + g * n0111 \
             - d * n0111

    dn1111 = n0111 * g * T_μ1 \
             + n1101 * g * T_μ2 \
             + n1011 * T_π1 \
             + n1110 * T_π2 \
             + g * n1111 \
             - d * n1111

    return list([dn00_1_00, dn0000_1_, dn0100, dn0001, dn1000, dn0010, dn1001, dn0110, dn1011, dn1110,
                 dn1100, dn0011, dn1101, dn0111, dn1111])

