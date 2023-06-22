# ======================================================================================================================
# Script Name: Optimisation.py
# File Name: Three_Phase_Mathematical_Model
# Date: 21-06-2021
# Creator: K.N. Lambrecht
# Purpose: This script is used for solution of the system of ODEs and the optimisation of parameters
# ======================================================================================================================
from Main import *
from Functions import *
from Optimisation_Variables import *
import time
from matplotlib import pyplot as plt

optimise_index_1 = np.arange(0, len(exp_vector1))

a1 = least_squares(objective, gs_bio, bounds=bounds, jac='3-point', x_scale='jac', loss='soft_l1', tr_solver='lsmr',
                   args=(p, x1, ts, time_array, regression_list, exp_vector1, stdev_vector1, optimise_index_1),
                   verbose=2)


new_p = vector_to_namespace(a1.x, regression_list)
dict1, dict2 = new_p.__dict__, fixed_values.__dict__
dict3 = {**dict1, **dict2}
p.param_all = types.SimpleNamespace(**dict3)

sol2 = speed_up(p, x1)
solution2 = vector_to_namespace(sol2, p.svar_list)
v2 = simple_calc(ts, solution2, p, time_array)
