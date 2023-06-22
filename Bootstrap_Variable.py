# ======================================================================================================================
# Script Name: Bootstrap_variables.py
# File Name: Three_Phase_Mathematical_Model
# Date: 21-06-2021
# Creator: K.N. Lambrecht
# Purpose: This script is used to declare number of itterations and generate randomised data
# ======================================================================================================================

from Main import *
from Functions import *
from Optimisation_Variables import *
import numpy as np

num_iterations = 5
bootstrap_indices = np.zeros((num_iterations, len(exp_vector1)), dtype=int)

for i in range(0, num_iterations):
    new_exp, new_sd = [], []
    bstp_index = bill(interest_list, time_array)

    bootstrap_indices[i, :] = bstp_index
