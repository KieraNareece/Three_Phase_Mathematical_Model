# ======================================================================================================================
# Script Name: Bootstrap.py
# File Name: Three_Phase_Mathematical_Model
# Date: 21-06-2021
# Creator: K.N. Lambrecht
# Purpose: This script is used to perform the likelihood for the parameter sets
# ======================================================================================================================

import time
from mpire import WorkerPool
from matplotlib import pyplot as plt
import numpy as np
from Functions import *
from Main import *
from Bootstrap_Variable import bootstrap_indices, num_iterations


if __name__ == '__main__':

    start = time.time()
    with WorkerPool(n_jobs=8) as pool:
        results = pool.map(regression, bootstrap_indices, iterable_len=num_iterations, chunk_size=1)

    stop = time.time()
    results_array = np.array(results)

    print(len(results_array), results_array, (stop-start))
    np.savetxt('D:\\Stellenbosch\\2. PhD\\2022\\0. Experimental\\Data Sets\\Arraysbootstrap_CS20_5itt.txt', results_array)
