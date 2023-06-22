# ======================================================================================================================
# Script Name: Main.py
# File Name: Three_Phase_Mathematical_Model
# Date: 21-06-2021
# Creator: K.N. Lambrecht
# Purpose: This script contains the main parameters, and performs basic data handling
# ======================================================================================================================

import pandas as pd
from Functions import *
import types
import numpy as np

# Define time and parameters
days = 9
minutes_fermentation = days*24*60 + 60.      # Calculates number of minutes that the fermentation runs for
ts = np.arange(0, minutes_fermentation)      # Define time space
V_MS = 20                                    # L, volume of tank
V_GR = 0.25*V_MS                             # Approximate volume of grapes
V_CL = 0.1*V_MS                              # Approximate volume of grapes

# Implement Simple Name Space for the parameters of the system
# The units and the definitions are defined below

p = types.SimpleNamespace()     # Line included to time the optimisation

# This parameter space contains the variables which will not be changed by the regression
p.param = types.SimpleNamespace(K_S=9.84, K_S2=149.073, K_N=0.01, K_E=1, Y_LX=0.00001, Y_SX=0.04, Y_CX=0.29, Y_EX=0.28,
                                Y_GX=0.14, Y_SCX=0.11, B=0.0024, Y_ES=3.44, Y_GS=0.003,  Y_SCS=0.00001)

# The following list contain the names of the parameters for indexing
p.param_names = ['K_S', 'K_S2', 'K_N', 'K_E', 'Y_LX', 'Y_SX', 'Y_CX', 'Y_EX', 'Y_GX', 'Y_SCX', 'B', 'Y_ES', 'Y_GS',
                 'Y_SCS']

# This parameter space contains variables which will change with the regression

p.param_biomass = types.SimpleNamespace(k_pp=5.500e-06, PI_LSPP=0.4, PI_GRA=0.76, PI_LSA=5, PI_GRT=0.36,
                                        PI_LST=1, PI_GRTA=0.36, PI_LSTA=1, PI_GRMA=2.14, PI_LSMA=1, a_CD=0.03,
                                        b_CD=0.02, a_TPI1=0.03, b_TPI1=0.02, c_TPI1=0.02, a_TPI2=0.01, b_TPI2=0.02,
                                        kla_ml_pp=0.01, kla_mg_a_i=0.0004, kla_ml_a=0.08, kla_mg_t_i=0.002,
                                        kla_ml_t=0.004, kla_mg_ta_i=0.0001, kla_ml_ta=0.0008, kla_mg_ma_i=0.00007,
                                        kla_ml_ma=0.0004, ea=0.002, et=0.00084, eta=0.0001, ema=0.000002,
                                        kla_max_a=0.021, kla_max_t=0.001, kla_max_ta=0.011, kla_max_ma=0.011)


p.param_names_biomass = ['k_pp', 'PI_LSPP', 'PI_GRA', 'PI_LSA', 'PI_GRT',
                         'PI_LST', 'PI_GRTA', 'PI_LSTA', 'PI_GRMA', 'PI_LSMA', 'a_CD', 'b_CD', 'a_TPI1', 'b_TPI1',
                         'c_TPI1', 'a_TPI2', 'b_TPI2', 'kla_ml_pp', 'kla_mg_a_i', 'kla_ml_a', 'kla_mg_t_i', 'kla_ml_t',
                         'kla_mg_ta_i', 'kla_ml_ta', 'kla_mg_ma_i', 'kla_ml_ma', 'ea', 'et', 'eta', 'ema', 'kla_max_a',
                         'kla_max_t', 'kla_max_ta', 'kla_max_ma']

p.param_reg_biomass = types.SimpleNamespace(K_S=9.84, K_S2=149.073, K_N=0.01, K_E=1, Y_LX=0.00001, Y_SX=0.04,
                                            Y_CX=0.29, Y_EX=0.28, Y_GX=0.14, Y_SCX=0.11, B=0.0024, Y_ES=3.44,
                                            Y_GS=0.003,  Y_SCS=0.00001, k_d=9.66e-06, u_max=0.000797, Y_NX=0.3)


p.param_reg_names_biomass = ['K_S', 'K_S2', 'K_N', 'K_E', 'Y_LX', 'Y_SX', 'Y_CX', 'Y_EX', 'Y_GX', 'Y_SCX', 'B', 'Y_ES',
                             'Y_GS', 'Y_SCS', 'k_d', 'u_max', 'Y_NX']


p.param_CA_CT_CPP = types.SimpleNamespace(K_S=9.84, K_S2=149.073, K_N=0.01, K_E=1, Y_LX=0.00001, Y_SX=0.04, Y_CX=0.29,
                                          Y_EX=0.28, Y_GX=0.14, Y_SCX=0.11, B=0.0024, Y_ES=3.44, Y_GS=0.003,
                                          Y_SCS=0.00001, k_d=9.66e-06, u_max=0.00079, Y_NX=0.307,  PI_GRTA=0.36,
                                          PI_LSTA=1, PI_GRMA=2.14, PI_LSMA=1, a_TPI2=0.01, b_TPI2=0.02,
                                          kla_mg_ta_i=0.0001, kla_ml_ta=0.0008, kla_mg_ma_i=0.00007, kla_ml_ma=0.0004,
                                          eta=0.0001, ema=0.000002, kla_max_ta=0.011, kla_max_ma=0.011)

# The following list contain the names of the parameters for indexing
p.param_names_CA_CT_CPP = ['K_S', 'K_S2', 'K_N', 'K_E', 'Y_LX', 'Y_SX', 'Y_CX', 'Y_EX', 'Y_GX', 'Y_SCX', 'B', 'Y_ES',
                           'Y_GS', 'Y_SCS', 'k_d', 'u_max', 'Y_NX', 'PI_GRTA', 'PI_LSTA', 'PI_GRMA', 'PI_LSMA',
                           'a_TPI2', 'b_TPI2', 'kla_mg_ta_i', 'kla_ml_ta', 'kla_mg_ma_i', 'kla_ml_ma', 'eta', 'ema',
                           'kla_max_ta', 'kla_max_ma']

p.param_reg_CA_CT_CPP = types.SimpleNamespace(k_pp=5.500e-06, PI_LSPP=0.4, PI_GRA=0.76, PI_LSA=5,  a_CD=0.03, b_CD=0.02,
                                              a_TPI1=0.03, b_TPI1=0.02, kla_ml_pp=0.01, kla_mg_a_i=0.0004,
                                              kla_ml_a=0.08, ea=0.002,  kla_max_a=0.021, PI_GRT=0.36, PI_LST=1,
                                              c_TPI1=0.02, kla_mg_t_i=0.002, kla_ml_t=0.004, et=0.00084,
                                              kla_max_t=0.001)

p.param_reg_names_CA_CT_CPP = ['k_pp', 'PI_LSPP', 'PI_GRA', 'PI_LSA', 'a_CD', 'b_CD', 'a_TPI1', 'b_TPI1', 'kla_ml_pp',
                               'kla_mg_a_i', 'kla_ml_a', 'ea', 'kla_max_a', 'PI_GRT', 'PI_LST', 'c_TPI1', 'kla_mg_t_i',
                               'kla_ml_t', 'et',  'kla_max_t']

p.param_CTA = types.SimpleNamespace(K_S=9.84, K_S2=149.073, K_N=0.01, K_E=1, Y_LX=0.00001, Y_SX=0.04, Y_CX=0.29,
                                    Y_EX=0.28, Y_GX=0.14, Y_SCX=0.11, B=0.0024, Y_ES=3.44, Y_GS=0.003,  Y_SCS=0.00001,
                                    k_d=9.66e-06, u_max=0.00079, Y_NX=0.307, k_pp=5.500e-06, PI_LSPP=0.4, PI_GRA=0.76,
                                    PI_LSA=5, PI_GRT=0.36, PI_LST=1,  PI_GRMA=2.14, PI_LSMA=1, a_CD=0.03, b_CD=0.02,
                                    a_TPI1=0.03, b_TPI1=0.02, c_TPI1=0.02, a_TPI2=0.01, b_TPI2=0.02, kla_ml_pp=0.01,
                                    kla_mg_a_i=0.0004, kla_ml_a=0.08, kla_mg_t_i=0.002, kla_ml_t=0.004,
                                    kla_mg_ma_i=0.00007, kla_ml_ma=0.0004, ea=0.002, et=0.00084,  ema=0.000002,
                                    kla_max_a=0.021, kla_max_t=0.001, kla_max_ma=0.011)

# The following list contain the names of the parameters for indexing
p.param_names_CTA = ['K_S', 'K_S2', 'K_N', 'K_E', 'Y_LX', 'Y_SX', 'Y_CX', 'Y_EX', 'Y_GX', 'Y_SCX', 'B', 'Y_ES', 'Y_GS',
                     'Y_SCS', 'k_d', 'u_max', 'Y_NX', 'k_pp', 'PI_LSPP', 'PI_GRA', 'PI_LSA', 'PI_GRT', 'PI_LST',
                     'PI_GRMA', 'PI_LSMA', 'a_CD', 'b_CD', 'a_TPI1', 'b_TPI1', 'c_TPI1', 'a_TPI2', 'b_TPI2',
                     'kla_ml_pp', 'kla_mg_a_i', 'kla_ml_a', 'kla_mg_t_i', 'kla_ml_t',  'kla_mg_ma_i', 'kla_ml_ma',
                     'ea', 'et',  'ema', 'kla_max_a', 'kla_max_t',  'kla_max_ma']

p.param_reg_CTA = types.SimpleNamespace(PI_GRTA=0.36, PI_LSTA=1, kla_mg_ta_i=0.0001, kla_ml_ta=0.0008, eta=0.0001,
                                        kla_max_ta=0.011)

p.param_reg_names_CTA = ['PI_GRTA', 'PI_LSTA', 'kla_mg_ta_i', 'kla_ml_ta', 'et', 'kla_max_ta']

p.param_CMA = types.SimpleNamespace(K_S=9.84, K_S2=149.073, K_N=0.01, K_E=1, Y_LX=0.00001, Y_SX=0.04, Y_CX=0.29,
                                    Y_EX=0.28, Y_GX=0.14, Y_SCX=0.11, B=0.0024, Y_ES=3.44, Y_GS=0.003,  Y_SCS=0.00001,
                                    k_d=9.66e-06, u_max=0.00079, Y_NX=0.307, k_pp=5.500e-06, PI_LSPP=0.4, PI_GRA=0.76,
                                    PI_LSA=5, PI_GRT=0.36, PI_LST=1, PI_GRTA=0.36, PI_LSTA=1,  a_CD=0.03, b_CD=0.02,
                                    a_TPI1=0.03, b_TPI1=0.02, c_TPI1=0.02, a_TPI2=0.01, b_TPI2=0.02, kla_ml_pp=0.01,
                                    kla_mg_a_i=0.0004, kla_ml_a=0.08, kla_mg_t_i=0.002, kla_ml_t=0.004,
                                    kla_mg_ta_i=0.0001, kla_ml_ta=0.0008,  ea=0.002, et=0.00084, eta=0.0001,
                                    kla_max_a=0.021, kla_max_t=0.001, kla_max_ta=0.011)

# The following list contain the names of the parameters for indexing
p.param_names_CMA = ['K_S', 'K_S2', 'K_N', 'K_E', 'Y_LX', 'Y_SX', 'Y_CX', 'Y_EX', 'Y_GX', 'Y_SCX', 'B', 'Y_ES', 'Y_GS',
                     'Y_SCS', 'k_d', 'u_max', 'Y_NX', 'k_pp', 'PI_LSPP', 'PI_GRA', 'PI_LSA', 'PI_GRT', 'PI_LST',
                     'PI_GRTA', 'PI_LSTA', 'a_CD', 'b_CD', 'a_TPI1', 'b_TPI1', 'c_TPI1', 'a_TPI2', 'b_TPI2',
                     'kla_ml_pp', 'kla_mg_a_i', 'kla_ml_a', 'kla_mg_t_i', 'kla_ml_t', 'kla_mg_ta_i', 'kla_ml_ta', 'ea',
                     'et', 'eta', 'kla_max_a', 'kla_max_t', 'kla_max_ta']

p.param_CMA = types.SimpleNamespace(PI_GRMA=2.14, PI_LSMA=1, kla_mg_ma_i=0.00007, kla_ml_ma=0.0004, ema=0.000002,
                                    kla_max_ma=0.011)

p.param_CMA = ['PI_GRMA', 'PI_LSMA', 'kla_mg_ma_i', 'kla_ml_ma', 'ema',  'kla_max_ma']

# The lower and upper bounds for parameter selection are declared and transformed in vectors
lp = types.SimpleNamespace(k_d=8.5e-06, u_max=0.0001, Y_NX=0.1, k_pp=3.500e-06, PI_LSPP=0, PI_GRA=0, PI_LSA=0, PI_GRT=0,
                           PI_LST=0, PI_GRTA=0, PI_LSTA=0, PI_GRMA=0, PI_LSMA=0, a_CD=0.01, b_CD=0.01, a_TPI1=0.01,
                           b_TPI1=0.01, c_TPI1=0.01, a_TPI2=0.01, b_TPI2=0.01, kla_ml_pp=0.001, kla_mg_a_i=0.00041,
                           kla_ml_a=0.01, kla_mg_t_i=0.0001, kla_ml_t=0.001, kla_mg_ta_i=0.00001, kla_ml_ta=0.0001,
                           kla_mg_ma_i=0.00001, kla_ml_ma=0.0001, ea=0.0001, et=0.0001, eta=0.00001, ema=0.000001,
                           kla_max_a=0.001, kla_max_t=0.0001, kla_max_ma=0.001, kla_max_ta=0.001)

up = types.SimpleNamespace(k_d=1.5e-05, u_max=0.005, Y_NX=0.1, k_pp=1.000e-05, PI_LSPP=10, PI_GRA=10, PI_LSA=10,
                           PI_GRT=10, PI_LST=10, PI_GRTA=10, PI_LSTA=10, PI_GRMA=10, PI_LSMA=10, a_CD=0.05, b_CD=0.05,
                           a_TPI1=0.05, b_TPI1=0.05, c_TPI1=0.05, a_TPI2=0.02, b_TPI2=0.05, kla_ml_pp=0.1,
                           kla_mg_a_i=0.001, kla_ml_a=0.1, kla_mg_t_i=0.01, kla_ml_t=0.01, kla_mg_ta_i=0.001,
                           kla_ml_ta=0.001, kla_mg_ma_i=0.0001, kla_ml_ma=0.001, ea=0.01, et=0.01, eta=0.001,
                           ema=0.00001, kla_max_a=0.1, kla_max_t=0.01, kla_max_ta=0.1, kla_max_ma=0.1)

lp_biomass = types.SimpleNamespace(K_S=0, K_S2=0, K_N=0, K_E=0, Y_LX=0, Y_SX=0, Y_CX=0, Y_EX=0, Y_GX=0, Y_SCX=0, B=0,
                                   Y_ES=0, Y_GS=0,  Y_SCS=0, k_d=0, u_max=0, Y_NX=0)

up_biomass = types.SimpleNamespace(K_S=20, K_S2=200, K_N=10, K_E=10, Y_LX=5, Y_SX=5, Y_CX=5, Y_EX=5, Y_GX=5, Y_SCX=5,
                                   B=1, Y_ES=5, Y_GS=5,  Y_SCS=1, k_d=1, u_max=1, Y_NX=20)

lp_CA_CT_CPP = types.SimpleNamespace(k_pp=1.500e-06, PI_LSPP=0, PI_GRA=0, PI_LSA=0, a_CD=0.01, b_CD=0.01, a_TPI1=0.01,
                                     b_TPI1=0.01, kla_ml_pp=0.001, kla_mg_a_i=0.00041, kla_ml_a=0.01, ea=0.0001,
                                     kla_max_a=0.001, PI_GRT=0, PI_LST=0, c_TPI1=0.01, kla_mg_t_i=0.0001,
                                     kla_ml_t=0.001, et=0.0001,  kla_max_t=0.0001)

up_CA_CT_CPP = types.SimpleNamespace(k_pp=1.000e-03, PI_LSPP=10, PI_GRA=10, PI_LSA=10, a_CD=1, b_CD=1,
                                     a_TPI1=1, b_TPI1=1, kla_ml_pp=0.1, kla_mg_a_i=0.1, kla_ml_a=0.1, ea=1,
                                     kla_max_a=0.1, PI_GRT=10, PI_LST=10, c_TPI1=1, kla_mg_t_i=0.1, kla_ml_t=0.1,
                                     et=1, kla_max_t=0.1)

lp_CTA = types.SimpleNamespace(PI_GRTA=0, PI_LSTA=0, kla_mg_ta_i=0.00001, kla_ml_ta=0.0001, eta=0.00001,
                               kla_max_ta=0.001)

up_CTA = types.SimpleNamespace(PI_GRTA=10, PI_LSTA=10,  kla_mg_ta_i=0.1, kla_ml_ta=0.1, eta=1, kla_max_ta=0.1)

lp_CMA = types.SimpleNamespace(PI_GRMA=0, PI_LSMA=0, kla_mg_ma_i=0.00001, kla_ml_ma=0.0001, ema=0.000001,
                               kla_max_ma=0.001)

up_CMA = types.SimpleNamespace(PI_GRMA=10, PI_LSMA=10, kla_mg_ma_i=0.1, kla_ml_ma=0.1, ema=0.1, kla_max_ma=0.1)


# The following lists are used to index the variables of interest when the model is solved
p.svar_list = ['X', 'N', 'S', 'CO2', 'E', 'GL', 'SC', 'M_LS', 'CPP', 'MPP_LS', 'CA', 'CA_GR', 'MA_LS', 'CT', 'CT_GR',
               'MT_LS', 'CTA', 'CTA_GR', 'MTA_LS', 'CMA', 'CMA_GR', 'MMA_LS']

# The following lists are used to sort the data frames and build dictionaries surrounding the experimental data
p.scenario = ['All', 'Temperature', 'Specific']
p.cultivars = ['CS', 'SH', 'ME']
p.temperature = ['20', '25', '28']

p.cap_components = ['Anthocyanins', 'Tannins', 'Total_Phenolic_Index', 'Tartaric_Acid', 'Malic_Acid']
p.lees_components = ['Anthocyanins', 'Polymeric_Pigments', 'Tannins', 'Tartaric_Acid', 'Malic_Acid']
p.must_components = ['Anthocyanins', 'Polymeric_Pigments', 'Tannins', 'Total_Phenolic_Index', 'Tartaric_Acid',
                     'Malic_Acid', 'Colour_Density', 'Succinic_Acid', 'Biomass', 'Sugar', 'Ethanol', 'pH', 'CO2',
                     'Nitrogen']

p.repeat = ['Repeat_1', 'Repeat_2', 'Repeat_3']
p.phase_name = ['Must', 'Cap', 'Lees', 'Grape']
p.interest_full = ['X', 'N', 'S', 'CO2', 'E', 'SC', 'CPP', 'CA', 'CT', 'CTA', 'CMA', 'cd', 'tpi', 'pH_must']
p.exp_components_full = ['Biomass', 'Nitrogen', 'Sugar', 'CO2', 'Ethanol', 'Succinic_Acid', 'Polymeric_Pigments',
                         'Anthocyanins', 'Tannins', 'Tartaric_Acid', 'Malic_Acid', 'Colour_Density',
                         'Total_Phenolic_Index', 'pH']

p.interest_biomass = ['X', 'N', 'S', 'CO2', 'E', 'SC']
p.exp_components_biomass = ['Biomass', 'Nitrogen', 'Sugar', 'CO2', 'Ethanol', 'Succinic_Acid']

p.interest_CA_CT_CPP = ['CPP', 'CA', 'CT', 'cd', 'tpi']
p.exp_components_CA_CT_CPP = ['Polymeric_Pigments', 'Anthocyanins', 'Tannins']

p.interest_CTA = ['CTA']
p.exp_components_CTA = ['Tartaric_Acid']

p.interest_CMA = ['CMA']
p.exp_components_CMA = ['Tartaric_Acid']

p.performance_list = ['X', 'N', 'S', 'CO2', 'E', 'SC', 'CPP', 'MPP_LS', 'CA', 'CA_GR', 'MA_LS', 'CT', 'CT_GR', 'MT_LS',
                      'CTA', 'CTA_GR', 'MTA_LS', 'CMA', 'CMA_GR', 'MMA_LS', 'cd', 'tpi', 'pH_must']

p.performance_list_case_1 = ['Must_CS_20_Biomass', 'Must_CS_20_Nitrogen', 'Must_CS_20_Sugar', 'Must_CS_20_CO2',
                             'Must_CS_20_Ethanol', 'Must_CS_20_Succinic_Acid', 'Must_CS_20_Polymeric_Pigments',
                             'Lees_CS_20_Polymeric_Pigments', 'Must_CS_20_Anthocyanins', 'Cap_CS_20_Anthocyanins',
                             'Lees_CS_20_Anthocyanins', 'Must_CS_20_Tannins', 'Cap_CS_20_Tannins', 'Lees_CS_20_Tannins',
                             'Must_CS_20_Tartaric_Acid', 'Cap_CS_20_Tartaric_Acid', 'Lees_CS_20_Tartaric_Acid',
                             'Must_CS_20_Malic_Acid', 'Cap_CS_20_Malic_Acid', 'Lees_CS_20_Malic_Acid',
                             'Must_CS_20_Colour_Density', 'Must_CS_20_Total_Phenolic_Index', 'Must_CS_20_pH']

p.performance_list_case_2 = ['Must_CS_25_Biomass', 'Must_CS_25_Nitrogen', 'Must_CS_25_Sugar', 'Must_CS_25_CO2',
                             'Must_CS_25_Ethanol', 'Must_CS_25_Succinic_Acid', 'Must_CS_25_Polymeric_Pigments',
                             'Lees_CS_25_Polymeric_Pigments', 'Must_CS_25_Anthocyanins', 'Cap_CS_25_Anthocyanins',
                             'Lees_CS_25_Anthocaynins', 'Must_CS_25_Tannins', 'Cap_CS_25_Tannins', 'Lees_CS_25_Tannins',
                             'Must_CS_25_Tartaric_Acid', 'Cap_CS_25_Tartaric_Acid', 'Lees_CS_25_Tartaric_Acid',
                             'Must_CS_25_Malic_Acid', 'Cap_CS_25_Malic_Acid', 'Lees_CS_25_Malicic_Acid',
                             'Must_CS_25_Colour_Density', 'Must_CS_25_Total_Phenolic_Index', 'Must_CS_25_pH']

p.performance_list_case_3 = ['Must_CS_28_Biomass', 'Must_CS_28_Nitrogen', 'Must_CS_28_Sugar', 'Must_CS_28_CO2',
                             'Must_CS_28_Ethanol', 'Must_CS_28_Succinic_Acid', 'Must_CS_28_Polymeric_Pigments',
                             'Lees_CS_28_Polymeric_Pigments', 'Must_CS_28_Anthocyanins', 'Cap_CS_28_Anthocyanins',
                             'Lees_CS_28_Anthocaynins', 'Must_CS_28_Tannins', 'Cap_CS_28_Tannins', 'Lees_CS_28_Tannins',
                             'Must_CS_28_Tartaric_Acid', 'Cap_CS_28_Tartaric_Acid', 'Lees_CS_28_Tartaric_Acid',
                             'Must_CS_28_Malic_Acid', 'Cap_CS_28_Malic_Acid', 'Lees_CS_28_Malicic_Acid',
                             'Must_CS_28_Colour_Density', 'Must_CS_28_Total_Phenolic_Index', 'Must_CS_28_pH']

p.performance_list_case_4 = ['Must_SH_20_Biomass', 'Must_SH_20_Nitrogen', 'Must_SH_20_Sugar', 'Must_SH_20_CO2',
                             'Must_SH_20_Ethanol', 'Must_SH_20_Succinic_Acid', 'Must_SH_20_Polymeric_Pigments',
                             'Lees_SH_20_Polymeric_Pigments', 'Must_SH_20_Anthocyanins', 'Cap_SH_20_Anthocyanins',
                             'Lees_SH_20_Anthocaynins', 'Must_SH_20_Tannins', 'Cap_SH_20_Tannins', 'Lees_SH_20_Tannins',
                             'Must_SH_20_Tartaric_Acid', 'Cap_SH_20_Tartaric_Acid', 'Lees_SH_20_Tartaric_Acid',
                             'Must_SH_20_Malic_Acid', 'Cap_SH_20_Malic_Acid', 'Lees_SH_20_Malicic_Acid',
                             'Must_SH_20_Colour_Density', 'Must_SH_20_Total_Phenolic_Index', 'Must_SH_20_pH']

p.performance_list_case_5 = ['Must_SH_25_Biomass', 'Must_SH_25_Nitrogen', 'Must_SH_25_Sugar', 'Must_SH_25_CO2',
                             'Must_SH_25_Ethanol', 'Must_SH_25_Succinic_Acid', 'Must_SH_25_Polymeric_Pigments',
                             'Lees_SH_25_Polymeric_Pigments', 'Must_SH_25_Anthocyanins', 'Cap_SH_25_Anthocyanins',
                             'Lees_SH_25_Anthocaynins', 'Must_SH_25_Tannins', 'Cap_SH_25_Tannins', 'Lees_SH_25_Tannins',
                             'Must_SH_25_Tartaric_Acid', 'Cap_SH_25_Tartaric_Acid', 'Lees_SH_25_Tartaric_Acid',
                             'Must_SH_25_Malic_Acid', 'Cap_SH_25_Malic_Acid', 'Lees_SH_25_Malicic_Acid',
                             'Must_SH_25_Colour_Density', 'Must_SH_25_Total_Phenolic_Index', 'Must_SH_25_pH']

p.performance_list_case_6 = ['Must_SH_28_Biomass', 'Must_SH_28_Nitrogen', 'Must_SH_28_Sugar', 'Must_SH_28_CO2',
                             'Must_SH_28_Ethanol', 'Must_SH_28_Succinic_Acid', 'Must_SH_28_Polymeric_Pigments',
                             'Lees_SH_28_Polymeric_Pigments', 'Must_SH_28_Anthocyanins', 'Cap_SH_28_Anthocyanins',
                             'Lees_SH_28_Anthocaynins', 'Must_SH_28_Tannins', 'Cap_SH_28_Tannins', 'Lees_SH_28_Tannins',
                             'Must_SH_28_Tartaric_Acid', 'Cap_SH_28_Tartaric_Acid', 'Lees_SH_28_Tartaric_Acid',
                             'Must_SH_28_Malic_Acid', 'Cap_SH_28_Malic_Acid', 'Lees_SH_28_Malicic_Acid',
                             'Must_SH_28_Colour_Density', 'Must_SH_28_Total_Phenolic_Index', 'Must_SH_28_pH']

p.performance_list_case_7 = ['Must_ME_20_Biomass', 'Must_ME_20_Nitrogen', 'Must_ME_20_Sugar', 'Must_ME_20_CO2',
                             'Must_ME_20_Ethanol', 'Must_ME_20_Succinic_Acid', 'Must_ME_20_Polymeric_Pigments',
                             'Lees_ME_20_Polymeric_Pigments', 'Must_ME_20_Anthocyanins', 'Cap_ME_20_Anthocyanins',
                             'Lees_ME_20_Anthocaynins', 'Must_ME_20_Tannins', 'Cap_ME_20_Tannins', 'Lees_ME_20_Tannins',
                             'Must_ME_20_Tartaric_Acid', 'Cap_ME_20_Tartaric_Acid', 'Lees_ME_20_Tartaric_Acid',
                             'Must_ME_20_Malic_Acid', 'Cap_ME_20_Malic_Acid', 'Lees_ME_20_Malicic_Acid',
                             'Must_ME_20_Colour_Density', 'Must_ME_20_Total_Phenolic_Index', 'Must_ME_20_pH']

p.performance_list_case_7 = ['Must_ME_25_Biomass', 'Must_ME_25_Nitrogen', 'Must_ME_25_Sugar', 'Must_ME_25_CO2',
                             'Must_ME_25_Ethanol', 'Must_ME_25_Succinic_Acid', 'Must_ME_25_Polymeric_Pigments',
                             'Lees_ME_25_Polymeric_Pigments', 'Must_ME_25_Anthocyanins', 'Cap_ME_25_Anthocyanins',
                             'Lees_ME_25_Anthocaynins', 'Must_ME_25_Tannins', 'Cap_ME_25_Tannins', 'Lees_ME_25_Tannins',
                             'Must_ME_25_Tartaric_Acid', 'Cap_ME_25_Tartaric_Acid', 'Lees_ME_25_Tartaric_Acid',
                             'Must_ME_25_Malic_Acid', 'Cap_ME_25_Malic_Acid', 'Lees_ME_25_Malicic_Acid',
                             'Must_ME_25_Colour_Density', 'Must_ME_25_Total_Phenolic_Index', 'Must_ME_25_pH']

p.performance_list_case_7 = ['Must_ME_28_Biomass', 'Must_ME_28_Nitrogen', 'Must_ME_28_Sugar', 'Must_ME_28_CO2',
                             'Must_ME_28_Ethanol', 'Must_ME_28_Succinic_Acid', 'Must_ME_28_Polymeric_Pigments',
                             'Lees_ME_28_Polymeric_Pigments', 'Must_ME_28_Anthocyanins', 'Cap_ME_28_Anthocyanins',
                             'Lees_ME_28_Anthocaynins', 'Must_ME_28_Tannins', 'Cap_ME_28_Tannins', 'Lees_ME_28_Tannins',
                             'Must_ME_28_Tartaric_Acid', 'Cap_ME_28_Tartaric_Acid', 'Lees_ME_28_Tartaric_Acid',
                             'Must_ME_28_Malic_Acid', 'Cap_ME_28_Malic_Acid', 'Lees_ME_28_Malicic_Acid',
                             'Must_ME_28_Colour_Density', 'Must_ME_28_Total_Phenolic_Index', 'Must_ME_28_pH']

# Import the experimental data
df_must = pd.read_csv('D:\\Stellenbosch\\2. PhD\\2022\\0. Experimental\\Data Sets\\wine_data_N.csv')
df_lees = pd.read_csv('D:\\Stellenbosch\\2. PhD\\2022\\0. Experimental\\Data Sets\\lees_data.csv')
df_cap = pd.read_csv('D:\\Stellenbosch\\2. PhD\\2022\\0. Experimental\\Data Sets\\cap_data_normalised.csv')
df_grape = pd.read_csv('D:\\Stellenbosch\\2. PhD\\2022\\0. Experimental\\Data Sets\\grape_data.csv')

# Functions used to sort the data frames and obtain the averages and standard deviations
df = types.SimpleNamespace(Must=df_must, Lees=df_lees, Cap=df_cap, Grape=df_grape)
data, list_values = experimental_data(df, p.phase_name)
sd, average = standard_deviations(df, p.phase_name, p)

# This declares the time points for plotting experimental data
time_points = (df_must['Hour_Must']).to_numpy()
time_array = time_points*60
time_array[10] = 13019

# Case determines the type of optimisation. The first element is the cultivar, the second is the temperature of interest
# the options for the first element are 'CS', 'SH', 'ME', or 'n/a'; the options ofr the second element are '20', '25',
# '28' or 'n/a'
case1, case2, case3, case4, case5, case6, case7, case8, case9, case10, case11, case12 = ['CS', '20'], ['CS', '25'], ['CS', '28'], ['SH', '20'], ['SH', '25'], ['SH', '28'], ['ME', '20'], ['ME', '25'], ['ME', '28'], ['n/a', '20'], ['n/a', '25'], ['n/a', '28']
