# ======================================================================================================================
# Script Name: Optimisation_Variables.py
# File Name: Three_Phase_Mathematical_Model
# Date: 21-06-2021
# Creator: K.N. Lambrecht
# Purpose: This script is used for solution of the system of ODEs and the optimisation of parameters
# ======================================================================================================================
import types

from Main import *
from Functions import *

# Declare variables of interest
interest_list = p.interest_biomass
experimental_list = p.exp_components_biomass

regression_value = p.param_reg_biomass
regression_list = p.param_reg_names_biomass
fixed_values = p.param_biomass
fixed_list = p.param_names_biomass
dict1, dict2 = regression_value.__dict__, fixed_values.__dict__
dict_bio = {**dict1, **dict2}
all_param_1 = types.SimpleNamespace(**dict_bio)

lower_bio = namespace_to_vector(lp_biomass, regression_list)
upper_bio = namespace_to_vector(up_biomass, regression_list)
bounds = (lower_bio, upper_bio)

gs_bio = (namespace_to_vector(regression_value, regression_list))
# The following is a stoichiometric matrix for the three-phase fermentation vessel

p.SM = types.SimpleNamespace(X=np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]),
                             N=np.array([0, -p.param_all.Y_NX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             S=np.array([-p.param_all.Y_ES*p.param_all.B, -p.param_all.Y_ES, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             CO2=np.array([p.param_all.B, p.param_all.Y_CX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             E=np.array([p.param_all.B, p.param_all.Y_EX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             GL=np.array([p.param_all.Y_GS, p.param_all.Y_GX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             SC=np.array([p.param_all.Y_SCS, p.param_all.Y_SCX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             CPP=np.array([0, 0, -1/V_MS, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             MPP_LS=np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             CA=np.array([0, 0, 0, -1, 1/V_MS, -1/V_MS, 0, 0, 0, 0, 0, 0, 0]),
                             CA_GR=np.array([0, 0, 0, 0, -1/V_GR, 0, 0, 0, 0, 0, 0, 0, 0]),
                             MA_LS=np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]),
                             CT=np.array([0, 0, 0, 0, 0, 0, 1/V_MS, -1/V_MS, 0, 0, 0, 0, 0]),
                             CT_GR=np.array([0, 0, 0, 0, 0, 0, -1/V_GR, 0, 0, 0, 0, 0, 0]),
                             MT_LS=np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]),
                             CTA=np.array([0, 0, 0, 0, 0, 0, 0, 0, 1/V_MS, -1/V_MS, 0, 0, 0]),
                             CTA_GR=np.array([0, 0, 0, 0, 0, 0, 0, 0, -1/V_GR, 0, 0, 0, 0]),
                             MTA_LS=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]),
                             CMA=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/V_MS, -1/V_MS, 0]),
                             CMA_GR=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/V_GR, 0, 0]),
                             MMA_LS=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]),
                             M_LS=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, p.param_all.Y_LX*V_MS]))

# The following generates two vectors containing the experimental values and standard deviations for each point. Based
# on the case
exp_vector1, stdev_vector1 = residuals(case1, average, time_array, sd, interest_list, experimental_list)
exp_vector2, stdev_vector2 = residuals(case2, average, time_array, sd, interest_list2, experimental_list2)
exp_vector3, stdev_vector3 = residuals(case3, average, time_array, sd, interest_list3, experimental_list3)
exp_vector4, stdev_vector4 = residuals(case4, average, time_array, sd, interest_list4, experimental_list4)
exp_vector5, stdev_vector5 = residuals(case5, average, time_array, sd, interest_list, experimental_list)
exp_vector6, stdev_vector6 = residuals(case6, average, time_array, sd, interest_list, experimental_list)
exp_vector7, stdev_vector7 = residuals(case7, average, time_array, sd, interest_list, experimental_list)
exp_vector8, stdev_vector8 = residuals(case8, average, time_array, sd, interest_list, experimental_list)
exp_vector9, stdev_vector9 = residuals(case9, average, time_array, sd, interest_list, experimental_list)
exp_vector10, stdev_vector10 = residuals(case10, average, time_array, sd, interest_list, experimental_list)
exp_vector11, stdev_vector11 = residuals(case11, average, time_array, sd, interest_list, experimental_list)
exp_vector12, stdev_vector12 = residuals(case12, average, time_array, sd, interest_list, experimental_list)

# The following name space contains the initial values for the system of Ordinary Differential Equations
x1 = types.SimpleNamespace(X=0.25, N=average.Must_CS_20_Nitrogen[0], S=average.Must_CS_20_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_CS_20_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_CS_20_Anthocyanins[0], CA_GR=1009.6, MA_LS=0,
                           CT=average.Must_CS_20_Tannins[0], CT_GR=5384.49, MT_LS=0,
                           CTA=average.Must_CS_20_Tartaric_Acid[0], CTA_GR=3, MTA_LS=0,
                           CMA=average.Must_CS_20_Malic_Acid[0], CMA_GR=2.11, MMA_LS=0)

x2 = types.SimpleNamespace(X=0.25, N=average.Must_CS_25_Nitrogen[0], S=average.Must_CS_25_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_CS_25_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_CS_25_Anthocyanins[0], CA_GR=data.Grape.CS_Repeat_3_Anthocyanins_Grape,
                           MA_LS=0, CT=average.Must_CS_25_Tannins[0], CT_GR=data.Grape.CS_Repeat_3_Tannins_Grape,
                           MT_LS=0, CTA=average.Must_CS_25_Tartaric_Acid[0],
                           CTA_GR=data.Grape.CS_Repeat_3_Tartaric_Acid_Grape,
                           MTA_LS=0, CMA=average.Must_CS_25_Malic_Acid[0],
                           CMA_GR=data.Grape.CS_Repeat_3_Malic_Acid_Grape, MMA_LS=0)

x3 = types.SimpleNamespace(X=0.25, N=average.Must_CS_28_Nitrogen[0], S=average.Must_CS_28_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_CS_28_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_CS_28_Anthocyanins[0], CA_GR=data.Grape.CS_Repeat_3_Anthocyanins_Grape,
                           MA_LS=0, CT=average.Must_CS_28_Tannins[0], CT_GR=data.Grape.CS_Repeat_3_Tannins_Grape,
                           MT_LS=0, CTA=average.Must_CS_28_Tartaric_Acid[0],
                           CTA_GR=data.Grape.CS_Repeat_3_Tartaric_Acid_Grape, MTA_LS=0,
                           CMA=average.Must_CS_28_Malic_Acid[0], CMA_GR=data.Grape.CS_Repeat_3_Malic_Acid_Grape,
                           MMA_LS=0)

x4 = types.SimpleNamespace(X=0.25, N=average.Must_SH_20_Nitrogen[0], S=average.Must_SH_20_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_SH_20_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_SH_20_Anthocyanins[0], CA_GR=data.Grape.SH_Repeat_1_Anthocyanins_Grape,
                           MA_LS=0, CT=average.Must_SH_20_Tannins[0], CT_GR=data.Grape.SH_Repeat_1_Tannins_Grape,
                           MT_LS=0, CTA=average.Must_SH_20_Tartaric_Acid[0],
                           CTA_GR=data.Grape.SH_Repeat_3_Tartaric_Acid_Grape,
                           MTA_LS=0, CMA=average.Must_SH_20_Malic_Acid[0],
                           CMA_GR=data.Grape.SH_Repeat_1_Malic_Acid_Grape, MMA_LS=0)

x5 = types.SimpleNamespace(X=0.25, N=average.Must_SH_25_Nitrogen[0], S=average.Must_SH_25_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_SH_25_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_SH_25_Anthocyanins[0], CA_GR=data.Grape.SH_Repeat_1_Anthocyanins_Grape,
                           MA_LS=0, CT=average.Must_SH_25_Tannins[0], CT_GR=data.Grape.SH_Repeat_1_Tannins_Grape,
                           MT_LS=0, CTA=average.Must_SH_25_Tartaric_Acid[0],
                           CTA_GR=data.Grape.SH_Repeat_3_Tartaric_Acid_Grape,
                           MTA_LS=0, CMA=average.Must_SH_25_Malic_Acid[0],
                           CMA_GR=data.Grape.SH_Repeat_1_Malic_Acid_Grape, MMA_LS=0)

x6 = types.SimpleNamespace(X=0.25, N=average.Must_SH_28_Nitrogen[0], S=average.Must_SH_28_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_SH_28_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_SH_28_Anthocyanins[0], CA_GR=data.Grape.SH_Repeat_1_Anthocyanins_Grape,
                           MA_LS=0, CT=average.Must_SH_28_Tannins[0], CT_GR=data.Grape.SH_Repeat_1_Tannins_Grape,
                           MT_LS=0, CTA=average.Must_SH_28_Tartaric_Acid[0],
                           CTA_GR=data.Grape.SH_Repeat_3_Tartaric_Acid_Grape,
                           MTA_LS=0, CMA=average.Must_SH_28_Malic_Acid[0],
                           CMA_GR=data.Grape.SH_Repeat_1_Malic_Acid_Grape, MMA_LS=0)

x7 = types.SimpleNamespace(X=0.25, N=average.Must_ME_20_Nitrogen[0], S=average.Must_ME_20_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_ME_20_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_ME_20_Anthocyanins[0], CA_GR=data.Grape.ME_Repeat_2_Anthocyanins_Grape,
                           MA_LS=0, CT=average.Must_ME_20_Tannins[0], CT_GR=data.Grape.ME_Repeat_2_Tannins_Grape,
                           MT_LS=0, CTA=average.Must_ME_20_Tartaric_Acid[0],
                           CTA_GR=data.Grape.ME_Repeat_2_Tartaric_Acid_Grape, MTA_LS=0,
                           CMA=average.Must_ME_20_Malic_Acid[0],
                           CMA_GR=data.Grape.ME_Repeat_2_Malic_Acid_Grape, MMA_LS=0)

x8 = types.SimpleNamespace(X=0.25, N=average.Must_ME_25_Nitrogen[0], S=average.Must_ME_25_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_ME_25_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_ME_25_Anthocyanins[0], CA_GR=data.Grape.ME_Repeat_2_Anthocyanins_Grape,
                           MA_LS=0, CT=average.Must_ME_25_Tannins[0], CT_GR=data.Grape.ME_Repeat_2_Tannins_Grape,
                           MT_LS=0, CTA=average.Must_ME_25_Tartaric_Acid[0],
                           CTA_GR=data.Grape.ME_Repeat_2_Tartaric_Acid_Grape, MTA_LS=0,
                           CMA=average.Must_ME_25_Malic_Acid[0], CMA_GR=data.Grape.ME_Repeat_2_Malic_Acid_Grape,
                           MMA_LS=0)

x9 = types.SimpleNamespace(X=0.25, N=average.Must_ME_28_Nitrogen[0], S=average.Must_ME_28_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_ME_28_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_ME_28_Anthocyanins[0], CA_GR=data.Grape.ME_Repeat_2_Anthocyanins_Grape,
                           MA_LS=0, CT=average.Must_ME_28_Tannins[0], CT_GR=data.Grape.ME_Repeat_2_Tannins_Grape,
                           MT_LS=0, CTA=average.Must_ME_28_Tartaric_Acid[0],
                           CTA_GR=data.Grape.ME_Repeat_2_Tartaric_Acid_Grape, MTA_LS=0,
                           CMA=average.Must_ME_28_Malic_Acid[0], CMA_GR=data.Grape.ME_Repeat_2_Malic_Acid_Grape,
                           MMA_LS=0)

