import numpy as np
import matplotlib.pyplot as plt


def Succes_detector(iterate_variables):
    Succes = [True]*4
    if Bearing_check(...,...,...)[0] == False:
        Succes[0] = False

    if Pull_through_check(...,...,...)[0] == False:
        Succes[1] = False

    if Thermal_check(...,...,...)[0] == False:
        Succes[2] = False


# --------- Calculation of reaction forces 


# --------- Lug stresses
#input: (n_lugs,h,t,w,D_1,F_vec,M_vec,Delta_yeld,Tau_strength)
#output: (True/False)

# --------- Fastener CG
#input: (w_plate,h_plate,x_off,z_off,D_2,Mat_type)
#output: (x_cg,z_cg,Fast_list)


# --------- List of materials
# Lug_materials_list = [[Delta_yield , Delta_bear, Tau_yield, Mat_type, Mat_name],[],[]]
# Fastener_materials_list = [[Delta_yield , Delta_bear, Tau_yield, Mat_type, Mat_name],[],[]]


# --------- Initial values

# Lug_material_prop_list = Materials_list[i]
# Fastener_material_prop_list = Materials_list[i]

# iterate_variables = [ Material_prop_list , Fastener_material_prop_list ,  t_1 , t_2 ,  h , w , D_1 , D_2 , w_plate , h_plate , x_off , z_off ]




min_thickness_1 = 3*10**(-3)
max_thickness_1 = 30*10**(-3)

t_1 = min_thickness_1

while t_1 < max_thickness_1:

    t_1 += (max_thickness_1-min_thickness_1)/100


#input: (w_plate,h_plate,x_off,z_off,D_2,Mat_type)
#output: (x_cg,z_cg,Fast_list)