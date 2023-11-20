import numpy as np
import matplotlib.pyplot as plt

import 
import 
import 

def Succes_detector(iterate_variables):


    Succes_list = [True]*4
    if Bearing_check(...,...,...)[0] == False:
        Succes_list[0] = False
        print("Fails by bearing check")

    if Pull_through_check(...,...,...)[0] == False:
        Succes_list[1] = False
        print("Fails by pull through check")

    if Thermal_check(...,...,...)[0] == False:
        Succes_list[2] = False
        print("Fails by thermal stress check")
    return Succes_list



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

#Material = [Youngus_Modulus,Yield_Strength,density, density]
Aluminium_7075 = [76*10**9, 225*10**6 , 2.8*10**3]


#-------- Step size configuration
min_thickness_1 = 3*10**(-3)
max_thickness_1 = 30*10**(-3)
thickness_1_steps = 100
thickness_1_stepsize = (max_thickness_1-min_thickness_1)/thickness_1_steps

min_thickness_2 = 3*10**(-3)
max_thickness_2 = 30*10**(-3)
thickness_2_steps = 100
thickness_2_stepsize = (max_thickness_2-min_thickness_2)/thickness_2_steps

min_D_1 =
max_D_1 =
D_1_steps =
D_1_step_size =

min_D_2 =
max_D_2 =
D_2_steps =
D_2_step_size =

min_w_plate =
max_w_plate =
w_steps = 
w_step_size = 

min_h_plate =
max_h_plate =
h_steps = 
h_step_size = 


iterate_variables = [Material_properties]
Valid_configurations =[]

if Success_dector(iteration_variables).count(True) == 3:
    Valid_configurations.append(iteration_variables)
else:
    if Succes_detector(iterate_variables)[1] == False:
        

        iterate_variables[3] = 
        




t_1 = min_thickness_1

while t_1 < max_thickness_1:

    t_1 += (max_thickness_1-min_thickness_1)/100


#input: (w_plate,h_plate,x_off,z_off,D_2,Mat_type)
#output: (x_cg,z_cg,Fast_list)