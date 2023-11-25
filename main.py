from calendar import c
import numpy as np
import matplotlib.pyplot as plt

import 
import 
import 

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
Aluminium_7075 = ["metal",76*10**9, 225*10**6 , 2.8*10**3]


#-------- Step size configuration
min_t_1 = 3*10**(-3)
max_t_1 = 30*10**(-3)
t_1_steps = 100
t_1_stepsize = (max_t_1-min_t_1)/t_1_steps

min_D_1 = 
max_D_1 =
D_1_steps =
D_1_step_size =

min_w_1 =
max_w_1 =
w_1_steps = 
w_1_stepsize = 

chosen_h_1 = 0.05

min_t_2 = 3*10**(-3)
max_t_2 = 30*10**(-3)
t_2_steps = 100
t_2_stepsize = (max_t_1-min_t_1)/t_1_steps


min_D_2 =
max_D_2 =
D_2_steps =
D_2_step_size =

min_w_2 =
max_w_2 =
w_2_steps = 
w_2_step_size = 

min_h_2 =
max_h_2 =
h_steps = 
h_step_size = 

chosen_n_r = 
chosen_n_c = 





def Succes_detector(iterate_variables):

    Succes_list = [True]*4

    if Lug_check(...,...,...)[0] == False:
        Succes_list[0] = False
        print("Fails by bearing check")
    
    if Bearing_check_fastener(...,...,...)[1] == False:
        Succes_list[1] = False
        print("Fails by bearing check")

    if Pull_through_check(...,...,...)[2] == False:
        Succes_list[2] = False
        print("Fails by pull through check")

    if Thermal_check(...,...,...)[3] == False:
        Succes_list[3] = False
        print("Fails by thermal stress check")
    return Succes_list



#iterate_variables = [Aluminium_7075, max_thickness_1,  max_w_1, max_D_1 , max_thickness_2 , max_D_2 , max_w_plate , max_h_plate ]

#---------------- FLANGE DESIGN ----------------------------
flange_variables = [Aluminium_7075, max_t_1,  max_w_1, chosen_h_1, max_D_1]

# Minimal thickness

while flange_variables[2] <= min_w_1:
    flange_variables[4] = max_D_1
    while flange_variables[4] <= min_D_1:
        #FUNCTION TO CALCULATE THICKNESS FOR AXIAL
        #FUNCTION TO CALCULATE THICKNESS FOR TRANSVERSAL
        flange_variables[1] = max_t_1
        
        while abs((Lug_check(...,...,...)[1])/(flange_variables[0][2]) - 1) <= 0.01 and (Lug_check(...,...,...)[1])/(flange_variables[0][2]) <= 1:
            flange_variables[1]*(Lug_check(...,...,...)[1])/(flange_variables[0][2])
        
        flange_variables[3] -= D_1_step_size

        if Lug_check(...,...,...)[0] == True:
            Valid_flange.append(flange_variables)

    flange_variables[2] -= w_1_stepsize

#---------------- PLATE DESIGN -------------------------

plate_variables = [Aluminium_7075, max_thickness_2,  max_w_2, max_h_2, max_D_2, chosen_n_r, chosen_n_c]

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