#---- CALCULATE DIMENSIONS of LUG ----

import numpy as np
import math
from Functions import Curve_generator as Cgen
from Functions import Reaction_force_calculator as RecF

Spacecraft_mass = 416
Solar_array_mass = 10
Thruster_Fz = 456
a_x = 0
a_y = 0
a_z = Thruster_Fz/Spacecraft_mass
alpha_x = 0
alpha_y = 0
alpha_z = 0


def Trans_Factor (Aav, Abr):
    x = Aav/Abr
    K_ty = Cgen.Closest_data_point("Functions\Reader Graphs\D-15","3",x)

    return K_ty

def Shear_out_Factor(w, D, t):
    e = (w/2)
    x = e/(D)

    if x >= 4:
        x = 4

    if t/D >= 0.6:
        Curve_choice = "0.60"

    if t/D >= 0.195 and t/D  <0.6:
        Curve_choice = '{:.2f}'.format(round(t/D,1))
    
    if t/D > 0.06 and t/D  < 0.195:
        Curve_choice = '{:.2f}'.format(round(t/D,2))
    
    if t/D <= 0.06:
        Curve_choice = "0.06"


    K_bry = Cgen.Closest_data_point("Functions\Reader Graphs\D-14", str(Curve_choice),x)

    #print("e/(D): ", x,"Curve_choice: ",Curve_choice,"Kbry: ",K_bry)
    if K_bry <= 0:
        K_bry = 0

    return K_bry

def TensionyieldFactor(w, D,Material):
    x = w/D
    K_t = Cgen.Closest_data_point("Functions\Reader Graphs\D-12", Material,x)

    return K_t




##-------------------------SAFETY FACTOR STUFF (SANTIAGO'S CODE)
# #Determine t1, w, D1

# import scipy.io

# F1 = (z*L3)/(2*h1)
# Pvec = np.array([0, y+F1, z])
# Pmag = np.linalg.norm(Pvec)


# def MS(t1, D, w):
#     sigma_yield = 123 #add yield stress of the material
#     A1 = ((D/2-cos(pi / 4) * D / 2) + (w - D) / 2) * t1
#     A2 = (w - D) * t1 / 2
#     A3 = (w - D) * t1 / 2
#     A4 = ((D/2-cos(pi / 4) * D / 2) + (w - D) / 2) * t1

#     #To avoid if denominator is 0
#     if A1 ==0 or A2==0 or A3 == 0 or A4 ==0:
#         return 0
    
#     Aav =  6 / (3 / A1 + 1 / A2 + 1 / A3 + 1 / A4)
#     Abr = D*t1

#     #Decomposing P into the transverse load and the axial load

#     axial = Pvec[1]
#     transverse  = Pvec[2]
#     #transverse shear factor
#     kty = TransFactor(Aav, Abr)
#     Pty = kty*Abr*sigma_yield
#     Rtr = transverse/Pty
#     #axial shearfactor
#     #will be tension and also shear
#     ksy = ShearyieldFactor(w, D)
#     Aeff =  2*((w-D)/2+D/2)*t1
#     Psy = ksy * Aeff * sigma_yield
#     kt = TensionyieldFactor(w, D)
#     Psyt = kt*Aeff*sigma_yield
#     Ra = axial/min(Psy, Psyt)

#     MS = (1/((Ra**1.6)+(Rtr**1.6))**0.625)-1
#     return MS

#-------------------------------------------------------------------------

#JUTTA'S CODE

#Material list
Aluminium_6061_T6= ["Metal","Aluminium",270*10**6,2.7*10**3]
Steel_8630 = ["Metal","Steel",550*10**6,7.85*10**3]
Mat_list = [Aluminium_6061_T6, Steel_8630 ]

#-----------------------LOAD CALCULATION ----------------

# WE SHOULD USE THE REACTION FORCE CALCULATOR FUNCTION HERE BY THE WAY TO FIND F_R_1 USING TWO FORCE MEMBER ASSUMPTION
#The F_r_1 used here is not the correct one, it assumes we are taking all the load on one lug 

#Equations for defining reaction forces
#F_y_1 + F_y_2 = 0
#F_z_1 + F_z_2 = a_z * Solar_array_mass
#(F_z_1 * d_y) + (F_y_1 * d_z) + (F_z_2 * d_y) + (F_y_2 * d_z) = 0

#Assuming 2 force members
#F_r_1*cos(theta) + F_r_1*cos(theta) = 0
#F_r_1*sin(theta) + F_r_2*sin(theta)= a_z * Solar_array_mass
#(F_z_1 * d_y) + (F_y_1 * d_z) + (F_z_2 * d_y) + (F_y_2 * d_z) = 0

F_r_1 = (Solar_array_mass * a_z)/2
Solar_boom_ang = math.radians(45)
P_ax = F_r_1 * np.cos(Solar_boom_ang)
P_trans = F_r_1 * np.sin(Solar_boom_ang)

print("Axial force: ", P_ax)
print("Tranversal force", P_trans)

#-----------------------SET STEPSIZE (Yan stuff) -------------------
min_t = 0.4*10**(-3) #[m]; min. size to 3D print aluminium is 0.381 mm
max_t = 20*10**(-3) #[m]
t_steps = 10
t_stepsize = (max_t-min_t)/t_steps

min_w = 2.1 *10**(-3) #[m]
max_w = 150 *10**(-3) #[m]
w_steps = 10
w_stepsize = (max_w-min_w)/w_steps

min_D = 2 *10**(-3) #[m]; min. hole size to 3D print metals
max_D = max_w * 0.96
D_steps = 10
D_stepsize = (max_D-min_D)/D_steps

#-----------------------ITERATIVE DESIGN CALCULATION (Jutta, Yan stuff) ----------------

#----Set initial values for the iteration
min_mass = 100

#iterate over materials: 0 = Aluminium, 1 = Steel
for Material in Mat_list:
    Sigma_y = Material[2]
    rho = Material[3]
    type_identifier = Material[1]

    D = min_D
    #iterate over D
    while D <= max_D:

        w = min_w
        #iterate over w
        while w <= max_w:

            if D > w * 0.95:
                w += w_stepsize
                continue
            
            t = min_t
            
            #iterate over t
            while t <= max_t:
                #find K's stuff
                Kt = TensionyieldFactor(w, D,Material[1])
                Kbry = Shear_out_Factor(w, D, t)
                
                #calculate A's
                
                sigma_yield = Sigma_y #add yield stress of the material
                A1 = ((D/2-np.cos(np.pi / 4) * D / 2) + (w - D) / 2) * t
                A2 = (w - D) * t / 2
                A3 = (w - D) * t / 2
                A4 = ((D/2-np.cos(np.pi / 4) * D / 2) + (w - D) / 2) * t
                Aav =  6 / (3 / A1 + 1 / A2 + 1 / A3 + 1 / A4)
                Abr = D*t
                At = w * t

                A_pin = (D/2)**2 * math.pi

                Kty = Trans_Factor (Aav, Abr)

                #calculate axial & transverse loads (NOT incl. safety margin)
                Pu = Kt * Sigma_y * At
                Pbry = Kbry * Sigma_y * Abr
                Pty = Kty * Abr * Sigma_y
                Pmin = min(Pu, Pbry)
                
                

                #calculate Ra and Rtr

                if Pmin <= 0:
                    print("No axial force allowed, Pmin = ", Pmin)
                    t += t_stepsize
                    continue

                Ra = P_ax / Pmin
                Rtr = P_trans / Pty

                #calculate safety margin
                MS = 1 / ((Ra**1.6 + Rtr**1.6)**0.625) -1

                #calculate shear stress on pin (pin on 1 lug carries half of the reaction force on the lug configuration)
                tau_pin = (F_r_1 * 0.5) / A_pin
                tau_max = 0.5*Sigma_y #Tresca failure criterion

                #check if the maximum loads are smaller than the ones we are facing, and if pin doesn't yield
                if P_ax <= Pmin/ MS and P_trans <= Pty/ MS and tau_pin < tau_max:
                    
                    #calculate meaningfull mass (i.e. only considering the circular part (tho this can be negative!!))
                    mass = rho * 0.5 * math.pi * ((0.5*w)**2 - (0.5*D)**2) * t 

                    #if new mass smaller than previous mass: store
                    if mass < min_mass:
                        min_mass = mass
                        Best_config = [Material[0],mass,t,D,w,MS]

                        print("\nBetter configuration found")
                        print("Flange allowed axial force:", Pu)
                        print("Kbry:", Kbry)
                        print("kt:", Kt)
                        print("Margin of safety used:", MS)
                        print("Allowed axial force",Pmin/MS)
                        print("Allowed transversal force",Pty/MS)
                        print("Mass: ", mass,"\n")
                        print("Pin stress: ", tau_pin)

                
                t += t_stepsize
            w += w_stepsize
        D += D_stepsize

#give final choice for material, D, w, t and resulting mass
print(Best_config)

# print("material: ",materialopt)
# print("diameter: D = ",Dopt*0.001,"mm; thickness: t = ",topt*0.001, "mm; width: w = ",wopt*0.001,"mm")
# print("used safety factor: MS = ",MSopt)
# print("resulting in a mass of: m = ",minmass,"kg")


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
#Aluminium_7075 = ["metal",76*10**9, 225*10**6 , 2.8*10**3]


#-------- Step size configuration
# min_t_1 = 3*10**(-3)
# max_t_1 = 30*10**(-3)
# t_1_steps = 100
# t_1_stepsize = (max_t_1-min_t_1)/t_1_steps

# min_D_1 = 
# max_D_1 =
# D_1_steps =
# D_1_step_size =

# min_w_1 =
# max_w_1 =
# w_1_steps = 
# w_1_stepsize = 

# chosen_h_1 = 0.05

# min_t_2 = 3*10**(-3)
# max_t_2 = 30*10**(-3)
# t_2_steps = 100
# t_2_stepsize = (max_t_1-min_t_1)/t_1_steps


# min_D_2 =
# max_D_2 =
# D_2_steps =
# D_2_step_size =

# min_w_2 =
# max_w_2 =
# w_2_steps = 
# w_2_step_size = 

# min_h_2 =
# max_h_2 =
# h_steps = 
# h_step_size = 

# chosen_n_r = 
# chosen_n_c = 





# def Succes_detector(iterate_variables):

#     Succes_list = [True]*4

#     if Lug_check(...,...,...)[0] == False:
#         Succes_list[0] = False
#         print("Fails by bearing check")
    
#     if Bearing_check_fastener(...,...,...)[1] == False:
#         Succes_list[1] = False
#         print("Fails by bearing check")

#     if Pull_through_check(...,...,...)[2] == False:
#         Succes_list[2] = False
#         print("Fails by pull through check")

#     if Thermal_check(...,...,...)[3] == False:
#         Succes_list[3] = False
#         print("Fails by thermal stress check")
#     return Succes_list



# #iterate_variables = [Aluminium_7075, max_thickness_1,  max_w_1, max_D_1 , max_thickness_2 , max_D_2 , max_w_plate , max_h_plate ]

# #---------------- FLANGE DESIGN ----------------------------
# flange_variables = [Aluminium_7075, max_t_1,  max_w_1, chosen_h_1, max_D_1]

# # Minimal thickness

# while flange_variables[2] <= min_w_1:
#     flange_variables[4] = max_D_1
#     while flange_variables[4] <= min_D_1:
#         #FUNCTION TO CALCULATE THICKNESS FOR AXIAL
#         #FUNCTION TO CALCULATE THICKNESS FOR TRANSVERSAL
#         flange_variables[1] = max_t_1
        
#         while abs((Lug_check(...,...,...)[1])/(flange_variables[0][2]) - 1) <= 0.01 and (Lug_check(...,...,...)[1])/(flange_variables[0][2]) <= 1:
#             flange_variables[1]*(Lug_check(...,...,...)[1])/(flange_variables[0][2])
        
#         flange_variables[3] -= D_1_step_size

#         if Lug_check(...,...,...)[0] == True:
#             Valid_flange.append(flange_variables)

#     flange_variables[2] -= w_1_stepsize

# #---------------- PLATE DESIGN -------------------------

# plate_variables = [Aluminium_7075, max_thickness_2,  max_w_2, max_h_2, max_D_2, chosen_n_r, chosen_n_c]

# if Success_dector(iteration_variables).count(True) == 3:
#     Valid_configurations.append(iteration_variables)
# else:
#     if Succes_detector(iterate_variables)[1] == False:
        

#         iterate_variables[3] = 
        




# t_1 = min_thickness_1

# while t_1 < max_thickness_1:

#     t_1 += (max_thickness_1-min_thickness_1)/100


#input: (w_plate,h_plate,x_off,z_off,D_2,Mat_type)
#output: (x_cg,z_cg,Fast_list)
