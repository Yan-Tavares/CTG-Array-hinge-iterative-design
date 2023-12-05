#---- CALCULATE DIMENSIONS of LUG ----

import numpy as np
import math
from Functions import Flange_k_factors as K_factors
from Functions import Thickness_Iterator as TI
from Functions import FastenerCG as Fcg
from Functions import Fastener_design as Fd

#------------ Functions

def Flanges_inertia(h,t,W):
    I_xx = ((t*W**3)/12)*2
    I_zz = ((W*t**3)/12 + t*W*(h/2+t/2)**2)*2

    return I_xx,I_zz

def Dual_flange_bending_stress(P_x,P_z,h,t,W,Flange_L,x,z):
    I_xx,I_zz = Flanges_inertia(h,t,W)
    Sigma_y = (P_x * Flange_L)*(x)/I_zz + (P_z* Flange_L)*(z)/I_xx

    return Sigma_y

def Single_flange_bending_stress(P_x,P_z,t,W,Flange_L,x,z):
    I_xx = ((t*W**3)/12)
    I_zz = ((W*t**3)/12)
    Sigma_y = (P_x * Flange_L)*(x)/I_zz + (P_z* Flange_L)*(z)/I_xx

    return Sigma_y

#------------ Highest loads
SF_x = 4.5
SF_y = 4.5
SF_z = 4.5

P_x = 254 * SF_x
P_y = 704 * SF_y
P_z = 496 * SF_z
M_y = 0

print("-----------------------------------------")
print("Maximimum reaction loads with safety factor")
print("-----------------------------------------\n")
print(f"{'R_x:':<15}{P_x:<8.0f}{'[N]':<15}{'SF:':<4}{SF_x:<}")
print(f"{'R_y:':<15}{P_y:<8.0f}{'[N]':<15}{'SF:':<4}{SF_y:<}")
print(f"{'R_z:':<15}{P_z:<8.0f}{'[N]':<15}{'SF:':<4}{SF_z:<}","\n")


#-------- Material Data
Aluminium_6061_T6= ["metal","Aluminium",270*10**6,2.7*10**3]
Steel_8630 = ["metal","Steel",550*10**6,7.85*10**3]

Mat_list_flanges = [Aluminium_6061_T6, Steel_8630]
Mat_list_fasteners = [Aluminium_6061_T6, Steel_8630]

#-----------------------SET STEPSIZE -------------------

min_t = 0.4*10**(-3) #[m]; min. size to 3D print aluminium is 0.381 mm
max_t = 20*10**(-3) #[m]
t_steps = 20
t_stepsize = (max_t-min_t)/t_steps

min_w = 4 *10**(-3) #[m]
max_w = 150 *10**(-3) #[m]
w_steps = 20
w_stepsize = (max_w-min_w)/w_steps

min_D = 2 *10**(-3) #[m]; min. hole size to 3D print metals
max_D = max_w - 3*10**(-3)
D_steps = 40
D_stepsize = (max_D-min_D)/D_steps


min_D_2 = 2 *10**(-3) #[m]; min. hole size to 3D print metals
max_D_2 = 150 * 10**(-3)
D_2_steps = 200
D_2_stepsize = (max_D-min_D)/D_steps

min_F_num = 3
max_F_num = 15

#-----------------------ITERATIVE DESIGN CALCULATION ----------------

#---- PRE DEFINED VALUES -----
h = 0.0125
min_mass = 10 #[kg]


print("------------------------------------------------------------")
print("--------------- Flange Iteration Process ------------------- \n")

#iterate over materials: 0 = Aluminium, 1 = Steel
for Material in Mat_list_flanges:
    type_identifier = Material[1]
    Sigma_yield = Material[2]
    rho = Material[3]

    D = min_D
    #iterate over D
    while D <= max_D:

        w = min_w
        #iterate over w
        while w <= max_w:

            if D > w -3*10**(-3):
                w += w_stepsize
                continue
            
            t = min_t
            #iterate over t
            while t <= max_t:
                #find K's stuff

                #---- Calculate Abr and Aav
                A1 = ((D/2-np.cos(np.pi / 4) * D / 2) + (w - D) / 2) * t
                A2 = (w - D) * t / 2
                A3 = A2
                A4 = A1
                Aav =  6 / (3 / A1 + 1 / A2 + 1 / A3 + 1 / A4)
                Abr = D*t
                At = w * t
                A_pin = (D/2)**2 * math.pi
                #-----------------------------------

                #---- Stress K factors
                Kty = K_factors.Trans_Factor (Aav, Abr)
                Kt = K_factors.TensionyieldFactor(w, D,Material[1])
                Kbry = K_factors.Shear_out_Factor(w, D, t)
                #-----------------------------------

                #---- Calculate allowed stresses on the ring region
                Sigma_bearing = ((w-D)/2)/D * Sigma_yield

                Pu = Kt * Sigma_yield * At
                Pbry = Kbry * Sigma_bearing * Abr
                Pty = Kty * Abr * Sigma_yield
                Pmin = min(Pu, Pbry)
                
                if Pmin <= 0: #Happens when no axial force is allowed
                    t += t_stepsize
                    continue
                #-----------------------------------
                
                #---- Calculate bending nominal stresses at flange/plate intersection
                Flange_L = w + w/2 #Venant Principle to let stresses even out
                Shoulder_fillet = t/2
                Sigma_y_ben_flange = Dual_flange_bending_stress(P_x,P_z,h,t,w,Flange_L,(h/2 +t),(w/2)) #Single_flange_bending_stress(P_x/2,P_z/2,t,w,Flange_L,t/2,w/2)

                Sigma_y_axial_flange = (P_y/2)/(w*t)
                Sigma_y_total_flange = Sigma_y_ben_flange + Sigma_y_axial_flange
                Kt_shoulder_fillet = K_factors.Edge_fillet_factor(3,Shoulder_fillet,t)

                Sigma_shoulder_fillet = Sigma_y_total_flange * Kt_shoulder_fillet


                #---- Calculate safety margin
                Ra = P_y / Pmin
                Rtr = P_z / Pty
                MS = 1 / ((Ra**1.6 + Rtr**1.6)**0.625)# -1

                if MS >= 1:
                    MS = 1
                #-----------------------------------
                
                #calculate shear stress on pin (pin on 1 lug carries half of the reaction force on the lug configuration)
                tau_pin = (P_y * 0.5) / A_pin
                tau_max = 0.5*Sigma_yield #Tresca failure criterion

                #Calculate

                #check if the maximum loads are smaller than the ones we are facing, and if pin doesn't yield
                if P_y <= Pmin * MS and P_z <= Pty * MS and tau_pin < tau_max and Sigma_shoulder_fillet < Sigma_yield:
                    
                    #calculate meaningfull mass (i.e. only considering the circular part (tho this can be negative!!))
                    half_ring_mass = rho * 0.5 * math.pi * ((0.5*w)**2 - (0.5*D)**2) * t
                    
                    Flange_mass = ((w * Flange_L * t -  0.5 * t * math.pi*(0.5*D)**2) * rho + half_ring_mass)* 2 


                    #if new mass smaller than previous mass: store
                    if Flange_mass < min_mass:
                        min_mass = Flange_mass
                        Best_config_flange = [Material[1],Flange_mass,t,D,w,MS]
                        
                        print("-----------------------------------------")
                        print("Better configuration for Flange found")
                        print("-----------------------------------------\n")
                        print("---- Efficiency and Concentration Factors\n")

                        print(f"{' Kbry_Ring Shear Bearing efficiency factor:':<60}{Kbry:<12.3f}")
                        print(f"{' Kt_Ring Tension efficiency factor : ':<60}{Kt:<12.3f}")
                        print(f"{' Kt Shoulder fillet: ':<60}{Kt_shoulder_fillet:<12.3f}","\n")

                        print("---- Ring Allowed Loads\n")

                        print(f"{' Allowed axial force: ':<60}{Pmin * MS:<12.1f}{'[N]':<}")
                        print(f"{' Allowed transversal force: ':<60}{Pty * MS:<12.1f}{'[N]':<}")
                        print(f"{' MS for oblique loading: ':<60}{MS:<12.3f}","\n")

                        print("---- Flange stresses\n")
                
                        print(f"{' Nominal root max bending stress: ':<60}{Sigma_y_ben_flange *10**(-6):<12.1f}{'[MPa]':<}")
                        print(f"{' Nominal root axial stress: ':<60}{Sigma_y_axial_flange *10**(-6):<12.1f}{'[MPa]':<}")
                        print(f"{' Total root nominal stress: ':<60}{Sigma_y_total_flange *10**(-6):<12.1f}{'[MPa]':<}")
                        print(f"{' Stress at shoulder fillet: ':<60}{Sigma_shoulder_fillet*10**(-6):<12.3f}{'[MPa]':<}")
                        print(f"{' Pin stress: ':<60}{tau_pin*10**(-6):<12.1f}{'[MPa]':<}","\n")

                        print("---- Flange Properties\n")
                        #print(f"{' Dimension ':<60}{'Value ':<12}{'Unity':<}")
                        print(f"{' Flange mass: ':<60}{Flange_mass* 10**3:<12.3f}{'[g]':<}")
                        print(f"{' W: ':<60}{w*10**3:<12.3f}{'[mm]':<}")
                        print(f"{' D: ':<60}{D*10**3:<12.3f}{'[mm]':<}")
                        print(f"{' t: ':<60}{t*10**3:<12.3f}{'[mm]':<}")
                        print(f"{' h: ':<60}{h*10**3:<12.3f}{'[mm]':<}")
                        print(f"{' L: ':<60}{(Flange_L+w/2)*10**3:<12.3f}{'[mm]':<}")
                        print(f"{' Shoulder fillet: ':<60}{Shoulder_fillet*10**3:<12.3f}{'[mm]':<}")

                t += t_stepsize
            w += w_stepsize
        D += D_stepsize

print("------------------------------------------------------------")
print("----------------- Plate Iteration Process -------------------\n")

w = Best_config_flange[4]
t_1 = Best_config_flange[2]

for Material in Mat_list_fasteners:
    Sigma_y = Material[2]
    rho = Material[3]

    D_2 = min_D_2
    while D_2 <= max_D_2:
        F_num = min_F_num

        while F_num <= max_F_num:
            x_i, z_i, x_cg, z_cg, area, sumx, width = Fcg.fastener_CG(h,t_1,D_2,w, Material[0],F_num)
            
            print(sumx,width)

            F_num += 1

        D_2 += D_2_stepsize








#-----------------------SET STEPSIZE -------------------

# for Material in Mat_list_flanges:
#     Sigma_y = Material[2]
#     rho = Material[3]
    
#     D_2 = min_D_2

#     while D_2 <= max_D_2:
#         F_num = min_F_num

#         while F_num <= max_F_num:
#             t_2 = TI.Min_thickness_finder(h, t_1, D_2, "metal", w, F_num, P_x, P_z, M_y)

#             Fas_and_plate_mass = 

#             if 
#             F_num +=1

#         D_2 += D_2_stepsize