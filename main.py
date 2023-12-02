#---- CALCULATE DIMENSIONS of LUG ----

import numpy as np
import math
from Functions import Flange_k_factors as K_factors
from Functions import Thickness_Iterator as TI

#------------ Functions




#------------ Highest loads
SF_x = 4.5
SF_y = 4.5
SF_z = 4.5

P_x = 254 * SF_x
P_y = 704 * SF_y
P_z = 496 * SF_z

print("-----------------------------------------")
print("Maximimum reaction loads with safety factor")
print("-----------------------------------------\n")
print(f"{'R_x:':<15}{P_y:<8.0f}{'[N]':<15}{'SF:':<4}{SF_x:<}")
print(f"{'R_y:':<15}{P_y:<8.0f}{'[N]':<15}{'SF:':<4}{SF_y:<}")
print(f"{'R_z:':<15}{P_y:<8.0f}{'[N]':<15}{'SF:':<4}{SF_z:<}","\n")


#-------- Material Data
Aluminium_6061_T6= ["Metal","Aluminium",270*10**6,2.7*10**3]
Steel_8630 = ["Metal","Steel",550*10**6,7.85*10**3]

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
max_D = max_w -3*10**(-3)
D_steps = 40
D_stepsize = (max_D-min_D)/D_steps


#-----------------------ITERATIVE DESIGN CALCULATION ----------------

#----Set initial values for the iteration
min_mass = 100

print("------------------------------------------------------------")
print("--------------- Flange Iteration Process ------------------- \n")

#iterate over materials: 0 = Aluminium, 1 = Steel
for Material in Mat_list_flanges:
    Sigma_y = Material[2]
    rho = Material[3]
    type_identifier = Material[1]

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
                Kt = K_factors.TensionyieldFactor(w, D,Material[1])
                Kbry = K_factors.Shear_out_Factor(w, D, t)
                
                sigma_yield = Sigma_y #add yield stress of the material
                A1 = ((D/2-np.cos(np.pi / 4) * D / 2) + (w - D) / 2) * t
                A2 = (w - D) * t / 2
                A3 = A2
                A4 = A1
                Aav =  6 / (3 / A1 + 1 / A2 + 1 / A3 + 1 / A4)
                Abr = D*t
                At = w * t

                A_pin = (D/2)**2 * math.pi

                Kty = K_factors.Trans_Factor (Aav, Abr)

                #calculate axial & transverse loads (NOT incl. safety margin)
                Pu = Kt * Sigma_y * At
                Pbry = Kbry * Sigma_y * Abr
                Pty = Kty * Abr * Sigma_y
                Pmin = min(Pu, Pbry)

                #calculate Ra and Rtr

                if Pmin <= 0:
                    #print("No axial force allowed, Pmin = ", Pmin)
                    t += t_stepsize
                    continue

                Ra = P_y / Pmin
                Rtr = P_z / Pty

                #calculate safety margin
                MS = 1 / ((Ra**1.6 + Rtr**1.6)**0.625)# -1

                if MS >= 1:
                    MS = 1

                #calculate shear stress on pin (pin on 1 lug carries half of the reaction force on the lug configuration)
                tau_pin = (P_y * 0.5) / A_pin
                tau_max = 0.5*Sigma_y #Tresca failure criterion

                #check if the maximum loads are smaller than the ones we are facing, and if pin doesn't yield
                if P_y <= Pmin * MS and P_z <= Pty * MS and tau_pin < tau_max:
                    
                    #calculate meaningfull mass (i.e. only considering the circular part (tho this can be negative!!))
                    half_ring_mass = rho * 0.5 * math.pi * ((0.5*w)**2 - (0.5*D)**2) * t
                    Flange_L = w + w/2 #Venant Principle
                    Flange_mass = (w * Flange_L * t -  0.5 * t * math.pi*(0.5*D)**2) * rho * 2

                    #if new mass smaller than previous mass: store
                    if Flange_mass < min_mass:
                        min_mass = Flange_mass
                        Best_config_flange = [Material[1],Flange_mass,t,D,w,MS]
                        
                        print("-----------------------------------------")
                        print("Better configuration for Flange found")
                        print("-----------------------------------------\n")
                        print(f"{'e/D:':<35}{w/(2*D):<12.3f}{'[none]':<}")
                        print(f"{'Kbry:':<35}{Kbry:<12.3f}{'[none]':<}")
                        print(f"{'kt: ':<35}{Kt:<12.3f}{'[none]':<}")
                        print(f"{'MS for oblique loading: ':<35}{MS:<12.3f}{'[none]':<}")
                        print(f"{'Allowed axial force: ':<35}{Pmin * MS:<12.1f}{'[N]':<}")
                        print(f"{'Allowed transversal force: ':<35}{Pty * MS:<12.1f}{'[N]':<}")
                        print(f"{'Ring mass: ':<35}{Flange_mass * 10**3:<12.2f}{'[g]':<}")
                        print(f"{'Pin stress: ':<35}{tau_pin*10**(-6):<12.1f}{'[GPa]':<}","\n")

                t += t_stepsize
            w += w_stepsize
        D += D_stepsize

print("------------------------------------------------------------")
print("----------------- Plate Iteration Process -------------------\n")


#-----------------------PRE DEFINED VALUES -------------------
h = 0.05
w = Best_config_flange[4]
t_1 = Best_config_flange[2]

for Material in Mat_list_fasteners:
    Sigma_y = Material[2]
    rho = Material[3]
    type_identifier = Material[1]

Min_thickness_finder(h, t_1, d2, "metal", w, fnum, P_x, P_z, My)

#-----------------------SET STEPSIZE -------------------

for Material in Mat_list_flanges:
    Sigma_y = Material[2]
    rho = Material[3]
    type_identifier = Material[1]

    D = min_D
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
                Kt = K_factors.TensionyieldFactor(w, D,Material[1])
                Kbry = K_factors.Shear_out_Factor(w, D, t)
                
                sigma_yield = Sigma_y #add yield stress of the material
                A1 = ((D/2-np.cos(np.pi / 4) * D / 2) + (w - D) / 2) * t
                A2 = (w - D) * t / 2
                A3 = A2
                A4 = A1
                Aav =  6 / (3 / A1 + 1 / A2 + 1 / A3 + 1 / A4)
                Abr = D*t
                At = w * t

                A_pin = (D/2)**2 * math.pi

                Kty = K_factors.Trans_Factor (Aav, Abr)

                #calculate axial & transverse loads (NOT incl. safety margin)
                Pu = Kt * Sigma_y * At
                Pbry = Kbry * Sigma_y * Abr
                Pty = Kty * Abr * Sigma_y
                Pmin = min(Pu, Pbry)

                #calculate Ra and Rtr

                if Pmin <= 0:
                    #print("No axial force allowed, Pmin = ", Pmin)
                    t += t_stepsize
                    continue

                Ra = P_y / Pmin
                Rtr = P_z / Pty

                #calculate safety margin
                MS = 1 / ((Ra**1.6 + Rtr**1.6)**0.625)# -1

                if MS >= 1:
                    MS = 1

                #calculate shear stress on pin (pin on 1 lug carries half of the reaction force on the lug configuration)
                tau_pin = (P_y * 0.5) / A_pin
                tau_max = 0.5*Sigma_y #Tresca failure criterion

                #check if the maximum loads are smaller than the ones we are facing, and if pin doesn't yield
                if P_y <= Pmin * MS and P_z <= Pty * MS and tau_pin < tau_max:
                    
                    #calculate meaningfull mass (i.e. only considering the circular part (tho this can be negative!!))
                    half_ring_mass = rho * 0.5 * math.pi * ((0.5*w)**2 - (0.5*D)**2) * t
                    Flange_L = w + w/2 #Venant Principle
                    Flange_mass = (w * Flange_L * t -  0.5 * t * math.pi*(0.5*D)**2) * rho * 2

                    #if new mass smaller than previous mass: store
                    if Flange_mass < min_mass:
                        min_mass = Flange_mass
                        Best_config = [Material[1],Flange_mass,t,D,w,MS]
                        
                        print("-----------------------------------------")
                        print("Better configuration for Flange found")
                        print("-----------------------------------------\n")
                        print(f"{'e/D:':<35}{w/(2*D):<12.3f}{'[none]':<}")
                        print(f"{'Kbry:':<35}{Kbry:<12.3f}{'[none]':<}")
                        print(f"{'kt: ':<35}{Kt:<12.3f}{'[none]':<}")
                        print(f"{'MS for oblique loading: ':<35}{MS:<12.3f}{'[none]':<}")
                        print(f"{'Allowed axial force: ':<35}{Pmin * MS:<12.1f}{'[N]':<}")
                        print(f"{'Allowed transversal force: ':<35}{Pty * MS:<12.1f}{'[N]':<}")
                        print(f"{'Ring mass: ':<35}{Flange_mass * 10**3:<12.2f}{'[g]':<}")
                        print(f"{'Pin stress: ':<35}{tau_pin*10**(-6):<12.1f}{'[GPa]':<}","\n")

                t += t_stepsize
            w += w_stepsize
        D += D_stepsize