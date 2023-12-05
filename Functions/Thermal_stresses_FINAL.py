import numpy as np
from math import *
def sigma_thermal(t2, d2):
    class variables:
        E_fastener = 190*10**9 #Steel 8630
        E_wall = 68.9*10**9 #Aluminium T-6
        Tcoeff_fastener = 17.4*10**-6 #steel 8630 
        Tcoeff_wall = 23.6 *10**-6 #aluminium t-6 
        Tcoeff_lug = 13 *10**-6 # steel 4130
        E_lug = 205*10**9

    var = variables()

    class hole:
        D2 = d2  #diameter of the fastener

    hole = hole()

    class temp:
        min_lug = -139.5 
        max_lug = 158.8
        min_wall= -139.5 
        max_wall = 158.8 
        tref = 288.15

    temp = temp()

    class tempdiffwall: #temperature difference in the spacecraft wall
        minwall = temp.min_wall-temp.tref
        maxwall = temp.max_wall -temp.tref

    class tempdifflug: #temperature difference in the lug
        minlug = temp.min_lug-temp.tref
        maxlug = temp.max_lug -temp.tref

    diff_wall = tempdiffwall()
    diff_lug = tempdifflug()
    

    diff_lug_vec = [diff_lug.maxlug, diff_lug.minlug]
    diff_wall_vec = [diff_wall.maxwall, diff_wall.maxwall]

    
    D_fi = hole.D2 #ad D2 here

    D_f0 = hole.D2 * 1.2   #add the diameter of the rivet here
   
    #Wall geometry and characteristics
    
    t_wall = 3*10**-2
    E_wall = var.E_wall
    delta_a1 = (4*t_wall)/(E_wall*pi*(D_f0**2-D_fi**2))

    t_lug = t2
    E_lug = var.E_lug
    delta_a2  =(4*t_lug)/(E_lug*pi*(D_f0**2-D_fi**2))

     #using table 7.1 
    L_head = 0.4 * D_fi
    L_engaged = 0.4* D_fi
    L_sha = t_wall+t_lug
    L_n = 0.4*D_fi

    A_thread = pi * (D_fi**2)
    #Using equation 7.5.5 to determine the compliance of the fastener in order to determine the force ratio
    delta_fastener = (1/var.E_fastener) * ((L_head/A_thread)+(L_engaged/A_thread)) + (L_sha/(A_thread*var.E_fastener))+(L_n/(var.E_fastener*A_thread))

    delta_b = delta_fastener

    #phi ratios for the lug with the fastener and wall with fastener

    phi_wall = delta_a1/(delta_a1 + delta_b)
    
    phi_lug = delta_a2/(delta_a2 + delta_b)
    
    #stiffness area

    A_stiff = 0.25*pi*D_fi**2

    #there are going to be two thermal expansions one will be the wall and the other one the lug

    thermal_wall = [] #these are the thermal stresses at the wall
    thermal_lug = [] #these are the thermal stresses at the lug

    for difftemp in diff_wall_vec:
        expansion_wall = var.E_fastener * (-var.Tcoeff_fastener + var.Tcoeff_wall)* difftemp* A_stiff*(1-phi_wall)
        #print(expansion_wall)
        if expansion_wall<0:
            expansion_wall = 0
        
        thermal_wall.append(expansion_wall/(hole.D2*t_wall)) #thermal stress wall

    for difftemp in diff_lug_vec:
        expansion_lug = var.E_fastener * (-var.Tcoeff_fastener + var.Tcoeff_lug)* difftemp* A_stiff*(1-phi_wall)
        if expansion_lug<0:
            expansion_lug = 0
        thermal_lug.append(expansion_lug/(hole.D2*t_lug)) #thermal stress lug


    thermal_stresswall = max(thermal_wall) 
    thermal_stresslug = max(thermal_lug) 
    print(thermal_stresswall*10**-6) 
    print(thermal_stresslug*10**-6) 
    print(thermal_wall) 
    return thermal_stresslug #, thermal_stresswall 
        
                            

