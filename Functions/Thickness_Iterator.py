import math as m
from Functions import Thermal_stresses_FINAL as TS

### ITERATE ON THICKNESS. ### 

h = 0.5
t1 = 0.003
d2 = 0.001
material = "metal"
W = 0.1
fnum = 4

Fcgx = 4.5* 254
Fcgz = 4.5* 496
My = 4.5* 177

def Min_thickness_finder(d2, fnum, F_inxp, F_inzp,F_inMyp, W):
    pi_list = []
    sigma_list = []
    stepsize = 0.0005

    for i in range(fnum): 
        pi_list.append(m.sqrt((F_inxp[i])**2+ (F_inzp[i])**2 + (F_inMyp[i])**2)) 

    t2 = 0.001
    

    for i in range(fnum):
        result = round(pi_list[i]/(d2*t2), 4) 
        sigma_list.append(result)

    sstrength = (240*10**6) ### Material Shear Strength for aluminum 

    Bearing_stress_list = sigma_list
    biggest_t2 = t2

    for i in range(fnum):
        t2 = 0.001 

        if sigma_list[i] >= sstrength:
            while sigma_list[i] >= sstrength:
                t2 = t2 + stepsize # stepsize
                sigma_list[i] = pi_list[i] / (d2*t2) 

        else:
            while sigma_list[i] < sstrength and t2 >= 0: 
                if t2 < 0.001 + stepsize: 
                    break

                if t2 - stepsize > 0:
                    t2 = t2 - stepsize # stepsize 
                
                sigma_list[i] = pi_list[i] / (d2*t2)

        if t2 > biggest_t2:
            biggest_t2 = t2
            Bearing_stress_list = sigma_list

    # # print(biggest_t2)
    # # print(sigma_list)
    # sigma_list2 = sigma_list
    # a = sigma_list2[0]

    # for i in range(len(sigma_list)):
    #     sigma_list[i] += TS.sigma_thermal(biggest_t2, d2) ### NEW STRESS LIST OBTAINED. 

    # for i in range(fnum):
    #     t2 = 0.001 

    #     if sigma_list[i] >= sstrength:
    #         while sigma_list[i] >= sstrength:
    #             t2 = t2 + stepsize # stepsize
    #             sigma_list[i] = pi_list[i] / (d2*t2) 

    #     else:
    #         while sigma_list[i] < sstrength and t2 >= 0: 
    #             if t2 < 0.001 + stepsize: 
    #                 break

    #             if t2 - stepsize > 0:
    #                 t2 = t2 - stepsize # stepsize 
                
    #             sigma_list[i] = pi_list[i] / (d2*t2)

    #     if t2 > biggest_t2:
    #         biggest_t2 = t2
    #         Bearing_stress = max(sigma_list)
    

    # print(biggest_t2)
    # print(sigma_list)
    b = sigma_list[0]
    # print(b-a)
    
    return biggest_t2 , Bearing_stress_list
    

# (Min_thickness_finder(h, t1, d2, material, fnum, Fcgx, Fcgz, My,W))

# For that you should start with no termal stress and find the minimal thickness.

# Then you calculate the thermal stress for that thickness and then add to the previously calculated stress, 
# then find a new thickness.
# 
# With the new thickness you calculate a new thermal stress, and then again add that termal 
# stress to the bearing stress you just had. Repeat that a sufficient amount
# of times such that the new thickness considering the new thermal stress does not change more than 5% or something