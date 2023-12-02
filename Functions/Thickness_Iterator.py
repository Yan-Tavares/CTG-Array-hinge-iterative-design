import In_Plane_Force_calc as B
import FastenerCG as F
import math as m

### ITERATE ON THICKNESS. ### 

pi_list = []
sigma_list = []

stepsize = 0.0005

F_inxp, F_inzp, F_inMyp, fnump, d2p = B.Bearing_check_II(0.5, 0.003, 0.001, "metal", 0.1, 4) 

for i in range(fnump): 
    pi_list.append(m.sqrt((F_inxp)**2+(F_inzp)**2+(F_inMyp[i])**2)) 

t2 = 0.001

for i in range(fnump):
    result = round(pi_list[i]/(d2p*t2), 4) 
    sigma_list.append(result)

sstrength = (386*10**6) ### Material Shear Strength for aluminum

for i in range(fnump):
    t2 = 0.001

    if sigma_list[i] >= sstrength:
        while sigma_list[i] >= sstrength:
            t2 = t2 + stepsize # stepsize
            sigma_list[i] = pi_list[i] / (d2p*t2) 

    else:
        while sigma_list[i] < sstrength and t2 >= 0: 
            if t2 < 0.001 + stepsize: 
                break

            if t2 - stepsize > 0:
                t2 = t2 - stepsize # stepsize 
            
            sigma_list[i] = pi_list[i] / (d2p*t2)

print(t2)

print(round(max(sigma_list), 4)) 
