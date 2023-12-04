import math as m
from Functions import In_Plane_Force_calc as In_plane_F

### ITERATE ON THICKNESS. ### 

# h = 0.5
# t1 = 0.003
# d2 = 0.001
# material = "metal"
# w = 0.1
# fnum = 4

# Fcgx = 4.5* 254
# Fcgz = 4.5* 496
# My = 4.5* 177

def Min_thickness_finder(h, t1, d2, material, fnum, Fcgx, Fcgz, My,W):
    pi_list = []
    sigma_list = []
    stepsize = 0.0005

    F_inxp, F_inzp, F_inMyp, fnump, d2p = In_plane_F.Fastener_loading(h, t1, d2, material, fnum, Fcgx, Fcgz, My,W) 

    for i in range(fnump): 
        pi_list.append(m.sqrt((F_inxp)**2+(F_inzp)**2+(F_inMyp[i])**2)) 

     
    t2 = 0.001
    biggest_t2 = t2

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

        if t2 > biggest_t2:
            biggest_t2 = t2

    return biggest_t2
