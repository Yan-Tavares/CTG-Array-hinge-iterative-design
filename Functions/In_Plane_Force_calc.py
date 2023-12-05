import math as m
from Functions import FastenerCG as F

def Fastener_loading(d2, fnum, Fcgx, Fcgz, My , x_i, z_i, x_cg, z_cg):

    area = m.pi()*(d2/2)**2

    F_inx = Fcgx / fnum 
    F_inz = Fcgz / fnum 

    F_inMy = []
    ri_list = []

    for i in range(fnum): 
        if i <= (len(x_i) - 1): 
            r = m.sqrt((x_cg-x_i[i])**2 + (z_cg-z_i[0])**2) 
            ri_list.append(round(r, 4)) 
        else:
            r = m.sqrt((x_cg-x_i[int(i-fnum/2)])**2 + (z_cg-z_i[1])**2) 
            ri_list.append(round(r, 4)) 

    denominator = sum(area * r**2 for r in ri_list) 

    for i in range(fnum): 
        F_inMy.append(round(My * area * ri_list[i] / denominator, 4)) 

    F_inMy_total = round(sum(F_inMy), 4) 

    return F_inx, F_inz, F_inMy



