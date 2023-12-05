import math as m
import numpy as np

def calculate_Fpi( Fy, nr ):#nr is number of fasteners
    Fpi = Fy / nr
    return Fpi

def In_plane_loading(d2, Fcgx, Fcgz, My, Mz , x_i, z_i, x_cg, z_cg,fnum):
    area = np.pi*(d2/2)**2

    F_xn = Fcgx / fnum 
    F_zn = Fcgz / fnum 

    F_inMy = []
    F_iMz = []
    ri_list = []

    for i in range(fnum): 
        if i <= (len(x_i) - 1):
            r = m.sqrt((x_cg-x_i[i])**2 + (z_cg-z_i[0])**2) 
            ri_list.append(round(r, 4)) 
        else:
            r = m.sqrt((x_cg-x_i[int(i-fnum/2)])**2 + (z_cg-z_i[1])**2) 
            ri_list.append(round(r, 4)) 

    denominator = sum(area * r**2 for r in ri_list) 

    F_xi = []
    F_zi = []

    for i in range(fnum): 
        F_inMy.append(round(My * area * ri_list[i] / denominator, 4))
        F_iMz.append(round(Mz * area * ri_list[i] / denominator, 4)) 
        F_xi.append(F_xn)
        F_zi.append(F_zn)

    return F_xi, F_zi, F_inMy , F_iMz


def calculate_Fp_Mz(Mz,D_2,x_cg ,z_cg,x_i, z_i,fnum):  #this needs fixing for the Ai and ri  
    nr = fnum
    A = np.pi*(D_2/2)**2
    Mz_List = []  
    r_list = []
    summ = 0

    for i in range(nr):
        r = m.sqrt((x_i[i]-x_cg)**2+(z_i[i-2]-z_cg)**2)
        r_list.append(r)
        summ = summ + (A*(r**2))
        
    for i in range(nr):
        Mz_List.append(Mz*A*r_list[i]/summ)
    
    return Mz_List


def calculate_out_of_plane_load (Mz_List, Fy, nr):
    Fpi = Fy / nr
    Fy_i = [Mz + Fpi for Mz in Mz_List]
    return Fy_i

def fastener_max_stress(D2, F_x_i, F_y_i, F_z_i):
    A = np.pi*(D2/2)**2
    #In plane stresses

    sig_y = max(F_y_i)/(A)
    tau_yz = max(F_z_i)/(A)
    tau_yx = max(F_x_i)/(A)

    sigma_von_mises = m.sqrt((sig_y)**2 + 3*tau_yx**2 + 3*tau_yz**2)
    return sigma_von_mises

def push_pull_through(Df0,D2,t_2,Fy_i,t_3):
    Fy_i_max = max(Fy_i)
    
    A_B = m.pi*((Df0/2)**2-(D2/2)**2)
    A_BP = 2*m.pi*D2*t_2
    A_SC = 2*m.pi*D2*t_3

    #Von Misses yield: yield if stress_BP/SC > tension yield
    sigma_BP = m.sqrt((Fy_i_max/A_B)**2 + 3*(Fy_i_max/A_BP)**2)*10**(-6)
    sigma_SC = m.sqrt((Fy_i_max/A_B)**2 + 3*(Fy_i_max/A_SC)**2)*10**(-6)
    
    return sigma_BP , sigma_SC