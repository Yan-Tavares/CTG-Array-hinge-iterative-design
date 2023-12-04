import math as m

def calculate_Fpi( Fy, nr ):#nr is number of fasteners
    Fpi = Fy / nr
    return Fpi

def calculate_Fp_Mz(Mz,x_cg ,z_cg,x_i, z_i,A):  #this needs fixing for the Ai and ri  
    nr = len(x_i)
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


def calculate_out_of_plane_load (Mz_List, Fpi):
    Fyi = [Mz + Fpi for Mz in Mz_List]
    return Fyi

def fastener_max_stress(D2, F_x_i, F_y_i, F_z_i):
    A = m.pi*(D2/2)**2
    #In plane stresses

    sig_y = F_y_i/(A)
    tau_yz = F_z_i/(A)
    tau_yx = F_x_i/(A)

    sigma_von_mises = m.sqrt((sig_y)**2 + 3*tau_yx**2 + 3*tau_yz**2)
    return sigma_von_mises

def push_pull_through(Df0,D2,t2,Fyi,t_3):
    A_B = m.pi*((Df0/2)**2-(D2/2)**2)
    A_BP = 2*m.pi*D2*t2
    A_SC = 2*m.pi*D2*t_3
    #Von Misses yield: yield if stress_BP/SC > tension yield
    sigma_BP = m.sqrt((Fyi/A_B)**2+3*(Fyi/A_BP)**2)*10**(-6)
    sigma_SC = m.sqrt((Fyi/A_B)**2+3*(Fyi/A_SC)**2)*10**(-6)
    
    return sigma_BP , sigma_SC