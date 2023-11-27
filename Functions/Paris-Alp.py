
import math as m


x_cg = 0
z_cg = 0

width = 0
h = 0
w = 0 # Width of the Lug
d_1 = 0 # Diameter of the Lug
d_2 = 2
mat = "composite"
nr = 4 # Fastener per row
nc = 2 # Column number
edge = 1.5*d_2
sum = 0

A = (m.pi) * ((d_2)**2)/4
x_i = [] 
z_i = []

for i in range(nr):
    x_i.append(1)
for i in range(nc):
   z_i.append(1)

for i in range(nr):
    if i == 0 or i == nr:
       x_i[i] = edge + d_2
    if mat == "metal":
       x_i[i] = x_i[i-1] + 2*d_2 ### Says 2 OR 3 on BS.

    elif mat == "composite":
       x_i[i] = x_i[i-1] + 4*d_2 ### Says 4 OR 5 on BS.

for i in range(nc):
    if i == 0 or i == nc:
       z_i[i] = edge + d_2
    if mat == "metal":
       z_i[i] = z_i[i-1] + 2*d_2 ### Says 2 OR 3 on BS.

    elif mat == "composite":
       z_i[i] = z_i[i-1] + 4*d_2 ### Says 4 OR 5 on BS. 

for i in range(nr):
   sum = sum + x_i[i]
x_cg = sum/nr

sum = 0

for i in range(nc): 
   sum = sum + z_i[i] 
z_cg = sum/nc

def fastener_CG():
    x_cg = 0
    z_cg = 0

    width = 0
    h = 0
    w = 0  # Width of the Lug
    d_1 = 0  # Diameter of the Lug
    d_2 = 2
    mat = "composite"
    nr = 4  # Fastener per row
    nc = 2  # Column number
    edge = 1.5 * d_2
    sum_x = 0
    sum_z = 0

    area = (m.pi) * ((d_2) ** 2) / 4
    x_i = []
    z_i = []

    for i in range(nr):
        x_i.append(1)
    for i in range(nc):
        z_i.append(1)

    for i in range(nr):
        if i == 0 or i == nr:
            x_i[i] = edge + d_2
        if mat == "metal":
            x_i[i] = x_i[i - 1] + 2 * d_2  # Says 2 OR 3 on BS.
        elif mat == "composite":
            x_i[i] = x_i[i - 1] + 4 * d_2  # Says 4 OR 5 on BS.

    for i in range(nc):
        if i == 0 or i == nc:
            z_i[i] = edge + d_2
        if mat == "metal":
            z_i[i] = z_i[i - 1] + 2 * d_2  # Says 2 OR 3 on BS.
        elif mat == "composite":
            z_i[i] = z_i[i - 1] + 4 * d_2  # Says 4 OR 5 on BS.

    for i in range(nr):
        sum_x = sum_x + x_i[i]
    x_cg = sum_x / nr

    for i in range(nc):
        sum_z = sum_z + z_i[i]
    z_cg = sum_z / nc

    return x_i, z_i, x_cg, z_cg

# Assign the return values to variables
x_i, z_i, x_cg, z_cg = fastener_CG()

#First need to determine out of plane loads
Fy = 169 #needs to be calculated (for now random value)
n_f = 4 #needs to be calculated (for now random value)
Mz = 150 #needs to be calculated (for now random value)
D2 = 24 #needs to be calculated (for now random value)
h = 7 #needs to be calculated (for now random value)
t1 = 10 #needs to be calculated (for now random value)
#Geometry


# First step calculate F_pi
def calculate_Fpi( Fy, nr ):#nr is number of fasteners
    Fpi = Fy / nr
    return Fpi

def calculate_Fp_Mz( Mz,x_cg ,z_cg,x_i, z_i,A):  #this needs fixing for the Ai and ri  
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

Fpi = calculate_Fpi(Fy, nr)
 
Mz_List = calculate_Fp_Mz(Mz,x_cg ,z_cg,x_i, z_i,A)
Fyi = calculate_out_of_plane_load(Mz_List, Fpi)

# Assuming multiple fasteners with same values, put them in a list
calculate_out_of_plane_load_values = Fyi * nr

#now the fasteners need to be ranked by their out-of-plane load
max_Fyi = max(calculate_out_of_plane_load_values)

print("Maximum out of plane load:", max_Fyi)

# Sorting out of plane loads in decending order

sorted_calculate_out_of_plane_load = sorted(enumerate(calculate_out_of_plane_load_values, start=1),
                                            key=lambda x: x[1], reverse=True)

print("Ranking by out of plane load in descending order:")
for idx, (item, strength) in enumerate(sorted_calculate_out_of_plane_load, start=1):
    print(f"{idx}. Fastener {item}: {strength}")
