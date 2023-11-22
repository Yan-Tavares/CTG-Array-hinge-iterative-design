import math
#First need to determine out of plane loads
Fy = 169 #needs to be calculated (for now random value)
n_f = 4 #needs to be calculated (for now random value)
Mz = 150 #needs to be calculated (for now random value)
D2 = 24 #needs to be calculated (for now random value)
h = 7 #needs to be calculated (for now random value)
t1 = 10 #needs to be calculated (for now random value)
#Geometry
A = math.pi*(D2 / 2) ** 2
r = math.sqrt((1.5 * D2) ** 2 + (h / 2 + t1 + 1.5 * D2) ** 2)
rz = 1.5*D2
rx = h/2 + t1 + 1.5 * D2

# First step calculate F_pi
def calculate_Fpi( Fy, n_f ):#nf is number of fasteners
    Fpi = Fy / n_f
    return Fpi

def calculate_Fp_Mz( Mz, A , r):  #this needs fixing for the Ai and ri 
   Fp_Mz = (Mz * A * r)/( A * (r ** 2))
   return Fp_Mz

def calculate_out_of_plane_load (Fp_Mz, Fpi):
    Fyi = Fp_Mz + Fpi
    return Fyi


Fpi = calculate_Fpi(Fy, n_f)

Fp_Mz = calculate_Fp_Mz(Mz, A, r)

Fyi = calculate_out_of_plane_load(Fp_Mz, Fpi)

# Assuming multiple fasteners with same values, put them in a list
calculate_out_of_plane_load_values = [Fyi] * n_f

#now the fasteners need to be ranked by their out-of-plane load
max_Fyi = max(calculate_out_of_plane_load_values)

print("Maximum out of plane load:", max_Fyi)

# Sorting out of plane loads in decending order

sorted_calculate_out_of_plane_load = sorted(enumerate(calculate_out_of_plane_load_values, start=1),
                                            key=lambda x: x[1], reverse=True)

print("Ranking by out of plane load in descending order:")
for idx, (item, strength) in enumerate(sorted_calculate_out_of_plane_load, start=1):
    print(f"{idx}. Fastener {item}: {strength}")

