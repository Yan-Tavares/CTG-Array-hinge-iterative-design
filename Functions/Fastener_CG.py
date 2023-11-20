import math as m

x_cg = 0
z_cg = 0

width = 0
h = 0
w = 0 # Width of the Lug
d_1 = 0 # Diameter of the Lug
d_2 = 2
mat = "composite"
nf = 4 # Fastener per row
nc = 2 # Column number
edge = 1.5*d_2
sum = 0

area = (m.pi) * ((d_2)**2)/4
x_i = [] 
z_i = []

for i in range(nf):
    x_i.append(1)
for i in range(nc):
   z_i.append(1)


for i in range(nf):
    if i == 0 or i == nf:
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


for i in range(nf):
   sum = sum + x_i[i]
x_cg = sum/nf

sum = 0

for i in range(nc): 
   sum = sum + z_i[i] 
z_cg = sum/nc

print(x_i, z_i, x_cg, z_cg) 