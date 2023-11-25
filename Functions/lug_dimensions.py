#---- CALCULATE DIMENSIONS of LUG ----

import numpy as np
import math

#define functions for Kt, Kty, Kbry

#define functions for At, Abr, Aav

#set constants
material = [Aluminium, Steel]
sigma_y = [30000000,250000000] #yield strengths of aluminium and steel, [Pa], min value found
rho = [2900,7900] #densities of aluminium and steel, [kg/m³], max value found
Pax =  / 2    #[N], load P/2 taken by each of the two lugs
Ptransv = / 2 #[N], load P/2 taken by each of the two lugs

#set initial values for t, D, w
t = 0.00001 #[m]
D = 0.001   #[m]
w = 0.001   #[m]
MS = 1

#iterate over materials: 0 = Aluminium, 1 = Steel
for material in range(0,1):
    #iterate over D
    while D < 0.9: #in meter
        D = D + 0.001
        #iterate over w
        while w < 1:
            w = w + 0.001
            while t < 0.1:
                t = t + 0.00001
                #find K's stuff
                Kt =
                Kbry =
                Kty =
                #calculate A's
                At =
                Abr =
                Aav = 
                #calculate axial & transverse loads (incl. safety margin)
                Pu = Kt * MS * sigma_y[material] * At
                Pbry = Kbry * MS * sigma_y[material] * Abr
                Pty = Kty * Abr * MS * sigma_y[material]
                Pmin = min(Pu, Pbry)
                #calculate Ra and Rtr
                Ra = Pax / Pmin
                Rtr = Ptransv / Pty
                #calculate safety margin
                MS = 1 / ((Ra**1.6 + Rtr**1.6)**0.625) -1
                #calculate meaningfull mass (i.e. only considering the circular part
                #of the lug, since the rectangular part is not being sized
                mass = rho[material] * 0.5 * math.pi * ((0.5*w)**2 - (0.5*d)**2)) * t
                #if new mass smaller than previous mass: store
                if mass < minmass:
                    minmass = mass
                    topt = t
                    Dopt = D
                    wopt = w
                    materialopt = material
                    MSopt = MS
                
#give final choice for material, D, w, t and resulting mass
print("material: ",materialopt)
print("diameter: D = ",Dopt*0.001,"mm; thickness: t = ",topt*0.001, "mm; width: w = ",wopt*0.001,"mm")
print("used safety factor: MS = ",MSopt)
print("resulting in a mass of: m = ",minmass,"kg")
