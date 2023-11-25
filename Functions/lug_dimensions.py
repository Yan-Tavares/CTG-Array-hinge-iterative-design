#---- CALCULATE DIMENSIONS of LUG ----

import numpy as np
import math
import Curve_generator as Cgen
import Reaction_force_calculator as RecF

Spacecraft_mass




def Trans_Factor (Aav, Abr):
    x = Aav/Abr
    K_ty_curve = Cgen.Polynomial_fit("Reader Graphs\D-15","3")
    K_ty = Cgen.Value_from_poly_fit(K_ty_curve,x)

    return K_ty

def Shear_out_Factor(w, D, t):
    e = (w/2)
    x = e/(D)

    Curve_choice = round(t/D,0)

    if Curve_choice > 0.6:
        Curve_choice = 0.6
    
    if Curve_choice < 0.06:
        Curve_choice = 0.06

    K_bry_curve = Cgen.Polynomial_fit("Reader Graphs\D-14", str(Curve_choice))
    K_bry= Cgen.Value_from_poly_fit(K_bry_curve,x)

    return K_bry

def TensionyieldFactor(w, D,Material):
    x = w/D
    K_t_curve = Cgen.Polynomial_fit("Reader Graphs\D-12", Material)
    K_t= Cgen.Value_from_poly_fit(K_t_curve,x)

    return K_t

#SANTIAGO'S CODE

#Determine t1, w, D1

import scipy.io

F1 = (z*L3)/(2*h1)
Pvec = np.array([0, y+F1, z])
Pmag = np.linalg.norm(Pvec)


def MS(t1, D, w):
    sigma_yield = 123 #add yield stress of the material
    A1 = ((D/2-cos(pi / 4) * D / 2) + (w - D) / 2) * t1
    A2 = (w - D) * t1 / 2
    A3 = (w - D) * t1 / 2
    A4 = ((D/2-cos(pi / 4) * D / 2) + (w - D) / 2) * t1

    #To avoid if denominator is 0
    if A1 ==0 or A2==0 or A3 == 0 or A4 ==0:
        return 0
    
    Aav =  6 / (3 / A1 + 1 / A2 + 1 / A3 + 1 / A4)
    Abr = D*t1

    #Decomposing P into the transverse load and the axial load

    axial = Pvec[1]
    transverse  = Pvec[2]
    #transverse shear factor
    kty = TransFactor(Aav, Abr)
    Pty = kty*Abr*sigma_yield
    Rtr = transverse/Pty
    #axial shearfactor
    #will be tension and also shear
    ksy = ShearyieldFactor(w, D)
    Aeff =  2*((w-D)/2+D/2)*t1
    Psy = ksy * Aeff * sigma_yield
    kt = TensionyieldFactor(w, D)
    Psyt = kt*Aeff*sigma_yield
    Ra = axial/min(Psy, Psyt)

    MS = (1/((Ra**1.6)+(Rtr**1.6))**0.625)-1
    return MS

#JUTTA'S CODE

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
