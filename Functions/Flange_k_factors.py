from Functions import Curve_generator as Cgen

def Trans_Factor (Aav, Abr):
    x = Aav/Abr
    K_ty = Cgen.Closest_data_point("Functions\Reader Graphs\D-15","3",x)

    return K_ty

def Shear_out_Factor(w, D, t):
    e = (w/2)
    x = e/(D)

    if x >= 4:
        x = 4

    if t/D >= 0.6:
        Curve_choice = "0.60"

    if t/D >= 0.195 and t/D  <0.6:
        Curve_choice = '{:.2f}'.format(round(t/D,1))
    
    if t/D > 0.06 and t/D  < 0.195:
        Curve_choice = '{:.2f}'.format(round(t/D,2))
    
    if t/D <= 0.06:
        Curve_choice = "0.06"


    K_bry = Cgen.Closest_data_point("Functions\Reader Graphs\D-14", str(Curve_choice),x)

    if K_bry <= 0:
        K_bry = 0

    return K_bry

def TensionyieldFactor(w, D,Material):
    x = w/D
    K_t = Cgen.Closest_data_point("Functions\Reader Graphs\D-12", Material,x)

    return K_t