import math as m


def fastener_CG(h, t1, d2, mat, fnum): # h is the input of D4.2; t1 is the output of D4.3; d2 is the fastener diameter; mat should be either
                                   # "metal" or "composite";  fnum is the number of total fasteners

    edge = 1.5 * d2
    mat = str(mat)
    sum_x = 0
    sum_z = 0

    

    # Calculate the number of columns:
    # Assume 2 rows 

    
    #if mat == "composite": 
    #    Shear_Strength = 260*(10**6) # CARBON FIBER http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
    #elif mat == "metal":
    #    Shear_Strength = 207*(10**6) # ALUMINIUM

    #Fx = 100 # Consistent with the drawing on Overleaf, fig 1.3
    #S = m.pi*(d2/2)**2 
    #S_total = Fx/Shear_Strength 
    #fastener_num = S_total/S 

    #if fastener_num < 1:
    #    fastener_num = 1

    col = fnum/ 2 

    #while col % 2 != 0:
    #    col = col + 1
    
    area = (m.pi) * ((d2) ** 2) / 4 
    x_i = [] 
    z_i = [] 
    
    for i in range(2):
        z_i.append(1)
        
    # for i in range(m.ceil(col)):
    #     x_i.append(1)

    if col % 2 == 0:
        symm = True
    else: 
        symm = False    

    sumy = 0

    for i in range(2): #  Using 2 rows seems reasonable
        if i == 0:
            if mat == "metal":
                z_i[i] = edge + d2/2
                sumy = sumy + edge + 2*d2
            elif mat == "composite":
                z_i[i] = edge + d2/2
                sumy = sumy + edge + 4*d2
            continue

        if i == 1:
                z_i[i] = sumy + edge + d2/2
                sumy = sumy + edge + d2
                continue
        
        if i != 0 and i % 2 != 0:
            if mat == "metal": 
                z_i[i] = z_i[i - 1] + 2 * d2  # Says 2 OR 3 on BS. 
            elif mat == "composite": 
                z_i[i] = z_i[i - 1] + 4 * d2  # Says 4 OR 5 on BS. 

    sumx = 0

    if symm:
        side_col = int(col/2) 

        for i in range(int(side_col)): 
            x_i.append(1)
            if i == 0: 
                if mat == "metal": 
                    x_i[i] = edge + d2/2
                    sumx = sumx + edge + 1*d2

                elif mat == "composite":
                    x_i[i] = edge + d2/2
                    sumx = sumx + edge + 4*d2
                continue
            if i == col/2-1:
                x_i[i] = sumx + 1.5*d2
                sumx = sumx + 3.5*d2
                continue
            if mat == "metal":
                x_i[i] = x_i[i - 1] + 2 * d2  # Says 2 OR 3 on BS.
                sumx = sumx + 2*d2
            elif mat == "composite":
                x_i[i] = x_i[i - 1] + 4 * d2  # Says 4 OR 5 on BS.
                sumx = sumx + 4*d2
        sumx = sumx + 2*t1 + h

        for i in range(side_col):
            x_i.append(1)
            if i == 0:
                if mat == "metal": 
                    x_i[i+side_col] = sumx + edge + d2/2
                    sumx = sumx + edge + 1*d2

                elif mat == "composite": 
                    x_i[i + side_col] = sumx + edge + d2/2
                    sumx = sumx + edge + 4*d2
                continue
            if i == col/2-1:
                x_i[i+side_col] = sumx  + 1.5*d2
                sumx = sumx + 3.5*d2
                continue
            if mat == "metal": 
                x_i[i+side_col] = x_i[i - 1 + side_col] + 2 * d2  # Says 2 OR 3 on BS. 
                sumx = sumx + 2*d2 
            elif mat == "composite": 
                x_i[i+side_col] = x_i[i - 1 + side_col] + 4 * d2  # Says 4 OR 5 on BS. 
                sumx = sumx + d2 

    else: 
        if fnum == 3 or fnum == 4 or fnum == 7 or fnum == 8 or fnum == 11 or fnum == 12 or fnum == 15 or fnum == 16: 
            side_col1 = int(m.ceil(col/2)) ### CHANGED
            side_col2 = int(m.ceil(col/2))
            col = side_col1+side_col2
            
        elif fnum == 5 or fnum == 6 or fnum == 9 or fnum == 10 or fnum == 13 or fnum == 14: 
            side_col1 = int(m.ceil(col/2))-1 ### CHANGED
            side_col2 = int(m.ceil(col/2))
            col = side_col1+side_col2

        for i in range(int(side_col1)):
            x_i.append(1)
            if i == 0:
                if mat == "metal":
                    x_i[i] = edge + d2/2
                    sumx = sumx + edge + 1*d2

                elif mat == "composite":
                    x_i[i] = edge + d2/2
                    sumx = sumx + edge + 4*d2
                continue
            if i == side_col1 - 1 :
                x_i[i] = sumx + 1.5*d2
                sumx = sumx + 3.5*d2
                continue
            if mat == "metal":
                x_i[i] = x_i[i - 1] + 2 * d2  # Says 2 OR 3 on BS.
                sumx = sumx + 2* d2
            elif mat == "composite":
                x_i[i] = x_i[i - 1] + 4 * d2  # Says 4 OR 5 on BS.
                sumx = sumx + d2
        sumx = sumx + 2*t1 + h

        for i in range(side_col2):
            x_i.append(1)
            if i == 0:
                if mat == "metal": 
                    x_i[i+side_col1] = sumx + edge + d2/2
                    sumx = sumx + edge + 1*d2

                elif mat == "composite": 
                    x_i[i + side_col1] = sumx + edge + d2/2
                    sumx = sumx + edge + 4*d2
                continue
            if i == side_col2 - 1:
                x_i[i+side_col1] = sumx + 1.5*d2
                sumx = sumx +3.5*d2
                continue
            if mat == "metal": 
                x_i[i+side_col1] = x_i[i - 1 + side_col1] + 2 * d2  # Says 2 OR 3 on BS. 
                sumx = sumx + d2 
            elif mat == "composite": 
                x_i[i+side_col1] = x_i[i - 1 + side_col1] + 4 * d2  # Says 4 OR 5 on BS. 
                sumx = sumx + d2 


    for i in range(int(col)):
        sum_x = sum_x + x_i[i]
    x_cg = sum_x / col

    for i in range(2):
        sum_z = sum_z + z_i[i]
    z_cg = sum_z / 2
       
    for i in range(2):
        if i == 0 or i == 2:
            z_i[i] = edge + d2 
        if mat == "metal":
            z_i[i] = z_i[i-1] + 2*d2 ### Says 2 OR 3 on BS.
    
        elif mat == "composite":
            z_i[i] = z_i[i-1] + 4*d2 ### Says 4 OR 5 on BS. 

    

    #print("x_i locations for fasteners are", x_i)
    #print("z_i locations for fasteners are", z_i)
    #print("x_cg is", x_cg)
    #print("z_cg is", z_cg) 
    #print("Number of fasteners:", fnum)
    #print("Number of columns:", col)
    #print("total length:", sumx)
    return x_i, z_i, x_cg, z_cg, area, sumx ### sumx is the length of the plate
    

#fastener_CG(6, 2, 1, "metal", 7)
       # h, t1, d2, mat, fnum 