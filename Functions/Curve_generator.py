import numpy as np

def Polynomial_fit(Graph_folder_directory,Curve_identifier):
    X_values = []
    Y_values = []
    with open (Graph_folder_directory +'\\' +Curve_identifier + ".txt",'r') as file:
        for l in enumerate(file):
            values = l[1]
            data_points = values.split(",")
            X_values.append(float(data_points[0]))
            Y_values.append(float(data_points[1]))
 
    Polynomial_constants = np.polyfit(X_values, Y_values, 5)
    
    return Polynomial_constants

def Value_from_poly_fit(Polynomial_constants,X):
    Y_value = 0
    for i in range (len(Polynomial_constants)):
        Y_value += Polynomial_constants[i]*X**(len(Polynomial_constants)-1-i)

    return Y_value

Function_constants = Polynomial_fit("Reader Graphs\D-14","0.6")
Y_value = Value_from_poly_fit(Function_constants,3)

print(Y_value)
