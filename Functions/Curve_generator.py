import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def Polynomial_fit(Graph_folder_directory,Curve_identifier):
    X_values = []
    Y_values = []

    with open (Graph_folder_directory +'\\' +Curve_identifier + ".txt",'r') as file:
        for l in enumerate(file):
            values = l[1]
            data_points = values.split(",")
            X_values.append(float(data_points[0]))
            Y_values.append(float(data_points[1]))

    def Second_degree_poly(x,A,B,C):
        y = A + B*x**1 + C*x**2
        return y

    Polynomial_constants,covariance = curve_fit(Second_degree_poly, X_values, Y_values)

    return Polynomial_constants


def Value_from_poly_fit(Polynomial_constants,X):
    Y_value = 0
    for i in range (len(Polynomial_constants)):
        Y_value += Polynomial_constants[i]*X**(len(Polynomial_constants))

    return Y_value



def Closest_data_point(Graph_folder_directory,Curve_identifier,X):
    with open (Graph_folder_directory +'\\' +Curve_identifier + ".txt",'r') as file:
        min_diff = X * 10**5

        for l in enumerate(file):
            values = l[1]
            data_points = values.split(",")

            if abs(X - float(data_points[0])) <= min_diff:
                min_diff = abs(X - float(data_points[0]))
                Y_best = float(data_points[1])

    return Y_best

# Function_constants = Polynomial_fit("Reader Graphs\D-14","0.6")
# Y_value = Value_from_poly_fit(Function_constants,3)

# print(Y_value)
