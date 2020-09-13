import matplotlib.pyplot as plt
import numpy as np

####################################################################
# 5 functions for analysis and their derivative functions
def function_1(x0):
    return np.sin(x0 ** 2)

def der_function_1(x0):
    return 2 * x0 * np.cos(x0 ** 2)    

def function_2(x0):
    return np.cos(np.sin(x0))   

def der_function_2(x0):
    return -np.sin(np.sin(x0)) * np.cos(x0)  

def function_3(x0):
    return np.exp(np.sin(np.cos(x0)))

def der_function_3(x0):
    return -np.exp(np.sin(np.cos(x0))) * np.sin(x0) * np.cos(np.cos(x0))

def function_4(x0):
    return np.log(x0 + 3)

def der_function_4(x0):
    return 1 / (x0 + 3)

def function_5(x0):
    return np.sqrt(x0 + 3)

def der_function_5(x0):
    return 1 / (2 * np.sqrt(x0 + 3))

####################################################################
# 5 functions for derivation

def derivative_1(func, x0, h):
    return (func(x0 + h) - func(x0)) / h

def derivative_2(func, x0, h):
    return (func(x0) - func(x0 - h)) / h   

def derivative_3(func, x0, h):
    return (func(x0 + h) - func(x0 - h)) / (2 * h)

def derivative_4(func, x0, h):
    return (4 / 3) * derivative_3(func, x0, h) - (1 / 3) * derivative_3(func, x0, 2 * h)

def derivative_5(func, x0, h):
    return (3 / 2) * derivative_3(func, x0, h) - (3 / 5) * derivative_3(func, x0, 2 * h) + (1 / 10) * derivative_3(func, x0, 3 * h)

####################################################################
# data from the problem statement
x0 = 5
numOfPoints = 21
points = np.linspace(1, numOfPoints, numOfPoints)
intervals = 2 / (2 ** points)

####################################################################
# finding array of errors
def errors(an_der, comp_der, func):
    err_arr = np.array([])
    for interval in intervals:
        err_arr = np.append(err_arr, abs(an_der(x0) - comp_der(func, x0, interval)))
    return err_arr

####################################################################
# function for easy making diagram
def makeplot(der_function, function, numfoo):
    fig, ax = plt.subplots(figsize=(16, 9)) 
    ax.set_title("Error versus interval func" + numfoo)
    ax.set_xlabel("Step value")                              
    ax.set_ylabel("Error value")   
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.plot(intervals, errors(der_function, derivative_1, function), 'ro',  color="blue", label="der_1", linestyle='-')
    ax.plot(intervals, errors(der_function, derivative_2, function), 'ro',  color="red", label="der_2", linestyle='-')
    ax.plot(intervals, errors(der_function, derivative_3, function), 'ro',  color="green", label="der_3", linestyle='-')
    ax.plot(intervals, errors(der_function, derivative_4, function), 'ro',  color="gold", label="der_4", linestyle='-')
    ax.plot(intervals, errors(der_function, derivative_5, function), 'ro',  color="black", label="der_5", linestyle='-')
    ax.legend()
    fig.savefig('.vscode/Lab1/' + numfoo + '.png', dpi=100) 

####################################################################
# making diagrams
makeplot(der_function_1, function_1, '1')
makeplot(der_function_2, function_2, '2')
makeplot(der_function_3, function_3, '3')
makeplot(der_function_4, function_4, '4')
makeplot(der_function_5, function_5, '5')
