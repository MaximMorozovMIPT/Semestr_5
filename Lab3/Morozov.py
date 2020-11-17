import numpy as np 
import math as m 

def func1(x): 
    return x * np.exp(-x**2)
def dfunc1(x):
     return (1 - 2 * x**2) * np.exp(-x**2)
def d2func1(x):  #Need it to find maximum value of func1 when dfunc1/dx = dfunc1 = 0
    return 2 * x * (2 * x**2 - 3) * np.exp(-x**2)

epsilon1 = 1e-10 #For accurate solution of equasions
epsilon2 = 1e-4  #For problem's condition

def newton_method(func, dfunc, x, eps):
    iteration = 0
    while 1:
        iteration += 1
        dx = -func(x) / dfunc(x)
        x += dx
        if abs(func(x)) < eps:
            return iteration, x

def simple_iteration(func, dfunc, x, eps):
    iteration = 0
    while 1:
        iteration += 1
        dx = -func(x) / 2
        x += dx
        if abs(func(x)) < eps:
            return iteration, x

def simple_iteration_rev(func, dfunc, x, eps):
    iteration = 0
    while 1:
        iteration += 1
        dx = func(x) / 2
        x += dx
        if abs(func(x)) < eps:
            return iteration, x

f = open('Results.txt', 'w')
f.write("Solutions of equasion y = x * exp(-x^2) by Newton and simple methods\n")
it1 , solution1 = newton_method(func1, dfunc1, 0.4, epsilon1)
f.write("Number of iterations in Newton method = " + str(it1) + " and if f(x) = 0 then x = " + str(solution1) + "\n")

it2 , solution2 = simple_iteration(func1, dfunc1, 0.4, epsilon1)
f.write("Number of iterations in simple method = " + str(it2) + " and if f(x) = 0 then x = " + str(solution2) + "\n\n")

f.write("Solution of problem in book Aristova \n")

it3 , solution3 = newton_method(dfunc1, d2func1, 0.2, epsilon1) #Find when df/dx = 0
fmax = func1(solution3)
f.write("Maximum value of f(x) = " + str(fmax) + " and x = " + str(solution3) + "\n")
def func2(x): #Func for function width at half-height
    return x * np.exp(-x**2) - fmax / 2

f.write("Find left and right values x of f(x) = fmax / 2 \n\n")
f.write("With Newton method \n")

it4 , solution4 = newton_method(func2, dfunc1, solution3 - 0.2, epsilon2)
f.write("Number of iterations = " + str(it4) + " and if f(x) = fmax / 2 then x left = " + str(solution4) + "\n")

it5 , solution5 = newton_method(func2, dfunc1, solution3 + 0.2, epsilon2)
f.write("Number of iterations = "+ str(it5) + " and if f(x) = fmax / 2 then x right = "+ str(solution5) + "\n")
f.write("Delta x = " + str(abs(solution5 - solution4)) + "\n\n")

f.write("With simple method \n")

it6 , solution6 = simple_iteration(func2, dfunc1, solution3 - 0.1, epsilon2)
f.write("Number of iterations = " + str(it6) + " and if f(x) = fmax / 2 then x left = "+ str(solution6) + "\n")

it7 , solution7 = simple_iteration_rev(func2, dfunc1, solution3 + 0.1, epsilon2)
f.write("Number of iterations = " + str(it7) + " and if f(x) = fmax / 2 then x right = "+ str(solution7) + "\n")
f.write("Delta x = " + str(abs(solution6 - solution7)) + "\n\n")
