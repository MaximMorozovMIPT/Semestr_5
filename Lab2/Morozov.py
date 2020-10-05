import numpy as np
import sys
from math import sqrt

np.set_printoptions(threshold=sys.maxsize)

####################################################################################################
# Number of unknown variables in the system of linear equations
numOfUnknown = 100

def CreateMatrix(type):
    ####################################################################################################
    # Matrix creation
    matrix = np.zeros((numOfUnknown, numOfUnknown + 1 - type))

    #First equation
    for i in range (numOfUnknown):
        matrix[0][i] = 1

    # 2-99 equations
    for i in range (1, numOfUnknown - 1):
        for j in range (numOfUnknown):
            if i == j:
                matrix[i][j] = 10
            if i == j + 1 or i == j - 1:
                matrix[i][j] = 1

    # 100 equation
    matrix[numOfUnknown - 1, numOfUnknown - 2] = 1
    matrix[numOfUnknown - 1, numOfUnknown - 1] = 1

    # equation values
    if type == 0:
        for i in range (numOfUnknown):
            matrix[i][numOfUnknown] = 100 - i
    return matrix

####################################################################################################
def PrintMatrix(matrix):
    f1 = open('matrix.txt', 'w')
    for k in range (numOfUnknown):
        count = 0
        for i in range (numOfUnknown):
            if i != numOfUnknown - 1:
                if matrix[k][i] == 0:
                    continue
                if count != 0:
                    f1.write('+ ')
                f1.write(str(int(matrix[k][i])) + '*x' + str(i + 1) + ' ')
                count += 1
            else:
                if matrix[k][i] == 0:
                    f1.write(' = ' + str(int(matrix[k][numOfUnknown])) + '\n')
                else:
                    f1.write(str(int(matrix[k][i])) + '*x' + str(i + 1))
                    f1.write(' = ' + str(int(matrix[k][numOfUnknown])) + '\n')

####################################################################################################
def GaussMethod(matrix, numOfUnknown):
    solution = np.zeros(numOfUnknown)
    # Gaussian forward
    for i in range(numOfUnknown):        
        for j in range(i + 1, numOfUnknown):
            divRation = matrix[j][i] / matrix[i][i] # division ration for every equasion
            for k in range(numOfUnknown+1):
                matrix[j][k] = matrix[j][k] - divRation * matrix[i][k] # Step for vanishing first before diag element

    ####################################################################################################
    # Gaussian reverse
    solution[numOfUnknown-1] = matrix[numOfUnknown-1][numOfUnknown] / matrix[numOfUnknown-1][numOfUnknown-1] # Solution for 
                                                                                                             # last equation
    for i in range(numOfUnknown-2, -1, -1):
        solution[i] = matrix[i][numOfUnknown]
        for j in range(i + 1, numOfUnknown):
            solution[i] = solution[i] - matrix[i][j] * solution[j] # Substituting known variables from [numOfUnknown - 1] to [i + 1]
        solution[i] = solution[i] / matrix[i][i]                   # Find i variable

    f1 = open('solutionGauss.txt', 'w')
    for i in range (numOfUnknown):
        f1.write('x' + str(i + 1) + ' = ' + str(round(solution[i], 6)) + '\n')
    return solution

####################################################################################################
def SeidelMethod(matrix, numOfUnknown):
    solution = np.zeros((2, numOfUnknown))
    i = 0
    while True:
        # Iterational step
        for j in range(0, numOfUnknown):
            if i == 0 :
                solution[i][j] = matrix[j][numOfUnknown] / matrix[j][j] # First approach
            else:
                sum = 0
                # Part of step for [j] variable when all previous are updated
                for k in range(0, numOfUnknown):
                    if k == j:
                        continue
                    else:
                        #Find actual summ at [j] step for [j + 1] var
                        if i % 2 == 1:
                            sum += solution[0][k] * matrix[j][k]
                        else:
                            sum += solution[1][k] * matrix[j][k]

                solution[i % 2][j] = (matrix[j][numOfUnknown] - sum) / matrix[j][j]

        ####################################################################################################
        # Iteration exit condition 
        # If difference in solutions of every variable between [n] and [n + 1] iteration is too small            
        if i > 0: 
            numOfExactSolutions = 0
            for j in range(0, numOfUnknown):
                diff = solution[0][j] - solution[1][j]
                if abs(diff) < 0.0001:
                    numOfExactSolutions += 1
            if numOfExactSolutions == numOfUnknown:
                f1 = open('solutionSeidel.txt', 'w')
                for l in range (numOfUnknown):
                    f1.write('x' + str(l + 1) + ' = ' + str(round(solution[i % 2][l], 6)) + '\n')
                return solution[i % 2]
        i += 1

####################################################################################################
def countNorm(matr):
    sum = np.zeros(numOfUnknown)
    for i in range (numOfUnknown):
        for j in range (numOfUnknown):
            sum[i] += matr[i][j]
    resNorm = np.amax(sum)
    return resNorm

def countMu(matr):
    return countNorm(matr) * countNorm(np.linalg.inv(matr))

####################################################################################################  
# Run methods
matrix1 = CreateMatrix(0)
matrix2 = CreateMatrix(0)
PrintMatrix(matrix1)

gauss = GaussMethod(matrix1, numOfUnknown)
seidel = SeidelMethod(matrix2, numOfUnknown)

matr = CreateMatrix(1)
f1 = open('additional.txt', 'w')
f1.write('Condition number mu = ' + str(countMu(matr)) + '\n')

residual = gauss - seidel
for i in range (numOfUnknown):
    residual[i] = abs(residual[i])
f1.write('Residual value = ' + str(np.amax(residual))) 
