import numpy as np
import sys
from math import sqrt

np.set_printoptions(threshold=sys.maxsize)

####################################################################################################
# Number of unknown variables in the system of linear equations
numOfUnknown = 100

####################################################################################################
# Matrix creation
matrix = np.zeros((numOfUnknown, numOfUnknown + 1))

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
for i in range (numOfUnknown):
    matrix[i][numOfUnknown] = 100 - i

np.savetxt('matrix.txt', matrix, fmt='%.0f')

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

    np.savetxt('solutionGauss.txt', solution, fmt='%.5f')

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
                np.savetxt('solutionSeidel.txt', solution[i % 2], fmt='%.5f')   
                break
        i += 1
     
# Run methods
GaussMethod(matrix, numOfUnknown)
SeidelMethod(matrix, numOfUnknown)
