import scipy.special
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as npl
import pandas as pd

def Lagrange(x,y):
    functions = []
    n = len(x)
    for i in range(0,n):
        denominator = 1
        for j in range(0,n):
            if i != j:
                denominator = denominator*(x[i] - x[j])
        newX = list(x)
        newX.pop(i)
        coefficients = multiplyPolynomial(newX)
        for j in range(0,n):
            coefficients[j] /= denominator
        functions.append(coefficients)
    finalCoefficients = []
    for i in range(0, n):
        sum = 0
        for j in range(0, n):
            sum += functions[j][i]*y[j]
        finalCoefficients.append(sum)
    return finalCoefficients

def multiplyPolynomial(x):
    n = len(x)
    if n <= 2:
        coef = []
        if n == 2:
            coef.append(1)
            coef.append(-1*x[0] + -1*x[1])
            coef.append((-1*x[0])*(-1*x[1]))
        else:
            coef.append(1)
            coef.append(x[0]*-1)
        return coef
    else:
        pol1 = multiplyPolynomial(x[0 : n // 2])
        pol2 = multiplyPolynomial(x[n // 2 : n])
        #matrix = []

        dict = {}
        for i in range(0, len(pol1)):
            row = []
            for j in range(0, len(pol2)):
                value = pol1[i]*pol2[j]
                if i+j in dict:
                    dict[i+j] += value
                else:
                    dict[i+j] = value
        coef = []
        for d in dict:
            coef.append(dict[d])
        return coef


#-----------------------------------------------

def splineInterpolation(x,y):
    matrix = []
    n = len(x)
    m = 2 + 2 + (n - 2) *4
    h = []
    z = 0
    s = 0

    for i in range(0,n):
        row = []
        row2 = []
        if i != n - 1:
            h.append( x[i + 1] - x[i])
        if i != n -1:
            for j in range(0, m):
                if j == s:
                    row.append(1)
                else:
                    row.append(0)
            matrix.append(row)
            s += 4
        if i != 0:
            for j in range(0, z):
              row2.append(0)
            k = 0
            for j in range(z, z + 4):
                row2.append(pow(h[i - 1],k))
                k += 1
            for j in range(z + 4, m):
                row2.append(0)
            matrix.append(row2)
            z += 4

    w = 0

    for i in range(0, n):
        row3 = []
        row4 = []
        if i == n - 2:
            for j in range(0,m):
                if j ==  2:
                    row3.append(1)
                else:
                    row3.append(0)
            matrix.append(row3)
        elif i == n - 1:
            for j in range(0, m):
                if j == m - 2:
                    row3.append(2)
                elif j == m - 1:
                    row3.append(h[n-2]*6)
                else:
                    row3.append(0)
            matrix.append(row3)
        else:
            for j in range(0, m):
                if j == w + 1:
                    row3.append(1)
                elif j == w + 2:
                    row3.append(2*h[i])
                elif j == w + 3:
                    row3.append(3*pow(h[i],2))
                elif j == w + 5:
                    row3.append(-1)
                else:
                    row3.append(0)
            matrix.append(row3)

            for j in range(0, m):
                if j == w + 2:
                    row4.append(2)
                elif j == w + 3:
                    row4.append(6*h[i])
                elif j == w + 6:
                    row4.append(-2)
                else:
                    row4.append(0)
            matrix.append(row4)
            w += 4
    b = []
    for i in range(0,n):
        if i == 0:
            b.append(y[i])
        elif i == n - 1:
            b.append(y[i])
        else:
            b.append(y[i])
            b.append(y[i])

    for i in range(0, (n - 2)*2 + 2):
        b.append(0)
    x = npl.solve(matrix,b)
    return x


def polynomial(coefficients,x):
    sum = 0
    n = len(coefficients)
    for i in range(0, n):
        sum += pow(x,n - i - 1)*coefficients[i]
    return sum

def polynomialSpline(coefficients,x,x2):
    sum = 0
    n = len(coefficients)
    for i in range(0, n):
        sum += pow(x - x2,i)*coefficients[i]
    return sum

data = pd.read_csv("MountEverest.csv")

#data = [[1,3,5],
 #       [6,-2,4]]

#x = splineInterpolation(data["D"],data["W"])

x = []
y = []
i = 0
for a in data["D"]:
    if i % 50 == 0:
        x.append(a)
    i += 1
i = 0
for b in data["W"]:
    if i % 50 == 0:
        y.append(b)
    i += 1

coefME = splineInterpolation(x,y)

j = 0
for i in range(0,len(x)):
    if i != len(x) - 1:
        argsME = np.arange(x[i], x[i + 1], 0.05)
        plt.plot(argsME, polynomialSpline(coefME[j:j + 4],argsME,x[i]).astype(np.double))
        j += 4
plt.title('Interpolacja drogi na Mount Everest')
plt.show()
