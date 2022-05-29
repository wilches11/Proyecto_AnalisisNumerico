import numpy as np
import sympy as sym 

def eliminacion(a):
    n = len(a)
    for k in range(0,n-1):
        for i in range(k+1,n):
            l = a[i][k]
            p = a[k][k]
            multi = a[i][k]/a[k][k]
            for j in range(k,n+1):
                a[i][j] = a[i][j] - (multi*a[k][j])

    return a



def susReg(ab):
    n = len(ab)-1
    x = [0 for i in range(n)]
    k = ab[n-1][n+1]/ab[n][n]
    x.insert(n,ab[n][n+1]/ab[n][n])
    for i in range(n-1,-1,-1):
        sum = 0
        for p in range(i+1,n+1):
            k= x[p]
            o = ab[i][p]
            sum = sum + (ab[i][p]*x[p])
        x[i] =(ab[i][n+1]-sum)/ab[i][i]
    return x[::-1]

