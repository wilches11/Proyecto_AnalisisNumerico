import numpy as np
import sympy as sym 

def gausPP(a):
    a = pivoteo(a,len(a),0)
    return __susReg__(__eliminacion__(a))

def pivoteo(a,n,k):
    mayor = a[k][k]
    filaMayor = k
    for s in range(k+1,n):
        if np.abs(a[s][k]) > mayor:
            mayor = np.abs(a[s][k])
            filaMayor = s
    if mayor == 0:
        return 'El sistema no tiene solución única'
    else:
        if filaMayor != k:
            a[k],a[filaMayor] =  a[filaMayor],a[k]
    return a

def __eliminacion__(a):
    n = len(a)
    for k in range(0,n-1):
        for i in range(k+1,n):
            l = a[i][k]
            p = a[k][k]
            multi = a[i][k]/a[k][k]
            for j in range(k,n+1):
                a[i][j] = a[i][j] - (multi*a[k][j])

    return a



def __susReg__(ab):
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

j = [ 
    [-7,2,-3,4,-12],
    [5,-1,14,-1,13],
    [1,9,-7,5,31],
    [-12,13,-8,-4,-32]
] 
sym.pprint(sym.Matrix((gausPP(j))) ) 

