from ast import Return
from turtle import shape
import numpy as np

def PT(A,B):
    A = np.array(A,dtype=float)
    B = np.array(B,dtype=float)
    f,c = np.shape(A)
    B = np.split(B,np.size(B))
    Af = np.append(A,B,axis=1)
    marcas = np.arange(c)
    marcas = np.split(marcas,np.size(marcas))
    for i in range(f):
        m = A[i:,i:]
        fp,cp = np.where(A == np.max(np.abs(m)))
        fp = fp[0]
        cp = cp[0]
        tmp = np.copy(Af[i])
        Af[i] = Af[fp]
        Af[fp] = tmp
        tmp = np.copy(Af[:,i])
        t1 = np.copy(marcas[i])
        marcas[i] = marcas[cp]
        marcas[cp] = t1
        Af[:,i] = Af[:,cp]
        Af[:,cp] = tmp
        pivote = Af[i,i]
        temp = Af[i]/pivote
        Af[i] = temp
        for j in range(i+1,f):
            p = Af[j][i]
            Af[j] = Af[j]-Af[i]*p
    for i in range(f-1):
        for j in range(f-1):
            p = Af[f-j-i-2][f-i-1]
            Af[f-j-2-i] = Af[f-j-2-i]-Af[f-i-1]*p
    respuesta = np.arange(c,dtype=float)
    for i in range(c):
        t = marcas[i][0]
        respuesta[marcas[i][0]] = (Af[i][c])
    return respuesta

A = [[1,2,3,4],[1,5,9,23],[4,5,4,29],[4,6,8,12]]
B = [4,10,2,11]
print(PT(A,B))