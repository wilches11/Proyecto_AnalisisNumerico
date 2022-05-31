import numpy as np
from symtable import Symbol
import sympy as sym 
from numpy.linalg import inv


def ReglaF(xi,xs,tol,n,f):
    fxi = f(xi)
    fxs = f(xs)
    if fxi == 0: 
        return [xi,0,0] 
    elif fxs == 0:
         return [xs,0,0]
    elif fxi*fxs < 0: 
        xm = xi-((fxi*(xs-xi)))/(fxs-fxi)
        fxm = f(xm)
        i = 1 
        error = tol +1 
        while error > tol and fxm != 0 and i < n:
            if fxi* fxm < 0:
                xs = xm
                fxs = fxm
            else:
                xi = xm
                fxi = fxm
            xaux = xm 
            xm = xi-((fxi*(xs-xi)))/(fxs-fxi)
            fxm = f(xm)
            error = np.abs(xm-xaux)
            rel = np.abs(error/xm)
            i += 1
        if fxm == 0:
            return [xm,0,0] 
        elif error<tol:
            return [xm,error,rel]
        else:
            return ["Fracasó en "+str(n)+" iteraciones"]
    else:
        return ["Intervalo inadecaudo"]


def bisec(xi,xs,tol,n,f):
    fxi = f(xi)
    fxs = f(xs)
    if fxi == 0: 
        return [xi,0,0] 
    elif fxs == 0:
         return [xs,0,0]
    elif fxi*fxs < 0: 
        xm = (xs+xi)/2
        fxm = f(xm)
        i = 1 
        error = tol +1 
        while error > tol and fxm != 0 and i < n:
            if fxi* fxm < 0:
                xs = xm
                fxs = fxm
            else:
                xi = xm
                fxi = fxm
            xaux = xm 
            xm = (xi+xs)/2
            fxm = f(xm)
            error = np.abs(xm-xaux)
            errorRel = np.abs(error/xm)
            i += 1
        if fxm == 0:
            return xm 
        elif error<tol:
            return [xm,error,errorRel]
        else:
            return ["Fracasó en "+str(n)+" iteraciones"]
    else:
        return ["Intervalo inadecaudo"]

def secante(x0,x1,tol,n,f):
    fx0 = f(x0)
    if fx0 == 0:
        return [x0,0,0]
    else:
        fx1 = f(x1)
        i = 0 
        error = tol +1 
        den = fx1-fx0
        while error>tol and fx1 != 0 and den != 0 and i <n:
            x2 = x1 - ((fx1*(x1-x0))/den)
            error = np.abs(x2-x1)
            Er = np.abs(error/x2)
            x0 = x1
            fx0 = fx1
            x1= x2
            fx1 = f(x1)
            den = fx1-fx0
            i += 1 
        if fx1 == 0:
            return [x1,0,0] 
        elif error<tol:
             return [x1,error, Er]
        elif den == 0:
            return ["Posible raíz multiple"]
        else:
            return  [str(x1)+"Fracasó en "+str(n)+" iteraciones"]

def RaicesMul(x0,tol,n,f,derivada1,derivdad2):
    fx0 = f(x0)
    if fx0 == 0:
        return [x0,0,0]
    else:
        fx0 = f(x0)
        fpx0 = derivada1(x0)
        fppx0 = derivdad2(x0)
        i = 0 
        error = tol +1 
        den = (fpx0**2)-(fx0*fppx0)
        while error>tol and fx0 != 0 and den != 0 and i <n:
            x1 = x0 - ((fx0*fpx0)/den)           
            fx0 = f(x1)
            fpx0 = derivada1(x1)
            fppx0 = derivdad2(x1)
            den = (fpx0**2)-(fx0*fppx0)
            error = np.abs(x1-x0)
            Er = np.abs(error/x1)
            x0 = x1
            i += 1 
        if fx0 == 0:
            return [x0,0,0] 
        elif error<tol:
             return [x0,error,Er]
        else:
            return  [x0,"Fracasó en "+str(n)+" iteraciones"]

def Newton(f,x0,tol,Nmax):
    x = sym.Symbol("x")
    xant = x0
    fant = f.subs(x,xant)
    E = 1000
    cont = 0

    while E > tol and cont < Nmax :
        df = sym.diff(f,x)
        temp = df.subs(x,xant)
        xact = xant-fant/temp
        fact = f.subs(x,xact)
        E = abs(xact-xant)
        Er = abs(E/xact)
        cont = cont + 1
        xant = xact
        fant = fact
    return [float(xact),cont,float(E),float(Er)]

def gaus(a):
    return __susReg__(__eliminacion__(a))

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

def jacobi(A,b,tol,numIter):
  n = np.size(A,0)
  L = - np.tril(A, -1)
  U = - np.triu(A,1)
  D = A+L+U
  x0 = np.zeros([n,1])
  Tj = np.matmul(inv(D),(L+U))
  autoval, autovec = np.linalg.eig(Tj)
  autoval= abs(autoval)

  for lam in autoval:
    if lam >= 1:
      return ["El método no converge con los parametros ingresados."]

  C = np.matmul(inv(D),b)
  xn = (np.matmul(Tj,x0))+C
  error = np.array(abs(xn - (np.dot(Tj,xn)+C)))
  error = np.amax(error)
  iter = 0
  while ((error > tol) and (iter<numIter)):
    nuevo = np.matmul(Tj,xn)+C
    error = np.array(abs(nuevo-xn))
    error = np.amax(error)
    xn = nuevo
    iter = iter +1
  return xn.tolist()

def GaussSeidel(A,b,tol,numIter):
  n = np.size(A,0)
  L = - np.tril(A, -1)
  U = - np.triu(A,1)
  D = A+L+U
  x0 = np.zeros([n,1])
  Tg = np.matmul(inv(D-L),U)
  autoval, autovec = np.linalg.eig(Tg)
  autoval = abs(autoval)

  for lam in autoval:
    if lam >= 1:
      return ["El método no pudo converger de acuerdo a los parametros ingresados"]

  C = np.matmul(inv(D-L),b)
  xn = (np.matmul(Tg,x0))+C
  error = np.array(abs(xn - (np.dot(Tg,xn)+C)))
  error = np.amax(error)
  iter = 0
  while ((error > tol) and (iter<numIter)):
    nuevo = np.matmul(Tg,xn)+C
    error = np.array(abs(nuevo-xn))
    error = np.amax(error)
    xn = nuevo
    iter = iter +1
  return xn.tolist()

def Splins(x,y):
    A = np.array([])
    t = np.size(x)-1
    for i in range(t):
        a = list(map(lambda x : x*0,range(i*3)))
        n = x[i]
        temp = list([(n)**2,n,1])
        b = list(map(lambda x : x*0,range((t-i-1)*3)))
        a.extend(temp)
        a.extend(b)
        a = np.array(a)
        if(i==0):
            A = a
        else:
            A = np.vstack((A,a))
    for i in range(t):
        a = list(map(lambda x : x*0,range(i*3)))
        n = x[i+1]
        temp = list([(n)**2,n,1])
        b = list(map(lambda x : x*0,range((t-i-1)*3)))
        a.extend(temp)
        a.extend(b)
        a = np.array(a)
        A = np.vstack((A,a))
    for i in range(t-1):
        a = list(map(lambda x : x*0,range(i*3)))
        n = x[i]
        temp = list([(n)**2,n,1])
        temp = list([2*n,1,0,-2*n,-1,0])
        b = list(map(lambda x : x*0,range((t-i-2)*3)))
        a.extend(temp)
        a.extend(b)
        a = np.array(a)
        A = np.vstack((A,a))
    a = [1]
    b = list(map(lambda x : x*0,range(1,(t*3))))
    a.extend(b)
    a = np.array(a)
    A = np.vstack((A,a))
    B = y[:np.size(y)-1]
    B.extend(y[1:])
    temp = list(map(lambda x: x*0,range(t)))
    B.extend(temp)
    X = np.linalg.solve(A,B)
    return np.split(X,t)

def diferenciasDivididas (puntos):
  x = sym.Symbol("x")
  n = np.size(puntos,0)
  X = puntos[:,0]
  Y = puntos[:,1]
  p = 0
  tabla = np.zeros([n, n])
  tabla[:,0] = Y
  for j in range(1,n):
    for i in range(n-j):
      tabla[i][j] = (tabla[i+1][j-1] - tabla[i][j-1]) / (X[i+j]-X[i])
  b = np.array(tabla[0,:])
  mult = 1
  for i in range(n):
    mult =1
    for j in range(i):
      mult = mult * (x-X[j])
    p = p + b[i]*(mult)
  p = sym.simplify(sym.expand(p))
  return p

def lagrange(puntos):
  x = sym.Symbol("x")
  n = np.size(puntos,0)
  p = 0
  X = puntos[:,0]
  Y = puntos[:,1]
  for k in range(n):
    L=1
    for i in range(n):
      if i != k:
        L = L*((x-X[i])/(X[k]-X[i]))
    p = p +L*(Y[k])
  p = sym.simplify(sym.expand(p))
  return p

def PuntoFijo(g,x0,tol,Nmax):
    x = sym.Symbol("x")
    #Inicialización
    xant = x0
    E = 1000
    cont = 0

    #Ciclo
    while E > tol and cont < Nmax :
        xact = g.subs(x,xant)
        E = abs(xact-xant)
        Er = abs(E/xact)
        cont = cont + 1
        xant = xact
    return [xact,cont,E,Er]

def factorizacionLU(A):
    n, m = A.shape
    P = np.identity(n)
    L = np.identity(n)
    U = A.copy()
    PF = np.identity(n)
    LF = np.zeros((n,n))
    for k in range(0, n - 1):
        index = np.argmax(abs(U[k:,k]))
        index = index + k 
        if index != k:
            P = np.identity(n)
            P[[index,k],k:n] = P[[k,index],k:n]
            U[[index,k],k:n] = U[[k,index],k:n] 
            PF = np.dot(P,PF)
            LF = np.dot(P,LF)
        L = np.identity(n)
        for j in range(k+1,n):
            L[j,k] = -(U[j,k] / U[k,k])
            LF[j,k] = (U[j,k] / U[k,k])
        U = np.dot(L,U)
    np.fill_diagonal(LF, 1)
    return  LF, U

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