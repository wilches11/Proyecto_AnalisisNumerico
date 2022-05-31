from symtable import Symbol
import sympy

x = sympy.Symbol("x")
f = x**2 - 2

def Newton(f,x0,tol,Nmax):

    xant = x0
    fant = f.subs(x,xant)
    E = 1000
    cont = 0

    while E > tol and cont < Nmax :
        df = sympy.diff(f,x)
        temp = df.subs(x,xant)
        xact = xant-fant/temp
        fact = f.subs(x,xact)
        E = abs(xact-xant)
        Er = abs(E/xact)
        cont = cont + 1
        xant = xact
        fant = fact
    return [float(xact),cont,float(E),float(Er)]


#print(Newton(f,-10,0.001,5))
