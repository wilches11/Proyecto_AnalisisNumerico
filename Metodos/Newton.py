from symtable import Symbol
import sympy

x = sympy.Symbol("x")
f = (x-4)**2

def Newton(f,x0,tol,Nmax):

    xant = 0
    fant = f.subs(x,xant)
    E = 1000
    cont = 0

    while E > tol and cont < Nmax :
        xact = xant-fant/(sympy.diff(f,x).subs(x,xant))
        fact = f.subs(x,xant)
        E = abs(xact-xant)
        cont = cont + 1
        xant = xact
        fant = fact
    return [xact,cont,E]

print(Newton(f,-10,0.001,100))