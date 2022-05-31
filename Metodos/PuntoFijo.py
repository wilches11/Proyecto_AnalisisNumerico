from symtable import Symbol
import sympy

x = sympy.Symbol("x")
f = (x)**2 - 100

def PuntoFijo(g,x0,tol,Nmax):

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

print(PuntoFijo_f(f,-12,0.01,100))
