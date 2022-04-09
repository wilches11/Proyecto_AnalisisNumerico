from symtable import Symbol
import sympy

x = sympy.Symbol("x")
f = (x)**2 - 100

def Cre_o_Decre(f,x0):
    return sympy.diff(f,x).subs(x,x0) > 0

def PuntoFijo_f(f,x0,tol,Nmax):

    if Cre_o_Decre(f,x0):
        f = f + 1/(x-1)
        g = 1/f + 1
    else :
        f = -f + 1/(x-1)
        g = 1/f + 1

    #Inicialización
    xant = x0
    E = 1000
    cont = 0

    #Ciclo
    while E > tol and cont < Nmax :
        xact = g.subs(x,xant)
        E = abs(xact-xant)
        cont = cont + 1
        xant = xact
    return [xact,cont,E]

def PuntoFijo(g,x0,tol,Nmax):

    #Inicialización
    xant = x0
    E = 1000
    cont = 0

    #Ciclo
    while E > tol and cont < Nmax :
        xact = g.subs(x,xant)
        E = abs(xact-xant)
        cont = cont + 1
        xant = xact
    return [xact,cont,E]

print(PuntoFijo_f(f,-12,0.01,100))