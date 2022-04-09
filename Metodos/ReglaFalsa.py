from symtable import Symbol
import sympy

x = sympy.Symbol("x")
f = x**2 - 10

# No permitir que f(a) = f(b)

def ReglaFalsa(f,a,b,tol,Nmax):

    #Inicialización
    fa = f.subs(x,a)
    fb = f.subs(x,b)
    if fa == fb:
        fa = fa + 1
    pm = (fb*a-fa*b)/(fb-fa)
    fpm = f.subs(x,pm)
    E = 1000
    cont = 1

    #Ciclo
    while E>tol and cont<Nmax :
        if fa*fb <0 :
            b = pm
        else :
            a = pm
        p0 = pm
        fa = f.subs(x,a)
        fb = f.subs(x,b)
        if fa == fb:
            fa = fa + 1
        pm = (fb*a-fa*b)/(fb-fa)
        E = abs(pm-p0)
        cont = cont + 1
    return [pm,cont,E]

print(ReglaFalsa(f,-10,10,0.1,100))