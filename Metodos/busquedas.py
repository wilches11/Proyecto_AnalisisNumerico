def busquedas(x0,delta,tol,n,f):
    fx0 = f(x0)
    if fx0 == 0: 
        return x0
    else:
        x1 = x1+delta
        i = 1
        fx1 = f(x1)
        while fx0*fx1 >0 and i<n:
            x0 = x1
            f0 = fx1
            x1 = x0+delta
            fx1 = f(x1)
            i +=1
        if fx1 == 0:
            return x1
        elif fx0*fx1<0:
            return f"hay una raíz entre{x0} y {x1}"
        else:
            return "Fracasó en "+str(n)+" iteraciones"
