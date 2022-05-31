import numpy as np

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

print(Splins([1,2,3,4],[3,4,6,7]))