from pywebio.input import *
from pywebio.output import *
from  metodos import bisec,ReglaF,secante,RaicesMul,Newton,gaus,gausPP,jacobi,GaussSeidel,Splins,diferenciasDivididas,lagrange,PuntoFijo,factorizacionLU,PT
import pandas as pd
import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import io
import pywebio

def restart():
    clear()
    main()

def main():


    # Opciones de metodos 
    indice = ['Bisección', 'Regla Falsa','Punto Fijo','Newton','Secante','Raices Multiples'
    ,'Gaus','Gaus con Pivoteo parcial','Gaus con Pivoteo Total','Factorización LU','Jacobi','Gauss Sediel','Splines','Diferencias Divididas',
    'Interpolación de Lagrange']


    #Despliega una barra donde se puede seleccionar alguno de los métodos diponibles
    metodo = select("Ingrese método que desea utilizar", indice)
    
    try:
        if metodo == 'Bisección':
            # info necesaria
            data = input_group("Ingrese la información de la función",[
            input("Ingrese función tipo python",name ='fun',required=True),
            input("X inicial",type=FLOAT,name ='xi',required=True),
            input("X final",type=FLOAT,name ='xf',required=True),
            input("Tolerancia",type=FLOAT,name ='tol',required=True),
            input("Numero de iteraciones",type=FLOAT,name ='n',required=True)      
                ])
            # Pregunta el tipo de error
            error =  select("Ingrese método que desea utilizar", ['Error absoluto', 'Error relativo'])
            f = eval('lambda x:'+data['fun'])
            
            #función para hallar la raíz    
            resultado = bisec(data['xi'],data['xf'],data['tol'],data['n'],f)
            
            if error == 'Error absoluto':
                # Output de los resultados
                if len(resultado) == 1:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Resultado', resultado[0]]
                    ])
                    put_button('home',onclick=restart)

                else:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  # equal to ['text', put_text('<hr/>')]
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Raíz', float(resultado[0])],
                    ['Error absoluto', float( resultado[1])]
                    ])
                    put_button('home',onclick=restart)
            else: 
                if len(resultado) == 1:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Resultado', resultado[0]]
                    ])
                    put_button('home',onclick=restart)

                else:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  # equal to ['text', put_text('<hr/>')]
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Raíz', float(resultado[0])],
                    ['Error relativo', float( resultado[2])]
                    ])
                    put_button('home',onclick=restart)

        elif metodo == 'Regla Falsa':
            data = input_group("Ingrese la información de la función",[
            input("Ingrese función tipo python",name ='fun',required=True),
            input("X inicial",type=FLOAT,name ='xi',required=True),
            input("X final",type=FLOAT,name ='xf',required=True),
            input("Tolerancia",type=FLOAT,name ='tol',required=True),
            input("Numero de iteraciones",type=FLOAT,name ='n',required=True)      
                ])

            error =  select("Ingrese método que desea utilizar", ['Error absoluto', 'Error relativo'])
            f = eval('lambda x:'+data['fun'])
            
            
            resultado = ReglaF(data['xi'],data['xf'],data['tol'],data['n'],f)
            #put_markdown(str(resultado))

            if error == 'Error absoluto':
                # Output de los resultados
                if len(resultado) == 1:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Resultado', resultado[0]]
                    ])
                    put_button('home',onclick=restart)

                else:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  # equal to ['text', put_text('<hr/>')]
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Raíz', float(resultado[0])],
                    ['Error absoluto', float( resultado[1])]
                    ])
                    put_button('home',onclick=restart)
            else: 
                if len(resultado) == 1:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Resultado', resultado[0]]
                    ])
                    put_button('home',onclick=restart)

                else:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  # equal to ['text', put_text('<hr/>')]
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Raíz', float(resultado[0])],
                    ['Error relativo', float( resultado[2])]
                    ])
                    put_button('home',onclick=restart)

        elif metodo == 'Punto Fijo':
            data = input_group("Ingrese la información de la función",[
            input("Ingrese función tipo python",name ='fun',required=True),
            input("X inicial",type=FLOAT,name ='xi',required=True),
            input("Tolerancia",type=FLOAT,name ='tol',required=True),
            input("Numero de iteraciones",type=FLOAT,name ='n',required=True)      
                ])

            error =  select("Ingrese método que desea utilizar", ['Error absoluto', 'Error relativo'])
            x = sym.Symbol("x")
            f = eval(data['fun'])   
            resultado = PuntoFijo(f,data['xi'],data['tol'],data['n'])
            if error == 'Error absoluto':            
                put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],   
                    ['x', float(resultado[0])],
                    ['iteración', resultado[1]],
                    ['Error absoluto', float(resultado[2])],
                    ])
                put_button('home',onclick=restart)
            else:
                put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],   
                    ['x', float(resultado[0])],
                    ['iteración', resultado[1]],
                    ['Error Relativo', float(resultado[3])],
                    ])
                put_button('home',onclick=restart)

        elif metodo == 'Newton':
            data = input_group("Ingrese la información de la función",[
            input("Ingrese función tipo python",name ='fun',required=True),
            input("X inicial",type=FLOAT,name ='xi',required=True),
            input("Tolerancia",type=FLOAT,name ='tol',required=True),
            input("Numero de iteraciones",type=FLOAT,name ='n',required=True)      
                ])

            error =  select("Ingrese método que desea utilizar", ['Error absoluto', 'Error relativo'])
            x = sym.Symbol("x")
            f = eval(data['fun'])   
            
            resultado = Newton(f,data['xi'],data['tol'],data['n'])
            if error == 'Error absoluto':
                put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],   
                    ['x', float(resultado[0])],
                    ['iteración', resultado[1]],
                    ['Error absoluto',  float(resultado[2])],
                    ])
                put_button('home',onclick=restart)
            else:
                put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],   
                    ['x', float(resultado[0])],
                    ['iteración', resultado[1]],
                    ['Error Relativo',  float(resultado[3])],
                    ])
                put_button('home',onclick=restart)
        
        elif metodo == 'Secante':
            # info necesaria
            data = input_group("Ingrese la información de la función",[
            input("Ingrese función tipo python",name ='fun',required=True),
            input("X inicial",type=FLOAT,name ='xi',required=True),
            input("X final",type=FLOAT,name ='xf',required=True),
            input("Tolerancia",type=FLOAT,name ='tol',required=True),
            input("Numero de iteraciones",type=FLOAT,name ='n',required=True)      
                ])
            # Pregunta el tipo de error
            error =  select("Ingrese método que desea utilizar", ['Error absoluto', 'Error relativo'])
            f = eval('lambda x:'+data['fun'])
            
            #función para hallar la raíz    
            resultado = secante(data['xi'],data['xf'],data['tol'],data['n'],f)


            if error == 'Error absoluto':
                # Output de los resultados
                if len(resultado) == 1:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Resultado', resultado[0]]
                    ])
                    put_button('home',onclick=restart)

                else:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  # equal to ['text', put_text('<hr/>')]
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Raíz', float(resultado[0])],
                    ['Error absoluto', float( resultado[1])]
                    ])
                    put_button('home',onclick=restart)
            else: 
                if len(resultado) == 1:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Resultado', resultado[0]]
                    ])
                    put_button('home',onclick=restart)

                else:
                    put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],  # equal to ['text', put_text('<hr/>')]
                    [put_html('X<sub>f</sub>'), data['xf']],  
                    ['Raíz', float(resultado[0])],
                    ['Error relativo', float( resultado[2])]
                    ])
                    put_button('home',onclick=restart)


        elif metodo == 'Raices Multiples':
            data = input_group("Ingrese la información de la función",[
            input("Ingrese función tipo python",name ='fun',required=True),
            input("X inicial",type=FLOAT,name ='xi',required=True),
            input("Tolerancia",type=FLOAT,name ='tol',required=True),
            input("Numero de iteraciones",type=FLOAT,name ='n',required=True)      
                ])

            error =  select("Ingrese método que desea utilizar", ['Error absoluto', 'Error relativo'])
            x = sym.Symbol("x")
            f = eval(data['fun'])   
            fp = sym.diff(f)
            fpp = sym.diff(fp)

            
            resultado = RaicesMul(data['xi'],data['tol'],data['n'],sym.lambdify(x,f),sym.lambdify(x,fp),sym.lambdify(x,fpp))
            
            
            if error == 'Error absoluto':
            
                put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],   
                    ['x', resultado[0]],
                    ['Error absoluto',  float(resultado[1])],
                    ])
                put_button('home',onclick=restart)
            else: 
                put_table([
                    ['Función', data['fun']],
                    [put_html('X<sub>i</sub>'), data['xi']],   
                    ['x', resultado[0]],
                    ['Error Relativo',  float(resultado[2])],
                    ])
                put_button('home',onclick=restart)

        elif metodo == 'Gaus':
            ecuaciones= []
            n = input('Número de ecuaciones',type=FLOAT,required=True)


            columnas = ['x'+str(l+1) for l in range(n)]

            for i in range(n):
                texto = "Ingrese el numero que acompaña acada variable de la ecuacion "
                ecuaciones.append(input_group(texto+str(i+1),[input('x'+str(i+1),name ='x'+str(i+1),required=True,type=FLOAT) for i in range(n)]+[input('Igual a',name='b',required=True,type=FLOAT)]))

            A = []
            b = []

            for i in ecuaciones:
                aux = []
                for j in columnas:
                    aux.append(float(i[j]))
                aux.append(float(i['b']))
                A.append(aux)





            put_text('Sistema')

            df = pd.DataFrame(A,columns=columnas+['b'],index=['ecuación '+str(m) for m in range(n)])
            put_html(df.to_html(border=0))

            put_text('Solución')

            df = pd.DataFrame(gaus(A),columns=['Resultado'],index=columnas)
            put_html(df.to_html(border=0))
            put_button('home',onclick=restart)

        elif metodo =='Gaus con Pivoteo parcial':
            ecuaciones= []
            n = input('Número de ecuaciones',type=FLOAT,required=True)


            columnas = ['x'+str(l+1) for l in range(n)]

            for i in range(n):
                texto = "Ingrese el numero que acompaña acada variable de la ecuacion "
                ecuaciones.append(input_group(texto+str(i+1),[input('x'+str(i+1),name ='x'+str(i+1),required=True,type=FLOAT) for i in range(n)]+[input('Igual a',name='b',required=True,type=FLOAT)]))

            A,b = [],[]
            
            for i in ecuaciones:
                aux = []
                for j in columnas:
                    aux.append(float(i[j]))
                aux.append(float(i['b']))
                A.append(aux)

            put_text('Sistema')

            df = pd.DataFrame(A,columns=columnas+['b'],index=['ecuación '+str(m) for m in range(n)])
            put_html(df.to_html(border=0))

            put_text('Solución')

            df = pd.DataFrame(gausPP(A),columns=['Resultado'],index=columnas[::-1])
            put_html(df.to_html(border=0))
            put_button('home',onclick=restart)
        
        elif metodo == 'Jacobi':
            ecuaciones= []
            
            data = input_group("Ingrese la información de la función",[
            input("Tolerancia",type=FLOAT,name ='tol',required=True),
            input("Numero de iteraciones",type=FLOAT,name ='n',required=True)      
                ])

            tol = data['tol']
            numIter = tol = data['n']

            n = input('Número de ecuaciones',type=FLOAT,required=True)




            columnas = ['x'+str(l+1) for l in range(n)]

            for i in range(n):
                texto = "Ingrese el numero que acompaña acada variable de la ecuacion "
                ecuaciones.append(input_group(texto+str(i+1),[input('x'+str(i+1),name ='x'+str(i+1),required=True,type=FLOAT) for i in range(n)]+[input('Igual a',name='b',required=True,type=FLOAT)]))

            A,b = [],[]
            
            for i in ecuaciones:
                aux = []
                for j in columnas:
                    aux.append(float(i[j]))
                b.append([float(i['b'])])
                A.append(aux)

            put_text('Sistema')
            df = pd.DataFrame(A,columns=columnas,index=['ecuación '+str(m) for m in range(n)])
            df['b'] = [k[0] for k in b]
            put_html(df.to_html(border=0))

            put_text('Solución')

            resultado = jacobi(np.array(A),np.array(b),tol,numIter)

            if len(resultado) == 1:
                put_text(resultado[0])
                put_button('home',onclick=restart)
            else:
                df = pd.DataFrame(resultado,columns=['Resultado'],index=columnas[::-1])
                put_html(df.to_html(border=0))
                put_button('home',onclick=restart)
        
        elif metodo == 'Gauss Sediel':
            ecuaciones= []
            
            data = input_group("Ingrese la información de la función",[
            input("Tolerancia",type=FLOAT,name ='tol',required=True),
            input("Numero de iteraciones",type=FLOAT,name ='n',required=True)      
                ])

            tol = data['tol']
            numIter = tol = data['n']

            n = input('Número de ecuaciones',type=FLOAT,required=True)




            columnas = ['x'+str(l+1) for l in range(n)]

            for i in range(n):
                texto = "Ingrese el numero que acompaña acada variable de la ecuacion "
                ecuaciones.append(input_group(texto+str(i+1),[input('x'+str(i+1),name ='x'+str(i+1),required=True,type=FLOAT) for i in range(n)]+[input('Igual a',name='b',required=True,type=FLOAT)]))

            A,b = [],[]
            
            for i in ecuaciones:
                aux = []
                for j in columnas:
                    aux.append(float(i[j]))
                b.append([float(i['b'])])
                A.append(aux)

            put_text('Sistema')
            df = pd.DataFrame(A,columns=columnas,index=['ecuación '+str(m) for m in range(n)])
            df['b'] = [k[0] for k in b]
            put_html(df.to_html(border=0))

            put_text('Solución')

            resultado = GaussSeidel(np.array(A),np.array(b),tol,numIter)

            if len(resultado) == 1:
                put_text(resultado[0])
                put_button('home',onclick=restart)
            else:
                df = pd.DataFrame(resultado,columns=['Resultado'],index=columnas[::-1])
                put_html(df.to_html(border=0))
                put_button('home',onclick=restart)

        elif metodo == 'Splines':
            data = input_group("Ingrese la información de la función",[
            input("Vector de x",help_text='[x1,x2,y3,...]',name ='x',required=True),
            input("Vector de y",help_text='[y1,y2,y3,...]',name ='y',required=True)      
                ])
            
            x = eval(data['x'])
            y = eval(data['y'])
        
            put_text('Hubo un error vuelve a intentarlo')
            put_button('home',onclick=restart)
            if len(x) != len(y):
                put_text('Error en las dimensiones,intentalo de nuevo')
                put_button('home',onclick=restart)
            else:
                k = Splins([1,2,3,4],[3,4,6,7])
                df = pd.DataFrame(k,index=[i for i in range(len(x))],columns=['a','b','c'])
                put_text('Solución')
                put_html(df.to_html(border=0))
                put_html('En la siguiente tabla se puede observar los coeficientes de cada uno de las funciones donde a es el coeficiente de x<sup>2</sup>,b el ceficiente de x y un coeficiente libre ')
                put_button('home',onclick=restart)
            
        elif metodo == 'Diferencias Divididas':
            n = input('Ingrese número de puntos',type=FLOAT)
            puntos = []
            x = []
            y = [] 
            for i in range(n):
                data = input_group("Ingrese el Punto",[
                    input("X",help_text='x',name ='x',type=FLOAT,required=True),
                    input("Y",help_text='y',name ='y',type=FLOAT,required=True)      
                    ])
                x.append(float(data['x']))
                y.append(float(data['y']))
                puntos.append([data['x'],data['y']])
            g = diferenciasDivididas(np.array(puntos))        
            gp = sym.lambdify(sym.Symbol('x'),g)
            fig, ax = plt.subplots()
            ax.plot(x,y,'or')
            x.sort()
            t = np.linspace(x[0],x[-1],15*len(x))
            ax.plot(t,gp(t),'b') 
            buf = io.BytesIO()  
            fig.savefig(buf)
            pywebio.output.put_image(buf.getvalue())

            put_text('La función interpolante es:')
            put_text(g)
            put_button('home',onclick=restart)
        
        elif metodo == 'Interpolación de Lagrange':
            n = input('Ingrese número de puntos',type=FLOAT,required=True)
            puntos = []
            x = []
            y = [] 
            for i in range(n):
                data = input_group("Ingrese el Punto",[
                    input("X",help_text='x',name ='x',type=FLOAT,required=True),
                    input("Y",help_text='y',name ='y',type=FLOAT,required=True)      
                    ])
                x.append(float(data['x']))
                y.append(float(data['y']))
                puntos.append([data['x'],data['y']])
            g = lagrange(np.array(puntos))        
            gp = sym.lambdify(sym.Symbol('x'),g)
            fig, ax = plt.subplots()
            ax.plot(x,y,'or')
            x.sort()
            t = np.linspace(x[0],x[-1],15*len(x))
            ax.plot(t,gp(t),'b') 
            buf = io.BytesIO()  
            fig.savefig(buf)
            pywebio.output.put_image(buf.getvalue())

            put_text('La función interpolante es:')
            put_text(g)
            put_button('home',onclick=restart)

        elif metodo == 'Factorización LU':
            ecuaciones= []
            n = input('Número filas',type=FLOAT,required=True)

            columnas = ['x'+str(l+1) for l in range(n)]
            

            for i in range(n):
                texto = "Ingrese el valor de cada columana para la fila"
                ecuaciones.append(input_group(texto+str(i+1),[input('x'+str(i+1),name ='x'+str(i+1),required=True,type=FLOAT) for i in range(n)]))

            A = []

            for i in ecuaciones:
                aux = []
                for j in columnas:
                    aux.append(float(i[j]))
                A.append(aux)
            L,U=factorizacionLU(np.array(A))
            put_text('Matriz inicial')
            df = pd.DataFrame(A)
            put_html(df.to_html(border=0))
            
            put_text('Matriz L')
            df = pd.DataFrame(L)
            put_html(df.to_html(border=0))

            put_text('Matriz U')
            df = pd.DataFrame(U)
            put_html(df.to_html(border=0))
            put_button('home',onclick=restart)
    
        elif metodo =='Gaus con Pivoteo Total':
            ecuaciones= []
            n = input('Número de ecuaciones',type=FLOAT,required=True)


            columnas = ['x'+str(l+1) for l in range(n)]

            for i in range(n):
                texto = "Ingrese el numero que acompaña acada variable de la ecuacion "
                ecuaciones.append(input_group(texto+str(i+1),[input('x'+str(i+1),name ='x'+str(i+1),required=True,type=FLOAT) for i in range(n)]+[input('Igual a',name='b',required=True,type=FLOAT)]))

            A,b = [],[]
            
            for i in ecuaciones:
                aux = []
                for j in columnas:
                    aux.append(float(i[j]))
                b.append(float(i['b']))
                A.append(aux)

            put_text('Sistema')

            df = pd.DataFrame(A,columns=columnas,index=['ecuación '+str(m) for m in range(n)])
            df['b'] = b
            put_html(df.to_html(border=0))

            put_text('Solución')

            df = pd.DataFrame(PT(A,b),columns=['Resultado'],index=columnas)
            put_html(df.to_html(border=0))
            put_button('home',onclick=restart)
       
    except:

        put_text('Hubo un error vuelve a intentarlo')
        put_image('https://www.ighniz.com/wp-content/uploads/2018/05/wtf-meme-face-12-512x350.jpg')
        put_button('home',onclick=restart)

# para poner la página online
# pywebio.start_server(main, port=8080, debug=True,remote_access=True)
main()

    
    