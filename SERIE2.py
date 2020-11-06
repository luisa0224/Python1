#TAREA EXAMEN 2
#Problema 1.
#f(0.43), si f(0) = 1, f(0.25) = 1.64872, f(0.5) = 2.71828, f(0.75) = 4.48169

import numpy as np
import pandas as pd
import sympy as sp
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial

#%%
def  lgrange(x,i,xm):
    """
    Evalúa el iésimo polinomio de interpolación de Lagrange a partir de
    valores de x y el punto a estimar.
    """
    n=len(xm)-1
    y=1
    for j in range(n+1):
        if i!=j:
            y*=(x-xm[j])/(xm[i]-xm[j])
    return y

def  interpolation(x, xm, ym):
    """
    Aproxima el punto deseado (x) en función del polinomio de Lagrange. 
    """
    n=len(xm)-1
    lagrpoly=[lgrange(x,i,xm) for i in range(n+1)]
    y = np.dot(ym, lagrpoly)
    return y

def coefficients(xm, i):
    """
    Extrae coeficientes y proporciona el polinomio simplificado en cada 
    iteración.
    i es el término a simplificar (L0, L1, etc).
    """
    x = sp.Symbol('x')
    n=len(xm)-1
    xc=1
    for j in range(n+1):
        if i!=j:
            xc*=((x-xm[j])/(xm[i]-xm[j]))
    xc=sp.simplify(xc)
    return xc

def simp_coef(x,y):
    """
    Simplifica e itera sobre toda la longitud del vector x para dar los
    coeficientes individuales multiplicados por f(x).
    """
    cf=[]
    for j in range(0,len(x)):
        c=coefficients(xm=x, i=j)
        c=c*y[j]
        cf.append(c)
    return cf
#%%

#%%
x=[0,0.25,0.5,0.75]
y=[1, 1.64872, 2.71828, 4.48169]

p_1=interpolation(x=0.43, xm=x, ym=y)
p_c1=simp_coef(x, y)
poly=lagrange(x, y)#Estas dos funciones del módulo scipy y polynomial
c1=Polynomial(poly).coef#Ya dan automáticamente los valores. Se colocan como
#comprobación.

x=[0,0.25,0.5]
y=[1, 1.64872, 2.71828]

p_2=interpolation(x=0.43, xm=x, ym=y)
p_c2=simp_coef(x,y)
poly=lagrange(x, y)
c2=Polynomial(poly).coef

x=[0.25,0.5]
y=[1.64872, 2.71828]

p_3=interpolation(x=0.43, xm=x, ym=y)
p_c3=simp_coef(x,y)
poly=lagrange(x, y)
c3=Polynomial(poly).coef
#%%

#Problema 2: f(0.25), si f(−1) = 0.86199480, f(−0.5) = 0.95802009, 
#f(0) = 1.0986123, f(0.5) = 1.2943767

#%%
x=[-1,-0.5, 0, 0.5]
y=[0.86199480, 0.95802009, 1.0986123, 1.2943767]

p_4=interpolation(x=0.25, xm=x, ym=y)
p_c4=simp_coef(x, y)
poly=lagrange(x, y)
c4=Polynomial(poly).coef

x=[-0.5, 0, 0.5]
y=[0.95802009, 1.0986123, 1.2943767]

p_5=interpolation(x=0.25, xm=x, ym=y)
p_c5=simp_coef(x, y)
poly=lagrange(x, y)
c5=Polynomial(poly).coef

x=[0, 0.5]
y=[1.0986123, 1.2943767]

p_6=interpolation(x=0.25, xm=x, ym=y)
p_c6=simp_coef(x, y)
poly=lagrange(x, y)
c6=Polynomial(poly).coef
#%%

#Problema 3: Población de EU
#%%
a=[1960,1970,1980,1990,2000,2010]
b=[179323000, 203302000, 226542000, 249633000, 281422000, 308746000]

poly=lagrange(a,b)
c7=Polynomial(poly).coef

p_7=interpolation(x=1950, xm=a, ym=b)
p_8=interpolation(x=1975, xm=a, ym=b)
p_9=interpolation(x=2014, xm=a, ym=b)
p_10=interpolation(x=2020, xm=a, ym=b)

er_1=(abs(150697360-p_7)/150697360)*100 
er_2=(abs(317298000-p_9)/317298000)*100
er_3=er_1-((abs(er_1-er_2))/(2014-1950))*(1975-1950)
er_4=er_2-((abs(er_1-er_2))/(2014-1950))*(2020-2014)

a=[1970,1990,2010]
b=[203302000, 249633000, 308746000]

poly=lagrange(a,b)
c8=Polynomial(poly).coef

p_11=interpolation(x=1950, xm=a, ym=b)
p_12=interpolation(x=1975, xm=a, ym=b)
p_13=interpolation(x=2014, xm=a, ym=b)
p_14=interpolation(x=2020, xm=a, ym=b)

er_5=(abs(150697360-p_11)/150697360)*100 
er_6=(abs(317298000-p_13)/317298000)*100
er_7=er_5-((abs(er_1-er_2))/(2014-1950))*(1975-1950)
er_8=er_6-((abs(er_1-er_2))/(2014-1950))*(2020-2014)
#%%

#%%
#Problema 4: Polinomios de interpolación de Newton
def inter_newton(x, y, xi):
    #Número de datos (lognitud de vector)
    n = len(x)
    #Inicio de diferencia dividida
    fdd = [[None for x in range(n)] for x in range(n)]
    #Y en diferentes puntos
    yint = [None for x in range(n)]
    #Error para el valor
    ea = [None for x in range(n)]
    #Encontrando diferencias divididas
    for i in range(n):
        fdd[i][0] = y[i]
    for j in range(1,n):
        for i in range(n-j):
            fdd[i][j] = (fdd[i+1][j-1] - fdd[i][j-1])/(x[i+j]-x[i])
    #Imprimiendo valores
    fdd_table = pd.DataFrame(fdd)
    print(fdd_table)
    #Interpolar x
    xterm = 1
    yint[0] = fdd[0][0]
    for order in range(1, n):
        xterm = xterm * (xi - x[order-1])
        yint2 = yint[order-1] + fdd[0][order]*xterm
        ea[order-1] = yint2 - yint[order-1]
        yint[order] = yint2
    #Poner en dataframe bonito de pandas
    return map(lambda yy, ee : [yy, ee], yint, ea)
#%%

#%%
x=[0,0.25,0.5,0.75]
y=[1, 1.64872, 2.71828, 4.48169]

n1=inter_newton(x=x, y=y, xi=0.43) #imprime diferencias divididas
df=pd.DataFrame(n1, columns=['f(x)','error']) #Imprime aproximación por grado
#%%

#%%
#Problema 5: 
x=[426.690, 483.297, 497.805, 568.920]
y=[2468, 2482, 2483, 2498]

n2=inter_newton(x=x, y=y, xi=430)
df1=pd.DataFrame(n2, columns=['f(x)','error'])
n3=inter_newton(x=x, y=y, xi=520.4)
df2=pd.DataFrame(n3, columns=['f(x)','error'])
#%%


