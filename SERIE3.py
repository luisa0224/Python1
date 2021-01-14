#Serie 3: Fernández Chirino Luisa y Montealegre Hernández Brandon
import math
from sympy import tan, pi
from sympy import*
x = symbols('x')
init_printing(use_unicode=True)

#Evaluación de polinomio realizada con SymPy.
#Definición de funciones de métodos.
def pfalsa(f,a,b,epsilon,it):
    z=lambdify(x,f,"math")
    if z(a)*z(b)<0:
      for i in range (it):
        c=a-z(a)*(b-a)/(z(b)-z(a))
        if z(a)*z(c)<0:
          b=c
        else:
          a=c
        if abs(z(c))<epsilon:
          break
    else:
      ans=('No existe la solución en ese intervalo')
    ans="Número de iteración: "+str(i)+"\nAproximación al cero: "+str(c)+"\n"
    return ans 

def msecante(f,a,b,epsilon,it):
    z=lambdify(x,f,"math")
    if z(a)*z(b)>=0:
        print("Falla del método.")
        return None
    p0=a
    p1=b
    q0=z(a)
    q1=z(b)
    for i in range(it):
      pn=p1-q1*((p1-p0)/(q1-q0))
      if abs(pn-p1)<epsilon:
          break
      p0=p1
      p1=pn
      q0=q1
      q1=z(pn)
    ans="Número de iteración: "+str(i)+"\nAproximación al cero: "+str(pn)+"\n"
    return ans

def newton(f,a,epsilon,it):
    z=lambdify(x,f,"math")
    fp=diff(f,x)
    d=lambdify(x,fp,"math")
    p0=a
    for i in range(it):
        if abs(z(p0))<epsilon:
            break
        elif d(p0)==0:
            break
        p0=p0-(z(p0)/d(p0))
    ans="Número de iteración: "+str(i)+"\nAproximación al cero: "+str(p0)+"\n"
    return ans

def cnewton(f,a,epsilon,it):
    z=lambdify(x,f,"math")
    fp=diff(f,x)
    d=lambdify(x,fp,"math")
    p0=a
    for i in range(it):
        if abs(z(p0))<epsilon:
            break
        elif d(p0)==0:
            break
        p0=p0-(z(p0)/d(p0))
    return p0

def biseccion(f,a,b,epsilon,it):
    z=lambdify(x,f,"math")
    for i in range(it):
      c=(a+b)/2
      if z(a)*z(c)<0:
          b=c
      elif z(c)*z(b)<0:
          a=c
      elif abs(z(c))<epsilon:
          break
      else: 
          break
          print("No hay solución en el intervalo.")
    ans="Número de iteración: "+str(i)+"\nAproximación al cero: "+str(c)+"\n"
    return ans

def Muller(f,it,epsilon): 
    x=symbols('x')
    ap=f
    y=[]
    while len(y)<3:
        r1=cnewton(f=ap, a=0, epsilon=epsilon, it=it)
        y.append(r1)
        factor=(-r1+x)
        ap=(div(f=ap, g=factor))[0]
    z=lambdify(x,ap,"math")
    p0=y[0]
    p1=y[1]
    p2=y[2]
    i=2
    h1=(p1-p0) 
    h2=(p2-p1)
    d1=(z(p1)-z(p0))/h1
    d2=(z(p2)-z(p1))/h2
    dd=(d2-d1)/(h2+h1)
    while (i<=it): 
        b=d2+h2*dd
        D=((b**2-4*z(p2)*dd))**(1/2)
        if (abs(b-D)<abs(b+dd)): 
            E=b+D
        else:
            E=b-D
        h=(-2*z(p2))/E
        p=p2+h
        if (abs(h)>epsilon):
            print ("Las raíz imaginaria es: "+str(p)+"."+"\nPolinomio: "+str(ap)+"\nEl resto de las raíces: "+str(y[0:])+"\nIteraciones: "+str(i))
            break
        p0=p1
        p1=p2
        p2=p
        d1=(z(p1)-z(p0))/h1
        d2=(z(p2)-z(p1))/h2
        dd=(d1-d2)/(h2+h1)
        i+=1

print("Problema 1. ")

#1.  Método de posición falsa, secante y Newton
expr=(230*x**4)+(18*x**3)+(9*x**2)-(221*x)-9

print("Método 1. Posición falsa.\n")
#1.1 Método de posición falsa
ap1=pfalsa(f=expr,a=-1,b=0,epsilon=1e-5, it=100)
print(ap1)

ap2=pfalsa(f=expr,a=0,b=1,epsilon=1e-5, it=100)
print(ap2)

print("Método 2. Secante.\n")
#1.2 Método de la secante
ap3=msecante(f=expr,a=-1,b=0,epsilon=1e-5, it=100)
print(ap3)

ap4=msecante(f=expr,a=0,b=1,epsilon=1e-5, it=100) 
print(ap4)

print("Método 3. Newton-Raphson.\n")
#1.3 Método de Newton
ap5=newton(f=expr,a=-0.5,epsilon=1e-5, it=100)
print(ap5)

ap6=newton(f=expr,a=0.50,epsilon=1e-5, it=100) 
print(ap6) #Tomando como punto de inicio 0.5 se aproxima la primera raíz.

ap6_1=newton(f=expr,a=0.60,epsilon=1e-5, it=100) 
print(ap6_1) 

print("Problema 2.")
#2.  Método de bisección, posición falsa y secante
expr=tan(pi*x)-6

print("Valor real: 0.4474315432887465700492218303236554454681266795729210493941276720")

print("Método 1. Bisección.\n")
#2.1 Método de bisección.
ap7=biseccion(f=expr, a=0, b=0.48, epsilon=1e-5, it=10)
print(ap7)

print("Método 2. Posición falsa.\n")
#2.2 Método de posición falsa.
ap8=pfalsa(f=expr, a=0, b=0.48, epsilon=1e-5, it=10)
print(ap8)

print("Método 3. Secante.\n")
#2.3 Método de la secante
ap9=msecante(f=expr, a=0, b=0.48, epsilon=1e-5, it=10)
print(ap9)
print("Este método no converge para esta ecuación.")

print("Problema 3.")
#3. Newton y Müller
expr=(x**5)+(11*x**4)-(21*x**3)-(10*x**2)-(21*x)-5

print("Primera raíz real por método de Newton.\n")
ap10=cnewton(f=expr, a=0, epsilon=1e-5, it=100)
print(ap10)

factor1=(-ap10+x)
expr2=(div(f=expr, g=factor1))[0]
print(expr2)

print("Segunda raíz real por método de Newton.\n")
ap11=cnewton(f=expr2, a=0, epsilon=1e-5, it=100)
print(ap11)

factor2=(-ap11+x)
expr3=(div(f=expr2, g=factor2))[0]
print(expr3)

print("Tercera raíz real por método de Newton.\n")
ap12=cnewton(f=expr3, a=0, epsilon=1e-5, it=100)
print(ap12)

factor3=(-ap12+x)
expr4=(div(f=expr3, g=factor3))[0]
print(expr4)

print("Aproximación a cero complejo con método de Newton-Muller.\n")
ap13=Muller(f=expr, it=1000, epsilon=1e-5)

print("Problema 4.\n")
expr=-1+(3*(0.25-(0.09*asin(x/0.3)-(x*((0.09-x**2)**(1/2))))))

ap14=pfalsa(f=expr, a=0.3, b=-0.3, epsilon=3e-10, it=1000)
print(ap14)

print("El cero real está en -0.283357 m (Wolfram). Aproximación dentro de +/- 0.3cm")