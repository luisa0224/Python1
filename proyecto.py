import math
import sympy
import numpy as np
from pylab import*
from sympy import*
x=symbols('x')
y=symbols('y')
C=symbols('C')
D=symbols('D')
F=symbols('F')
K=symbols('K')
L=symbols('L')
M=symbols('M')
N=symbols('N')
init_printing(use_unicode=True)

#Proceso semi-automático para resolución de EDO de coeficientes indeterminados
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

def pprint(A):
    n = len(A)
    for i in range(0, n):
        line = ""
        for j in range(0, n + 1):
            line += str(A[i][j]) + "\t"
            if j == n - 1:
                line += "| "
        print(line)
    print("")


def gauss(A):
    n = len(A)
    for i in range(0, n):
        maxEl = abs(A[i][i])
        maxRow = i
        for k in range(i + 1, n):
            if abs(A[k][i]) > maxEl:
                maxEl = abs(A[k][i])
                maxRow = k
        for k in range(i, n + 1):
            tmp = A[maxRow][k]
            A[maxRow][k] = A[i][k]
            A[i][k] = tmp
        for k in range(i + 1, n):
            c = -A[k][i] / A[i][i]
            for j in range(i, n + 1):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]
    x = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        x[i] = A[i][n] / A[i][i]
        for k in range(i - 1, -1, -1):
            A[k][n] -= A[k][i] * x[i]
    return x
    pprint(A)
    x = gauss(A)
  
#Ejemplo 1: y''-1y'-2y=E**(2*x)-2*(x**2+x-1)
#Ejemplo 2: y'''-3y''+3y'-1y=2*cos(x)+(1+10*x**3)
print("Acercamiento numérico para resolución de ecuaciones diferenciales no homogéneas de segundo y tercer orden.")
p0=input("Introducir ecuación: ")
p0_1=p0.split("=")
print("\nSe iguala a 0 para obtener el polinomio característico:\n")
if len(p0_1[0])<12:
    pc=str(p0_1[0][0]+"**2"+p0_1[0][3:5]+"*y"+p0_1[0][-3:-1])
    print(pc+"=0\n")
    print("Se resuelve el polinomio ya sea factorizando:\n")
    p0_2=parse_expr(pc)
    p1=factor(p0_2)
    print(p1)
    p2=solve(p1,y)
    print(p2)
    print("\nO por método de Newton (20 iteraciones con una tolerancia de E-5):\n")
    px1=p0_2.subs(y,x)
    px=cnewton(px1,0,0.00001,20)
    px2=cnewton(px1,1,0.00001,20)
    print(px2,px)
    print("\nSe sustituye en Yh.\n")
    yh=C*E**(p2[0]*x)+D*E**(p2[1]*x)
    print("Yh=C1*e**("+str(p2[0])+"x)+C2*e**("+str(p2[1])+"x)\n")
    print("Tomamos la composición de variables independientes en las cuales estaba anteriormente igualada la ecuación homogénea y formamos Yp.\n")
    p0_3=p0_1[1].split("-"or"+",1)
    yp1=parse_expr(p0_3[0])
    yp1_i=K*x*yp1
    print("Yp1=K*x*"+str(p0_3[0])+"\n")
    print("Se aplican dos derivaciones:\n")
    yp1_1=diff(yp1_i,x)
    yp1_2=diff(yp1_1,x)
    print(yp1_1)
    print(yp1_2)
    c2=1
    c1=int(p0_1[0][3:5])
    c0=int(p0_1[0][-3:-1])
    a1=yp1_i*c0
    a2=yp1_1*c1
    a3=yp1_2*c2
    k1=a1+a2+a3
    print("Se suman los 3 arreglos multiplicados por su respectiva constante y se iguala con "+str(p0_3[0])+":\n")
    print(str(k1)+"="+str(p0_3[0]))
    k_s=k1-parse_expr(p0_3[0])
    yp1_3=float(solve(k_s,K)[0])
    print(yp1_3)
    yp_e=yp1_3*yp1
    print("Yp1="+str(yp_e)+"\n")
    print("Se construye Yp 2 con la otra parte de la ecuación: ")
    yp2_i=parse_expr("-"+p0_3[1])
    yp2=K+L*x+M*x**2
    yp2_1=diff(yp2,x)
    yp2_2=diff(yp2_1,x)
    print("\nSe realizan las derivadas: \n")
    print(yp2)
    print(yp2_1)
    print(yp2_2)
    print("\nSe suman y arreglan las ecuaciones. Se iguala con "+str(yp2_i)+"\n")
    b1=yp2*c0
    b2=yp2_1*c1
    b3=yp2_2*c2
    k2=b1+b2+b3
    e1=[0,0]
    e2=[0]
    e1.append(sympy.poly(k2.coeff(x,2)).coeffs()[0])
    e2.append(sympy.poly(k2.coeff(x,1)).coeffs()[0])
    e2.append(sympy.poly(k2.coeff(x,1)).coeffs()[1])
    e3=sympy.poly(k2.coeff(x,0)).coeffs()[0:4]
    m2=[e1,e2,e3]
    m1=[[yp2_i.coeff(x,2)],[yp2_i.coeff(x,1)],[yp2_i.coeff(x,0)]]                                           
    A = np.hstack((m2, m1))
    print("\nSe escribe como matriz aumentada igualando a "+str(yp2_i)+" :\n")
    pprint(A)
    print("Se resuelve la matriz por eliminación gaussiana:\n")
    g1=gauss(A)
    print(g1)
    yp2=yp2.subs(K,g1[0])
    yp2=yp2.subs(L,g1[1])
    yp2=yp2.subs(M,g1[2])
    print("\nQueda Yp2 simplificada: "+str(yp2)+"\n")
    yx=yh+yp_e+yp2
    print("\nSe forma la función: "+str(yx))
else:
    pc=str(p0_1[0][0]+"**3"+p0_1[0][4:6]+"*y**2"+p0_1[0][9:11]+"*y"+p0_1[0][-3:-1])
    print(pc+"=0\n")
    print("Se resuelve el polinomio ya sea factorizando:\n")
    p0_2=parse_expr(pc)
    p1=factor(p0_2)
    print(p1)
    p2=solve(p1,y)
    if len(p2)<3:
        for i in range(0,2):
            p2.append(p2[0])
    print(p2)
    print("\nO por método de Newton (100 iteraciones con una tolerancia de E-5):\n")
    px1=p0_2.subs(y,x)
    print("Primera raíz real por método de Newton:\n")
    ap10=cnewton(f=px1, a=0, epsilon=1e-5, it=100)
    print(ap10)
    factor1=(-ap10+x)
    expr2=(div(f=px1, g=factor1))[0]
    print("\nSegunda raíz real por método de Newton.\n")
    ap11=cnewton(f=expr2, a=0, epsilon=1e-5, it=100)
    print(ap11)
    factor2=(-ap11+x)
    expr3=(div(f=expr2, g=factor2))[0]
    print("\nTercera raíz real por método de Newton.\n")
    ap12=cnewton(f=expr3, a=0, epsilon=1e-5, it=100)
    print(ap12)
    factor3=(-ap12+x)
    expr4=(div(f=expr3, g=factor3))[0]
    
    print("\nSe sustituye en Yh (Yh=C1e**(l1x)+C2e**(l2x)+C3e**(l2x))\n")
    yh=E**(x)*(C+D*x+F*x**2)
    print("Al ser iguales, se factoriza: ")
    print("Yh=e**(x)*(C1*"+str(p2[0])+"+C2*"+str(p2[1])+"*x+C3*"+str(p2[2])+"*x**2)\n")
    print("Tomamos la composición de variables independientes en las cuales estaba anteriormente igualada la ecuación homogénea y formamos Yp.\n")
    p0_3=p0_1[1].split("+",1)
    yp1=parse_expr(p0_3[1])
    yp1_i=K+L*x+M*x**2+N*x**3
    print("Yp1="+str(yp1_i)+"\n")
    print("Se aplican tres derivaciones:\n")
    yp1_1=diff(yp1_i,x)
    yp1_2=diff(yp1_1,x)
    yp1_3=diff(yp1_2,x)
    print(yp1_1)
    print(yp1_2)
    print(yp1_3)
    c3=1
    c2=int(p0_1[0][4:6])
    c1=int(p0_1[0][9:11])
    c0=int(p0_1[0][-3:-1])
    a0=yp1_i*c0
    a1=yp1_1*c1
    a2=yp1_2*c2
    a3=yp1_3*c3
    k1=a0+a1+a2+a3
    
    print("Se suman los 4 arreglos multiplicados por su respectiva constante y se iguala con "+str(p0_3[1])+":\n")

    e0=[0,0,0]
    e1=[0,0]
    e2=[0]
    
    e0.append(sympy.poly(k1.coeff(x,3)).coeffs()[0])
    e1.append(sympy.poly(k1.coeff(x,2)).coeffs()[0])
    e1.append(sympy.poly(k1.coeff(x,2)).coeffs()[1])
    e2.append(sympy.poly(k1.coeff(x,1)).coeffs()[0])
    e2.append(sympy.poly(k1.coeff(x,1)).coeffs()[1])
    e2.append(sympy.poly(k1.coeff(x,1)).coeffs()[2])
    e3=sympy.poly(k1.coeff(x,0)).coeffs()[0:5]
    
    m2=[e0,e1,e2,e3]
    m1=[[yp1.coeff(x,3)],[yp1.coeff(x,2)],[yp1.coeff(x,1)],[yp1.coeff(x,0)]]
                                               
    A = np.hstack((m2, m1))
    print("\nSe escribe como matriz aumentada igualando a "+str(yp1)+" :\n")
    pprint(A)
    print("Se resuelve la matriz por eliminación gaussiana:\n")
    g1=gauss(A)
    print(g1)
    
    yp1_i=yp1_i.subs(K,g1[0])
    yp1_i=yp1_i.subs(L,g1[1])
    yp1_i=yp1_i.subs(M,g1[2])
    yp1_i=yp1_i.subs(N,g1[3])
    
    print("\nQueda Yp1 simplificada: "+str(yp1_i)+"\n")
    
    print("Tomamos la composición de variables independientes en las cuales estaba anteriormente igualada la ecuación homogénea y formamos Yp.\n")
    yp2=parse_expr(p0_3[1])
    yp2_i=K*cos(x)+L*sin(x)
    print("Yp2="+str(yp2_i)+"\n")
    print("Se aplican tres derivaciones:\n")
    yp2_1=diff(yp2_i,x)
    yp2_2=diff(yp2_1,x)
    yp2_3=diff(yp2_2,x)
    print(yp2_1)
    print(yp2_2)
    print(yp2_3)
    
    b0=yp2_i*c0
    b1=yp2_1*c1
    b2=yp2_2*c2
    b3=yp2_3*c3
    k2=b0+b1+b2+b3
    print("\nSe suman los 4 arreglos multiplicados por su respectiva constante y se iguala con "+str(p0_3[0])+":\n")
    print(str(k2)+"="+str(p0_3[0]))
    k_s=k2-parse_expr(p0_3[0])
    yp2_4=collect(k_s.expand(trig=True), [sin(x), cos(x)], evaluate=False)
    y_cos=yp2_4[cos(x)]
    y_sin=yp2_4[sin(x)]
    yp2_5=solve((y_cos,y_sin),K,L)
    print(yp2_5)
    yp2f=yp2_5[K]*cos(x)+yp2_5[L]*sin(x)
    print("Yp2="+str(yp2f)+"\n")
    
    yx1=yh+yp2f+yp1_i
    
    print("\nSe forma la función: "+str(yx1))

print("\nDeterminación de constantes. Se establecen condiciones iniciales ejemplo.\n")
print("y(0)=5\ny(1)=1-e**(2)\ny(2)=-e**(4)-1")
y1=5
x1=0
y2=1
x2=1
y3=3
x3=2
ncon=input("\nSelecciona el número de constantes a determinar (2 o 3): ")
if ncon=="2":
    p1=yx.subs(x,x1)-y1
    p2=yx.subs(x,x2)-y2
    constants=solve((p1,p2),C,D)
    eqf=yx.subs(C,constants[C])
    eqf1=eqf.subs(D,constants[D])
    print(eqf1)
else:
    p1=yx1.subs(x,x1)-y1
    p2=yx1.subs(x,x2)-y2
    p3=yx1.subs(x,x3)-y3    
    constants=solve((p1,p2,p3),C,D,F)
    eqf=yx1.subs(C,constants[C])
    eqf1=eqf.subs(D,constants[D])
    eqf2=eqf1.subs(F,constants[F])
    print(eqf2)

print(constants)





