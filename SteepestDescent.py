import math
import pylab
import numpy
import matplotlib.pyplot as plt

alpha = 0.0
a = 0.0 
b = 0.0
eps = 0.0001
x1 = 9.0
x2 = 10.0
i = 0

def func(x, y):
    return (-x*x-13*x-y*y-4*y-37)

def f(alp, p1, p2):
    global x1
    global x2
    return (-(x1+alp*p1)*(x1+alp*p1)-13*(x1+alp*p1)-(x2+alp*p2)*(x2+alp*p2)-4*(x2+alp*p2)-37)

def swann(alphaa, p1, p2):
    global a
    global b
    global eps
    h = 5.0
    al = 0.0
    k = 0
    if (f(alphaa-h, p1, p2) >= f(alphaa, p1, p2)) and (f(alphaa+h, p1, p2) >= f(alphaa, p1, p2)):
        a = alphaa-h
        b = alphaa+h
    elif (f(alphaa, p1, p2) >= f(alphaa-h, p1, p2)) and (f(alphaa, p1, p2) >= f(alphaa+h, p1, p2)):
        alphaa = 2*(alphaa+1)
    
    if (f(alphaa-h, p1, p2) <= f(alphaa, p1, p2)) and (f(alphaa+h, p1, p2) >= f(alphaa, p1, p2)):
        alphaa = alphaa+h
        al = alphaa+(2**k)*h
    else:
        if (f(alphaa-h, p1, p2) >= f(alphaa, p1, p2)) and (f(alphaa+h, p1, p2) <= f(alphaa, p1, p2)):
            alphaa = alphaa-h
            h = -h
            al = alphaa+(2**k)*h
        else:
            a = alphaa-h
            b = alphaa+h

    while f(al, p1, p2) >= f(alphaa, p1, p2):
        k = k + 1
        alphaa = al
        al = alphaa + (2**k)*h

    if h>0:
        a = alphaa
        b = al
    else:
        a = al
        b = alphaa

def kvad_interpol(alphaa, p1, p2):
    global a
    global b
    global eps
    h = 0.5
    alpha1 = (a+b)/2
    alpha2 = alpha1+h
    alpha3 = 0.0
    q = 0.0

    if f(alpha1, p1, p2) > f(alpha2, p1, p2):
        alpha3=alpha1+2*h
    else:
        alpha3=alpha1-h

    while math.fabs(f(alpha1, p1, p2) - f(q, p1, p2)) >= eps:
        A=(alpha2-alpha3)*(alpha3-alpha1)*(alpha1-alpha2)
        a=((alpha3-alpha2)*f(alpha1, p1, p2)+(alpha1-alpha3)*f(alpha2, p1, p2)+(alpha2-alpha1)*f(alpha3, p1, p2))/A
        b=((alpha2*alpha2-alpha3*alpha3)*f(alpha1, p1, p2)+(alpha3*alpha3-alpha1*alpha1)*f(alpha2, p1, p2)+(alpha1*alpha1-alpha2*alpha2)*f(alpha3, p1, p2))/A
        q=(-b/(2*a))

        if (f(alpha1, p1, p2) < f(alpha2, p1, p2)) and (alpha3 < q) and (q < alpha1):
            alpha2=alpha1
            alpha1=q
        else:
            if (f(alpha1, p1, p2) < f(alpha2, p1, p2)) and (alpha1 < q) and (q < alpha2):
                alpha3=alpha1
                alpha1=q
            else:
                if (f(alpha1, p1, p2) > f(alpha2, p1, p2)) and (alpha2 < q) and (q < alpha3):
                    alpha3 = alpha2
                    alpha2 = q
                else:
                    alpha1 = alpha2
                    alpha2 = q

    alphaa=q
    return alphaa

def makeData():
    x = numpy.arange(-12, 10, 0.01)
    y = numpy.arange(-16, 15, 0.01)
    xgrid, ygrid = numpy.meshgrid(x, y)
 
    zgrid = xgrid**2 + 13*xgrid +ygrid**2 + 4*ygrid + 37;
    return xgrid, ygrid, zgrid

def Spusk():
    global i
    global eps
    global x1
    global x2
    p_x1=(-2*x1-13.0)/math.sqrt((-2*x1-13)*(-2*x1-13)+(-2*x2-4)*(-2*x2-4))
    p_x2=(-2*x2-4.0)/math.sqrt((-2*x1-13)*(-2*x1-13)+(-2*x2-4)*(-2*x2-4))

    alpha = 0.0

    swann(alpha, p_x1, p_x2)
    alpha=kvad_interpol(alpha, p_x1, p_x2)

    new_x1=x1+alpha*p_x1
    new_x2=x2+alpha*p_x2
    
    if math.fabs(func(new_x1, new_x2) - func(x1, x2)) >= eps:
        x1=new_x1
        x2=new_x2
        fout = open("data.txt", 'a')
        fout.write('%s %s %s\n' % (x1, x2, func(x1, x2)))
        fout.close()
        i = i + 1
        Spusk()
    else:
        x1=new_x1
        x2=new_x2

if __name__ == '__main__':
    fout = open("data.txt", 'w')
    fout.write('%s %s %s\n' % (x1, x2, func(x1, x2)))
    fout.close()
    Spusk()
    print("Spusk: x1=", x1, " x2=", x2, " f(x1,x2)=", func(x1, x2), " In iterations ", i, "\n")
    x, y, z = makeData()
    pylab.contour(x, y, z)
    data = numpy.loadtxt("data.txt")
    plt.plot(data[:,0], data[:,1])
    plt.plot(data[:,0], data[:,1], 'ro')
    pylab.show()

