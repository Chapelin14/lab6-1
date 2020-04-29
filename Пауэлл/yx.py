import pylab
import numpy
import matplotlib.pyplot as plt

def makeData():
    x = numpy.arange(-12, 12, 0.01)
    y = numpy.arange(-16, 15, 0.01)
    xgrid, ygrid = numpy.meshgrid(x, y)
 
    zgrid = xgrid**2 + 13*xgrid +ygrid**2 + 4*ygrid + 37;
    return xgrid, ygrid, zgrid

if __name__ == '__main__':
    x, y, z = makeData()
    pylab.contour(x, y, z)
    data = numpy.loadtxt("pt.dat")
    plt.plot(data[:,0], data[:,1])
    plt.plot(data[:,0], data[:,1], 'ro')
    pylab.show()