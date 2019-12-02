#---------------------------------------------------------------------------
# gauss.py
#
# (c) 30.1.2004 Ennes Sarradj
#
# calculation of points and weights for Gauss-Legendre numerical integration
#---------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  Imports:
#-------------------------------------------------------------------------------
from numpy import *
from numpy.linalg import eig
#from LinearAlgebra import *
#from MLab import *

Float = float
#-------------------------------------------------------------------------------
#  make a function remember its values
#  unknown author, stolen from the web
#-------------------------------------------------------------------------------
class MemoizedFunction:
    def __init__(self, function):
        self.function = function
        self.memo = {}

    def __call__(self, *args):
        try:
            return self.memo[args]
        except KeyError:
            # uncomputed value
            v = self.function(*args)
            self.memo[args] = v
            return v
        except TypeError:
            # args was not hashable; don't cache the result
            return apply(self.function, args)


def memo(function):
    return MemoizedFunction(function).__call__

#---------------------------------------------------------------------------
#  calculation of points and weights for 2D Gauss-Legendre integration
#  interval: (-1,1)
#---------------------------------------------------------------------------
def GaussPointsAndWeights(order):
    u=arange(order-1)+1
    u=u/((2*u)**2-1)**0.5
    (evalu,evec)=eig(diag(u,-1)+diag(u,1))
    return (evalu,2*evec.T[:,0]**2)
    
MGaussPointsAndWeights=memo(GaussPointsAndWeights)

#---------------------------------------------------------------------------
#  calculation of points and weights for 2D Gauss-Legendre integration
#  scaled interval: (a,b)
#---------------------------------------------------------------------------
def ScaledGauss(a,b,order):
    (points,weights)=MGaussPointsAndWeights(order)
    return (a+(b-a)*(points+1)*0.5,(b-a)*0.5*weights)

#---------------------------------------------------------------------------
#  calculation of points and weights for 2D Gauss-Legendre integration
#  with break points
#  br: ordered list of interval and break points
#  num: max total number of points
#---------------------------------------------------------------------------
def BreakGauss(br,num):
    ninterv=len(br)-1
    order=int(num/ninterv)
    ppoints=ones((ninterv,order),Float)
    wweights=ones((ninterv,order),Float)
    for i in range(ninterv):
        (ppoints[i],wweights[i])=ScaledGauss(br[i],br[i+1],order)
    xsort=argsort(ppoints.flat)
    return take(ppoints.flat,xsort),take(wweights.flat,xsort)

if __name__ == '__main__':
    print(GaussPointsAndWeights(3))
    print(ScaledGauss(-1,1,3))
    print(BreakGauss((-1,1,3),4))
