import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab
from RQAnalysis import*

sns.set()

def Rossler(x,y,z,a=0.2,b=0.2,c=6.3):
    x_dot = -y-z
    y_dot = x + a*y
    z_dot = b + x*z -c*z
    return(x_dot,y_dot,z_dot)

def Lorenz(X,t,sigma,beta,rho):
    x,y,z = X
    dx = -sigma*(x-y)
    dy = rho*x - y - x*z
    dz = -beta*z + x*y
    return(dx,dy,dz)

def logistic_map(x0, r, T):
    """
    Returns a time series of length T using the logistic map
    x_(n+1) = r*x_n(1-x_n) at parameter r and using the initial condition x0.

    INPUT: x0 - Initial condition, 0 <= x0 <= 1
            r - Bifurcation parameter, 0 <= r <= 4
            T - length of the desired time series
    TODO: Cythonize
    """
    #  Initialize the time series array
    timeSeries = np.empty(T)
    timeSeries[0] = x0
    for i in range(1, len(timeSeries)):
        xn = timeSeries[i-1]
        timeSeries[i] = r * xn * (1 - xn)

    return timeSeries

def Henon_Map(xn, yn, a, b):
    """
    Two-dimensional dissipative map.
    :param xn: First evolution point.
    :param yn: Second evolution point.
    :param a: First system parameter.
    :param b: Second system parameter.
    :return: Next iteration.
    """
    return yn + 1 - (a * xn**2), b*xn
