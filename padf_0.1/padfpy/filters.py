
import numpy as np


def make_gaussian_1d( nx, rad=None, cenx=None, norm=False, power=20 ):

    if rad is None: rad = nx/2
    if cenx is None: cenx = nx/2
    radsq = float(rad)**power
  
    x = np.arange( nx ) - nx/2
    g = np.exp ( -(x**power)/float(radsq) )
    if norm == True: g *= 1./np.sum(g)
    return g


def low_pass_sharp_filter( nx, wid ):

    pow = 12
    x = np.arange( nx/2 )
    ilow = np.where( x<wid)
    lp1d = np.zeros( nx/2 )
    lp1d[ilow] = 1.0
    
    return lp1d



def high_pass_circle_filter( nx, wid ):

    pow = 12
    x = np.arange( nx/2)
    hp1d = 1.0 - np.exp( - (x/float(wid))**pow )
    
    ilow = np.where( hp1d < 1e-2*np.max(hp1d) )
    hp1d[ilow] = 0.0
    print "DEBUG <filters.py> ilow ", len(ilow[0])

    return hp1d


def gaussian_filter( nx, gwid ):

     g = make_gaussian_1d( nx, power=2 )
     return g[nx/2:]
