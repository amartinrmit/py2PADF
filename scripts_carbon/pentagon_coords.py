import numpy as np



sides = 5
dth = 2.0* np.pi / float(sides)

for i in np.arange(sides):

    th = i * dth
    print 0.25*np.cos( th ) + 0.5, 0.25*np.sin( th ) + 0.5
    
