import numpy as np
import matplotlib.pyplot as plt

noise = 1.e-15
a = 2.0 * np.pi * np.random.rand()

n = 1000
r = np.random.rand(n)

p = np.zeros( (n,2) )
p[:,0] = r*np.cos(a) + np.random.rand() + noise*0.5*( 2.0*np.random.rand(n) - 1.0 )
p[:,1] = r*np.sin(a) + np.random.rand() + noise*0.5*( 2.0*np.random.rand(n) - 1.0 )
    
fig, ax = plt.subplots()
ax.plot( p[:,0], p[:,1], '.' )
ax.set_aspect('equal', adjustable='box')
plt.show()

np.savetxt( 'xy.dat', p )
