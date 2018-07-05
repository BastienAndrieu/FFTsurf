import numpy as np
import matplotlib.pyplot as plt
import sys


# arguments:
# 0: nom_script
# 1: n1
# 2: n2
args = sys.argv
if len(args) < 2:
    n1 = 10
else:
    n1 = int( args[1] )
if len(args) < 3:
    n2 = 50
else:
    n2 = int( args[2] )


dim = 2


p1 = np.random.rand(n1,dim)

s = 0.15
p = np.empty( (0,dim), float )
for i in range(n1):
    ni = np.random.randint(n2)
    pi = s * ( 2.0*np.random.rand(ni,dim) - 1.0 )
    for j in range(dim):
        pi[:,j] = pi[:,j] + p1[i,j]
    p = np.vstack( (p,pi) )

    
fig, ax = plt.subplots()
ax.plot( p[:,0], p[:,1], '.' )
ax.set_aspect('equal', adjustable='box')
plt.show()

np.savetxt( 'xy.dat', p )
