import numpy as np
import matplotlib.pyplot as plt
import sys

# arguments:
# 0: nom_script
# 1: dim
# 2: n1
# 3: n2

dim_def = 2
n1_def = 15
n2_def = 10


args = sys.argv
if len(args) < 2:
    dim = dim_def
else:
    dim = int(args[1])

if len(args) < 3:
    n1 = n1_def
else:
    n1 = int(args[2])

if len(args) < 4:
    n2 = n2_def
else:
    n2 = int(args[3])



p1 = np.random.rand(n1,dim)

s = 0.15
p = np.empty( (0,dim), float )
for i in range(n1):
    ni = np.random.randint(n2)
    pi = s * ( 2.0*np.random.rand(ni,dim) - 1.0 )
    for j in range(dim):
        pi[:,j] = pi[:,j] + p1[i,j]
    p = np.vstack( (p,pi) )

"""
if dim == 2:
    fig, ax = plt.subplots()
    ax.plot( p[:,0], p[:,1], '.' )
    ax.set_aspect('equal', adjustable='box')
    plt.show()
"""

hdr = str(p.shape[0])+' '+str(p.shape[1])
np.savetxt( 'points.dat', p, header=hdr, comments='' )
