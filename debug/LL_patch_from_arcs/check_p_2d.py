import numpy
import matplotlib.pyplot as plt

import sys
sys.path.append('/home/bastien/GitHub/Code/Python/')
from lib_compgeom import minimal_OBB


############################################################
def OBB_to_polygon(ctr, rng, axe):
    verts = [ctr - axe[:,0]*rng[0]*(-1)**i - axe[:,1]*rng[1]*(-1)**j
             for j in range(2)
             for i in range(2)
    ]
    return numpy.asarray([verts[i] for i in [0,1,3,2,0]])
############################################################
p = numpy.loadtxt('p_2d.dat')
print 'p_2d ='
print p

p_OBB = numpy.loadtxt('p_OBB.dat')

ctr, rng, axe = minimal_OBB(p)#, tol=1e-9)
print('ctr =')
print(ctr)
print('rng =')
print(rng)
print('area =')
print(rng[0]*rng[1])
print 'axe = '
print axe




fig, ax = plt.subplots()

ax.plot(p[:,0], p[:,1], 'k*')

b1 = OBB_to_polygon(p_OBB[:,0], p_OBB[:,1], p_OBB[:,2:4])
ax.plot(b1[:,0], b1[:,1], 'r')

b2 = OBB_to_polygon(ctr, rng, axe)
ax.plot(b2[:,0], b2[:,1], 'b--')

b3 = OBB_to_polygon(
    [-0.29461423,  0.33369714],
    [1.17943359, 1.17630811],
    numpy.array(
        [[-0.56029786, -0.8282912 ],
         [ 0.8282912,  -0.56029786]]
    )
)
ax.plot(b3[:,0], b3[:,1], 'g-')

ax.set_aspect('equal')
plt.show()
