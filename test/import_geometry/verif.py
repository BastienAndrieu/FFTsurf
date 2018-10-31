import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

xyz = np.loadtxt('xyz.dat', dtype=float)
tri = np.loadtxt('tri.dat', dtype=int) - 1



fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot_trisurf(xyz[:,0], xyz[:,1], xyz[:,2], triangles=tri, lw=1.0, antialiased=True)

f = open('paths.dat','r')
npaths = int(f.readline())
for ipath in range(npaths):
    ihype = int(f.readline())
    nv = int(f.readline())
    verts = np.zeros(nv, dtype=int)
    for ivert in range(nv):
        verts[ivert] = int(f.readline()) - 1
    ax.plot(xyz[verts,0], xyz[verts,1], xyz[verts,2], '.-')
    #print verts+1

plt.show()
