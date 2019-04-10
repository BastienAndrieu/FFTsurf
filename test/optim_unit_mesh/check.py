import numpy

import matplotlib.pyplot as plt

import sys
sys.path.append('/d/bandrieu/GitHub/Code/Python/')
import lib_halfedge as lhe

################################################
prefs = ['initial','optimized']#['initial','temporary','optimized']
fig, ax = plt.subplots(1,len(prefs))

for iax, pref in enumerate(prefs):
    # load mesh data  
    xyz = numpy.loadtxt(pref+'/xyz.dat', dtype=float)
    tri = numpy.loadtxt(pref+'/tri.dat', dtype=int) - 1
    twin = numpy.loadtxt(pref+'/twin.dat', dtype=int) - 1
    v2h = numpy.loadtxt(pref+'/v2h.dat', dtype=int) - 1
    #
    mesh = lhe.SurfaceMesh(verts=[], f2v=tri, twin=[])
    for i in range(len(xyz)):
        mesh.verts.append(lhe.Vertex(co=xyz[i], edge=v2h[i,[1,0]], index=i))

    for i in range(len(tri)):
        twins = []
        for j in range(3):
            if twin[i,2*j] > -1 and twin[i,2*j+1] > -1:
                twins.append([int(twin[i,2*j+1]), int(twin[i,2*j])])
            else:
                twins.append(None)
        mesh.twin.append(twins)
   
    lhe.plot_mesh(mesh, ax[iax],
                  faces=False,
                  halfedges=True,
                  vertices=False,
                  v2h=True,
                  v2f=False)
    ax[iax].set_aspect('equal')
    ax[iax].set_title(pref)
plt.show()

