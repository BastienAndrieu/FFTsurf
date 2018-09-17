import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon

args = sys.argv
if len(args) < 2:
    iface = 1
else:
    iface = int(args[1])

strnum = format(iface, "03")
trif = np.loadtxt("uv/tri_" + strnum + ".dat", dtype=int) - 1
uvf = np.loadtxt("uv/uv_" + strnum + ".dat", dtype=float)


idstyp = np.loadtxt("idstyp_smooth.dat", dtype=int)
uv = np.loadtxt("uv_smooth.dat", dtype=float)
"""
idstyp = np.loadtxt("idstyp.dat", dtype=int)
uv = np.loadtxt("uv.dat", dtype=float)
"""

tri = np.loadtxt("tri.dat", dtype=int) - 1


edges = np.loadtxt("../../brep/edges.dat", dtype=int)


fig, ax = plt.subplots()
mrg = 5.e-2
hw = 1.0 + mrg
cl = (1.000, 0.835, 0.505)

isonface = np.zeros(len(uv), dtype=bool)
x = np.zeros((len(uv),2))


for j in range(len(trif)):
    tj = Polygon((uvf[trif[j,0],0:2], uvf[trif[j,1],0:2], uvf[trif[j,2],0:2]),
                 fc=cl, ec=cl, lw=0.5)
    ax.add_artist(tj)

#ax.triplot(uvf[:,0], uvf[:,1], trif, 'k-', lw=0.5)

for i in range(len(idstyp)):
    if idstyp[i,1] == 2 and idstyp[i,0] == iface:
        ax.plot(uv[i,0], uv[i,1], '.', markersize=5, color='k')
        #ax.text(uv[i,0], uv[i,1], str(i+1))
        isonface[i] = True
        x[i] = uv[i,0:2]
    elif idstyp[i,1] == 1:
        iedge = idstyp[i,0] - 1
        for j in range(2):
            ih = 2*iedge + j
            if edges[ih,0] == iface:
                k = (j+1)%2
                ax.plot(uv[i,2*k], uv[i,2*k+1], '.', markersize=5, color='r')
                #ax.text(uv[i,2*k], uv[i,2*k+1], str(i+1))
                isonface[i] = True
                x[i] = uv[i,2*k:2*(k+1)]
"""
for j in range(len(tri)):
    if all(isonface[tri[j]]):
        tj = Polygon((x[tri[j,0],0:2], x[tri[j,1],0:2], x[tri[j,2],0:2]),
                 fc='none', ec='b', lw=0.5)
    ax.add_artist(tj)
"""

ax.set_xlim(-hw, hw)
ax.set_ylim(-hw, hw)
ax.set_aspect('equal')
plt.grid(False)
#plt.axis('off')
#ax.tick_params(labelbottom='off', labelleft='off')

print 'uv min =', np.amin(x), ', max =', np.amax(x)


plt.show()



        

