import numpy as np
import matplotlib.pyplot as plt

class polyline:
    def __init__( self, xy ):
        self.xy = xy

f = open('result.dat','r')
n = int(f.readline())
if n > 0:
    isegm = np.zeros((2,n), dtype=int)
    for i in range(n):
        isegm[:,i] = [ int(a)-1 for a in f.readline().split() ]
    lambd = np.zeros((2,n), dtype=float)
    for i in range(n):
        lambd[:,i] = [ float(a) for a in f.readline().split() ]
f.close()



clr = ['b','r']
fig, ax = plt.subplots()
ply = []

huge = 1e6
xymin = huge*np.ones(2)
xymax = -huge*np.ones(2)

for i in range(2):
    print 'polyline #',i
    xy = np.loadtxt( 'xy'+str(i+1)+'.dat', dtype=float, skiprows=1 )
    xymin[0] = min(xymin[0], np.amin(xy[:,0]))
    xymin[1] = min(xymin[1], np.amin(xy[:,1]))
    xymax[0] = max(xymin[0], np.amax(xy[:,0]))
    xymax[1] = max(xymin[1], np.amax(xy[:,1]))
    ply.append(polyline(xy))
    plt.plot(xy[:,0], xy[:,1], '.-', color=clr[i], markersize=5)
    for j in range(len(xy)):
        plt.text(xy[j,0], xy[j,1], str(j+1), color=clr[i])
    y = np.zeros((n,2))
    for j in range(n):
        z = (1.0 - lambd[i,j])*xy[isegm[i,j],:] + lambd[i,j]*xy[isegm[i,j]+1,:]
        y[j] = z
        plt.plot(z[0], z[1], 'kx', markersize=6)
    for j in range(n-1):
        for k in range(j+1,n):
            print j,k,np.sqrt(np.sum(np.power((y[j] - y[k]),2)))
"""
d = ply[0].xy[0] - ply[1].xy[3]
print np.hypot(d[0],d[1])
"""

mid = 0.5*(xymax + xymin)
rng = 0.5*(xymax - xymin)
rng = 1.05*np.maximum(rng, 0.4*np.amax(rng))

xymin = mid - rng
xymax = mid + rng

#ax.set_xlim(xymin[0], xymax[0])
#ax.set_ylim(xymin[1], xymax[1])
mrg = 1e-2
ax.set_xlim(-1.0 - mrg, 1.0 + mrg)
ax.set_ylim(-1.0 - mrg, 1.0 + mrg)
ax.set_aspect('equal', adjustable='box')
plt.show()
