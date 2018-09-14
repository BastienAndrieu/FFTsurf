import numpy as np
import matplotlib.pyplot as plt


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

print isegm


clr = ['b','r']
fig, ax = plt.subplots()
for i in range(2):
    print 'polyline #',i
    xy = np.loadtxt( 'xy'+str(i+1)+'.dat', dtype=float, skiprows=1 )
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
ax.set_aspect('equal', adjustable='box')
plt.show()
