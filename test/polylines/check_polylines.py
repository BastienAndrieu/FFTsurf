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


clr = ['r','b']
fig, ax = plt.subplots()
for i in range(2):
    xy = np.loadtxt( 'xy'+str(i+1)+'.dat', dtype=float, skiprows=1 )
    plt.plot(xy[:,0], xy[:,1], '-', color=clr[i], markersize=3)
    for j in range(n):
        z = (1.0 - lambd[i,j])*xy[isegm[i,j],:] + lambd[i,j]*xy[isegm[i,j]+1,:]
        plt.plot(z[0], z[1], 'ko', markersize=4)
ax.set_aspect('equal', adjustable='box')
plt.show()
