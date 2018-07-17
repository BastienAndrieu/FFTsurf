import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



n_pts = [400, 250]
nc = 10
tc = np.linspace(0.0, 1.0, nc)

noisex = 0.1
noisey = 0.5

clr = ['r','b']
fig, ax = plt.subplots()
for i in range(2):
    xc = 2.0*tc - 1.0 + noisex*(2.0*np.random.rand(nc) - 1.0)
    yc = noisey*(2.0*np.random.rand(nc) - 1.0)
    a = 2.0*np.pi*np.random.rand()
    c = np.cos(a)
    s = np.sin(a)
    pc = np.stack( (c*xc + s*yc, -s*xc + c*yc) )

    fx = interp1d(tc, pc[0,:], kind='cubic')
    fy = interp1d(tc, pc[1,:], kind='cubic')
    t = np.linspace(0.0, 1.0, n_pts[i])
    x = fx(t)
    y = fy(t)
    #plt.plot(pc[0,:], pc[1,:], 'o', color=clr[i])
    plt.plot(x, y, '.-', color=clr[i])
    np.savetxt('xy'+str(i+1)+'.dat', 
               np.transpose(np.stack((x,y))), 
               fmt='%.18e', 
               header=str(len(x))+' 2', 
               comments='')

ax.set_aspect('equal', adjustable='box')
plt.show()
