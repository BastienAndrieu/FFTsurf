import numpy as np
import matplotlib.pyplot as plt
from random import randint
import sys
sys.path.append('/stck/bandrieu/Bureau/Python/mylibs/')
#sys.path.append('/home/bastien/Bureau/Python/')
import my_colors as myc


def plot_face(ax, outer, inner, xy, he, fcl, ecl):
    outerxy = []
    e = outer
    while True:
        outerxy.append([xy[he[e,4],0], xy[he[e,4],1]])
        e = he[e,3]
        if e == outer:
            break
    p = np.array(outerxy)
    ax.fill(p[:,0], p[:,1], ec=ecl, fc=fcl)
    return


######## Read DCEL data ########
f = open('data.dat','r')

f.readline()
nv = int(f.readline())

f.readline()
ne = int(f.readline())

f.readline()
nhe = int(f.readline())

f.readline()
nf = int(f.readline())

f.readline()
xy = np.zeros((nv,2), dtype=float)
for i in range(nv):
    xy[i] = [float(a) for a in f.readline().split()]

f.readline()
e2h = np.zeros(ne, dtype=int)
for i in range(ne):
    e2h[i] = int(f.readline())-1

f.readline()
esmooth = np.zeros(ne, dtype=bool)
for i in range(ne):
    l = int(f.readline())
    if l == 0:
        esmooth[i] = False
    else:
        esmooth[i] = True

f.readline()
he = np.zeros((nhe,6), dtype=int)
for i in range(nhe):
    he[i] = [int(a)-1 for a in f.readline().split()]

f.readline()
face_outer = np.zeros(nf, dtype=int)
for i in range(nf):
    face_outer[i] = int(f.readline())-1

f.readline()
face_inner = []
for i in range(nf):
    ni = int(f.readline())
    if ni == 0:
        face_inner.append([])
    else:
        face_inner.append([int(a)-1 for a in f.readline().split()])

f.readline()
v2h = np.zeros(nv, dtype=int)
for i in range(nv):
    v2h[i] = int(f.readline())-1

f.close()


######## plot DCEL ########
cm = plt.get_cmap('Set2')#'rainbow')
cle = [cm(1.*i/ne) for i in range(ne)]
#clf = myc.random_pastels(nf)
cm = plt.get_cmap('Pastel2')
clf = [cm(1.*i/nf) for i in range(nf)]

fv = 0.07
fw = 0.03

fig, ax = plt.subplots()

xmin = np.amin(xy[:,0])
xmax = np.amax(xy[:,0])
ymin = np.amin(xy[:,1])
ymax = np.amax(xy[:,1])

rngx = xmax - xmin
rngy = ymax - ymin
mrgx = 0.1*rngx
mrgy = 0.1*rngy

for i in range(nf):
    plot_face(ax, face_outer[i], face_inner[i], xy, he, clf[i], 'none')


plt.plot(xy[:,0], xy[:,1], 'ko')

for i in range(nhe):
    a = xy[he[i,4]]
    b = xy[he[he[i,3],4]]
    v = b - a
    lv = np.hypot(v[0], v[1])
    w = np.array((-v[1], v[0]))
    ie = he[i,5]
    
    ax.arrow(a[0] + fv*v[0] + fw*w[0], a[1] + fv*v[1] + fw*w[1],
             (1.0 - 2.0*fv)*v[0], (1.0 - 2.0*fv)*v[1],
             fc=cle[ie], ec=cle[ie],
             head_length=fv*lv, head_width=0.8*fv*lv,
             length_includes_head=True, shape='right')

for i in range(ne):
    a = xy[he[e2h[i],4]]
    b = xy[he[he[e2h[i],3],4]]
    if esmooth[i]:
        lsty = '--'
    else:
        lsty = '-'
    plt.plot([a[0],b[0]], [a[1],b[1]], color=cle[i], ls=lsty)

ax.set_xlim(xmin - mrgx, xmax + mrgx)
ax.set_ylim(ymin - mrgy, ymax + mrgy)
ax.set_aspect('equal', adjustable='box')
plt.show()







