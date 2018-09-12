import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/stck/bandrieu/Bureau/Python/mylibs/')
#sys.path.append('/home/bastien/Bureau/Python/')
import my_colors as myc

#############################################################
class trimming_wire:
    def __init__( self, nod, arc ):
        self.nod = nod
        self.arc = arc
#############################################################
def solve_2x2(a11, a12, a21, a22, b):
    det = a11*a22 - a12*a21
    if abs(det) < 1.e-15:
        return 0, False
    else:
        return np.array( (b[0]*a22 - a12*b[1], a11*b[1] - b[0]*a21) ) / det, True
#############################################################


####################################################
# READ EMBEDDED GRAPH
####################################################
args = sys.argv
if len(args) < 2:
    num = 1
else:
    num = int(args[1])
    
f = open('graph_'+format(num,'03')+'.dat', 'r')

narc = int(f.readline())
arc_end = np.zeros((narc,2), dtype=int)
for iarc in range(narc):
    arc_end[iarc] = [int(a)-1 for a in f.readline().split()]
arc_ang = np.zeros((narc,2), dtype=float)
for iarc in range(narc):
    arc_ang[iarc] = [float(a) for a in f.readline().split()]

nnod = int(f.readline())
xy = np.zeros((nnod,2), dtype=float)
for inod in range(nnod):
    xy[inod] = [float(a) for a in f.readline().split()]

f.close()



f = open('wires_'+format(num,'03')+'.dat', 'r')
nwires = int(f.readline())
Wires = []
for i in range(nwires):
    wirearc = [int(a)-1 for a in f.readline().split()]
    wirenod = [int(a)-1 for a in f.readline().split()]
    Wires.append(trimming_wire(wirenod, wirearc))
f.close()


####################################################
# PLOT RAW GRAPH AND WIRES
####################################################

t = np.linspace(0,1,20)
b = np.stack( [np.power(1.0-t,2) , 2.0*(1.0-t)*t , np.power(t,2)] )

# Bezier control points of curved arc embeddings (cubic curves)
cp = np.zeros((narc,3,2), dtype=float)
t = np.zeros((2,2))
for iarc in range(narc):
    for iend in range(2):
        t[0,iend] = np.cos(arc_ang[iarc,iend])
        t[1,iend] = np.sin(arc_ang[iarc,iend])
        dx = xy[arc_end[iarc,1]] - xy[arc_end[iarc,0]]
        s, nonsingular = solve_2x2(t[0,0], t[0,1], t[1,0], t[1,1], dx)
        cp[iarc,0] = xy[arc_end[iarc,0]]
        if nonsingular:
            cp[iarc,1] = xy[arc_end[iarc,0]] + s[0]*t[:,0]
        else:
            cp[iarc,1] = 0.5*(xy[arc_end[iarc,0]] + xy[arc_end[iarc,1]])
        cp[iarc,2] = xy[arc_end[iarc,1]]

# xy bounding box
mn = np.zeros(2)
mx = np.zeros(2)
for idim in range(2):
    mn[idim] = np.amin(cp[:,:,idim])
    mx[idim] = np.amax(cp[:,:,idim])
rng = mx - mn
mrg = 0.1*rng
mx = mx + mrg
mn = mn - mrg


# plot parameters
lwa = 1.0               # arc linewidth
lq = 0.033*np.amax(rng) # arrow length
wq = 0.45*lq            # arrow width
la = lq



fig, ax = plt.subplots(1,2, figsize=(12,6))

# raw graph
cm = plt.get_cmap('Set2')#'rainbow')
cl = [cm(1.*i/narc) for i in range(narc)]

for iarc in range(narc):
    if True:
        xya = np.matmul(np.transpose(cp[iarc,:,:]), b)
        ax[0].plot(xya[0], xya[1], '-', color=cl[iarc], lw=lwa)
        xym = 0.25*cp[iarc,0] + 0.5*cp[iarc,1] + 0.25*cp[iarc,2]
        dxym = cp[iarc,2] - cp[iarc,0]
        ldxym = np.hypot(dxym[0], dxym[1])
        if ldxym > 1.e-15:
            dxym = lq*dxym/np.hypot(dxym[0], dxym[1])
            ax[0].arrow(xym[0], xym[1], dxym[0], dxym[1],
                        fc=cl[iarc], ec=cl[iarc],
                        head_length=lq, head_width=2.0*wq, length_includes_head=True, shape="right")
        for i,inod in enumerate(arc_end[iarc]):
            ax[0].arrow(xy[inod,0], xy[inod,1], lq*np.cos(arc_ang[iarc,i]), lq*np.sin(arc_ang[iarc,i]),
                        fc=cl[iarc], ec=cl[iarc],
                        head_length=0.3*lq, head_width=wq, length_includes_head=True)
for inod in range(nnod):
    if True:
        ax[0].plot(xy[inod,0], xy[inod,1], 'k.', markersize=10)
        ax[0].text(xy[inod,0], xy[inod,1], str(inod))
ax[0].set_title('Raw graph', size = 12, weight = 'bold')



cm = plt.get_cmap('Set2')#'rainbow')
#cl = myc.random_pastels(len(Wires))
cl = [cm(1.*i/len(Wires)) for i in range(len(Wires))]

for il, wire in enumerate(Wires):
    for iarc in wire.arc:
        xya = np.matmul(np.transpose(cp[iarc,:,:]), b)
        ax[1].plot(xya[0], xya[1], '-', color=cl[il], lw=lwa)
        xym = 0.25*cp[iarc,0] + 0.5*cp[iarc,1] + 0.25*cp[iarc,2]
        dxym = cp[iarc,2] - cp[iarc,0]
        dxym = lq*dxym/np.hypot(dxym[0], dxym[1])
        ax[1].arrow(xym[0], xym[1], dxym[0], dxym[1],
                    fc=cl[il], ec=cl[il],
                    head_length=lq, head_width=2.0*wq, length_includes_head=True, shape="right")
    ax[1].plot(xy[wire.nod,0], xy[wire.nod,1], 'k.', markersize=10)
ax[1].set_title('Wires', size = 12, weight = 'bold')


for a in ax:
    a.set_xlim(mn[0], mx[0])
    a.set_ylim(mn[1], mx[1])
    a.set_xticks([])
    a.set_yticks([])
    a.set_aspect('equal', adjustable='box')

#plt.tight_layout()
plt.show()



