import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/stck/bandrieu/Bureau/Python/mylibs/')
#sys.path.append('/home/bastien/Bureau/Python/')
import my_colors as myc



#############################################################
class trimming_loop:
    def __init__( self, nod, arc ):
        self.nod = nod
        self.arc = arc
#############################################################
def diff_angle(a1, a2):
    c1 = np.cos(a1)
    s1 = np.sin(a1)
    c2 = np.cos(a2)
    s2 = np.sin(a2)
    return np.arctan2(s1*c2 - c1*s2, c1*c2 + s1*s2)
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
    fname = '/d/bandrieu/stck/Bureau/FFTsurf2/test/graph/data_test_loop.dat'
else:
    num = int(args[1])
    fname = 'test'+format(num,'02')+'.dat'
    
f = open(fname, 'r')

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


####################################################
# COMPUTE INITIAL NODE -> ARC INCIDENCE
# (distinguish outgoing/ingoing arcs)
####################################################
nod_out = [None]*nnod
nod_in  = [None]*nnod
for inod in range(nnod):
    nod_out[inod] = []
    nod_in[inod]  = []
nout = np.zeros(nnod, dtype=int)
nin  = np.zeros(nnod, dtype=int)
for iarc in range(narc):
    inod = arc_end[iarc][0]
    nout[inod] += 1
    nod_out[inod].append(iarc)
    jnod = arc_end[iarc][1]
    nin[jnod] += 1
    nod_in[jnod].append(iarc)

explored = np.zeros(narc, dtype=bool)
####################################################
# DISABLE INVALID NODES/ARCS
validarc = np.ones(narc, dtype=bool)
for inod in range(nnod):
    if nout[inod] < 1 or nin[inod] < 1:
        for iarc in nod_out[inod]:
            print 'arc #', iarc, arc_end[iarc],' is invalid'
            explored[iarc] = True
            validarc[iarc] = False
            nod_out[inod].remove(iarc)
            nout[inod] -= 1
            jnod = arc_end[iarc][1]
            nod_in[jnod].remove(iarc)
            nin[jnod] -= 1
        for iarc in nod_in[inod]:
            print 'arc #', iarc, arc_end[iarc],' is invalid'
            explored[iarc] = True
            validarc[iarc] = False
            nod_in[inod].remove(iarc)
            nin[inod] -= 1
            jnod = arc_end[iarc][0]
            nod_out[jnod].remove(iarc)
            nout[jnod] -= 1
####################################################



####################################################
# MAKE LOOPS
####################################################
Loops = []
usedinaloop = np.zeros(narc, dtype=bool)
inod = 0
for inod in range(nnod):
    print 'node #',inod,', nod_out =',nod_out[inod], ', nod_in =',nod_in[inod]
    if nout[inod] < 1 or nin[inod] < 1:
        # such a node cannot be part of a new loop
        for iarc in nod_out[inod]:
            explored[iarc] = True
        for iarc in nod_in[inod]:
            explored[iarc] = True
        continue
    for karc in nod_out[inod]:
        loopnod = []
        loopnod.append(inod)
        looparc = []
        looparc.append(karc)
        jnod = arc_end[karc,1]
        valid = True
        while True: #jnod != inod:
            if nout[jnod] > 0:
                # get leftmost ingoing arc at current node (except current arc)
                mxangin = -np.pi
                for jarc in nod_in[jnod]:
                    if jarc != looparc[-1]:
                        angj = diff_angle(arc_ang[jarc,1]+np.pi, arc_ang[looparc[-1],1])
                        mxangin = max(mxangin, angj)
                # get leftmost outgoing arc at current node
                mxangout = -np.pi
                iarc = -1
                for jarc in nod_out[jnod]:
                    angj = diff_angle(arc_ang[jarc,0], arc_ang[looparc[-1],1])
                    if angj > mxangout:
                        mxangout = angj
                        iarc = jarc
                if mxangout < mxangin:
                    valid = False
                    break
            else:
                valid = False
                break
            if jnod == loopnod[0]:
                if iarc != looparc[0]:
                    valid = False
                break
            if usedinaloop[iarc]:
                print '!!! arc already used in a previous loop!'
            loopnod.append(jnod)
            looparc.append(iarc)
            jnod = arc_end[iarc,1]

        for iarc in looparc:
            explored[iarc] = True

        if valid:
            l = trimming_loop(loopnod, looparc)
            print '\t+1 loop, nodes=',loopnod,', arcs=',looparc
            Loops.append(l)
            # update node -> arc incidences
            for inod in loopnod:
                if nout[inod] == 1:
                    nout[inod] = 0
                    nod_out[inod] = []
                else:
                    for iarc in looparc:
                        if iarc in nod_out[inod]:
                            nout[inod] -= 1
                            nod_out[inod] = [a for a in nod_out[inod] if a != iarc]
            # update mask for used arcs
            for iarc in looparc:
                usedinaloop[iarc] = True
        else:
            print '\tinvalid loop, nodes=',loopnod,', arcs=',looparc
            print 'last node, out :',nod_out[loopnod[-1]],', in :',nod_in[loopnod[-1]]
            continue
            
print 'explored =',explored



####################################################
# PLOT RAW GRAPH AND LOOPS
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




fig, ax = plt.subplots(1,2, figsize=(12,6))

# raw graph
cm = plt.get_cmap('Set2')#'rainbow')
cl = [cm(1.*i/narc) for i in range(narc)]

for iarc in range(narc):
    if True:#validarc[iarc]:
        xya = np.matmul(np.transpose(cp[iarc,:,:]), b)
        ax[0].plot(xya[0], xya[1], '-', color=cl[iarc], lw=lwa)
        xym = 0.25*cp[iarc,0] + 0.5*cp[iarc,1] + 0.25*cp[iarc,2]
        dxym = cp[iarc,2] - cp[iarc,0]
        dxym = lq*dxym/np.hypot(dxym[0], dxym[1])
        ax[0].arrow(xym[0], xym[1], dxym[0], dxym[1],
                        fc=cl[iarc], ec=cl[iarc],
                        head_length=lq, head_width=wq, length_includes_head=True)
ax[0].plot(xy[:,0], xy[:,1], 'k.', markersize=10)
for inod in range(nnod):
    ax[0].text(xy[inod,0], xy[inod,1], str(inod))
ax[0].set_title('Raw graph', size = 12, weight = 'bold')


cm = plt.get_cmap('Set2')#'rainbow')
#cl = myc.random_pastels(len(Loops))
cl = [cm(1.*i/len(Loops)) for i in range(len(Loops))]

for il, loop in enumerate(Loops):
    for iarc in loop.arc:
        xya = np.matmul(np.transpose(cp[iarc,:,:]), b)
        ax[1].plot(xya[0], xya[1], '-', color=cl[il], lw=lwa)
        xym = 0.25*cp[iarc,0] + 0.5*cp[iarc,1] + 0.25*cp[iarc,2]
        dxym = cp[iarc,2] - cp[iarc,0]
        dxym = lq*dxym/np.hypot(dxym[0], dxym[1])
        ax[1].arrow(xym[0], xym[1], dxym[0], dxym[1],
                    fc=cl[il], ec=cl[il],
                    head_length=lq, head_width=wq, length_includes_head=True)
    ax[1].plot(xy[loop.nod,0], xy[loop.nod,1], 'k.', markersize=10)
ax[1].set_title('Loops', size = 12, weight = 'bold')

for a in ax:
    a.set_xlim(mn[0], mx[0])
    a.set_ylim(mn[1], mx[1])
    a.set_xticks([])
    a.set_yticks([])
    a.set_aspect('equal', adjustable='box')

#plt.tight_layout()
plt.show()




