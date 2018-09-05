import numpy as np
import matplotlib.pyplot as plt


#############################################################
class node:
    def __init__( self, level, parent, children ):
        self.level = level
        self.parent = parent
        self.children = children
#############################################################
class face:
    def __init__( self, outer, inner ):
        self.outer = outer
        self.inner = inner
#############################################################
def assign_xy(nodes, inod, width, x, y):
    m = len(nodes[inod].children)
    for i,n in enumerate(nodes[inod].children):
        y[n] = y[inod] - 1.0
        x[n] = x[inod]
        if m > 1:
            x[n] += -0.4*width + 0.8*width*float(i)/float(m-1)
        x, y = assign_xy(nodes, n, 0.5*width, x, y)
    return x, y
#############################################################
def make_face(nodes, inod, faces):
    if nodes[inod].level%2 == 0:
        f = face(outer=inod, inner=nodes[inod].children)
        faces.append(f)
    for jnod in nodes[inod].children:
        faces = make_face(nodes, jnod, faces)
    return faces
#############################################################

arc2nod = np.loadtxt('arcs.dat', dtype=int)
arc2nod = np.random.permutation(arc2nod) # randomly permute arcs

#arc2nod = np.array(([5,2],[4,3],[1,3],[6,4],[6,3]))
print arc2nod


narc = len(arc2nod)
nnod = np.amax(arc2nod)


child2parent = [None]*nnod
for inod in range(nnod):
    child2parent[inod] = []
    
for iarc in range(narc):
    inod = arc2nod[iarc,0] - 1
    jnod = arc2nod[iarc,1] - 1
    child2parent[inod].append(jnod)

maxlevel = 0
for inod in range(nnod):
    print 'node #',inod,', ascendants :',child2parent[inod],' (length = ',len(child2parent[inod]),')'
    maxlevel = max(maxlevel, len(child2parent[inod]))
print 'maxlevel = ',maxlevel,'\n'

    
nodes = [None]*nnod
for inod in range(nnod):
    nodes[inod] = node(0, [], [])

# build tree from top down level by level
for ilev in range(maxlevel+1):
    for inod in range(nnod):
        if len(child2parent[inod]) == ilev:
            for jnod in child2parent[inod]:
                if nodes[jnod].level == ilev - 1:
                    nodes[inod].parent = jnod
                    nodes[jnod].children.append(inod)
                    break
            nodes[inod].level = ilev


#############################################################
for i,n in enumerate(nodes):
    print 'node #',i
    print 'level =',n.level
    print 'parent =',n.parent
    print 'children =',n.children,'\n'
    
#############################################################

# get faces
faces = []
for inod in range(nnod):
    if nodes[inod].level == 0:
        faces = make_face(nodes, inod, faces)

for f in faces:
    print 'outer : ', f.outer
    print 'inner : ', f.inner, '\n'
        
#############################################################


cl = ['b','r']

l0 = []
for inod in range(nnod):
    if nodes[inod].level == 0:
        l0.append(inod)

x = np.zeros(nnod, dtype=float)
y = np.zeros(nnod, dtype=float)
w0 = float(maxlevel)
for i, inod in enumerate(l0):
    x[inod] = float(i)*w0/float(len(l0))
    x, y = assign_xy(nodes, inod, 0.5*w0, x, y)

fig, ax = plt.subplots()

for inod in range(nnod):
    plt.text(x[inod], y[inod], str(inod+1),
             color=cl[nodes[inod].level%2],
             ha="center", va="center",
             bbox=dict(boxstyle="round",
                       ec=cl[nodes[inod].level%2],
                       fc='w'))
    a = np.array((x[inod], y[inod]))
    for jnod in nodes[inod].children:
        b = np.array((x[jnod], y[jnod]))
        v = b - a
        #a = a + 0.1*v
        #b = a + 0.8*v
        plt.plot([a[0], b[0]], [a[1], b[1]], 'k-')


mn = np.array((np.amin(x), np.amin(y)))
mx = np.array((np.amax(x), np.amax(y)))
rng = mx - mn
mrg = 0.1*rng
mx = mx + mrg
mn = mn - mrg

ax.set_xlim(mn[0], mx[0])
ax.set_ylim(mn[1], mx[1])
ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect('equal', adjustable='box')
plt.axis('off')
plt.show()


