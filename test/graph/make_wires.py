import matplotlib.pyplot as plt
import numpy as np
from math import cos, sin, atan2, pi
import sys

#############################################################
#############################################################
class GraphNode:
    def __init__(self, co=(0,0), ingoing=None, outgoing=None):
        self.co = co
        if ingoing is None:
            self.ingoing = []
        else:
            self.ingoing = ingoing
        if outgoing is None:
            self.outgoing = []
        else:
            self.outgoing = outgoing
        return
#############################################################
class GraphArc:
    def __init__(self, orig=None, dest=None, angles=None):
        if orig is None:
            self.orig = []
        else:
            self.orig = orig
        if dest is None:
            self.dest = []
        else:
            self.dest = dest
        if angles is None:
            self.angles = []
        else:
            self.angles = angles
        return
    
    def curve(self, n=100, tension=1):
        return curved_arc([self.orig.co, self.dest.co], self.angles, n, tension)
    
    def plot(self, ax, n=100, tension=1, color='k', arrow_length=0):
        aabb = plot_curved_arc(ax, [self.orig.co, self.dest.co], self.angles, n, tension, color, arrow_length)
        return aabb
#############################################################
"""
class GraphCycle:
    def __init__(self, nodes=None, arcs=None):
        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes
        if arcs is None:
            self.arcs = []
        else:
            self.arcs = arcs
        return
"""
class GraphCycle:
    def __init__(self, points=None, angles=None):
        if points is None:
            self.points = []
        else:
            self.points = points
        if angles is None:
            self.angles = []
        else:
            self.angles = angles
        return
#############################################################
def remove_dangling_branches(nodes, arcs):
    while True:
        changes = False
        for node in nodes:
            remove = False
            if len(node.ingoing) < 1:
                remove = True
                for arc in node.outgoing:
                    arc.dest.ingoing.remove(arc)
                    arcs.remove(arc)
            if len(node.outgoing) < 1:
                remove = True
                for arc in node.ingoing:
                    arc.orig.outgoing.remove(arc)
                    arcs.remove(arc)
            if remove:
                changes = True
                nodes.remove(node)
        if not changes: return nodes, arcs
#############################################################
def plot_graph(ax, nodes, arcs, node_color='k', arc_color='k'):
    huge = 1e6
    aabb = [huge, -huge, huge, -huge]
    for arc in arcs:
        aabb_arc = arc.plot(ax, color=arc_color, arrow_length=0.033)
        for i in range(2):
            aabb[2*i] = min(aabb[2*i], aabb_arc[2*i])
            aabb[2*i+1] = max(aabb[2*i+1], aabb_arc[2*i+1])
    for node in nodes:
        ax.plot(node.co[0], node.co[1], 'o', c=node_color, mew=0)
    return aabb
#############################################################
def plot_curved_arc(ax, points, angles, n=100, tension=1, color='k', arrow_length=0):
    p = curved_arc(points, angles, n, tension)
    ax.plot(p[:,0], p[:,1], color=color, lw=1)
    if arrow_length > 0:
        #m = int(n/2)
        m = curve_midpoint_index(p)
        pm = p[m]
        dpm = p[m+1] - pm
        dpm = arrow_length*dpm/np.hypot(dpm[0], dpm[1])
        ax.arrow(
            x=pm[0], y=pm[1],
            dx=dpm[0], dy=dpm[1],
            color=color,
            head_length=arrow_length, 
            head_width=arrow_length, 
            length_includes_head=True
        )
    xymin = np.amin(p, axis=0)
    xymax = np.amax(p, axis=0)
    return (xymin[0], xymax[0], xymin[1], xymax[1])
#############################################################
def diff_angle(a1, a2):
    c1 = cos(a1)
    s1 = sin(a1)
    c2 = cos(a2)
    s2 = sin(a2)
    return atan2(s1*c2 - c1*s2, c1*c2 + s1*s2)
#############################################################
def curved_arc(points, angles, n=100, tension=1):
    l = np.hypot(points[1][0] - points[0][0], points[1][1] - points[0][1])
    p0 = np.array([points[0][0], points[0][1]])
    p1 = np.array([points[1][0], points[1][1]])
    m0 = tension*l*np.array([cos(angles[0]), sin(angles[0])])
    m1 = tension*l*np.array([cos(angles[1]), sin(angles[1])])
    t = np.linspace(0,1,n)
    h00 = 2*t**3 - 3*t**2 + 1
    h10 = t**3 - 2*t**2 + t
    h01 = -2*t**3 + 3*t**2
    h11 = t**3 - t**2
    c = np.zeros((n,2))
    for i in range(n):
        c[i] = h00[i]*p0 + h01[i]*p1 + h10[i]*m0 + h11[i]*m1
    return c
#############################################################
def curve_midpoint_index(points):
    n = len(points)
    s = np.zeros(n)
    for i in range(1,n):
        s[i] = s[i-1] + np.sqrt(np.sum((points[i] - points[i-1])**2))
    smid = 0.5*s[-1]
    for i in range(n-1):
        if s[i] <= smid and s[i+1] >= smid:
            return i
#############################################################
def extended_axlim((xmin, xmax), (ymin, ymax), mrg=0.05):
    xrng = xmax - xmin
    yrng = ymax - ymin
    xmin -= mrg*xrng
    xmax += mrg*xrng
    ymin -= mrg*yrng
    ymax += mrg*yrng
    return xmin, xmax, ymin, ymax
#############################################################


####################################################
# READ EMBEDDED GRAPH
####################################################
args = sys.argv
if len(args) < 2:
    fname = '/stck/bandrieu/Bureau/FFTsurf2/test/graph/data_test_loop.dat'
else:
    num = int(args[1])
    fname = 'test'+format(num,'02')+'.dat'
    
f = open(fname, 'r')

narc = int(f.readline())
arc_end = []
for iarc in range(narc):
    arc_end.append([int(a)-1 for a in f.readline().split()])
arc_ang = []
for iarc in range(narc):
    arc_ang.append([float(a) for a in f.readline().split()])

nnod = int(f.readline())
nodes = []
for inod in range(nnod):
    nodes.append(
        GraphNode(
            co=[float(a) for a in f.readline().split()],
            ingoing=[],
            outgoing=[]
        )
    )
f.close()

arcs = []
for iarc in range(narc):
    arcs.append(
        GraphArc(
            orig=nodes[arc_end[iarc][0]],
            dest=nodes[arc_end[iarc][1]],
            angles=arc_ang[iarc]
        )
    )

####################################################
# COMPUTE INITIAL NODE -> ARC INCIDENCE
# (distinguish outgoing/ingoing arcs)
####################################################
for arc in arcs:
    arc.orig.outgoing.append(arc)
    arc.dest.ingoing.append(arc)

####################################################
# PLOT INITIAL GRAPH
####################################################
mrg = 5e-2

fig, axes = plt.subplots(1,2)
aabb = plot_graph(axes[0], nodes, arcs)
xmin, xmax, ymin, ymax = extended_axlim(aabb[0:2], aabb[2:4], mrg)

    
####################################################
# MAKE WIRES(CYCLES)
####################################################
cycles = []

while len(arcs) > 0:
    nodes, arcs = remove_dangling_branches(nodes, arcs)
    print('%d cycle(s), %d arc(s) left' % (len(cycles), len(arcs)))

    # choose an arbitrary arc from the graph
    arc = arcs[0]

    # initiate a new (temp) cycle starting from that arc
    cycle_arcs = [arc]
    cycle_nodes = [arc.orig, arc.dest]

    # set current arc's destination as current node
    node = arc.dest
    while True:
        # get leftmost ingoing arc at current node (relative to current arc)
        maxangle_in = -pi
        for inarc in node.ingoing:
            if inarc == arc: continue
            delta_angle = diff_angle(inarc.angles[1] + pi, arc.angles[1])
            maxangle_in = max(maxangle_in, delta_angle)
        
        # get leftmost outgoing arc at current node (relative to current arc)
        maxangle_out = -2*pi # must be < -pi
        leftmost_arc = None
        for outarc in node.outgoing:
            delta_angle = diff_angle(outarc.angles[0], arc.angles[1])
            if delta_angle > maxangle_out :
                leftmost_arc = outarc
                maxangle_out = max(maxangle_out, delta_angle)
        if leftmost_arc is None:
            exit('leftmost_arc is None')
        if maxangle_out < maxangle_in : break # this is an invalid cycle
        
        arc = leftmost_arc
        cycle_arcs.append(arc)

        node = arc.dest
        if node == cycle_nodes[0]:
            print('+1 cycle (used %d arcs)' % len(cycle_arcs))
            cycles.append(
                GraphCycle(
                    points=[n.co for n in cycle_nodes],
                    angles=[a.angles for a in cycle_arcs]
                )
            )
            for arc in cycle_arcs:
                arc.orig.outgoing.remove(arc)
                arc.dest.ingoing.remove(arc)
                arcs.remove(arc)
            break
        else:
            cycle_nodes.append(node)

####################################################
# PLOT CYCLES
####################################################
cmap = plt.get_cmap('Set2')#'rainbow')
cycle_color = [cmap(1.*i/len(cycles)) for i in range(len(cycles))]

for ic, cycle in enumerate(cycles):
    x = [p[0] for p in cycle.points]
    y = [p[1] for p in cycle.points]
    m = len(cycle.points)
    for i in range(m):
        xmin = min(xmin, x[i])
        xmax = max(xmax, x[i])
        ymin = min(ymin, y[i])
        ymax = max(ymax, y[i])
        j = (i+1)%m
        axes[1].plot(x, y, 'o', c=cycle_color[ic], mew=0)
        aabb_arc = plot_curved_arc(
            axes[1], 
            points=[
                [x[i], y[i]],
                [x[j], y[j]]
            ],
            angles=cycle.angles[i], 
            color=cycle_color[ic],
            arrow_length=0.033
        )
for ax in axes:
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
plt.show()