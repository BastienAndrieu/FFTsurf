import matplotlib.pyplot as plt
import numpy as np
from math import cos, sin, atan2, pi
import sys
sys.path.append('/d/GitHub/FFTsurf/test/graph/')
from lib_curved_graph import *


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
def diff_angle(a1, a2):
    c1 = cos(a1)
    s1 = sin(a1)
    c2 = cos(a2)
    s2 = sin(a2)
    return atan2(s1*c2 - c1*s2, c1*c2 + s1*s2)
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
args = sys.argv
if len(args) < 2:
    #fname = '/stck/bandrieu/Bureau/FFTsurf2/test/graph/data_test_loop.dat'
    fname = 'test_demo.dat'
else:
    num = int(args[1])
    if num < 0:
        #fname = '/stck/bandrieu/Bureau/FFTsurf2/test/graph/data_test_loop.dat'
        fname = 'test_demo.dat'
    else:
        fname = 'test'+format(num,'02')+'.dat'

if len(args) < 3:
    show_steps = False
else:
    show_steps = (int(args[2]) != 0)
print 'SHOW STEPS?', show_steps
    
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
perm = np.random.permutation(narc)
for jarc in range(narc):
    iarc = perm[jarc]
    arcs.append(
        GraphArc(
            orig=nodes[arc_end[iarc][0]],
            dest=nodes[arc_end[iarc][1]],
            angles=arc_ang[iarc]
        )
    )
####################################################


####################################################
# COMPUTE INITIAL NODE -> ARC INCIDENCE
# (distinguish outgoing/ingoing arcs)
for arc in arcs:
    arc.orig.outgoing.append(arc)
    arc.dest.ingoing.append(arc)
####################################################


####################################################
# PLOT INITIAL GRAPH
mrg = 5e-2

fig, axes = plt.subplots(1,2)
aabb = plot_graph(axes[0], nodes, arcs)
xmin, xmax, ymin, ymax = extended_axlim(aabb[0:2], aabb[2:4], mrg)
if show_steps:
    plt.close()
####################################################

    
####################################################
# MAKE WIRES(CYCLES)
cycles = []

last_arcs = []

        
while len(arcs) > 0:
    if show_steps:
        figsteps, axessteps = plt.subplots(1, 3, figsize=(14.4, 4.8))
        plot_graph(axessteps[0], nodes, arcs)
        axessteps[0].set_title('graph')
        
    nodes, arcs = remove_dangling_branches(nodes, arcs)
    
    if show_steps:
        plot_graph(axessteps[1], nodes, arcs)
        axessteps[1].set_title('rm dangling branches')
        
    print('%d cycle(s), %d arc(s) left' % (len(cycles), len(arcs)))

    # choose an arbitrary arc from the graph
    for arc in arcs:
        if arc not in last_arcs: break
    #arc = arcs[0]
    print 'START FROM ARC #', arcs.index(arc)

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
        if maxangle_out < maxangle_in :
            last_arcs.append(cycle_arcs[0])
            #if show_steps: plt.close()
            if show_steps:
                plot_graph(axessteps[2], nodes, arcs)
                axessteps[2].set_title('wrong cycle')
                points=[n.co for n in cycle_nodes]
                angles=[a.angles for a in cycle_arcs]
                x = [p[0] for p in points]
                y = [p[1] for p in points]
                m = len(points)
                for i in range(m-1):
                    axessteps[2].plot(x, y, 'o', c='r', mew=0)
                    plot_curved_arc(
                        axessteps[2], 
                        points=[
                            [x[i], y[i]],
                            [x[i+1], y[i+1]]
                        ],
                        angles=angles[i], 
                        color='r',
                        arrow_length=0.033
                    )
                for ax in axessteps:
                    ax.set_xlim(xmin, xmax)
                    ax.set_ylim(ymin, ymax)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.set_aspect('equal')
                plt.show()
                plt.close()
            break # this is an invalid cycle
        
        arc = leftmost_arc
        cycle_arcs.append(arc)

        node = arc.dest
        if node == cycle_nodes[0]:
            print('+1 cycle (used %d arcs)' % len(cycle_arcs))
            last_arcs = []
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
            if show_steps:
                plot_graph(axessteps[2], nodes, arcs)
                axessteps[2].set_title('extracted cycle')
                cycle = cycles[-1]
                x = [p[0] for p in cycle.points]
                y = [p[1] for p in cycle.points]
                m = len(cycle.points)
                for i in range(m):
                    j = (i+1)%m
                    axessteps[2].plot(x, y, 'o', c='g', mew=0)
                    plot_curved_arc(
                        axessteps[2], 
                        points=[
                            [x[i], y[i]],
                            [x[j], y[j]]
                        ],
                        angles=cycle.angles[i], 
                        color='g',
                        arrow_length=0.033
                    )
                for ax in axessteps:
                    ax.set_xlim(xmin, xmax)
                    ax.set_ylim(ymin, ymax)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.set_aspect('equal')
                plt.show()
                plt.close()
            break
        else:
            cycle_nodes.append(node)
####################################################


####################################################
# PLOT CYCLES
if not show_steps:
    cmap = plt.get_cmap('Set2')#'rainbow')
    cycle_color = [cmap(1.*i/len(cycles)) for i in range(len(cycles))]

    for ic, cycle in enumerate(cycles):
        x = [p[0] for p in cycle.points]
        y = [p[1] for p in cycle.points]
        m = len(cycle.points)
        for i in range(m):
            """xmin = min(xmin, x[i])
            xmax = max(xmax, x[i])
            ymin = min(ymin, y[i])
            ymax = max(ymax, y[i])"""
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
####################################################

####################################################
# MAKE FACES
faces = make_faces(cycles, n=100, tension=1)
face_color = [cmap(1.*i/len(faces)) for i in range(len(faces))]

fig, ax = plt.subplots()

for iface, face in enumerate(faces):
    for cycle in [face.outer] + face.inner:
        x = [p[0] for p in cycle.points]
        y = [p[1] for p in cycle.points]
        m = len(cycle.points)
        for i in range(m):
            j = (i+1)%m
            ax.plot(x, y, 'o', c=cycle_color[iface], mew=0)
            plot_curved_arc(
                ax, 
                points=[
                    [x[i], y[i]],
                    [x[j], y[j]]
                ],
                angles=cycle.angles[i], 
                color=cycle_color[iface],
                arrow_length=0.033
            )

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_aspect('equal')
plt.show()
####################################################


####################################################
# EXPORT FACES (for TikZ)
f = open('graph_faces_nwires.dat', 'w')
for face in faces:
    f.write('%d\n' % len(face.inner))
f.close()


for iface, face in enumerate(faces):
    f = open('graph_face_%d_outer.dat' % (iface+1), 'w')
    cycle = face.outer
    x = [p[0] for p in cycle.points]
    y = [p[1] for p in cycle.points]
    m = len(cycle.points)
    for i in range(m):
        j = (i+1)%m
        b = curved_arc_bezier_control_points(
            [(x[i], y[i]), (x[j], y[j])],
            cycle.angles[i],
            tension=1
        )
        for bx, by in b:
            f.write('%s, %s, ' % (bx, by))
        f.write('\n')
    f.close()

    for k, cycle in enumerate(face.inner):
        f = open('graph_face_%d_inner_%d.dat' % ((iface+1), (k+1)), 'w')
        x = [p[0] for p in cycle.points]
        y = [p[1] for p in cycle.points]
        m = len(cycle.points)
        for i in range(m):
            j = (i+1)%m
            b = curved_arc_bezier_control_points(
                [(x[i], y[i]), (x[j], y[j])],
                cycle.angles[i],
                tension=1
            )
            for bx, by in b:
                f.write('%s, %s, ' % (bx, by))
            f.write('\n')
        f.close()
####################################################



####################################################
# EXPORT FACES (full TikZ code)
f = open('graph_faces_tikzcode.tex', 'w')
# filled faces
for iface, face in enumerate(faces):
    f.write('\\path[fill=facecolor%d!30!white]\n' % (iface+1))
    for cycle in [face.outer] + face.inner:
        x = [p[0] for p in cycle.points]
        y = [p[1] for p in cycle.points]
        m = len(cycle.points)
        for i in range(m):
            j = (i+1)%m
            b = curved_arc_bezier_control_points(
                [(x[i], y[i]), (x[j], y[j])],
                cycle.angles[i],
                tension=1
            )
            if i == 0:
                f.write('(%s, %s) ' % (b[0][0], b[0][1]))
            f.write(
                '.. controls (%s, %s) and (%s, %s) .. (%s, %s)\n' %
                (b[1][0], b[1][1],
                 b[2][0], b[2][1],
                 b[3][0], b[3][1])
            )
            if i == m-1:
                f.write('-- cycle\n')
    f.write(';\n')
f.write('\n\n')

# halfedges
for iface, face in enumerate(faces):
    for icycle, cycle in enumerate([face.outer] + face.inner):
        if icycle == 0: # outer wire
            shade = 60
        else:           # inner wire
            shade = 80
        x = [p[0] for p in cycle.points]
        y = [p[1] for p in cycle.points]
        m = len(cycle.points)
        for i in range(m):
            j = (i+1)%m
            b = curved_arc_bezier_control_points(
                [(x[i], y[i]), (x[j], y[j])],
                cycle.angles[i],
                tension=1
            )
            f.write('\\draw[halfedge, facecolor%d!%d!black] ' % ((iface+1), shade))
            f.write(
                '(%s, %s) .. controls (%s, %s) and (%s, %s) .. (%s, %s);\n' %
                (b[0][0], b[0][1],
                 b[1][0], b[1][1],
                 b[2][0], b[2][1],
                 b[3][0], b[3][1])
            )
f.close()
####################################################
