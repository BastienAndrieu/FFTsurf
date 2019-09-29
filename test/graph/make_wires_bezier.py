import matplotlib.pyplot as plt
import numpy
from math import cos, sin, atan2, pi
import sys
sys.path.append('/d/bandrieu/GitHub/FFTsurf/test/graph/')
from lib_curved_graph import *
sys.path.append('/d/bandrieu/GitHub/Code/Python/')
import lib_bezier as lbez


#############################################################
def diff_angle(a1, a2):
    c1 = cos(a1)
    s1 = sin(a1)
    c2 = cos(a2)
    s2 = sin(a2)
    return atan2(s1*c2 - c1*s2, c1*c2 + s1*s2)
#############################################################

cmap = plt.get_cmap('Set2')#'rainbow')

########################################
args = sys.argv

if len(args) < 2:
    numtest = 1
else:
    numtest = int(args[1])

if len(args) < 3:
    show_steps = False
else:
    show_steps = (int(args[2]) != 0)
########################################

########################################
# READ INKSCAPE EXPORT
segments = []
f = open('inkscape_export_' + format(numtest,'02') + '.dat', 'r')
while True:
    line = f.readline()
    if ("" == line): break # end of file
    nbp = int(line)
    #
    bp = numpy.zeros((nbp,2))
    for i in range(nbp):
        bp[i] = [float(x) for x in f.readline().split()]
    segments.append(bp)
f.close()
########################################

########################################
# FLIP Y AND SHIFT/RESCALE
HUGE = 1e9
MRG = 5e-2

xymin = HUGE*numpy.ones(2)
xymax = -xymin
t = numpy.linspace(0,1,100)
for bp in segments:
    bp[:,1] = -bp[:,1]
    aabb = lbez.AABB_2d_bezier_curve(bp)
    xymin = numpy.minimum(xymin, aabb[0])
    xymax = numpy.maximum(xymax, aabb[1])

ctr = 0.5*(xymin + xymax)
rng = 0.5*numpy.amax(xymax - xymin)
for bp in segments:
    for j in range(2):
        bp[:,j] = (1 - MRG)*(bp[:,j] - ctr[j])/rng
########################################

"""
t = numpy.linspace(0,1,100)
j = int(len(t)/2)

fig, ax = plt.subplots()

for i, bp in enumerate(segments):
    xy = lbez.eval_bezier_curve(bp, t)
    ax.plot(xy[:,0], xy[:,1])
    ax.text(xy[j,0], xy[j,1], str(i))

ax.set_aspect('equal')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)

plt.show()
exit()
"""

########################################
# RANDOMIZE
perm = numpy.random.permutation(len(segments))
print 'perm = ', perm
segments = [segments[i] for i in perm]
########################################


########################################
# INTERPRET AS GRAPH
SPT = 1e-2
nodes = []
arcs = []
for bp in segments:
    #
    arcs.append(
        GraphArcB(
            bp=bp
        )
    )
    #
    for i, xy in enumerate(bp[[0,-1]]):
        alias = None
        for nod in nodes:
            if numpy.sum((nod.co - xy)**2) < SPT**2:
                alias = nod
                break
        if alias is None:
            nodes.append(
                GraphNode(
                    co=xy,
                    ingoing=[],
                    outgoing=[]
                )
            )
            alias = nodes[-1]
        if i == 0:
            arcs[-1].orig = alias
            alias.outgoing.append(arcs[-1])
        else:
            arcs[-1].dest = alias
            alias.ingoing.append(arcs[-1])
########################################




########################################
# PLOT

t = numpy.linspace(0,1,100)

fig, axes = plt.subplots(1,2)
"""
for bp in segments:
    xy = lbez.eval_bezier_curve(bp, t)
    axes[0].plot(xy[:,0], xy[:,1])
"""
plot_graph(axes[0], nodes, arcs)
"""
ax.set_aspect('equal')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
"""
#    plt.show()
if show_steps:
    plt.close()
########################################


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
    
    if len(arcs) < 1:
        print('Done!')
        if show_steps:
            for ax in axessteps:
                ax.set_xlim(-1,1)
                ax.set_ylim(-1,1)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_aspect('equal')
            plt.show()
        break

    # choose an arbitrary arc from the graph
    for arc in arcs:
        if arc not in last_arcs: break
    #arc = arcs[0]
    print 'START FROM ARC #', arcs.index(arc)

    # initiate a new (temp) cycle starting from that arc
    cycle_arcs = [arc]
    cycle_nodes = [arc.orig, arc.dest]
    #cycle_bp = [arc.bp]

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
            if show_steps:
                plot_graph(axessteps[2], nodes, arcs)
                axessteps[2].set_title('wrong cycle')
                for bp in [a.bp for a in cycle_arcs]:#cycle_bp:
                    plot_directed_bezier_curve(
                        axessteps[2], 
                        bp,
                        color='r',
                        arrow_length=0.033
                    )
                for ax in axessteps:
                    ax.set_xlim(-1,1)
                    ax.set_ylim(-1,1)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.set_aspect('equal')
                plt.show()
                plt.close()
            break # this is an invalid cycle
        
        arc = leftmost_arc
        cycle_arcs.append(arc)
        #cycle_bp.append(arc.bp)

        node = arc.dest
        if node == cycle_nodes[0]:
            print('+1 cycle (used %d arcs)' % len(cycle_arcs))
            last_arcs = []
            cycles.append(
                GraphCycleB(
                    segments=[a.bp for a in cycle_arcs]
                    #segments=cycle_bp
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
                for bp in cycle.segments:
                    axessteps[2].plot(bp[[0,-1],0], bp[[0,-1],1], 'o', c='g', mew=0)
                    plot_directed_bezier_curve(
                        axessteps[2], 
                        bp,
                        color='g',
                        arrow_length=0.033
                    )
                for ax in axessteps:
                    ax.set_xlim(-1,1)
                    ax.set_ylim(-1,1)
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
cycle_color = [cmap(1.*i/len(cycles)) for i in range(len(cycles))]

if not show_steps:
    for ic, cycle in enumerate(cycles):
        for bp in cycle.segments:
            axes[1].plot(bp[[0,-1],0], bp[[0,-1],1], 'o', c=cycle_color[ic], mew=0)
            plot_directed_bezier_curve(
                axes[1],
                bp,
                color=cycle_color[ic],
                arrow_length=0.033
            )
    for ax in axes:
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        ax.set_aspect('equal')
    plt.show()
####################################################




####################################################
# MAKE FACES
faces = make_faces(cycles, n=100)
face_color = [cmap(1.*i/len(faces)) for i in range(len(faces))]

fig, ax = plt.subplots()

for iface, face in enumerate(faces):
    for cycle in [face.outer] + face.inner:
        for bp in cycle.segments:
            ax.plot(bp[[0,-1],0], bp[[0,-1],1], 'o', c=cycle_color[iface], mew=0)
            plot_directed_bezier_curve(
                ax,
                bp,
                color=cycle_color[iface],
                arrow_length=0.033
            )
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_aspect('equal')
plt.show()
####################################################

