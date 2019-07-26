import sys
sys.path.append('/d/GitHub/FFTsurf/test/graph/')
from lib_curved_graph import *

from math import radians

#############################################################
#############################################################

nodes_xy = [
    (   -1,    -1), # 0
    (    1,    -1), # 1
    (    1,     1), # 2
    (   -1,     1), # 3
    (-0.25,    -1), # 4
    (    1,     0), # 5
    ( 0.65, -0.35), # 6
    ( 0.55, -0.85), # 7
    (   -1,  -0.5), # 8
    (    1,   0.5), # 9
    (-0.75,   0.8), # 10
    ( -0.2,   0.1), # 11
    (-0.65,  0.65), # 12
    (-0.35,  0.25), # 13
    (  0.2,   0.7), # 14
    (  0.7,  0.75), # 15
    ( 0.25, -0.65)  # 16
]

arc_end = [
    ( 3,  8), # 0
    ( 8,  0), # 1
    ( 0,  4), # 2
    ( 4,  1), # 3
    ( 1,  5), # 4
    ( 5,  9), # 5
    ( 9,  2), # 6
    ( 2,  3), # 7
    ( 6,  7), # 8
    ( 7, 16), #( 7,  6), # 9
    (10, 11), # 10
    (13, 12), #(12, 13), # 11
    (12, 13), #(13, 12), # 12
    (11, 10), # 13
    (14, 15), # 14
    (15, 14), # 15
    ( 5,  4), # 16
    ( 8,  9), # 17
    (16,  6)  # 18
]


arc_ang = [
    (270, 270), # 0
    (270, 270), # 1
    (  0,   0), # 2
    (  0,   0), # 3
    ( 90,  90), # 4
    ( 90,  90), # 5
    ( 90,  90), # 6
    (180, 180), # 7
    (350, 190), #(  0, 180), # 8
    (190, 120), #(180,   0), # 9
    ( 45, 195), #(225,  15), #(225,  45), # 10
    ( 45, 225), #( 45, 225), #( 45, 225), #(315, 135), # 11
    (270, 270),#( 90,  90),#(225,  45), #(135, 315), # 12
    (195,  45), #( 15, 225), #( 45, 225), # 13
    ( 90, 270), # 14
    (270,  90), # 15
    (185, 260), # 16
    (  0,   0), # 17
    (120, 350)  # 18    
]






####################################################
# EXPORT EMBEDDED GRAPH (for make_wire.py)
####################################################
f = open('test_demo.dat', 'w')

f.write('%d\n' % len(arc_end))
for a, b in arc_end:
    f.write('%d %d\n' % (a+1, b+1))
for a, b in arc_ang:
    f.write('%s %s\n' % (radians(float(a)), radians(float(b))))

f.write('%d\n' % len(nodes_xy))
for x, y in nodes_xy:
    f.write('%s %s\n' % (float(x), float(y)))

f.close()



####################################################
# 
####################################################
nodes = []
for xy in nodes_xy:
    nodes.append(
        GraphNode(
            co=[float(x) for x in xy],
            ingoing=[],
            outgoing=[]
        )
    )


narc = len(arc_end)
arcs = []
for iarc in range(narc):
    arcs.append(
        GraphArc(
            orig=nodes[arc_end[iarc][0]],
            dest=nodes[arc_end[iarc][1]],
            angles=[radians(float(a)) for a in arc_ang[iarc]]
        )
    )

####################################################
# EXPORT EMBEDDED GRAPH (for TikZ)
####################################################
f = open('graph_bezier_control_points.dat', 'w')

for arc in arcs:
    b = curved_arc_bezier_control_points(
        [arc.orig.co, arc.dest.co],
        arc.angles,
        tension=1
    )
    for x, y in b:
        f.write('%s, %s, ' % (x, y))
    f.write('\n')
f.close()

####################################################
# COMPUTE INITIAL NODE -> ARC INCIDENCE
# (distinguish outgoing/ingoing arcs)
####################################################
for arc in arcs:
    arc.orig.outgoing.append(arc)
    arc.dest.ingoing.append(arc)

####################################################
# PLOT EMBEDDED GRAPH
####################################################
mrg = 5e-2
lim = 1 + mrg

fig, ax = plt.subplots()
plot_graph(ax, nodes, arcs)
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect('equal')
plt.show()



"""
sys.path.append('/d/bandrieu/GitHub/Code/Python/')
from lib_bezier import eval_bezier_curve
fig, ax = plt.subplots()

for arc in arcs:
    b = curved_arc_bezier_control_points(
        [arc.orig.co, arc.dest.co],
        arc.angles,
        tension=1
    )
    p = eval_bezier_curve(b, np.linspace(0,1,100))
    ax.plot(p[:,0], p[:,1], 'k-')

ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect('equal')
plt.show()
"""







"""
fig, ax = plt.subplots()

for i, arc in enumerate(arc_end):
    print i, arc
    x = [float(nodes_xy[v][0]) for v in arc]
    y = [float(nodes_xy[v][1]) for v in arc]
    print nodes_xy[arc[0]], nodes_xy[arc[1]]
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    ax.arrow(
            x[0], y[0], dx, dy,
            color='k',
            length_includes_head=True
        )
    ax.text(0.5*(x[0]+x[1]), 0.5*(y[0]+y[1]), str(i), color='b')

for i, xy in enumerate(nodes_xy):
    ax.text(xy[0], xy[1], str(i), color='r')
    
lim = 1.05
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect('equal')
plt.show()
"""
