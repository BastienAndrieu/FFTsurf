import matplotlib.pyplot as plt
import numpy as np
from math import cos, sin, atan2, pi
import sys
sys.path.append('/d/bandrieu/GitHub/Code/Python/')
import lib_bezier as lbez

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
    #
    def curve(self, n=100, tension=1):
        return curved_arc([self.orig.co, self.dest.co], self.angles, n, tension)
    #
    def plot(self, ax, n=100, tension=1, color='k', arrow_length=0):
        aabb = plot_curved_arc(ax, [self.orig.co, self.dest.co], self.angles, n, tension, color, arrow_length)
        return aabb
#############################################################
class GraphArcB:
    def __init__(self, orig=None, dest=None, bp=None):
        self.orig = orig
        self.dest = dest
        self.bp = bp
        self.compute_angles()
        return
    #
    def compute_angles(self):
        if self.bp is not None:
            self.angles = [
                np.arctan2(v[1], v[0])
                for v in [self.bp[1] - self.bp[0], self.bp[-1] - self.bp[-2]]
            ]
        return
    #
    def curve(self, n=100, tension=1):
        t = np.linspace(0,1,n)
        return lbez.eval_bezier_curve(self.bp, t)
    #
    def plot(self, ax, n=100, tension=1, color='k', arrow_length=0):
        #xy = self.curve(n, tension)
        #plot_directed_curve(ax, xy, color, arrow_length)
        xy = plot_directed_bezier_curve(ax, self.bp, n, color, arrow_length)
        xymin = np.amin(xy, axis=0)
        xymax = np.amax(xy, axis=0)
        return (xymin[0], xymax[0], xymin[1], xymax[1])


#############################################################
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
class GraphCycleB:
    def __init__(self, segments=None):
        self.segments = segments
        return
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
def curved_arc_end_tangents(points, angles, tension=1):
    l = np.hypot(points[1][0] - points[0][0], points[1][1] - points[0][1])
    m0 = tension*l*np.array([cos(angles[0]), sin(angles[0])])
    m1 = tension*l*np.array([cos(angles[1]), sin(angles[1])])
    return m0, m1
#############################################################
def curved_arc(points, angles, n=100, tension=1):
    """
    Hermite cubic curve
    """
    p0 = np.array([points[0][0], points[0][1]])
    p1 = np.array([points[1][0], points[1][1]])
    m0, m1 = curved_arc_end_tangents(points, angles, tension)
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
def curved_arc_bezier_control_points(points, angles, tension=1):
    p0 = np.array([points[0][0], points[0][1]])
    p1 = np.array([points[1][0], points[1][1]])
    m0, m1 = curved_arc_end_tangents(points, angles, tension)
    return np.array([p0, p0+m0/3., p1-m1/3., p1])
#############################################################
def plot_directed_curve(ax, p, color='k', arrow_length=0):
    ax.plot(p[:,0], p[:,1], color=color, lw=1)
    if arrow_length > 0:
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
    return 
#############################################################
def plot_curved_arc(ax, points, angles, n=100, tension=1, color='k', arrow_length=0):
    p = curved_arc(points, angles, n, tension)
    plot_directed_curve(ax, p, color, arrow_length)
    xymin = np.amin(p, axis=0)
    xymax = np.amax(p, axis=0)
    return (xymin[0], xymax[0], xymin[1], xymax[1])
#############################################################
def plot_directed_bezier_curve(ax, bp, n=100, color='k', arrow_length=0):
    t = np.linspace(0,1,n)
    xy = lbez.eval_bezier_curve(bp, t)
    if ax is not None: plot_directed_curve(ax, xy, color, arrow_length)
    return xy
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
def make_faces(cycles, n=100, tension=1):
    sys.path.append('/d/bandrieu/GitHub/Code/Python/')
    from lib_compgeom import is_inside_polygon
    import lib_brep
    
    # make cycle polygons
    polys = []
    for cycle in cycles:
        poly = np.empty((0,2), dtype=float)
        if isinstance(cycle, GraphCycleB):
            for bp in cycle.segments:
                xy = plot_directed_bezier_curve(None, bp, n)
                poly = np.vstack([poly, xy[:-1]])
        else:
            m = len(cycle.points)
            x = [p[0] for p in cycle.points]
            y = [p[1] for p in cycle.points]
            for i in range(m):
                j = (i+1)%m
                curve = curved_arc(
                    points=[
                        [x[i], y[i]],
                        [x[j], y[j]]
                    ],
                    angles=cycle.angles[i],
                    n=n,
                    tension=tension
                )
                poly = np.vstack([poly, curve[:-1]])
        polys.append(poly)

    # get inclusion relationships
    ascendants = []
    for cycle in cycles:
        ascendants.append([])
    ncycles = len(cycles)
    for i in range(ncycles-1):
        for j in range(i+1,ncycles):
            j_inside_i = is_inside_polygon(
                point=polys[j][0],
                verts=polys[i],
                edges=None,
                BIG=100
            )
            if j_inside_i:
                ascendants[j].append(i)
            else:
                i_inside_j = is_inside_polygon(
                    point=polys[i][0],
                    verts=polys[j],
                    edges=None,
                    BIG=100
                )
                if i_inside_j:
                    ascendants[i].append(j)
        

    # make nesting tree
    level = [0]*ncycles
    maxlevel = 0
    children = []
    for i in range(ncycles):
        maxlevel = max(maxlevel, len(ascendants[i]))
        children.append([])
    for lv in range(maxlevel+1):
        for i in range(ncycles):
            if len(ascendants[i]) == lv:
                for j in ascendants[i]:
                    if len(ascendants[j]) == lv-1:
                        children[j].append(i)

    # make faces
    faces = []
    for i in range(ncycles):
        if len(ascendants[i])%2 == 0:
            faces.append(
                lib_brep.Face(
                    outer=cycles[i],
                    inner=[cycles[j] for j in children[i]],
                    index=len(faces)
                )
            )
    return faces
