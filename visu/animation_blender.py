import bpy
from mathutils import Vector
import numpy
import sys
sys.path.append('/stck/bandrieu/Bureau/Python/mylibs/')
import my_lib

########################################################################

def insert_keyframe(fcurves, frame, values):
    for fcu, val in zip(fcurves, values):
        fcu.keyframe_points.insert(frame, val, {'FAST'})
        
########################################################################

cas = 'jouke'
pth = '/d/bandrieu/GitHub/FFTsurf/cases/' + cas + '/output/'

reprises = numpy.loadtxt(pth + 'reprises.dat', dtype=int)
reprises = numpy.hstack([reprises, 62])

########################################################################










########################################################################
my_lib.blankScene()
scn = bpy.context.scene


xyzini = numpy.loadtxt(pth + 'pos_001.dat', dtype=float)
tri = numpy.loadtxt(pth + 'connect_01.dat', dtype=int)-1
ref = numpy.loadtxt(pth + 'faceref_01.dat', dtype=int)

verts = []
for x in xyzini:
    verts.append(list(x))
faces = []
for t in tri:
    faces.append([int(t[0]), int(t[1]), int(t[2])])

name = "mesh_rep_00"
#create mesh and object
msh = bpy.data.meshes.new(name)
obj = bpy.data.objects.new(name, msh)
#set mesh location
obj.location = [0,0,0]
bpy.context.scene.objects.link(obj)
#create mesh from python data
msh.from_pydata(verts,[],faces)
msh.update(calc_edges=True)
# show mesh as smooth
obj.data.use_auto_smooth = 1
mypolys = msh.polygons
for q in msh.polygons:
    q.use_smooth = True

####################################################
action = bpy.data.actions.new("MeshAnimation_00")
msh.animation_data_create()
msh.animation_data.action = action

data_path = "vertices[%d].co"
vec_z = Vector((0.0, 0.0, 1.0))

frames = 0, 25, 50, 75, 100
values = 0.0, 3.0, 5.0, 3.0, 0.0

"""
for v in msh.vertices:
    fcurves = [action.fcurves.new(data_path % v.index, i) for i in range(3)]
    co_rest = v.co

    for t, value in zip(frames, values):
        co_kf = co_rest + value * vec_z
        insert_keyframe(fcurves, t, co_kf)  
"""

fcurves = [[action.fcurves.new(data_path % j, i) for i in range(3)]
           for j in range(len(msh.vertices))]

for instant in range(1,reprises[1]):
    xyz = numpy.loadtxt(pth + 'pos_' + format(instant,'03') + '.dat', dtype=float)
    for iv, v in enumerate(msh.vertices):
        insert_keyframe(fcurves[iv], instant, Vector(xyz[iv]))

