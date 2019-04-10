import bpy
import bmesh


obj = bpy.context.active_object

me = obj.data
# Get a BMesh representation
bm = bmesh.new()
bm.from_mesh(me)

# Triangulate
bmesh.ops.triangulate(bm, faces=bm.faces[:], quad_method=0, ngon_method=0)

# Finish up, write the bmesh back to the mesh
bm.to_mesh(me)
bm.free()


###########################
# export mesh data
pth = '/d/bandrieu/GitHub/FFTsurf/test/optimmesh/'

f = open(pth + 'test_mesh_xyz.dat', 'w')
for v in me.vertices:
    for x in v.co:
        f.write('%.8s ' %(x))
    f.write('\n')
f.close()

f = open(pth + 'test_mesh_tri.dat', 'w')
for t in me.polygons:
    for v in t.vertices:
        f.write('%d ' %(v+1))
    f.write('\n')
f.close()
