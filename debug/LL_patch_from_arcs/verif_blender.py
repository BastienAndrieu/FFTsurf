ROOT = '/home/bastien/'#'/d/bandrieu/'#

import bpy
from mathutils import Vector, Matrix

import numpy
from numpy.polynomial.chebyshev import chebgrid2d, chebval2d

import sys
sys.path.append(ROOT + 'GitHub/Code/Python/')
import lib_blender_util as lbu
import lib_blender_edit as lbe
import lib_chebyshev as lcheb
from lib_compgeom import complete_orthonormal_matrix


################################################################
def read_Vector(f):
    return Vector([float(x) for x in f.readline().split()])
################################################################


################################################################
# CHECKER TEXTURE
imgchecker = bpy.data.images.load(filepath=ROOT + 'GitHub/These/memoire/figures/images/checker6.png')
texchecker = bpy.data.textures.new('texture_checker', 'IMAGE')
texchecker.image = imgchecker

################################################################


################################################################
# SET SCENE
scene = bpy.context.scene
lbu.clear_scene(meshes=True, lamps=False, cameras=False)

resx = 1024
resy = resx
lbu.set_scene(
    resolution_x=resx,
    resolution_y=resy,
    resolution_percentage=100,
    alpha_mode='TRANSPARENT',
    horizon_color=(1,1,1),
    light_samples=4,
    use_environment_light=True,
    environment_energy=0.3,
    environment_color='PLAIN'
)
scene.render.engine = 'BLENDER_RENDER'

lamp = bpy.data.objects['Lamp']
cam = scene.camera
################################################################


################################################################
# LOAD DATA
pthin = ROOT + 'GitHub/FFTsurf/debug/LL_patch_from_arcs/'

s = numpy.loadtxt(pthin + 's.dat')
Btmp = numpy.loadtxt(pthin + 'Btmp.dat')
p = numpy.loadtxt(pthin + 'p.dat')
p_2d = numpy.loadtxt(pthin + 'p_2d.dat')
p_OBB = numpy.loadtxt(pthin + 'p_OBB.dat')
B = numpy.loadtxt(pthin + 'B.dat')
lonlat_OBB = numpy.loadtxt(pthin + 'lonlat_OBB.dat')

nplanes = len(s)
################################################################


################################################################
# ADD SPHERE
bpy.ops.mesh.primitive_uv_sphere_add(
    size=1,
    segments=100,
    ring_count=100
)

sphere = bpy.data.objects['Sphere']
lbe.set_smooth(sphere)
"""sphere.draw_type = 'WIRE'
sphere.show_all_edges = True"""

mat = bpy.data.materials.new('mat_sphere')
mat.diffuse_color = (0.8,0.8,0.8)
mat.diffuse_intensity = 1
mat.specular_intensity = 0.346
mat.specular_hardness = 20
mat.emit = 0.4
mat.use_transparency = True
#mat.raytrace_transparency.fresnel = 2.
mat.alpha = 0.5

sphere.data.materials.append(mat)
sphere.show_transparent = True
################################################################

"""
################################################################
# ADD CORNERS
corners = lbu.pydata_to_mesh(
    verts=s,
    faces=[],
    edges=[(i, (i+1)%nplanes) for i in range(nplanes)],
    name='s'
)
################################################################
"""




################################################################
# ADD BASIS
euler_angles = Matrix(B).to_euler()
bpy.ops.object.empty_add(
    rotation=euler_angles,
    type='ARROWS'
)
basis = bpy.data.objects['Empty']
basis.name = 'B'
################################################################

"""
################################################################
# ADD LL PATCH
n = 40
u = lonlat_OBB[0][0] + lonlat_OBB[1][0]*numpy.linspace(-1, 1, n)
v = lonlat_OBB[0][1] + lonlat_OBB[1][1]*numpy.linspace(-1, 1, n)

x = numpy.outer(numpy.cos(u), numpy.cos(v))
y = numpy.outer(numpy.sin(u), numpy.cos(v))
z = numpy.outer(numpy.ones(n), numpy.sin(v))

verts, faces = lbu.tensor_product_mesh_vf(x, y, z)
LLpatch = lbu.pydata_to_mesh(
    verts=verts,
    faces=faces,
    edges=[],
    name='LLpatch'
)
LLpatch.rotation_euler = euler_angles
################################################################
"""

################################################################
# ADJUST CAMERA
bpy.ops.object.select_all(action='DESELECT')
sphere.select = True
bpy.ops.view3d.camera_to_view_selected()
cam.data.angle += numpy.radians(10.)
################################################################



################################################################
# ARCS
ncc = 10
ucc = numpy.linspace(0,1,ncc)
edges = [(i, i+1) for i in range(ncc-1)]

box_verts = [((-1)**(i+1), (-1)**(j+1), 0) for j in range(2) for i in range(2)]
box_edges = [(0, 1), (1, 3), (3, 2), (2, 0)]

f = open(pthin + 'arcs.dat', 'r')
vxyz = read_Vector(f)
nplanes = int(f.readline())
for i in range(nplanes):
    cc_center = read_Vector(f)
    normal = read_Vector(f)
    c1 = read_Vector(f)
    c2 = read_Vector(f)
    r1 = read_Vector(f)
    r1perp = read_Vector(f)
    print('r1perp = ', r1perp)
    angle = float(f.readline())
    minmax_lon = [float(x)/numpy.pi for x in f.readline().split()]
    minmax_lat = [float(x)/numpy.pi for x in f.readline().split()]

    ilayer = i + 1
    ui = angle*ucc
    arc = [cc_center - vxyz + r1*numpy.cos(ui[j]) + r1perp*numpy.sin(ui[j]) for j in range(ncc)]
    lon = [numpy.arctan2(xyz.dot(B[:,1]), xyz.dot(B[:,0]))/numpy.pi for xyz in arc]
    lat = [numpy.arcsin(xyz.dot(B[:,2])/xyz.length)/numpy.pi for xyz in arc]
    for j in range(ncc):
        print(ui[j])
        print(arc[j])
        print(j+1, lon[j], lat[j])
    print('\n\n')

    scene.layers[ilayer] = True
    scene.layers[ilayer-1] = False
    obj = lbu.pydata_to_mesh(
        verts=[(lon[j], lat[j], 0) for j in range(ncc)],
        faces=[],
        edges=edges,
        name='lonlat_arc_'+str(i)
    )

    box = lbu.pydata_to_mesh(
        verts=box_verts,
        faces=[],
        edges=box_edges,
        name='lonlat_box_'+str(i)
    )
    box.location.x = 0.5*(minmax_lon[0] + minmax_lon[1])
    box.location.y = 0.5*(minmax_lat[0] + minmax_lat[1])
    box.scale.x = 0.5*(minmax_lon[1] - minmax_lon[0])
    box.scale.y = 0.5*(minmax_lat[1] - minmax_lat[0])
################################################################
