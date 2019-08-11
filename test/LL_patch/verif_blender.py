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
pthin = ROOT + 'GitHub/FFTsurf/test/LL_patch/'
f = open(pthin + 'input.dat', 'r')
nplanes = int(f.readline())
vxyz = Vector([float(x) for x in f.readline().split()])
cc_center = []
for i in range(nplanes):
    cc_center.append(Vector([float(x) for x in f.readline().split()]))
normal = []
for i in range(nplanes):
    normal.append(Vector([float(x) for x in f.readline().split()]))
corner = []
for i in range(nplanes):
    corner.append(Vector([float(x) for x in f.readline().split()]))
f.close()
################################################################


################################################################
# ADD SPHERE
rho = 0.0
for xyz in corner:
    rho += (xyz - vxyz).length
rho /= float(nplanes)
print('rho = %s' % rho)

bpy.ops.mesh.primitive_uv_sphere_add(
    location=vxyz,
    size=rho,
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

B = numpy.loadtxt(pthin + 'B.dat')
sphere.rotation_euler = Matrix(B).to_euler()
################################################################


################################################################
# ADD ARCS
ncc = 100
ucc = numpy.linspace(0,1,ncc)
for i in range(nplanes):
    j = (i-1)%nplanes
    r1 = corner[j] - cc_center[i]
    r2 = corner[i] - cc_center[i]

    r1perp = r1.cross(normal[i])
    angle = numpy.arctan2(r2.dot(r1perp), r2.dot(r1))%(2*numpy.pi)

    xyz = [cc_center[i] + r1*numpy.cos(ucc[k]*angle) + r1perp*numpy.sin(ucc[k]*angle)
           for k in range(ncc)]

    arc = lbu.pydata_to_polyline(
        xyz,
        name='arc_'+str(i),
        thickness=1e-3
    )
################################################################

"""
################################################################
# ADD CHARACTERISTIC CIRCLES
for i in range(nplanes):
    ri = numpy.sqrt(rho**2 - (cc_center[i] - vxyz).length_squared)
    rot = Matrix(complete_orthonormal_matrix(normal[i], i=2)).transposed()
    bpy.ops.mesh.primitive_circle_add(
        location=cc_center[i],
        rotation=rot.to_euler(),
        radius=ri,
        vertices=100,
        fill_type='NOTHING'
    )
################################################################
"""

"""
################################################################
polygon = lbu.pydata_to_mesh(
    verts=corner,
    faces=[],
    edges=[(i, (i+1)%nplanes) for i in range(nplanes)],
    name='polygon'
)
################################################################
"""

################################################################
# LOAD AND ADD LL-PATCH
c = lcheb.read_polynomial2(pthin + 'surf.cheb')
m = 100
u = numpy.linspace(-1,1,m)

xyz = chebgrid2d(u, u, c)
mverts, mfaces = lbu.tensor_product_mesh_vf(xyz[0], xyz[1], xyz[2])
LLpatch = lbu.pydata_to_mesh(
    mverts,
    mfaces,
    name='LLpatch'
)
lbe.set_smooth(LLpatch)
#LLpatch.data.materials.append(mat)
################################################################


################################################################
# ADJUST CAMERA
bpy.ops.object.select_all(action='DESELECT')
sphere.select = True
bpy.ops.view3d.camera_to_view_selected()
cam.data.angle += numpy.radians(10.)
################################################################
