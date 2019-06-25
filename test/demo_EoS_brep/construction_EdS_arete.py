import bpy

import sys
sys.path.append('/d/bandrieu/GitHub/Code/Python/')
import lib_blender_util as lbu
import lib_blender_edit as lbe
import lib_color as lco

sys.path.append('/stck/bandrieu/Bureau/Python/mylibs/')
import my_lib as myl

import numpy
from numpy.polynomial.chebyshev import chebgrid2d





################################
## Set Scene
scene = bpy.context.scene
lbu.clear_scene(meshes=True, lamps=True, cameras=False)

lbu.set_scene(resolution_x=800,
              resolution_y=800,
              resolution_percentage=100,
              alpha_mode='SKY',
              horizon_color=[1,1,1],
              light_samples=4,#16,
              use_environment_light=True,
              environment_energy=0.3,
              environment_color='PLAIN')

## Set Lighting
lamp = lbu.add_point_light(name="lamp",
                           energy=1.2,
                           shadow_method='RAY_SHADOW',
                           shadow_ray_samples=16,
                           shadow_soft_size=2.0,
                           location=[3.75,1.65,3.20])
for i in range(4):
    lamp.layers[i] = True

## Set Camera
cam = scene.camera
for i in range(4):
    cam.layers[i] = True
cam.location = [2.102,1.798,1.104]
cam.rotation_euler = numpy.radians(numpy.array([66.7,0.778,132.2]))
cam.data.angle = numpy.radians(37.72)
################################

pthin = '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/'
pthout = '/d/bandrieu/GitHub/These/memoire/figures/data/BRep/'


clf = numpy.loadtxt(pthout + 'face_color.dat')
cle = lco.cc_hsv(clf, fs=0.9, fv=1.1)
clf = lco.cc_hsv(clf, fs=1.2, fv=1.0)
nclf = len(clf)


ifaceRL = [4,6]
iedge = 21

# materials
mat_surf = bpy.data.materials.new('mat_surf')
mat_surf.diffuse_color = (1,1,1)

mat_faceR = bpy.data.materials.new('mat_faceR')
mat_faceR.diffuse_color = cle[(ifaceRL[0]-1)%nclf]

mat_faceL = bpy.data.materials.new('mat_faceL')
mat_faceL.diffuse_color = cle[(ifaceRL[1]-1)%nclf]

#mat_edge = bpy.data.materials.new('mat_edge')
#mat_edge.diffuse_color = clf[(iedge-1)%nclf]

diffuse_intensity = 1
specular_intensity = 0.2
specular_hardness = 30

for mat in [mat_faceL, mat_faceR]:
    mat.use_transparency = False#
    mat.alpha = 0.5


# load initial interface mesh
bpy.ops.import_scene.obj(filepath=pthin+'mesh/mesh.obj',
                         axis_forward='Y', axis_up='Z')
surf = bpy.data.objects['mesh']
surf.name = 'surface'
# set smooth shading
lbe.set_smooth(surf)

# set materials
for slot in surf.material_slots:
    iface = int(slot.material.name[3:])
    mat = slot.material
    if iface == ifaceRL[0]:
        mat.diffuse_color = clf[(ifaceRL[0]-1)%nclf]
    elif iface == ifaceRL[1]:
        mat.diffuse_color = clf[(ifaceRL[1]-1)%nclf]
    else:
        mat.diffuse_color = (1,1,1)

# load R/L faces
ptheos = pthin + 'brepmesh_eos/'
facesRL = []
suffixe = 'RL'
for j, iface in enumerate(ifaceRL):
    strf = format(iface, '03')
    tri = numpy.loadtxt(ptheos + 'tri_' + strf + '.dat', dtype=int)-1
    xyz = numpy.loadtxt(ptheos + 'xyz_' + strf + '.dat', dtype=float)

    verts = [[x for x in p] for p in xyz]
    faces = [[int(v) for v in t] for t in tri]
    name = 'eos_face' + suffixe[j]

    obj = lbu.pydata_to_mesh(verts,
                             faces,
                             name=name )
    lbe.set_smooth(obj)

    # material
    if j == 0:
        obj.data.materials.append(mat_faceR)
    elif j == 1:
        obj.data.materials.append(mat_faceL)
   
    facesRL.append(obj)



# load edge
m = 200
n = 50
u = numpy.linspace(-1.0,1.0,m)
v = numpy.linspace(-1.0,1.0,n)
c = myl.readCoeffs(ptheos + 'c_' + format(iedge,'03') + '.cheb')
edge = myl.addTensorProductPatch(
    chebgrid2d(u, v, c[:,:,0]),
    chebgrid2d(u, v, c[:,:,1]),
    chebgrid2d(u, v, c[:,:,2]),
    name='edge',
    location=(0,0,0),
    smooth=True,
    color=cle[(iedge-1)%nclf],
    alpha=0
)


# ajust camera
bpy.ops.object.select_all(action='DESELECT')
edge.select = True
bpy.ops.view3d.camera_to_view_selected() # set camera to fit 'edge' in screen
cam.data.angle += numpy.pi/180.0 # increase FOV angle by 1 degree


# ajust all materials
for mat in bpy.data.materials:
    mat.diffuse_intensity = diffuse_intensity
    mat.specular_intensity = specular_intensity
    mat.specular_hardness = specular_hardness
    mat.use_transparency = False
