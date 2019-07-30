import bpy
import numpy

import sys
sys.path.append('/d/bandrieu/GitHub/Code/Python/')
import lib_blender_edit as lbe
import lib_blender_util as lbu
import lib_color as lco

pthin = '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_MAT/brepmesh_eos/'

#################################################
## Set Scene
scene = bpy.context.scene
lbu.clear_scene(meshes=True, lamps=False, cameras=False)

lbu.set_scene(
    resolution_x=1024,
    resolution_y=768,
    resolution_percentage=100,
    alpha_mode='SKY',
    horizon_color=[1,1,1],
    light_samples=4,
    use_environment_light=True,
    environment_energy=0.1,
    environment_color='PLAIN'
)
#################################################


#################################################
# CHECKER TEXTURE
imgchecker = bpy.data.images.load(filepath='/d/bandrieu/GitHub/These/memoire/figures/images/checker6.png')
texchecker = bpy.data.textures.new('texture_checker', 'IMAGE')
texchecker.image = imgchecker
#################################################


#################################################
# ADD FACES
nfaces = 10
for iface in range(nfaces):
    strf = format(iface+1, '02')
    tri = numpy.loadtxt(pthin + 'tri_' + strf + '.dat', dtype=int)-1
    xyz = numpy.loadtxt(pthin + 'xyz_' + strf + '.dat', dtype=float)
    uv = numpy.loadtxt(pthin + 'uv_' + strf + '.dat', dtype=float)

    verts = [[x for x in p] for p in xyz]
    faces = [[int(v) for v in t] for t in tri]

    obj = lbu.pydata_to_mesh(
        verts,
        faces,
        name='face_'+strf
    )

    lbe.set_smooth(obj)

    # material
    mat = bpy.data.materials.new('mat_face_'+strf)
    mat.diffuse_color = (1,1,1)#color[iface]
    mat.diffuse_intensity = 1
    mat.specular_intensity = 0

    slot = mat.texture_slots.add()
    slot.texture = texchecker
    slot.texture_coords = 'UV'
    slot.blend_type = 'MULTIPLY'
    slot.diffuse_color_factor = 0.1
    
    obj.data.materials.append(mat)

    # UV coordinates
    scene.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT')
    print('   unwrap UVs...')
    bpy.ops.uv.unwrap(
        method='ANGLE_BASED',
        fill_holes=True,
        correct_aspect=True,
        use_subsurf_data=False,
        margin=0.001
    )
    bpy.ops.object.mode_set(mode='OBJECT')
    uvlayer = obj.data.uv_layers.active
    print('   edit UVs...')
    for i, f in enumerate(obj.data.polygons):
        for j in range(f.loop_total):
            k = f.loop_start + j
            for l in range(2):
                uvlayer.data[k].uv[l] = 0.5*(uv[f.vertices[j],l] + 1.0)
#################################################

#################################################
# ADJUST CAMERA
bpy.ops.object.select_all(action='SELECT')
bpy.ops.view3d.camera_to_view_selected()
#################################################
