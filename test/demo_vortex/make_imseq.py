import bpy
import sys
sys.path.append('/stck/bandrieu/Bureau/Python/mylibs/')
import mycolors as myc
import my_lib1 as myl
import my_lib as mylb
import math
import numpy as np

args = sys.argv

if len(args) < 4:
    istart = 1
else:
    istart = int(args[3])

if len(args) < 5:
    iend = 1
else:
    iend = int(args[4])
    

palette = [
    [0.5502,    0.8571,    1.0000],
    [0.6342,    0.9873,    1.0000],
    [0.7074,    1.0000,    0.8390],
    [0.8960,    1.0000,    0.7081],
    [1.0000,    0.8842,    0.6134],
    [1.0000,    0.7593,    0.6060]]
# cc
palette = myc.cc_sv(palette,1.3,0.9)


pth = '/d/bandrieu/GitHub/FFTsurf/test/demo_vortex/'


edge_thickness = 7.5e-4 
vertex_radius = 2.0*edge_thickness


for instant in range(istart,iend+1):
    mylb.blankScene()
    s_instant = format(instant, '03')
    print ('instant : '+s_instant)

    # Add surface mesh
    for isurf in range(6):
        #print (pth + 'snap_'+s_instant+'c_'+format(isurf+1,'03')+'.cheb')
        c = mylb.readCoeffs(pth + 'data/snap_'+s_instant+'c_'+format(isurf+1,'03')+'.cheb')
        mylb.addPatch(c, 200, palette[isurf])
    
    # add BRep vertices
    verts = np.loadtxt(pth + 'data/verts_'+s_instant+'.dat')
    mylb.addVertices(verts, [0.0,0.0,0.0], vertex_radius, 0)

    # add BRep edges
    fe = open(pth + 'data/edges_'+s_instant+'.dat', 'r')
    ne = int(fe.readline())
    for ie in range(ne):
        npts = int(fe.readline())
        exyz = []
        for ip in range(npts):
            exyz.append([float(a) for a in fe.readline().split()])
        mylb.addPolyline(exyz, [0.0,0.0,0.0], edge_thickness, 0)
    fe.close()

    # Environment lighting
    scene = bpy.context.scene
    scene.world.light_settings.samples = 8
    scene.world.light_settings.use_environment_light = True
    scene.world.light_settings.environment_energy = 0.15
    scene.world.light_settings.environment_color = 'PLAIN'
    #scene.world.light_settings.use_ambient_occlusion = True
    #scene.world.light_settings.use_ao_blend_type = 'MULTIPLY'
    

    scene.camera.rotation_mode = 'XYZ'
    bpy.data.cameras["Camera"].lens = 52.0
    scene.camera.rotation_euler = (math.radians(107.245), math.radians(0.805), math.radians(-5.389))
    scene.camera.location = (0.52821, -0.48525, 0.28909)

    # Lamp1
    # Create new lamp datablock
    lamp_data = bpy.data.lamps.new(name="Lamp1", type='POINT')
    lamp_data.energy = 1.0
    lamp_data.shadow_method = 'RAY_SHADOW'#
    lamp_data.shadow_ray_samples = 12
    lamp_data.shadow_soft_size = 2.0

    # Create new object with our lamp datablock
    lamp_object = bpy.data.objects.new(name="Lamp1", object_data=lamp_data)

    # Link lamp object to the scene so it'll appear in this scene
    bpy.context.scene.objects.link(lamp_object)

    # Place lamp to a specified location
    #lamp_object.location = [-0.52197, 1.78084, 1.82178]
    lamp_object.location = [-0.47971, -1.17763, 0.38481]

    # And finally select it make active
    lamp_object.select = True
    bpy.context.scene.objects.active = lamp_object


    # Lamp2
    # Create new lamp datablock
    lamp_data2 = bpy.data.lamps.new(name="Lamp1", type='POINT')
    lamp_data2.energy = 0.7
    lamp_data2.shadow_method = 'NOSHADOW'#

    # Create new object with our lamp datablock
    lamp_object2 = bpy.data.objects.new(name="Lamp2", object_data=lamp_data2)

    # Link lamp object to the scene so it'll appear in this scene
    bpy.context.scene.objects.link(lamp_object2)

    # Place lamp to a specified location
    lamp_object2.location = [1.94414, 2.57809, 2.28225]

    # And finally select it make active
    lamp_object2.select = True
    bpy.context.scene.objects.active = lamp_object2
    

    # Render properties
    scene.render.resolution_x = 800
    scene.render.resolution_y = scene.render.resolution_x
    scene.render.resolution_percentage = 100

    bpy.data.worlds["World"].horizon_color = (1.0,1.0,1.0)

    scene.render.image_settings.file_format='JPEG'
    scene.render.image_settings.compression = 90
    scene.render.use_border = False
    bpy.data.scenes['Scene'].render.filepath = pth + 'imseq/im_' + s_instant
    bpy.ops.render.render( write_still=True ) 

