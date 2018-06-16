import bpy
import numpy as np
import sys
sys.path.append('/home/bastien/Bureau/Python/')
import my_lib1 as myl
import my_lib as mylb


uvxyz = np.loadtxt('uv_xyz.dat')
pth = '/home/bastien/GitHub/FFTsurf/test/dev_intersection_simple_surface/'

## set up a blank scene
mylb.blankScene()

## add both surfaces
n = 100
clr = [ (1.0, 0.3, 0.1), (0.1, 0.3, 1.0) ]
for isurf in range(2):
    c = myl.readCoeffs2( pth + 'c_' + str(isurf+1) + '.cheb' )
    mylb.addPatch( c, n, clr=clr[isurf], alp=1)

## add intersection points
mylb.addVertices(
    uvxyz[:,4:7],
    clr=(1,1,1),
    radius=1e-2,
    layer=0 )

## set camera
for ob in bpy.context.scene.objects:
    if ob.type == 'MESH':
        ob.select = True
    else:
        ob.select = False
        
bpy.ops.view3d.camera_to_view_selected()

