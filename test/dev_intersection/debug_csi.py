import bpy
import numpy as np
import sys
sys.path.append('/stck/bandrieu/Bureau/Python/mylibs/')
#sys.path.append('/home/bastien/Bureau/Python/')
import my_lib1 as myl
import my_lib as mylb

pth = '/d/bandrieu/GitHub/FFTsurf/test/dev_intersection/'
#pth = '/home/bastien/GitHub/FFTsurf/test/dev_intersection/'

###########################################################

bc = myl.readCoeffs1( pth + 'debugbsi_c.bern' )
bs = myl.readCoeffs2( pth + 'debugbsi_s.bern' )



## set up a blank scene
mylb.blankScene()

## add Bezier control curve
mylb.addPolyline(
    bc,
    clr=(0,0,0),
    thickness=5.e-4,
    layer=0 )

## add Bezier control surface
# vertices
m = bs.shape[0]
n = bs.shape[1]
verts = np.zeros( (m*n,3), dtype=float )
for j in range (n):
    for i in range(m):
        verts[m*j+i] = bs[i,j]

# faces
faces = []
for j in range(n-1):
    for i in range(m-1):
        faces.append( [m*j+i, m*j+i+1, m*j+i+m+1, m*j+i+m] )
        
# create mesh and object
mesh = bpy.data.meshes.new("surf")
object = bpy.data.objects.new("surf",mesh)
    
# set mesh location
object.location = [0.0,0.0,0.0]
bpy.context.scene.objects.link(object)

# create mesh from python data
mesh.from_pydata(verts,[],faces)

# material
mat = mylb.surfaceMaterial("surf_mat",color=(1,1,1),alpha=1)
mylb.setMaterial(object,mat)


## set camera
for ob in bpy.context.scene.objects:
    if ob.type == 'MESH':
        ob.select = True
    else:
        ob.select = False
        
bpy.ops.view3d.camera_to_view_selected()

