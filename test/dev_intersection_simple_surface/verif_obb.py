import bpy
import numpy as np
import sys
sys.path.append('/home/bastien/Bureau/Python/')
import my_lib1 as myl
import my_lib as mylb



###################################################

def addOBB( ctr, rng, axe, clr, alp=0.5, layer=0 ):
    faces = [ [0,2,3,1],
              [0,4,6,2],
              [0,1,5,4],
              [1,3,7,5],
              [2,6,7,3],
              [4,5,7,6] ]
    xyz = np.zeros( (8,3) )
    for k in range(2):
        for j in range(2):
            for i in range(2):
                xyz[4*k + 2*j + i] = ctr + np.matmul( -np.power( -1.0, np.array([i,j,k],dtype=float) ) * rng, np.transpose( axe ) )

    #create mesh and object
    mesh = bpy.data.meshes.new("obb")
    object = bpy.data.objects.new("obb",mesh)
    
    #set mesh location
    object.location = [0.0,0.0,0.0]
    bpy.context.scene.objects.link(object)
    
    #create mesh from python data
    mesh.from_pydata( xyz, [], faces )
    mesh.update( calc_edges=True )

    # material
    mat = mylb.surfaceMaterial("mat_obb",clr,alp)
    mylb.setMaterial(object,mat)
    
    return

###################################################




pth = '/home/bastien/GitHub/FFTsurf/test/dev_intersection_simple_surface/'

bc = myl.readCoeffs1( pth + 'region_c_bezier.bern' )
bs = myl.readCoeffs2( pth + 'region_s_bezier.bern' )



boxc = np.loadtxt( pth + 'xyzbox_c.dat' )
boxs = np.loadtxt( pth + 'xyzbox_s.dat' )






## set up a blank scene
mylb.blankScene()


## add Bezier control curve
mylb.addPolyline(
    bc,
    clr=(0,0,0),
    thickness=max(1.e-4,0.01*np.amin(boxc[1])),#1.e-3,
    layer=0 )

## add curve OBB
addOBB(
    boxc[0,:],
    boxc[1,:],
    np.transpose(boxc[2:5,:]),
    clr=(0,0,1),
    alp=0.3 )


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


## add surface OBB
addOBB(
    boxs[0,:],
    boxs[1,:],
    np.transpose(boxs[2:5,:]),
    clr=(1,0,0),
    alp=0.3 )

## set camera
for ob in bpy.context.scene.objects:
    if ob.type == 'MESH':
        ob.select = True
    else:
        ob.select = False
        
bpy.ops.view3d.camera_to_view_selected()

