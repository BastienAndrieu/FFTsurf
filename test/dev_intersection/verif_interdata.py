import bpy
import numpy as np
import sys
sys.path.append('/stck/bandrieu/Bureau/Python/mylibs/')
#sys.path.append('/home/bastien/Bureau/Python/')
import my_lib1 as myl
import my_lib as mylb
import my_colors as myc

##########################################################

pth = '/d/bandrieu/GitHub/FFTsurf/test/dev_intersection/'
#pth = '/home/bastien/GitHub/FFTsurf/test/dev_intersection/'

##########################################################

class intersection_curve:
    def __init__( self, endpoints, uvbox, uv, xyz ):
        self.endpoints = endpoints
        self.uvbox = uvbox
        self.uv = uv
        self.xyz = xyz

##########################################################


args = sys.argv
if len(args) < 4:
    nsurf = 2
else:
    nsurf = int(args[3])



points = np.loadtxt( pth + 'interdataglobal_points.dat' )

f = open( pth + 'interdataglobal_curves.dat', 'r' )
curves = []
nc = int(f.readline())
aabb = np.vstack((1.e6*np.ones(3), -1.e6*np.ones(3)))

for ic in range(nc):
    #endpoints = [ int(a) for a in f.readline().split() ]
    uvbox = np.zeros( (2,2,2) )
    for isurf in range(2):
        umin, umax, vmin, vmax = [ float(a) for a in f.readline().split() ]
        uvbox[0,0,isurf] = umin
        uvbox[1,0,isurf] = umax
        uvbox[0,1,isurf] = vmin
        uvbox[1,1,isurf] = vmax
    
    ns = int(f.readline())
    endpoints = np.zeros((ns,2))
    for ip in range(ns):
        endpoints[ip] = [ int(a)-1 for a in f.readline().split() ]
    
    n = int(f.readline())
    uv = np.zeros( (n,2,2) )
    xyz = np.zeros( (n,3) )
    for ip in range(n):
        u1, v1, u2, v2, x, y, z =  [ float(a) for a in f.readline().split() ]
        uv[ip,0,0] = u1
        uv[ip,0,1] = v1
        uv[ip,1,0] = u2
        uv[ip,1,1] = v2
        xyz[ip,0] = x
        xyz[ip,1] = y
        xyz[ip,2] = z
    for idim in range(3):
        aabb[0,idim] = min(aabb[0,idim], np.amin(xyz[:,idim]))
        aabb[1,idim] = max(aabb[1,idim], np.amax(xyz[:,idim]))
    c = intersection_curve( endpoints, uvbox, uv, xyz )
    curves.append( c )

f.close()
##########################################################

scn_diameter = np.sqrt( np.sum( np.power(aabb[1] - aabb[0], 2) ) )
print("scene diameter =",scn_diameter)

cls = myc.random_pastels(nsurf)
clc = [1.0,0.0,0.0]
thc = 2.e-3 * scn_diameter#2.5e-3
rad = 1.5*thc
ns = 200


## set up a blank scene
mylb.delCube()

## add surfaces
for isurf in range(nsurf):
    c = myl.readCoeffs2(pth + 'surfroot' + str(isurf+1) + '_x.cheb')
    mylb.addPatch(c, ns, cls[isurf])


## add intersection curves
for c in curves:
    mylb.addPolyline(
        c.xyz,
        clr=clc,
        thickness=thc,
        layer=0 )

## add intersection points
mylb.addVertices(points, clr=clc, radius=rad, layer=0)

## set camera
for ob in bpy.context.scene.objects:
    if ob.type == 'MESH':
        ob.select = True
    else:
        ob.select = False
        
bpy.ops.view3d.camera_to_view_selected()
#bpy.data.screens['Default'].spaces.lock_camera = True

scene = bpy.context.scene
# Environment lighting
scene.world.light_settings.use_environment_light = True
scene.world.light_settings.environment_energy = 0.25
scene.world.light_settings.environment_color = 'PLAIN'
scene.world.light_settings.samples = 10
scene.world.horizon_color = (1.0,1.0,1.0)

bpy.data.lamps['Lamp'].energy = 1.2
bpy.data.lamps['Lamp'].shadow_ray_samples = 8
bpy.data.lamps['Lamp'].shadow_soft_size = 2.0

# Render properties
scene.render.resolution_x = 1024
scene.render.resolution_y = 768
scene.render.resolution_percentage = 100
scene.render.use_border = False
