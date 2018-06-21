import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class intersection_curve:
    def __init__( self, endpoints, uvbox, uv, xyz ):
        self.endpoints = endpoints
        self.uvbox = uvbox
        self.uv = uv
        self.xyz = xyz


##########################################################


nv = 0

f = open( 'interdat_curves.dat', 'r' )
curves = []
nc = int(f.readline())

for ic in range(nc):
    endpoints = [ int(a) for a in f.readline().split() ]
    nv = max( nv, max(endpoints) )
    uvbox = np.zeros( (2,2,2) )
    for isurf in range(2):
        umin, umax, vmin, vmax = [ float(a) for a in f.readline().split() ]
        uvbox[0,0,isurf] = umin
        uvbox[1,0,isurf] = umax
        uvbox[0,1,isurf] = vmin
        uvbox[1,1,isurf] = vmax
    
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
    c = intersection_curve( endpoints, uvbox, uv, xyz )
    curves.append( c )

f.close()


uv = np.zeros( (nv,2,2) )
for ic in range(nc):
    for i, ip in enumerate(curves[ic].endpoints):
        j = i * (len(curves[ic].uv) - 1)
        uv[ip-1,:,:] = curves[ic].uv[j,:,:]


##########################################################

fig, axes = plt.subplots(nrows=1, ncols=2)
lw_tree = 0.5
lw_curve = 1.0
cl_tree = '0.7'#'k'
cl_curve = 'r'
cl_point = 'r'
cl_box = 'y'

for isurf in range(2):
    tree = np.loadtxt( 'tree_' + str(isurf+1) + '.dat' )
    if len(tree.shape) < 2:
        rect = patches.Rectangle( ( tree[0], tree[2] ), # lower-left corner
                                  tree[1] - tree[0],    # u-extent
                                  tree[3] - tree[2],    # v-extent
                                  fill = False,
                                  ec = cl_tree,
                                  lw = lw_tree )
        axes[isurf].add_patch( rect )
    else:
        for region in tree:
            rect = patches.Rectangle( ( region[0], region[2] ), # lower-left corner
                                      region[1] - region[0],    # u-extent
                                      region[3] - region[2],    # v-extent
                                      fill = False,
                                      ec = cl_tree,
                                      lw = lw_tree )
            axes[isurf].add_patch( rect )
    for ic in range(nc):
        rect = patches.Rectangle( curves[ic].uvbox[0,:,isurf],                               # lower-left corner
                                  curves[ic].uvbox[1,0,isurf] - curves[ic].uvbox[0,0,isurf], # u-extent
                                  curves[ic].uvbox[1,1,isurf] - curves[ic].uvbox[0,1,isurf], # v-extent                    
                                  fc = cl_box,
                                  ec ='none',
                                  alpha = 0.2 )
        axes[isurf].add_patch( rect )
        axes[isurf].plot( curves[ic].uv[:,isurf,0], curves[ic].uv[:,isurf,1], 
                          color = cl_curve, 
                          linewidth = lw_curve )
    
    for iv in range(nv):
        axes[isurf].plot( uv[iv,isurf,0], uv[iv,isurf,1], 
                          '.', 
                          color = cl_point )
        """
        axes[isurf].text( uv[iv,isurf,0], uv[iv,isurf,1], str(iv+1), 
                          color = cl_point, 
                          horizontalalignment = 'center', 
                          verticalalignment = 'bottom' )
        """
    axes[isurf].set_aspect('equal', adjustable='box')
    axes[isurf].set_xlim( -1, 1 )
    axes[isurf].set_ylim( -1, 1 )
    axes[isurf].set_title( 'Surface #' + str(isurf+1),
                           size = 12, weight = 'bold' )
    axes[isurf].set_xticks( np.linspace(-1.0,1.0,5) )
    axes[isurf].set_yticks( np.linspace(-1.0,1.0,5) )
    
plt.show()
