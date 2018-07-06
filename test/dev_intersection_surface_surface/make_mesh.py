import numpy as np

##################################################
class intersection_curve:
    def __init__( self, endpoints, uvbox, uv, xyz ):
        self.endpoints = endpoints
        self.uvbox = uvbox
        self.uv = uv
        self.xyz = xyz

##################################################

# Read intersection polylines
f = open( 'interdat_curves.dat', 'r' )
curves = []
nc = int(f.readline())
for ic in range(nc):
    endpoints = [ int(a) for a in f.readline().split() ]
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
##################################################

# Common 'info.dat' file
hmax = 0.0
hmin = 1.e9
for c in curves:
    n = len( c.xyz )
    for i in range(n-1):
        h = np.sum( np.power( c.xyz[i+1] - c.xyz[i], 2 ) )
        hmax = max( hmax, h )
        hmin = min( hmin, h )

mrg = 0.1
hmin = (1.0 - mrg) * np.sqrt(hmin)
hmax = (1.0 + mrg) * np.sqrt(hmax)
tol = 5.e-2

f = open( 'meshgen_info.dat', 'w' )
f.write( str(hmin) + '\n' + str(hmax) + '\n' + str(tol) + '\n' )
f.close()


nb = 10#50
t = np.linspace( -1.0, 1.0, nb )
EPSuv = 1.e-9
lim = 1.0 - EPSuv
# Generate mesh for each surface
for isurf in range(2):    
    # boundary points on intersection curves
    uvbc1 = np.empty( (0,2) )
    uvbc2 = np.empty( (0,2) )
    uvbc3 = np.empty( (0,2) )
    uvbc4 = np.empty( (0,2) )
    for c in curves:
        n = len( c.xyz )
        for ipt in range(n):
            if   c.uv[ipt,isurf,1] < -lim and abs(c.uv[ipt,isurf,0]) < lim:
                uvbc1 = np.vstack( ( uvbc1, c.uv[ipt,isurf] ) )
            elif c.uv[ipt,isurf,0] >  lim and abs(c.uv[ipt,isurf,1]) < lim:
                uvbc2 = np.vstack( ( uvbc2, c.uv[ipt,isurf] ) )
            elif c.uv[ipt,isurf,1] >  lim and abs(c.uv[ipt,isurf,0]) < lim:
                uvbc3 = np.vstack( ( uvbc3, c.uv[ipt,isurf] ) )
            elif c.uv[ipt,isurf,0] < -lim and abs(c.uv[ipt,isurf,1]) < lim:
                uvbc4 = np.vstack( ( uvbc4, c.uv[ipt,isurf] ) )
    
    # 'natural' boundary points
    uvb1 = np.transpose( np.vstack( ( t[0:-1:1]     , -np.ones(nb-1) ) ) )
    uvb2 = np.transpose( np.vstack( ( np.ones(nb-1) , t[0:-1:1]      ) ) )
    uvb3 = np.transpose( np.vstack( ( t[-1:0:-1]    , np.ones(nb-1)  ) ) )
    uvb4 = np.transpose( np.vstack( ( -np.ones(nb-1), t[-1:0:-1]     ) ) )
    
    # merge
    if len(uvbc1) > 0:
        uvb1 = np.vstack( ( uvb1, uvbc1 ) )
    if len(uvbc2) > 0:
        uvb2 = np.vstack( ( uvb2, uvbc2 ) )
    if len(uvbc3) > 0:
        uvb3 = np.vstack( ( uvb3, uvbc3 ) )
    if len(uvbc4) > 0:
        uvb4 = np.vstack( ( uvb4, uvbc4 ) )
    
    # sort
    uvb1.sort( axis=0 )
    uvb2.sort( axis=0 )
    uvb3[::-1].sort( axis=0 )
    uvb4[::-1].sort( axis=0 )
    
    uvb = np.vstack( ( uvb1, uvb2, uvb3, uvb4 ) )

##################################################
