import bpy
import numpy
from numpy.polynomial.chebyshev import chebgrid2d
import sys
sys.path.append('/d/bandrieu/GitHub/Code/Python/')
import lib_blender_util as lbu
import lib_blender_edit as lbe
import lib_chebyshev as lcheb
from lib_linalg import matmul

####################################
def add_surface_with_obbs(itest, isurf, compute_obbs=False):
    # add surfaces    
    c = lcheb.read_polynomial2(prefC + 'C%d_test%2.2d.txt' % (isurf, itest))
    degr = c.shape[0:2]
    print('degree = (%d, %d)' % (degr[0]-1, degr[1]-1))
    xyz = chebgrid2d(u, u, c)
    v, f = lbu.tensor_product_mesh_vf(xyz[0], xyz[1], xyz[2])
    obj = lbu.pydata_to_mesh(
        v,
        f,
        name='surf%d_test%d' % (isurf, itest)
    )
    lbe.set_smooth(obj)
    #
    # add bezier control polygon
    print('computing transformation matrices...')
    L = lcheb.C2Bmatrix_new(degr[0])
    if degr[1] == degr[0]:
        R = L.T
    else:
        R = lcheb.C2Bmatrix_new(degr[1]).T
    print('done')
    b = numpy.zeros((3,degr[0],degr[1]))
    for k in range(3):
        b[k] = matmul(matmul(L, c[:,:,k]), R)
    #
    v, f = lbu.tensor_product_mesh_vf(b[0], b[1], b[2])
    bp = lbu.pydata_to_mesh(
        v,
        f,
        name='bezier_polygon'
    )
    bp.draw_type = 'WIRE'
    #
    # add obbs
    vol = [0,0]
    for ibox in range(2):
        data = numpy.loadtxt(prefB + 'test%2.2dsurf%dbox%d.dat' % (itest, isurf, ibox+1))
        obb = lbu.obb_to_mesh(
            center=data[0],
            ranges=data[1],
            axes=data[2:5].T
        )
        vol[ibox] = numpy.prod(data[1])
        print('box %s: volume = %s' % (bases[ibox], vol[ibox]))
        obb.name = 'OBB_' + bases[ibox]
        obb.draw_type = 'WIRE'
    print('%s/%s volume ratio = %s' % (bases[1], bases[0], vol[1]/vol[0]))
    return obj
####################################

    
lbu.clear_scene(
    meshes=True,
    lamps=False,
    cameras=False
)


prefC = '/d/bandrieu/GitHub/FFTsurf/test/coeffstest/'
prefB = '/d/bandrieu/GitHub/FFTsurf/test/compare_bern_cheb/'
bases = ['Chebyshev', 'Bernstein']

m = 100
u = numpy.linspace(-1,1,m)

arg = k = int(sys.argv[3])

k = int(arg)
isurf = 2 - k%2
itest = 1 + int((k - isurf)/2)
print('test#%d, surf#%d' % (itest, isurf))
add_surface_with_obbs(itest, isurf)


