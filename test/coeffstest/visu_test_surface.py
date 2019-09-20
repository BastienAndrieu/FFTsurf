import bpy
import numpy
from numpy.polynomial.chebyshev import chebgrid2d
import sys
sys.path.append('/d/bandrieu/GitHub/Code/Python/')
import lib_blender_util as lbu
import lib_blender_edit as lbe
import lib_chebyshev as lcheb

####################################
def add_surface(itest, isurf):
    c = lcheb.read_polynomial2(pref + 'C%d_test%2.2d.txt' % (isurf, itest))
    #
    xyz = chebgrid2d(u, u, c)
    #
    v, f = lbu.tensor_product_mesh_vf(xyz[0], xyz[1], xyz[2])
    #
    obj = lbu.pydata_to_mesh(
        v,
        f,
        name='surf%d_test%d' % (isurf, itest)
    )
    lbe.set_smooth(obj)
    return obj

####################################

lbu.clear_scene(
    meshes=True,
    lamps=False,
    cameras=False
)


pref = '/d/bandrieu/GitHub/FFTsurf/test/coeffstest/'

m = 100
u = numpy.linspace(-1,1,m)

args = sys.argv[3:]
nargs = len(args)
if nargs == 1:
    k = int(args[0])
    isurf = 1 + k%2
    itest = 1 + int((k - isurf)/2)
    add_surface(itest, isurf)
elif nargs > 1:
    ntest = int(nargs/2)
    for jtest in range(ntest):
        itest = int(args[2*jtest])
        isurf = int(args[2*jtest+1])
        add_surface(itest, isurf)
