ROOT = '/d/bandrieu/'#'/home/bastien/'

import bpy

import numpy
from numpy import pi, sin, cos, arcsin
import sys
sys.path.append(ROOT + 'GitHub/Code/Python/')
import lib_blender_util as lbu
import lib_cadcheb as lcad


def add_primitive(xyz, name='primitive'):
    v, f = lbu.tensor_product_mesh_vf(xyz[0], xyz[1], xyz[2])
    obj = lbu.pydata_to_mesh(
        v,
        f,
        name=name
    )
    obj.show_wire = True
    obj.show_all_edges = True
    return obj


#############################################
n = 16
offset = 2.0
scale = 1e-2
tolchord = 1e-3
hmin = 8e-3
hmax = 1e-2
#############################################


m = n

## ********* Géométrie ********* ##
b0 = numpy.eye(3)
ex = b0[:,0]
ey = b0[:,1]
ez = b0[:,2]

angle = 2*pi/5

r1 = 10.0
r2 = 20.0
l1 = 100.0
l2 = 10.0

wp1 = 2.1*r1
wp2 = 2.1*r2

r3 = 50.0
r_CC1 = 10.0
r_CC2 = 5.0
r_f = 1.0

w = 1.5

alpha = numpy.arctan2(r2-r1,l1)
beta = numpy.radians(125.)
#############################################

BREP = lcad.BREP_t()


# 0 (fond avant, plan)
o = [-offset, 0.5*wp1, -0.5*wp1]
base = lcad.matrix_from_columns([-ey, ez, -ex])
xyz = lcad.rectangle_cgl(o, wp1, wp1, base, m, n)
BREP.patches.append(lcad.Patch_t(xyz))

# 1 (fond arriere, plan)
o = [l1+l2+offset, -0.5*wp2, -0.5*wp2]
base = lcad.matrix_from_columns([ey, ez, ex])
xyz = lcad.rectangle_cgl( o, wp2, wp2, base, m, n )
BREP.patches.append(lcad.Patch_t(xyz))

# **** CONE-CYL ****
nrepcc = 5#3#
anglecc = 2*pi/nrepcc

# 2 (cone)
c_o = [-offset*sin(alpha),0,0]
r1_o = r1 + offset*cos(alpha)
r2_o = r2 + offset*cos(alpha)
base = lcad.matrix_from_columns([ey, ez, ex])
xyz = lcad.cone_cgl(
    c_o,
    r1_o,
    r2_o,
    l1,
    base,
    0.5*(pi-anglecc),
    0.5*(pi+anglecc),
    m,
    n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[3,5]
    )
)

# 3 (raccord tore avant)
base = lcad.matrix_from_columns([ez, ey, -ex])
xyz = lcad.torus_cgl(
    [0,0,0],
    r1,
    offset,
    base,
    -0.5*anglecc,
    0.5*anglecc,
    alpha,
    pi/2,
    m,
    n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[0,2]
    )
)

# 4 (cylindre)
base = lcad.matrix_from_columns([ey, ez, ex])
xyz = lcad.cylinder_cgl(
    [l1,0,0],
    r2+offset,
    l2,
    base,
    0.5*(pi-anglecc),
    0.5*(pi+anglecc),
    m,
    n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[5,6]
    )
)

# 5 (raccord tore milieu)
base = lcad.matrix_from_columns([ez, ey, -ex])
xyz = lcad.torus_cgl(
    [l1,0,0],
    r2,
    offset,
    base,
    -0.5*anglecc,
    0.5*anglecc,
    0,
    alpha,
    m,
    n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[2,4]
    )
)

# 6 (raccord tore arriere)
base = lcad.matrix_from_columns([ez, ey, -ex])
xyz = lcad.torus_cgl(
    [l1+l2,0,0],
    r2,
    offset,
    base,
    -0.5*anglecc,
    0.5*anglecc,
    -pi/2,
    0,
    m,
    n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[1,4]
    )
)

# Periodicite ConeCyl
first_percc = 2
last_percc = 6
n_percc = last_percc - first_percc + 1
for k in range(1,nrepcc):
    a = k*anglecc
    for i in range(n_percc):
        #j = first_percc + k*n_percc + i
        iPa = first_percc + i
        xyz = BREP.patches[iPa].xyz.copy()
        adj = BREP.patches[iPa].adj.copy()

        y = -sin(a)*xyz[2] + cos(a)*xyz[1]
        z = cos(a)*xyz[2] + sin(a)*xyz[1]
        xyz[1] = y
        xyz[2] = z

        #l = find(adj >= first_percc & adj <= last_percc)
        #adj(l) = (k-1)*n_percc + adj(l);
        for l in range(len(adj)):
            if adj[l] >= first_percc and adj[l] <= last_percc: adj[l] += (k-1)*n_percc
        BREP.patches.append(
            lcad.Patch_t(
                xyz,
                adj
            )
        )

for k in range(nrepcc):
    for i in range(n_percc):
        j = first_percc + k*n_percc + i
        jm = first_percc + ((k-1)%nrepcc)*n_percc + i
        jp = first_percc + ((k+1)%nrepcc)*n_percc + i
        BREP.patches[j].adj += [jm,jp]


## Fins
npa = len(BREP.patches)
nrepfin = 5#9
anglefin = 2*pi/nrepfin

pC = numpy.array([80.0, 40.0])
pO1 = pC - [0,r_CC1]
pB = pO1 - r_CC1*numpy.array([sin(beta), cos(beta)])
pD = [95.0, 40,0]
pE = [100.0, 35.0]
pO2 = [95.0, 35.0]

pF = [100.0, r2*cos(arcsin(w/r2))]

# (Rectangle) 'front plane'
c = [pB[0], w, pB[1]]
c_o = c - offset*numpy.array([sin(beta), 0, cos(beta)])

lp = 26.0

base = lcad.rotate_basis(b0, 2, pi)
base = lcad.rotate_basis(base, 1, beta)
base = lcad.rotate_basis(base, 0, 0)

xyz = lcad.rectangle_cgl(c_o, lp, 2*w, base, m, n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa+1,npa+5]
    )
)

# (Cylindre) 'front round corner'
c = [pO1[0], -w, pO1[1]]

base = lcad.rotate_basis(b0, 2, pi)
base = lcad.rotate_basis(base, 1, beta)

xyz = lcad.cylinder_cgl(c,r_CC1+offset,-2*w,base,0,beta-pi,m,n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa+0,npa+2,npa+6]
    )
)

# (Rectangle) 'top plane'
c = numpy.array([pC[0], -w, pC[1]])
c_o = c + [0,0,offset]

base = lcad.rotate_basis(b0, 2, 0)

xyz = lcad.rectangle_cgl(c_o, pD[0]-pC[0], 2*w, base, m, n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa+1,npa+4,npa+7]
    )
)

# (Cylindre) 'rear round corner'
base = lcad.rotate_basis(b0, 1, 0)
c = [pO2[0], -w, pO2[1]]

xyz = lcad.cylinder_cgl(c,r_CC2+offset,2*w,base,0,0.5*pi,m,n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa+2,npa+4,npa+8]
    )
)

# (Rectangle) 'rear'
base = lcad.rotate_basis(b0,1,-pi/2)
base = lcad.rotate_basis(base,1,0)

c = numpy.array([pE[0], -w, pE[1]])
c_o = c + [offset, 0, 0]

xyz = lcad.rectangle_cgl(c_o, pE[1]-pF[1], 2*w, base, m, n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa+3,npa+9]
    )
)

# (Cylindre) 'front fillet'
c = [pB[0], w, pB[1]]

base = lcad.rotate_basis(b0, 2, pi)
base = lcad.rotate_basis(base, 1, beta)
base = lcad.rotate_basis(base, 0, 0)
v = base[:,0].copy()
base[:,0] = -base[:,1]
base[:,1] = v.copy()
base = base[:,[2,0,1]]

c = c - r_f*base[:,0]

xyz = lcad.cylinder_cgl(c,r_f+offset,lp,base,0,0.5*pi,m,n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa,npa+6,npa+10]
    )
)

# (Tore) 'top-front fillet'
base = base[:,[0,2,1]]
base[:,1] = -base[:,1]

c = [pB[0], w, pB[1]]
c = c - r_CC1*base[:,0]

xyz = lcad.torus_cgl(
    c,
    r_CC1-r_f,
    r_f+offset,
    base,
    0,
    pi-beta,
    0,
    pi/2,
    m,
    n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa+i-1 for i in [2,6,8,11]]
    )
)

# (Cylindre) 'top fillet'
c = [pC[0], w, pC[1]-r_f]

base = b0[:,[1,2,0]]

xyz = lcad.cylinder_cgl(c,r_f+offset,pD[0]-pC[0],base,0,pi/2,m,n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa+i-1 for i in [3,7,9,11]]
    )
)

# (Tore) 'top-rear fillet'
c = [pD[0], w, pD[1]-r_CC2]

base = b0[:,[2,0,1]]

xyz = lcad.torus_cgl(
    c,
    r_CC2-r_f,
    r_f+offset,
    base,
    0,
    pi/2,
    0
    ,pi/2,
    m,
    n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa+i-1 for i in [4,8,10,11]]
    )
)

# (Cylindre) 'rear fillet'
base = b0[:,[1,0,2]]
base[:,2] = -base[:,2]

c = [pE[0]-r_f, w ,pE[1]]

xyz = lcad.cylinder_cgl(c,r_f+offset,pE[1]-pF[1],base,0,pi/2,m,n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa+i-1 for i in [5,9,11]]
    )
)

# (plan)
mrg = 0#1.33 * r_f
pG = BREP.patches[npa+5].xyz[:,-1,-1]
pH = BREP.patches[npa+5].xyz[:,-1,0]
pI = BREP.patches[npa+6].xyz[:,-1,-1]
pL = BREP.patches[npa+9].xyz[:,0,-1]
pP = numpy.array([pL[0], pL[1], pI[2]])
pO = pG + (pH - pG) * (pI[2] - pG[2])/(pH[2] - pG[2])

l1 = pL[0] - pG[0] + 2*mrg
l2 = pI[2] - pG[2] + 2*mrg
o = [pL[0]+mrg, pG[1], pG[2]-mrg]
base = lcad.matrix_from_columns([-ex, ez, ey])

xyz = lcad.bilinear_cgl(pL, pG, pO, pP, m, n)
BREP.patches.append(
    lcad.Patch_t(
        xyz,
        adj=[npa+i-1 for i in [6,7,8,9,10]]
    )
)
                 
isideplanes = [npa + 10]

# Symmetry FIN
first_sym = npa+5
last_sym = npa+10
n_sym = last_sym - first_sym + 1
for i in range(n_sym):
    j = last_sym + i + 1
    iPa = first_sym + i
    xyz = lcad.symmetry(BREP.patches[iPa].xyz, iaxe=1)
    adj = BREP.patches[iPa].adj.copy()
    #l = find(adj >= first_sym & adj <= last_sym);
    #adj(l) = adj(l) + n_sym;
    #BREP.Pa(j).adj = adj;
    BREP.patches.append(
        lcad.Patch_t(
            xyz,
            adj
        )
    )
    
    if iPa == isideplanes[0]: isideplanes.append(j)

#  Periodicite FIN
first_perfin = first_sym - 5
last_perfin = last_sym + n_sym
n_perfin = last_perfin - first_perfin + 1
for k in range(1,nrepfin):
    a = k*anglefin
    for i in range(n_perfin):
        #j = first_perfin - 1 + (k-1)*n_perfin + i;
        iPa = first_perfin + i
        xyz = BREP.patches[iPa].xyz.copy()
        adj = BREP.patches[iPa].adj.copy()

        y = -sin(a)*xyz[2] + cos(a)*xyz[1]
        z = cos(a)*xyz[2] + sin(a)*xyz[1]
        xyz[1] = y
        xyz[2] = z
        
        for l in range(len(adj)):
            if adj[l] >= first_perfin and adj[l] <= last_perfin: adj[l] += (k-1)*n_perfin
        
        BREP.patches.append(
            lcad.Patch_t(
                xyz,
                adj
            )
        )
        #if ~isempty(intersect(first_perfin - 1 + i, isideplanes))
        #    isideplanes = [isideplanes, j];



# Tags
for Pa in BREP.patches:
    Pa.tag = 1

"""
# ********* Topologie *********
BREP.G1_continuity = zeros(length(BREP.Pa));
for i = 1:length(BREP.Pa)
    for j = BREP.Pa(i).adj
        BREP = add_edge(BREP, i, j, 1 );
    end
end
"""
        
# Scaling
for Pa in BREP.patches:
    Pa.xyz = scale*Pa.xyz

"""
# tangential intersection curves (--> smooth edges)
nfaces = length(BREP.Pa);
tgc_faces = [];
for i = 1:nfaces-1
    for j = i+1:nfaces
        if BREP.G1_continuity(i,j)
            tgc_faces = [tgc_faces;[i,j]];
        end
    end
end
"""

iplanes = [1,2,isideplanes];

curves = trace_tangential_intersection_curves( ...
        BREP, tgc_faces, iplanes, hmin, hmax, tolchord );











    
lbu.clear_scene(True, True, True)


for iPa, Pa in enumerate(BREP.patches):
    add_primitive(Pa.xyz, 'patch_'+format(iPa, '03'))
