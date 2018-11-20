import os
import sys
import numpy as np
from my_vtk import write_vtk_mesh#, write_vtk_displacement


# clear all visu data
#os.system('rm -rf rep*')


# arguments:
# 0: nom_script
# 1: cas
# 2: 1er instant
# 3: dernier instant
args = sys.argv
if len(args) < 2:
	args.append('compensateur2')
if len(args) < 3:
	args.append('0')
if len(args) < 4:
	args.append(args[2])
if len(args) < 5:
	args.append('0')

cas = args[1]
pth = '../cases/' + cas + '/output/'

rep = np.loadtxt(pth + 'reprises.dat', dtype=int)
rep = np.hstack([rep, 10000])


tri = np.loadtxt(pth + 'connect_01.dat', dtype=int)
ref = np.loadtxt(pth + 'faceref_01.dat', dtype=int)

irep = -1
# instant(s)
i1 = int(args[2])
i2 = int(args[3])

if not os.path.exists(cas):
        os.mkdir(cas)
if not os.path.exists(cas + '/vtk/'):
        os.mkdir(cas + '/vtk/')
if not os.path.exists(cas + '/img/'):  
        os.mkdir(cas + '/img/')
for instant in range(i1,i2+1):
    newrep = False
    while True:
        if instant >= rep[irep] and instant < rep[irep+1]:
            break
        irep += 1
        newrep = True
        srep = format(irep+1,'02')
        tri = np.loadtxt(pth + 'connect_' + srep + '.dat', dtype=int)
        ref = np.loadtxt(pth + 'faceref_' + srep + '.dat', dtype=int)
    sinst = format(instant, '03')
    print 'instant ' + sinst + ', reprise ' + str(irep)
    xyz = np.loadtxt(pth + 'pos_' + sinst + '.dat', dtype=float)
    """
    if newrep:
        os.system('mkdir rep_' + srep)
        xyz0 = xyz
        write_vtk_mesh('rep_' + srep + '/initmesh.vtk', xyz0, tri, ref)
    """
    #else:
    #write_vtk_displacement('rep_' + srep + '/disp_' + format(instant,'03'), xyz - xyz0)
    write_vtk_mesh(cas + '/vtk/instant_' + format(instant,'03') + '.vtk', xyz, tri, ref)
