import os
import sys
import numpy as np


#############################################################
def write_obj( xyz, faces, filename ):
    f = open( filename, "w" )
    f.write( filename + "\n\n\n" )
    
    # vertices
    for v in xyz:
        f.write( "v ")
        for x in v:
            f.write( str(x) + " " )
        f.write( "\n" )
    f.write( "\n" )

    # faces
    for fa in faces:
        f.write( "f ")
        for i in fa:
            f.write( str(i) + " " )
        f.write( "\n" )
    f.write( "\n" )

    f.close()
    return
#############################################################



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
if not os.path.exists(cas + '/obj/'):
        os.mkdir(cas + '/obj/')

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
    
write_obj(xyz, tri, cas + '/obj/instant_' + format(instant,'03') + '.obj')
