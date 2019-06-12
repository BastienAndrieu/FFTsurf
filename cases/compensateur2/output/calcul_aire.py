import sys
import numpy as np

float_fmt = "08.6e"

propergol = [4,5]

# arguments:
# 0: nom_script
# 2: 1er instant
# 3: dernier instant
# 4: nb symetries
args = sys.argv
if len(args) < 2:
	args.append('0')
if len(args) < 3:
	args.append('0')
if len(args) < 4:
	args.append('0')

tri = np.loadtxt('connect_01.dat', dtype=int)
ref = np.loadtxt('faceref_01.dat', dtype=int)

# instant(s)
i1 = int(args[1])
i2 = int(args[2])


rep = np.loadtxt('reprises.dat', dtype=int)
rep = np.hstack([rep, 10000])


f = open('../compensateur2.opt','r')
while True:
        line = f.readline()
        if line[0:8] == 'timestep':
                break
dt = float(f.readline())
f.close()

time = dt*np.arange(i1,i2+1,1)
aire = np.zeros(i2-i1+1)

irep = 0
for instant in range(i1,i2+1):
        while True:
                if instant >= rep[irep] and instant < rep[irep+1]:
                        break
                irep += 1
                tri = np.loadtxt('connect_'+format(irep+1,'02')+'.dat', dtype=int)
                ref = np.loadtxt('faceref_'+format(irep+1,'02')+'.dat', dtype=int)
        sinst = format(instant, '03')
        print 'instant ' + sinst
        xyz = np.loadtxt('pos_' + sinst + '.dat', dtype=float)
        for i,t in enumerate(tri):
                if ref[i] in propergol:
                        u = xyz[t[1]-1] - xyz[t[0]-1]
                        v = xyz[t[2]-1] - xyz[t[0]-1]
                        w = np.cross(u,v)
                        aire[instant-i1] += 0.5*np.sqrt(np.sum(np.power(w,2)))

data = np.transpose(np.vstack([time,4.0*aire]))
np.savetxt('aire.dat', data)
