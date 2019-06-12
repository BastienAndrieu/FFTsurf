import sys
import numpy as np

float_fmt = "08.6e"

# arguments:
# 0: nom_script
# 2: 1er instant
# 3: dernier instant
args = sys.argv
if len(args) < 2:
	args.append('0')
if len(args) < 3:
	args.append('0')

tri = np.loadtxt('connect.dat', dtype=int)

# instant(s)
i1 = int(args[1])
i2 = int(args[2])
for instant in range(i1,i2+1):
    sinst = format(instant, '03')
    print 'instant ' + sinst
    xyz = np.loadtxt('pos_' + sinst + '.dat', dtype=float)
    f = open('tecplot/instant_' + sinst + '.dat', 'w')
    f.write('VARIABLES = "X" "Y" "Z"\n')
    f.write('ZONE T="instant ' + sinst + '"\n')
    f.write('STRANDID=0, SOLUTIONTIME=0\n')
    f.write('Nodes=' + str(len(xyz)) + ', Elements=' + str(len(tri)) + '\n')
    f.write('ZONETYPE=FETriangle\n')
    f.write('DATAPACKING=POINT\n')
    for v in xyz:
        for x in v:
            f.write(format(x, float_fmt) + ' ')
        f.write('\n')
    for t in tri:
        for i in t:
            f.write(str(i) + ' ')
        f.write('\n')
    f.close()
