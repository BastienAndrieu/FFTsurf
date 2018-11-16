import os
import sys
from PIL import Image

# arguments:
# 0: nom_script
# 1: cas
# 2: fps
# 3: duree premiere image
# 4: duree derniere image

args = sys.argv
if len(args) < 2:
    args.append('compensateur2')
if len(args) < 3:
    args.append('20')

cas = args[1]
fps = int(args[2])

im0 = Image.open(cas + "/img/instant_000.jpeg")
width, height = im0.size

scl = 0.5
width *= scl
height *= scl

cmd = "-ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:vbitrate=32000000:mbd=2:keyint=132:vqblur=1.0:cmp=2:subcmp=2:dia=2:mv0:last_pred=3 -mf type=jpg:fps="+str(fps)+" -vf scale="+str(int(width))+":"+str(int(height))+" -nosound -o "

"""
os.system("ls " + cas + "/img/instant_* > liste")

flist = open("liste", "r")
files = []
for l in flist:
    files.append(l)
#print files
flist.close()
"""

files = os.popen("ls " + cas + "/img/*.jpeg").readlines()

if len(args) < 4:
    nrep0 = 0
else:
    nrep0 = int(round(float(args[3])*float(fps)))

if len(args) < 5:
    nrep1 = nrep0
else:
    nrep1 = int(round(float(args[4])*float(fps)))

f = open("liste", "w")
for i in range(nrep0):
    f.write(files[0])
for l in files:
    f.write(l)
for i in range(nrep1):
    f.write(files[-1])
f.close()

os.system("mencoder " + cmd + " /dev/null mf://@liste")
os.system("mencoder " + cmd + " " + cas +"/film.mp4 mf://@liste")
os.system("rm -f divx2pass.log liste")

