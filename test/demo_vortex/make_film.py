import os
import sys
from PIL import Image

# arguments:
# 0: nom_script
# 1: fps
# 2: nom fichier film

args = sys.argv
if len(args) < 2:
    args.append('25')

fps = int(args[1])

if len(args) < 3:
    filename = "film"
else:
    filename = args[2]

im0 = Image.open("imseq/im_001.jpg")
width, height = im0.size

scl = 1
width *= scl
height *= scl

cmd = "-ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:vbitrate=32000000:mbd=2:keyint=132:vqblur=1.0:cmp=2:subcmp=2:dia=2:mv0:last_pred=3 -mf type=jpg:fps="+str(fps)+" -vf scale="+str(int(width))+":"+str(int(height))+" -nosound -o "

files = os.popen("ls imseq/*.jpg").readlines()

f = open("liste", "w")
for l in files:
    f.write(l)
f.close()

os.system("mencoder " + cmd + " /dev/null mf://@liste")
os.system("mencoder " + cmd + filename + ".mp4 mf://@liste")
os.system("rm -f divx2pass.log liste")

