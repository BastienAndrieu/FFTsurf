from PIL import Image, ImageFont, ImageDraw, ImageChops, ImageEnhance
import os
import sys

args = sys.argv
if len(args) < 2:
    shading = "flat"#"shaded"#
else:
    shading = args[1]

nframes = 22

width = 1200
height = 675
left = 0
top = 76
box = (left, top, left+width, top+height)

###################################################################################
base = Image.open("base_"+shading+".png") # background (still) image
base = ImageEnhance.Color(base).enhance(1.2) # boost saturation
if shading == "shaded":
    base = ImageEnhance.Brightness(base).enhance(1.2) # boost brightness

for i in range(1,nframes+1):
    print "frame "+str(i)+"/"+str(nframes)
    frame = Image.open("im_%s.png" %i)
    image = ImageChops.multiply(frame,base)
    image = image.crop(box)
    """
    font = ImageFont.truetype("/usr/share/fonts/khmeros/KhmerOS.ttf", 36)
    draw = ImageDraw.Draw(image)
    j = max(0,i-2)
    draw.text((10,-10), "%s" %j, (0,0,0), font=font) 
    """
    image.save(shading+"_"+format(i,'02')+".png")

###################################################################################
fps = 6
cmd = "-ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:vbitrate=32000000:mbd=2:keyint=132:vqblur=1.0:cmp=2:subcmp=2:dia=2:mv0:last_pred=3 -mf type=png:fps="+str(fps)+" -vf scale="+str(width)+":"+str(height)+" -nosound -o"


os.system("ls "+shading+"_* > liste")

os.system("mencoder " + cmd + " /dev/null mf://@liste")
os.system("mencoder " + cmd + " " + shading + ".mp4 mf://@liste")
os.system("rm -f divx2pass.log liste")

