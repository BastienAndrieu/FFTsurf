from paraview.simple import *
import numpy as np

renderView = GetRenderView()

initmesh = OpenDataFile('/d/bandrieu/GitHub/FFTsurf/debug/test.vtk')

initmeshDisplay = Show(initmesh, renderView)
initmeshDisplay.ColorArrayName = [None, '']
initmeshDisplay.DiffuseColor = [1,1,1]
initmeshDisplay.EdgeColor = [0,0,0]
initmeshDisplay.LineWidth = 0.1
initmeshDisplay.Specular = 0.0
initmeshDisplay.Representation = 'Surface With Edges' # show mesh edges

renderView.CameraViewAngle = 40.0
renderView.CameraFocalPoint = [0,0,0.5]

renderView.Background = [1,1,1] # white background
renderView.AnnotationColor = [0,0,0] #?
renderView.KeyLightWarmth  = 0.5 # def. 0.6
renderView.FillLightWarmth = 0.5 # def. 0.4
renderView.BackLightWarmth = 0.5 # def. 0.5
renderView.HeadLightWarmth = 0.5 # def. 0.5

renderView.KeyLightIntensity = 0.8 # def. 0.75
renderView.FillLightKFRatio = 2.0 # def. 3.0
renderView.BackLightKBRatio = 2.0 # def. 3.5
renderView.HeadLightKHRatio = 3.0 # def. 3.0

renderView.PPI = 96 # def 96
renderView.ViewSize = [800, 600]
renderView.OrientationAxesVisibility = 0 # hide xyz axes
#renderView.UseLight = 0 # no shading

#SaveScreenshot('/d/bandrieu/Bureau/test4.png', magnification=2, quality=100, view=renderView)

camera = GetActiveCamera()
camera.Dolly(3)
#camera.Pitch(-20)
first = 0
last = 90
camera.Azimuth(first)
for i in range(last - first + 1):
    camera.Azimuth(1)
    SaveScreenshot('/stck/bandrieu/Bureau/Python/paraview/frame_%03d'%(i)+'.png', magnification=1, quality=100, view=renderView)
