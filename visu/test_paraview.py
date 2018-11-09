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
renderView.CameraPosition = [-0.3440828157422071, -0.2491378250921422, 1.6159702216458034]
renderView.CameraViewUp = [-0.23629658042452942, 0.8937815256214137, 0.3812066506995025]
renderView.CameraFocalPoint = [0.6491938550658729, 0.8015521676750476, -0.23179246211295618]

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

SaveScreenshot('/d/bandrieu/Bureau/test4.png', magnification=2, CompressionLevel=9, view=renderView)
SaveScreenshot('/d/bandrieu/Bureau/test5.jpg', magnification=2, Quality=100, view=renderView)
SaveScreenshot('/d/bandrieu/Bureau/test6.jpg', magnification=2, Quality=0, view=renderView)

