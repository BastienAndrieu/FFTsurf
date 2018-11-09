from paraview.simple import *

renderView = GetRenderView()

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
#renderView.ViewSize = [1600, 900]
renderView.OrientationAxesVisibility = 0 # hide xyz axes
#renderView.UseLight = 0 # no shading

Render(renderView)
