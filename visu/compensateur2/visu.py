from paraview.simple import *
import os
import math

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

renderView = GetRenderView()

liste = os.popen("ls /d/bandrieu/GitHub/FFTsurf/visu/compensateur2/vtk/*.vtk").readlines()
for i in range(len(liste)):
    liste[i] = liste[i][:-1]
ninstants = len(liste)
#instants = LegacyVTKReader(FileNames=liste)
instants = OpenDataFile(liste)

# get animation scene
anim = GetAnimationScene()

# update animation scene based on data timesteps
anim.UpdateAnimationUsingDataTimeSteps()


disp = Show(instants, renderView)
disp.ColorArrayName = [None, '']
disp.DiffuseColor = [1,1,1]
disp.EdgeColor = [0,0,0]
disp.LineWidth = 0.1
disp.Specular = 0.0
disp.Representation = 'Surface With Edges' # show mesh edges

renderView.CameraViewAngle = 25.0
renderView.CameraPosition = [-0.9982480908232798, -0.29763221088743713, 1.8403594976731341]
renderView.CameraViewUp = [-0.16299842403304848, 0.9405733152064597, 0.29791500882680017]
renderView.CameraFocalPoint = [3.5549028556856603, 1.8819705625484349, -2.5498851661619217]


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


if False:
    renderView.UseLight = 0 # no shading
    # create a new 'Mesh Quality'
    meshQuality1 = MeshQuality(Input=instants)
    meshQuality1.TriangleQualityMeasure = 'Shape'
    # show data in view
    meshQuality1Display = Show(meshQuality1, renderView)
    # get Range
    """
    temporalStatistics1 = TemporalStatistics(Input=meshQuality1)
    minq = min(temporalStatistics1.CellData.GetArray("Quality_minimum").GetRange())
    maxq = max(temporalStatistics1.CellData.GetArray("Quality_maximum").GetRange())
    minq = math.floor(minq)
    maxq = math.ceil(maxq)
    print minq, maxq
    """
    minq = 0.0
    maxq = 1.0
     # show mesh edges
    meshQuality1Display.Representation = 'Surface With Edges'
    # get color transfer function/color map for 'Quality'
    qualityLUT = GetColorTransferFunction('Quality')
    # get opacity transfer function/opacity map for 'Quality'
    qualityPWF = GetOpacityTransferFunction('Quality')
    # hide data in view
    Hide(instants, renderView)
    # hide color bar/color legend
    meshQuality1Display.SetScalarBarVisibility(renderView, True)
    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    qualityLUT.ApplyPreset('RdOr', True)
    # invert the transfer function
    #qualityLUT.InvertTransferFunction()
    qualityLUT.AutomaticRescaleRangeMode = 'Never'
    # Rescale transfer function
    qualityLUT.RescaleTransferFunction(minq, maxq)
    # Rescale transfer function
    qualityPWF.RescaleTransferFunction(minq, maxq)
    # get color legend/bar for qualityLUT in current view
    qualityLUTColorBar = GetScalarBar(qualityLUT, renderView)
    # Properties modified on qualityLUTColorBar
    qualityLUTColorBar.Orientation = 'Horizontal'
    qualityLUTColorBar.HorizontalTitle = 1
    qualityLUTColorBar.TitleShadow = 1
    qualityLUTColorBar.LabelShadow = 1
    qualityLUTColorBar.TitleBold = 1
    qualityLUTColorBar.TitleColor = [0,0,0]
    qualityLUTColorBar.LabelColor = [0,0,0]
    qualityLUTColorBar.Title = 'Quality'
    qualityLUTColorBar.ComponentTitle = ''
    qualityLUTColorBar.LabelFormat = '%-#6.2g'
    qualityLUTColorBar.RangeLabelFormat = '%-#6.2g'
    qualityLUTColorBar.DrawTickMarks = 0
    qualityLUTColorBar.DrawTickLabels = 0

    qualityLUTColorBar.TitleFontFamily = 'File'
    qualityLUTColorBar.TitleFontSize = 18
    qualityLUTColorBar.TitleFontFile = '/d/bandrieu/GitHub/FFTsurf/visu/fonts/texgyreheros-bold-webfont.ttf'

    # Properties modified on qualityLUTColorBar
    qualityLUTColorBar.LabelFontFamily = 'File'
    qualityLUTColorBar.LabelFontFile = '/d/bandrieu/GitHub/FFTsurf/visu/fonts/texgyreheros-regular-webfont.ttf'

# save animation
SaveAnimation('/d/bandrieu/GitHub/FFTsurf/visu/compensateur2/img/instant_.jpeg',
              renderView,
              ImageResolution=[1600, 1200],
              FrameWindow=[0, int(anim.EndTime)], 
              # JPEG options
              Quality=80,
              SuffixFormat='%03d')

