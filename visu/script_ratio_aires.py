from paraview.simple import *





renderView = GetRenderView()


# find source
instant_0 = FindSource('instant_0*')
# hide data in view
Hide(instant_0, renderView)

# create a new 'Point Data to Cell Data'
# get hTarget at cells
pointDatatoCellData1 = PointDatatoCellData(Input=instant_0)


# create a new 'Cell Size'
cellSize1 = CellSize(Input=pointDatatoCellData1)
# only compute cell areas
cellSize1.ComputeVertexCount = 0
cellSize1.ComputeLength = 0
cellSize1.ComputeVolume = 0

"""
# create a new 'Python Calculator'
# compute unscaled ideal cell area
pythonCalculator1 = PythonCalculator(Input=cellSize1)
pythonCalculator1.ArrayAssociation = 'Cell Data'
pythonCalculator1.Expression = 'hTarget**2'
pythonCalculator1.ArrayName = 'IdealAreaNN'

# create a new 'Python Calculator'
pythonCalculator2 = PythonCalculator(Input=pythonCalculator1)
pythonCalculator2.ArrayAssociation = 'Cell Data'
pythonCalculator2.Expression = 'IdealAreaNN*sum(Area)/sum(IdealAreaNN)'
pythonCalculator2.ArrayName = 'Actual/Ideal area ratio'
"""
pythonCalculator2 = PythonCalculator(Input=cellSize1)
pythonCalculator2.ArrayAssociation = 'Cell Data'
pythonCalculator2.Expression = 'Area*sum(hTarget**2)/(sum(Area)*hTarget**2)'
pythonCalculator2.ArrayName = 'Actual/Ideal area ratio'

# show data in view
pythonCalculator2Display = Show(pythonCalculator2, renderView)
# show color bar/color legend
pythonCalculator2Display.SetScalarBarVisibility(renderView, True)
pythonCalculator2Display.SetRepresentationType('Surface With Edges')
# update the view to ensure updated data information
renderView.Update()

