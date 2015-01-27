# in-situ script to create images

try: paraview.simple
except: from paraview.simple import *

from paraview import numpy_support
from paraview import coprocessing
from paraview import simple
simple._DisableFirstRenderCameraReset()

import numpy
import numpy.ma

import datetime
current_date = datetime.datetime.now()

import mpi4py
import mpi4py.MPI
import mpi4py.MPI as MPI

#
# -------------------------------------------------------------------------
# main configuration options
# -------------------------------------------------------------------------

# choose grid type:
# X_Y: Cartesian domain 
# X_Y_Z: Spherical domain, on globe
# LON_LAT: Spherical domain, Mercator projection
# primal: cell-centered variables
# dual: variables on vertices
# 1LAYER: One z-layer is visualized at a time
# NLAYER: image is a 3D view of all layers together.

#gridname = 'X_Y_NLAYER-primal'
#gridname = 'X_Y_NLAYER-dual'
#gridname = 'X_Y_Z_1LAYER-primal'
#gridname = 'X_Y_Z_1LAYER-dual'
#gridname = 'X_Y_Z_NLAYER-primal'
#gridname = 'X_Y_Z_NLAYER-dual'
gridname = 'LON_LAT_1LAYER-primal'
#gridname = 'LON_LAT_1LAYER-dual'
#gridname = 'LON_LAT_NLAYER-primal'
#gridname = 'LON_LAT_NLAYER-dual'

# choose image sets to create.
# format: (<variable>, <filename pattern>, [layer])
# the 'layer' argument will be used only when using a 1LAYER grid.
# in 'filename pattern', %t will be replaced by the timestep, and %l will be
# replaced by the layer number.
variable_action = [
                   ('temperature', 'temperature_k%l_t%t.png', 0),
                   ('temperature', 'temperature_k%l_t%t.png', 25),
                   ('salinity', 'salinity_k%l_t%t.png', 6),
                   ('salinity', 'salinity_k%l_t%t.png', 12),
                  ]

# show axes or not
display_axes = False

# the zoom level of the image.
# greater than 1.0 is zoom-in and less than 1.0 is zoom-out
image_zoom = 2.3

# how often to do in-situ
frequency = 1

# annotation flags
annotationShowHeading = True
annotationShowRightSide = True
annotationShowLeftSide = True
annotationShowSimulationTime = True
annotationShowDate = True

# strings for annotations
annotationHeading = 'Title goes here'
annotationRightSide = 'Right side annotation'
annotationLeftSide = 'Left side annotation'
annotationColor = [1.0, 1.0, 1.0]

# the format of how the date annotation is shown
# for information about how to customize it, see the python datetime module
dateFormat = '%m/%d/%Y'

# how the simulation time annotation is shown
timeFormat = 'Timestep: %0.0f'

# image size
imageSize = [1024, 512]

# save data to a vtk file
save_scripts = ['mpas_vtkwriter']
save_data = False

# -------------------------------------------------------------------------

# Basics
RenderView1 = None
ScalarBar1 = None
LookupTable1 = None
Calculator1 = None
Threshold1 = None

# Available filter representations
BasicRep = None
ThresholdRep = None

# Annotation representations
HeadingRep = None
RightSideRep = None
LeftSideRep = None
TimeRep = None
DateRep = None

# Live visualization inside ParaView
live_visu_active = False
pv_host = "localhost"
pv_port = 22222

use_layer = False
if '1LAYER' in gridname:
    use_layer = True

if save_data:
  save_modules = []
  for script in save_scripts:
    save_modules.append(__import__(script))

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      pv_output_ = coprocessor.CreateProducer( datadescription, gridname )
      CreateMPASViews(pv_output_, variable_action[0][0], datadescription)

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  freqs = {gridname: [frequency]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

coprocessor = CreateCoProcessor()

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
  "Callback to populate the request for current timestep"
  global coprocessor

  if datadescription.GetForceOutput() == True:
    for i in range(datadescription.GetNumberOfInputDescriptions()):
      datadescription.GetInputDescription(i).AllFieldsOn()
      datadescription.GetInputDescription(i).GenerateMeshOn()
    return

  coprocessor.LoadRequestedData(datadescription)

  if save_data:
    for m in save_modules:
      m.RequestDataDescription(datadescription)

    # manually activate all required grids
    datadescription.GetInputDescriptionByName(gridname).AllFieldsOn()
    datadescription.GetInputDescriptionByName(gridname).GenerateMeshOn()
    for m in save_modules:
      for name, dataset in m.datasets.items():
        grid = dataset['grid']
        datadescription.GetInputDescriptionByName(grid).AllFieldsOn()
        datadescription.GetInputDescriptionByName(grid).GenerateMeshOn()

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
  "Callback to do co-processing for current timestep"
  global coprocessor

  timestep = datadescription.GetTimeStep()

  # Update the coprocessor by providing it the newly generated simulation data
  # If the pipeline hasn't been setup yet, this will setup the pipeline.
  coprocessor.UpdateProducers(datadescription)

  # Write output data, if appropriate.
  coprocessor.WriteData(datadescription);

  # Write images and charts
  WriteMPASImages(datadescription, timestep)

  # Live Visualization, if enabled.
  coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)

  if save_data:
    for m in save_modules:
      m.DoCoProcessing(datadescription)

# ----------------------- Write images ----------------------------------------

def WriteMPASImages(datadescription, timestep):
  global LookupTable1, ScalarBar1, RenderView1
  global BasicRep
  global HeadingRep, RightSideRep, LeftSideRep, TimeRep, DateRep
  global Calculator1, Threshold1, ThresholdRep

  grid = datadescription.GetInputDescriptionByName(gridname).GetGrid()

  for view in servermanager.GetRenderViews():
    if (timestep % view.cpFrequency == 0 or 
        datadescription.GetForceOutput() == True):

      for varTuple in variable_action:

        varName = varTuple[0]
        fileName = varTuple[1]
        if use_layer:
          layer = varTuple[2]

        HeadingRep.Visibility = 0
        RightSideRep.Visibility = 0
        LeftSideRep.Visibility = 0
        TimeRep.Visibility = 0
        DateRep.Visibility = 0

        # Output file name
        fileName = fileName.replace('%t', str(timestep))
        if use_layer:
          fileName = fileName.replace('%l', str(layer))

        ScalarBar1.Title = varName

        # Scalar range for colormap
        # need to mask out land values. first convert to a numpy array, mask
        # the values, and find the min and max of the remaining values
        data = numpy_support.vtk_to_numpy(grid.GetCellData().GetArray(varName))
        data = numpy.ma.masked_equal(data, -1e34)
        mins = data.min(axis=0)
        maxes = data.max(axis=0)
        if use_layer:
          minval = mins[layer]
          maxval = maxes[layer]
          LookupTable1.VectorMode = 'Component'
          LookupTable1.VectorComponent = layer
          ScalarBar1.ComponentTitle = str(layer)
        else:
          minval = mins
          maxval = maxes

        # find the global min and max by doing an mpi allreduce
        comm = mpi4py.MPI.COMM_WORLD

        sendbuf = numpy.array(minval, 'f')
        rcvbuf = numpy.array(0.0, 'f')
        comm.Allreduce([sendbuf, MPI.FLOAT], [rcvbuf, MPI.FLOAT], op=MPI.MIN)
        minval = rcvbuf

        sendbuf = numpy.array(maxval, 'f')
        rcvbuf = numpy.array(0.0, 'f')
        comm.Allreduce([sendbuf, MPI.FLOAT], [rcvbuf, MPI.FLOAT], op=MPI.MAX)
        maxval = rcvbuf

        # set the scalar bar range to the global min and max
        LookupTable1.RGBPoints=[minval, 0.23, 0.299, 0.754, 
                                (minval+maxval)/2.0, 0.865, 0.865, 0.865, 
                                maxval, 0.706, 0.016, 0.15]

        # if using layers, update the threshold to remove the land cells
        if use_layer:
          # set color bar title here
          array_name = '%s_%i' % (varName, layer)
          Calculator1.Function = array_name
          Calculator1.ResultArrayName = array_name

          BasicRep.Visibility = 0
          ThresholdRep.Visibility = 1
          ThresholdRep.ColorAttributeType = 'CELL_DATA'
          ThresholdRep.SelectionCellFieldDataArrayName = array_name
          ThresholdRep.ColorArrayName = ('CELL_DATA', array_name)
        else:
          BasicRep.ColorArrayName = varName

        if annotationShowHeading:
          HeadingRep.Visibility = 1
        if annotationShowRightSide:
          RightSideRep.Visibility = 1
        if annotationShowLeftSide:
          LeftSideRep.Visibility = 1
        if annotationShowSimulationTime:
          TimeRep.Visibility = 1
        if annotationShowDate:
          DateRep.Visibility = 1

        # camera bounds
        view.SMProxy.ResetCamera()
        view.SMProxy.GetActiveCamera().Zoom(image_zoom)

        view.ViewTime = datadescription.GetTime()
        Render(view)
        WriteImage(fileName, view, Magnification=view.cpMagnification)

# ----------------------- View definition -----------------------

def CreateMPASViews(producer, varName, datadescription):
  global coprocessor
  global LookupTable1, ScalarBar1, RenderView1
  global BasicRep
  global XYChartView1, PlotOverLine1, PlotOverLineRep
  global HeadingRep, RightSideRep, LeftSideRep, TimeRep, DateRep
  global Calculator1, Threshold1, ThresholdRep

# ----------------------- Filter definitions -----------------------

  # create a version of the data with point data
  #CellDatatoPointData1 = CellDatatoPointData(guiName="CellDatatoPointData1",
  #                                           PieceInvariant=1, PassCellData=0)

# ----------------------- Create render view -----------------------

  RenderView1 = coprocessor.CreateView( CreateRenderView, " ", 1, 0, 1,
                                        imageSize[0], imageSize[1] )

  RenderView1.LightSpecularColor = [1.0, 1.0, 1.0]
  RenderView1.UseOutlineForLODRendering = 0
  RenderView1.KeyLightAzimuth = 10.0
  RenderView1.UseTexturedBackground = 0
  RenderView1.UseLight = 1
  RenderView1.CameraPosition = [2202.3693, -1405.9497, 1382320.5176]
  RenderView1.FillLightKFRatio = 3.0
  RenderView1.Background2 = [0.0, 0.0, 0.165]
  RenderView1.FillLightAzimuth = -10.0
  RenderView1.LODResolution = 0.2733
  RenderView1.BackgroundTexture = []
  RenderView1.InteractionMode = '2D'
  RenderView1.StencilCapable = 1
  RenderView1.LightIntensity = 1.0
  RenderView1.CameraFocalPoint = [2202.3693, -1405.9497, 0.0]
  RenderView1.ImageReductionFactor = 2
  RenderView1.CameraViewAngle = 30.0
  RenderView1.CameraParallelScale = 32655.4490
  RenderView1.EyeAngle = 2.0
  RenderView1.HeadLightKHRatio = 3.0
  RenderView1.StereoRender = 0
  RenderView1.KeyLightIntensity = 0.75
  RenderView1.BackLightAzimuth = 110.0
  RenderView1.OrientationAxesInteractivity = 0
  RenderView1.UseInteractiveRenderingForSceenshots = 0
  RenderView1.UseOffscreenRendering = 0
  RenderView1.Background = [0.3199, 0.3400, 0.4299]
  RenderView1.UseOffscreenRenderingForScreenshots = 0
  RenderView1.NonInteractiveRenderDelay = 2.0
  RenderView1.CenterOfRotation = [160000.0, 0.0, 0.0]
  RenderView1.CameraParallelProjection = 0
  RenderView1.CompressorConfig = 'vtkSquirtCompressor 0 3'
  RenderView1.HeadLightWarmth = 0.5
  RenderView1.MaximumNumberOfPeels = 4
  RenderView1.LightDiffuseColor = [1.0, 1.0, 1.0]
  RenderView1.StereoType = 'Red-Blue'
  RenderView1.DepthPeeling = 1
  RenderView1.BackLightKBRatio = 3.5
  RenderView1.StereoCapableWindow = 1
  RenderView1.CameraViewUp = [0.0, 1.0, 0.0]
  RenderView1.LightType = 'HeadLight'
  RenderView1.LightAmbientColor = [1.0, 1.0, 1.0]
  RenderView1.RemoteRenderThreshold = 0.5
  RenderView1.CacheKey = 0.0
  RenderView1.UseCache = 0
  RenderView1.KeyLightElevation = 50.0
  RenderView1.CenterAxesVisibility = 0
  RenderView1.MaintainLuminance = 0
  RenderView1.StillRenderImageReductionFactor = 1
  RenderView1.BackLightWarmth = 0.5
  RenderView1.FillLightElevation = -75.0
  RenderView1.MultiSamples = 0
  RenderView1.FillLightWarmth = 0.4
  RenderView1.AlphaBitPlanes = 1
  RenderView1.LightSwitch = 0
  RenderView1.OrientationAxesVisibility = 1
  RenderView1.CameraClippingRange = [1368497.3124, 1403055.3253]
  RenderView1.BackLightElevation = 0.0
  RenderView1.ViewTime = 0.0
  RenderView1.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
  RenderView1.LODThreshold = 5.1
  RenderView1.CollectGeometryThreshold = 100.0
  RenderView1.UseGradientBackground = 0
  RenderView1.KeyLightWarmth = 0.6
  RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]

  LookupTable1 = GetLookupTableForArray( varName, 1, RGBPoints=[0.0, 0.23,
      0.299, 0.754, 0.5, 0.865, 0.865, 0.865, 1.0, 0.706, 0.016, 0.15],
      VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0],
      ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

  ScalarBar1 = CreateScalarBar( Title=varName, TextPosition=1, 
      Position=[0.1, 0.25], Position2=[0.13, 0.5], TitleOpacity=1.0, 
      TitleShadow=0, AutomaticLabelFormat=1, TitleFontSize=10,
      TitleColor=[1.0, 1.0, 1.0], AspectRatio=20.0, NumberOfLabels=5,
      ComponentTitle='', Resizable=0, TitleFontFamily='Arial', Visibility=1,
      LabelFontSize=8, LabelFontFamily='Arial', TitleItalic=0, Selectable=0,
      LabelItalic=0, Enabled=0, LabelColor=[1.0, 1.0, 1.0], LabelBold=0,
      UseNonCompositedRenderer=1, LabelOpacity=1.0, TitleBold=0,
      LabelFormat='%-#6.3g', Orientation='Vertical', LabelShadow=0,
      LookupTable=LookupTable1, Repositionable=0 )

  if use_layer:
      LookupTable1.VectorMode = 'Component'
      LookupTable1.VectorComponent = 0
      ScalarBar1.ComponentTitle = '0'

  RenderView1.Representations.append(ScalarBar1)

# ----------------------- Producer render -----------------------

  SetActiveSource(producer)
  SetActiveView(RenderView1)

  BasicRep = Show()
  BasicRep.ColorArrayName = varName 
  BasicRep.LookupTable = LookupTable1
  if display_axes:
    BasicRep.CubeAxesVisibility = 1
  else:
    BasicRep.CubeAxesVisibility = 0

  if use_layer:
    # we need to find land cells at a given layer, and extract it out
    # first we need to use the calculator to extract the correct layer
    # values, then we threshold cells which have land value -1e34

    SetActiveSource(producer)

    array_name = '%s_%i' % (variable_action[0][0], variable_action[0][2])
    Calculator1 = Calculator()
    Calculator1.Function = array_name
    Calculator1.ResultArrayName = 'layer'
    Calculator1.AttributeMode = 'Cell Data'

    SetActiveSource(Calculator1)
    Threshold1 = Threshold()
    Threshold1.ThresholdRange = [-1.0e+30, float('inf')]

    SetActiveSource(Threshold1)
    SetActiveView(RenderView1)

    ThresholdRep = Show()
    ThresholdRep.ColorAttributeType = 'CELL_DATA'
    ThresholdRep.SelectionCellFieldDataArrayName = array_name
    ThresholdRep.ColorArrayName = ('CELL_DATA', array_name)
    ThresholdRep.LookupTable = LookupTable1
    if display_axes:
      ThresholdRep.CubeAxesVisibility = 1
    else:
      ThresholdRep.CubeAxesVisibility = 0
    BasicRep.Visibility = 0

# ----------------------- Annotations -----------------------

  SetActiveView(RenderView1)
  HeadingText = Text(Text=annotationHeading)
  HeadingRep = Show()
  HeadingRep.WindowLocation = 'UpperCenter'
  HeadingRep.Color = annotationColor
  HeadingRep.FontSize = 18
  HeadingRep.TextScaleMode = 'Viewport'

  RightSideText = Text(Text=annotationRightSide)
  RightSideRep = Show()
  RightSideRep.WindowLocation = 'UpperRightCorner'
  RightSideRep.Color = annotationColor
  RightSideRep.FontSize = 10
  RightSideRep.Orientation = 90
  RightSideRep.TextScaleMode = 'Viewport'

  # get the version of mpas
  #grid = datadescription.GetInputDescriptionByName(gridname).GetGrid()
  #versionArray = grid.GetFieldData().GetArray("version")
  #modver = versionArray.GetTuple1(0)
  #modext = versionArray.GetTuple1(1)
  #versionstring = "MPAS %d.%d" % (modver, modext)

  LeftSideText = Text(Text=annotationLeftSide)
  LeftSideRep = Show()
  LeftSideRep.WindowLocation = 'UpperLeftCorner'
  LeftSideRep.Color = annotationColor
  LeftSideRep.FontSize = 10
  LeftSideRep.Orientation = 90
  LeftSideRep.TextScaleMode = 'Viewport'

  AnnotateTimeFilter1 = AnnotateTimeFilter()
  TimeRep = Show()
  AnnotateTimeFilter1.Format = timeFormat
  TimeRep.WindowLocation = 'LowerLeftCorner'
  TimeRep.Color = annotationColor
  TimeRep.FontSize = 10
  TimeRep.TextScaleMode = 'Viewport'

  datestring = current_date.strftime(dateFormat)
  DateText = Text(Text=datestring)
  DateRep = Show()
  DateRep.WindowLocation = 'LowerRightCorner'
  DateRep.Color = annotationColor
  DateRep.FontSize = 10
  DateRep.TextScaleMode = 'Viewport'
