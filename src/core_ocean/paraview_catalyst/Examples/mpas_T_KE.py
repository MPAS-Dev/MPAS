
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 4.2.0 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 4.2.0

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1360, 600]
      renderView1.InteractionMode = '2D'
      renderView1.OrientationAxesVisibility = 0
      renderView1.CenterOfRotation = [0.07033830543157371, 0.006836779980497865, 0.0]
      renderView1.CameraPosition = [0.01894335010112129, 0.06757627264375982, 17.916893783590513]
      renderView1.CameraFocalPoint = [0.01894335010112129, 0.06757627264375982, 0.0]
      renderView1.CameraParallelScale = 1.4016805999214295
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='mpas_T_%t.png', freq=1, fittoscreen=0, magnification=1, width=1360, height=600)

      # Create a new 'Render View'
      renderView2 = CreateView('RenderView')
      renderView2.ViewSize = [1360, 600]
      renderView2.OrientationAxesVisibility = 0
      renderView2.CenterOfRotation = [-0.009500407345893302, 0.10769542137877294, 0.0]
      renderView2.CameraPosition = [-0.009500407345893302, 0.10769542137877294, 5.256500455626923]
      renderView2.CameraFocalPoint = [-0.009500407345893302, 0.10769542137877294, 0.0]
      renderView2.CameraParallelScale = 4.637233340272187
      renderView2.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView2,
          filename='mpas_KE_%t.png', freq=1, fittoscreen=0, magnification=1, width=1360, height=600)

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Partitioned Unstructured Grid Reader'
      # create a producer from a simulation input
      mpas_data_1pvtu = coprocessor.CreateProducer(datadescription, 'LON_LAT_1LAYER-primal')

      # create a new 'XML Partitioned Unstructured Grid Reader'
      # create a producer from a simulation input
      mpas_data_1pvtu_1 = coprocessor.CreateProducer(datadescription, 'LON_LAT_1LAYER-primal')

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'kineticEnergyCell'
      kineticEnergyCellLUT = GetColorTransferFunction('kineticEnergyCell')
      kineticEnergyCellLUT.RGBPoints = [0.0001, 0.0, 0.0, 0.0, 0.0021478304741305307, 0.0, 0.0, 0.501961, 0.04613175745603795, 0.0, 0.501961, 1.0, 1.0, 1.0, 1.0, 1.0]
      kineticEnergyCellLUT.UseLogScale = 1
      kineticEnergyCellLUT.LockScalarRange = 1
      kineticEnergyCellLUT.ColorSpace = 'RGB'
      kineticEnergyCellLUT.NanColor = [1.0, 1.0, 0.0]
      kineticEnergyCellLUT.ScalarRangeInitialized = 1.0
      kineticEnergyCellLUT.VectorMode = 'Component'

      # get opacity transfer function/opacity map for 'kineticEnergyCell'
      kineticEnergyCellPWF = GetOpacityTransferFunction('kineticEnergyCell')
      kineticEnergyCellPWF.Points = [0.0001, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
      kineticEnergyCellPWF.ScalarRangeInitialized = 1

      # get color transfer function/color map for 'temperature'
      temperatureLUT = GetColorTransferFunction('temperature')
      temperatureLUT.RGBPoints = [-1.642665147781372, 0.278431, 0.278431, 0.858824, 2.695169249296188, 0.0, 0.0, 0.360784, 7.002669140100478, 0.0, 1.0, 1.0, 11.37083804345131, 0.0, 0.501961, 0.0, 15.678337934255598, 1.0, 1.0, 0.0, 20.01617233133316, 1.0, 0.380392, 0.0, 24.35400672841072, 0.419608, 0.0, 0.0, 28.69184112548828, 0.878431, 0.301961, 0.301961]
      temperatureLUT.ColorSpace = 'RGB'
      temperatureLUT.NanColor = [1.0, 1.0, 0.0]
      temperatureLUT.ScalarRangeInitialized = 1.0
      temperatureLUT.VectorMode = 'Component'

      # get opacity transfer function/opacity map for 'temperature'
      temperaturePWF = GetOpacityTransferFunction('temperature')
      temperaturePWF.Points = [-1.642665147781372, 0.0, 0.5, 0.0, 28.69184112548828, 1.0, 0.5, 0.0]
      temperaturePWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from mpas_data_1pvtu
      mpas_data_1pvtuDisplay = Show(mpas_data_1pvtu, renderView1)
      # trace defaults for the display properties.
      mpas_data_1pvtuDisplay.ColorArrayName = ['CELLS', 'temperature']
      mpas_data_1pvtuDisplay.LookupTable = temperatureLUT
      mpas_data_1pvtuDisplay.ScalarOpacityUnitDistance = 0.28142083362291853

      # show color legend
      mpas_data_1pvtuDisplay.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for temperatureLUT in view renderView1
      temperatureLUTColorBar = GetScalarBar(temperatureLUT, renderView1)
      temperatureLUTColorBar.Position = [0.7657688540646425, 0.548812351543943]
      temperatureLUTColorBar.Position2 = [0.11999999999999988, 0.4299999999999998]
      temperatureLUTColorBar.Title = 'temperature'
      temperatureLUTColorBar.ComponentTitle = '0'

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView2'
      # ----------------------------------------------------------------

      # show data from mpas_data_1pvtu
      mpas_data_1pvtuDisplay_1 = Show(mpas_data_1pvtu, renderView2)
      # trace defaults for the display properties.
      mpas_data_1pvtuDisplay_1.ColorArrayName = ['CELLS', 'kineticEnergyCell']
      mpas_data_1pvtuDisplay_1.LookupTable = kineticEnergyCellLUT
      mpas_data_1pvtuDisplay_1.ScalarOpacityUnitDistance = 0.28142083362291853

      # show color legend
      mpas_data_1pvtuDisplay_1.SetScalarBarVisibility(renderView2, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for kineticEnergyCellLUT in view renderView2
      kineticEnergyCellLUTColorBar = GetScalarBar(kineticEnergyCellLUT, renderView2)
      kineticEnergyCellLUTColorBar.Position = [0.5777409860191316, 0.5391485809682804]
      kineticEnergyCellLUTColorBar.Position2 = [0.12, 0.43000000000000005]
      kineticEnergyCellLUTColorBar.Title = 'kineticEnergyCell'
      kineticEnergyCellLUTColorBar.ComponentTitle = '0'
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'LON_LAT_1LAYER-primal': [1, 1, 1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(False, 1)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
