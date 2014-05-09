
# write MPAS output to vtk format (*.pvtu file)

try: paraview.simple
except: from paraview.simple import *

################################################################################
# To use this file, copy it to a file named 'mpas.py' in the directory where
# the simulation will be run. Also copy or symlink 'mpas_common.py' beside it.
# Then edit the 'datasets' dictionary according to its comments.
################################################################################

from mpas_common import *



datasets = {
     'MPAS_OUTPUT': {
         # choose grid type:
         # X_Y: Cartesian domain 
         # X_Y_Z: Spherical domain, on globe
         # LON_LAT: Spherical domain, Mercator projection
         # primal: cell-centered variables
         # dual: variables on vertices
         # 1LAYER: One z-layer is visualized at a time
         # NLAYER: image is a 3D view of all layers together.
         #'grid': 'X_Y_NLAYER-primal',
         #'grid': 'X_Y_NLAYER-dual',
         #'grid': 'X_Y_Z_1LAYER-primal',
         #'grid': 'X_Y_Z_1LAYER-dual',
         #'grid': 'X_Y_Z_NLAYER-primal',
         #'grid': 'X_Y_Z_NLAYER-dual',
         'grid': 'LON_LAT_1LAYER-primal',
         #'grid': 'LON_LAT_1LAYER-dual',
         #'grid': 'LON_LAT_NLAYER-primal',
         #'grid': 'LON_LAT_NLAYER-dual',
 
         # fields to output
         'fields': ['salinity', 'temperature'],
 
         'writers': [
             {
                 # filename for output. %t will be replaced with timestep
                 'pattern': 'one_%t.pvtu',
 
                 # how often to produce output (1 is every timestep)
                 'frequency': 2,
 
                 # do not change
                 'source': 'simulation',
                 'function': XMLPUnstructuredGridWriter,
             }
         ]
     },

     'MPAS_OUTPUT2': {
         #'grid': 'X_Y_NLAYER-primal',
         #'grid': 'X_Y_NLAYER-dual',
         #'grid': 'X_Y_Z_1LAYER-primal',
         #'grid': 'X_Y_Z_1LAYER-dual',
         #'grid': 'X_Y_Z_NLAYER-primal',
         #'grid': 'X_Y_Z_NLAYER-dual',
         'grid': 'LON_LAT_1LAYER-primal',
         #'grid': 'LON_LAT_1LAYER-dual',
         #'grid': 'LON_LAT_NLAYER-primal',
         #'grid': 'LON_LAT_NLAYER-dual',

         # fields to output
         'fields': ['salinity', 'temperature'],
 
         'writers': [
             {
                 # filename for output. %t will be replaced with timestep
                 'pattern': 'two_%t.pvtu',

                 # how often to produce output (1 is every timestep)
                 'frequency': 1,

                 # do not change
                 'source': 'simulation',
                 'function': XMLPUnstructuredGridWriter,
             }
         ]
     },

}

coprocessor = MPASCreateCoProcessor(datasets)

# To use other scripts, you may import them here:
scripts = []
modules = []
for script in scripts:
    modules.append(importlib.import_module(script))

# These functions is required and is called from Catalyst without arguments.
# Instead, pass the datasets we want to export to MPASCreateCoProcessor.
def RequestDataDescription(datadescription):
    MPASRequestDataDescription(coprocessor, datadescription)
    for module in modules:
        module.RequestDataDescription(datadescription)

def DoCoProcessing(datadescription):
    MPASDoCoProcessing(coprocessor, datadescription)
    for module in modules:
        module.DoCoProcessing(datadescription)
