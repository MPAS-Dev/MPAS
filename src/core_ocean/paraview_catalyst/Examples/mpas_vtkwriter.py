
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
         # the grid to output
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
         'fields': ['salinity', 'temperature', 'kineticEnergyCell', 'relativeVorticityCell'],
 
         'writers': [
             {
                 # filename for output. %t will be replaced with timestep
                 'pattern': 'mpas_data_%t.pvtu',
 
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
