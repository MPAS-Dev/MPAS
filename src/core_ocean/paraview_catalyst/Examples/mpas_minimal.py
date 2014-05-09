try: paraview.simple
except: from paraview.simple import *

################################################################################
# To use this file, copy it to a file named 'mpas.py' in the directory where
# the simulation will be run. Also copy or symlink 'mpas_common.py' beside it.
# Then edit the 'datasets' dictionary according to its comments.
################################################################################

from mpas_common import *

datasets = { }

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
