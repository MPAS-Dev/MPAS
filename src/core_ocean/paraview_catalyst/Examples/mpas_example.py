try: paraview.simple
except: from paraview.simple import *

################################################################################
# To use this file, copy it to a file named 'mpas.py' in the directory where
# the simulation will be run. Also copy or symlink 'mpas_common.py' beside it.
# Then edit the 'datasets' dictionary according to its comments.
################################################################################

from mpas_common import *

# A dictionary of datasets to export is used to determine what is wanted when
# running a simulation. Each entry in the dictionary describes a grid to export
# along with information describing how it should be exported.
datasets = {
    # This is the name associated with the output. It is used in the default
    # filename patterns if another name is not given.
    'EXAMPLE': {
        ########################################################################
        # This entry is an example only and will be removed if it not the
        # only one in the map so that it may be used as a reference.
        # Renaming it will also preserve it.
        ########################################################################

        ########################################################################
        # Core information
        #-----------------------------------------------------------------------
        # REQUIRED: This is required and is the name of the grid to export from
        # the simulation. Each entry may only export a single grid, but the
        # same grid may be used multiple times in different outputs.
        'grid': 'LON_LAT_1LAYER-primal',
        # OPTIONAL: The list of fields to use from the generated data. If
        # empty, all fields will be used (default is empty). These are case
        # sensitive.
        'fields': ['salinity', 'temperature'],

        ########################################################################
        # Filters
        #-----------------------------------------------------------------------
        # Filters may be used to transform data and as input for writers or
        # other filters.
        'filters': [
            {
                # The source to use for this filter. The top-level source is
                # named "simulation". It is recommended to be explicit and
                # always specify the source.
                'source': 'simulation',
                # The name of the filter. This may be used as a 'source' for
                # other filters and writers.
                'name': 'myfilter',
                # If True, the filter will be rendered in the window. This is
                # useful for filling the background or loading a base model to
                # render on top of.
                'show': True,
                # The constructor for the filter. A VTK class will typically be
                # given here.
                'function': LegacyVTKReader,
                # Any keyword arguments which must be passed to the constructor.
                'kwargs': {
                    'FileNames': ['my_base_model.vtk']
                }
            }
        ],

        ########################################################################
        # Writers
        #-----------------------------------------------------------------------
        # These are added as writers and updated as necessary. Multiple writers
        # may be used. There must be at least one writer.
        'writers': [
            {
                # The source to use for this writer. The top-level source is
                # named "simulation". It is recommended to be explicit and
                # always specify the source.
                'source': 'simulation',
                # The function to construct an instance of the writer class.
                # A VTK class will typically be given here.
                'function': XMLUnstructuredGridWriter,
                # The filename pattern to use when writing grid files. Use '%t'
                # to use the timestep of the grid.
                'pattern': 'example_lonlat1_%t.pvtu',
                # How often to use this writer.
                'frequency': 5
            }
            # Further writers may go here.
        ]
    }
}

# Remove the example from the datasets in case it is wanted for documentation.
if len(datasets) > 1 and 'EXAMPLE' in datasets:
    del datasets['EXAMPLE']

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
