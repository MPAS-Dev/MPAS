try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing
try:
    from paraview import data_exploration as wx
except:
    wx = None

import math

DEFAULT_VIEW_OPTIONS = {
    'fit_to_screen': 1,
    'magnification': 1.0,
    'width': 960,
    'height': 540,
}
DEFAULT_VIEW_PROPERTIES = {
    'CacheKey': 0.0,
    'StereoType': 0,
    'UseLight': 1,
    'StereoRender': 0,
    'CameraPosition': [15000.0, 99592.9248046875, 647834.944767225],
    'StereoCapableWindow': 0,
    'CameraClippingRange': [441356.5953195527, 1361552.4689387335],
    'LightSwitch': 0,
    'ViewTime': 0.0,
    'Background': [0, 0, 0],
    'CameraFocalPoint': [15000.0, 99592.9248046875, -200000.0],
    'CameraParallelScale': 219435.83080920158,
    'CenterOfRotation': [15000.0, 99592.9248046875, -200000.0],
    'ViewSize': [1000, 500],
    'OrientationAxesVisibility': 0,
    'CenterAxesVisibility': 0
}

if wx is not None:
    def slice_writer(options, fng=None):
        if fng is None:
            fng = wx.FileNameGenerator(options.get('dir', '.'), options['pattern'])

        def create_slice_explorer():
            data = GetActiveSource()
            view = GetActiveView()
            return wx.SliceExplorer(fng, view, data,
                    options['colors'],
                    options.get('slices', 10),
                    options.get('normal', [0, 0, 1]),
                    options.get('viewup', [0, 1, 0]),
                    options.get('bound_range', [0, 1]),
                    options.get('scale_ratio', 2))
        return create_slice_explorer

    def rotate_writer(options, fng=None):
        if fng is None:
            fng = wx.FileNameGenerator(options.get('dir', '.'), options['pattern'])

        def create_rotate_explorer():
            view = GetActiveView()
            return wx.ThreeSixtyImageStackExporter(fng, view,
                    options.get('focal_point', [0, 0, 0]),
                    options.get('distance', 100),
                    options.get('axis', [0, 0, 1]),
                    options.get('step', [10, 15]))
        return create_rotate_explorer

    def contour_explorer(writer, options, fng=None):
        if fng is None:
            fng = wx.FileNameGenerator(options.get('dir', '.'), options['pattern'])

        def create_contour_explorer():
            data = GetActiveSource()
            view = GetActiveView()
            explorer = wx.ContourExplorer(fng, data, options['contours'], options['range'], options['steps'])
            proxy = explorer.getContour()
            rep = Show(proxy)

            rep.LookupTable = options['lut']
            rep.ColorArrayName = options['contours']

            internal = writer(options['inner'], fng)

            class ContourProxy(object):
                def __init__(self, explorer, loopee):
                    self.explorer = explorer
                    self.loopee = loopee

                def add_attribute(self, name, value):
                    setattr(self, name, value)

                def UpdatePipeline(self, time=0):
                    for steps in self.explorer:
                        self.loopee.UpdatePipeline(time)
                    self.explorer.reset()

            return ContourProxy(explorer, internal)
        return create_contour_explorer

def MPASCreateCoProcessor(datasets, options={}):
    freqs = {}
    for (name, dataset) in datasets.items():
        grid = dataset['grid']

        # Prepare the frequency map.
        if grid not in freqs:
            freqs[grid] = []

        # Any frequencies must be known here.
        if 'image_frequency' in dataset:
            freqs[grid].append(dataset['image_frequency'])
        for writer in dataset['writers']:
            freqs[grid].append(writer['frequency'])

    def mpas_create_pipeline(coprocessor, datadescription):
        class MPASPipeline(object):
            for (name, dataset) in datasets.items():
                grid = dataset['grid']

                image_pattern = dataset.get('image_pattern', '%s_%%t.png' % name)
                # Use pi if no images are wanted since it will never have a
                # zero modulo with an integer. Using zero causes
                # ZeroDivisionError exceptions in coprocessing.
                image_frequency = dataset.get('image_frequency', math.pi)

                view_options = DEFAULT_VIEW_OPTIONS.copy()
                if 'view' in options:
                    view_options.update(options['view'])

                view_props = DEFAULT_VIEW_PROPERTIES.copy()
                if 'view_properties' in options:
                    view_props.update(options['view_properties'])

                view = coprocessor.CreateView(CreateRenderView,
                        image_pattern,
                        image_frequency,
                        view_options['fit_to_screen'],
                        view_options['magnification'],
                        view_options['width'],
                        view_options['height'])
                for (k, v) in view_props.items():
                    setattr(view, k, v)

                filters = {}
                producer = coprocessor.CreateProducer(datadescription, grid)
                SetActiveSource(producer)
                filters['simulation'] = producer

                fields = dataset.get('fields', [])
                if fields:
                    filters['simulation_orig'] = producer
                    filters['simulation'] = PassArrays(CellDataArrays=fields)

                if 'filters' in dataset:
                    for filter_desc in dataset['filters']:
                        args = filter_desc.get('args', [])
                        kwargs = filter_desc.get('kwargs', {})
                        filt = filter_desc['function'](*args, **kwargs)
                        if 'source' in filter_desc:
                            source = filter_desc['source']
                            if source in filters:
                                SetActiveSource(filters[source])
                            else:
                                raise RuntimeError('Unknown source for filter: %s' % source)
                        if 'name' in filter_desc:
                            filters[filter_desc['name']] = filt
                        if 'show' in filter_desc:
                            if filter_desc['show']:
                                Show(filt)
                            else:
                                Hide(filt)
                        if 'properties' in filter_desc:
                            for (k, v) in filter_desc['properties'].items():
                                setattr(filt, k, v)

                for writer in dataset['writers']:
                    if 'source' in writer:
                        source = writer['source']
                        if source in filters:
                            SetActiveSource(filters[source])
                        else:
                            raise RuntimeError('Unknown source for web view: %s' % source)
                    pattern = writer.get('pattern', '')
                    writer_obj = coprocessor.CreateWriter(writer['function'], pattern, writer['frequency'])
                    if 'properties' in writer:
                        for (k, v) in writer['properties'].items():
                            setattr(writer_obj, k, v)

        return MPASPipeline()

    class MPASCoProcessor(coprocessing.CoProcessor):
        def CreatePipeline(self, datadescription):
            self.Pipeline = mpas_create_pipeline(self, datadescription)

    coprocessor = MPASCoProcessor()
    coprocessor.SetUpdateFrequencies(freqs)

    # Enable Live-Visualizaton with ParaView
    coprocessor.EnableLiveVisualization(False)

    return coprocessor

# ---------------------- Data Selection method ----------------------

def MPASRequestDataDescription(coprocessor, datadescription):
    "Callback to populate the request for current timestep"
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

def MPASDoCoProcessing(coprocessor, datadescription):
    "Callback to do co-processing for current timestep"
    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
