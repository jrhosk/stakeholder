import panel as pn
import os
import sys
from bokeh.models.annotations import Title
import numpy as np

from bokeh.layouts import column, row, Spacer
from bokeh.models import CustomJS, Slider, RadioButtonGroup, TextInput, Button, MultiChoice
from bokeh.models import BoxAnnotation, PreText, Range1d, LinearAxis, Span, HoverTool, DataTable, TableColumn
from bokeh.events import SelectionGeometry
from bokeh.plotting import ColumnDataSource, figure, show

from . import notebook_testing as nt
from . import tclean_options 

class TCleanPanel(tclean_options.TCleanOptionsBaseClass):
    
    def __init__(self, terminal=False, config_file=os.path.dirname(os.path.realpath(__file__)) + '/config/tclean.yaml'):
        self.terminal = terminal

        self.config_file = config_file
        self.standard = nt.Test_standard()
        self.standard.setUp()
        self.standard.config_file = config_file

        self.layout = pn.template.GoldenTemplate(
            title="Interactive Stakeholder Tests"
        )
        
        # Image

        self.file_widget = pn.widgets.FileSelector(os.getcwd(), name="File Selector")
        self.file_widget.param.watch(self._update_file, 'value')
        
        self.imsize_widget = pn.Param(
            self.standard.param.imsize, 
            widgets={
                'imsize': pn.widgets.LiteralInput
        })
        
        self.npixels_widget = pn.Param(
            self.standard.param.npixels, 
            widgets={
                'npixels': pn.widgets.IntInput
        })

        # Analysis

        self.cell_widget = pn.Param(
            self.standard.param.cell, 
            widgets={
                'cell': pn.widgets.TextInput
        })
        
        self.specmode_widget = pn.Param(
            self.standard.param.specmode, 
            widgets={
                'specmode': {'widget_type': pn.widgets.Select, 'options':['mfs', 'cube', 'cubedata']}
        })

        self.interpolation_widget = pn.Param(
            self.standard.param.interpolation, 
            widgets={
                'Interpolation': {'widget_type': pn.widgets.Select, 'options': ['nearest', 'linear', 'cubic']}
        })

        self.start_widget = pn.Param(
            self.standard.param.start, 
            widgets={
                'start': pn.widgets.TextInput
        })
        
        self.width_widget = pn.Param(
            self.standard.param.width, 
            widgets={
                'width': pn.widgets.TextInput
        })

        self.pblimit_widget = pn.Param(
            self.standard.param.pblimit, 
            widgets={
                'pblimit': pn.widgets.FloatInput
        })

        self.deconvolver_widget = pn.Param(
            self.standard.param.deconvolver, 
            widgets={
                'deconvolver': {'widget_type': pn.widgets.Select, 'options': ['hogbom', 'clark', 'multiscale', 'mtmfs', 'mem', 'clarkstokes']}
        })


        self.scales_widget = pn.Param(
            self.standard.param.scales, 
            widgets={
                'scales': {'widget_type': pn.widgets.LiteralInput, 'type':  list}
        })

        self.field_widget = pn.Param(
            self.standard.param.field, 
            widgets={
                'field': pn.widgets.TextInput
        })

        self.spw_widget = pn.Param(
            self.standard.param.spw, 
            widgets={
                'spw': {'widget_type': pn.widgets.LiteralInput, 'type':  list}
        })

        self.antenna_widget = pn.Param(
            self.standard.param.antenna, 
            widgets={
                'antenna': {'widget_type': pn.widgets.LiteralInput, 'type':  list}
        })

        self.scan_widget = pn.Param(
            self.standard.param.scan, 
            widgets={
                'scan': {'widget_type': pn.widgets.LiteralInput, 'type':  list}
        })

        self.intent_widget = pn.Param(
            self.standard.param.intent, 
            widgets={
                'intent': {'widget_type': pn.widgets.LiteralInput, 'type':  list}
        })

        self.datacolumn_widget = pn.Param(
            self.standard.param.datacolumn, 
            widgets={
                'datacolumn': pn.widgets.TextInput,
        })

        self.stokes_widget = pn.Param(
            self.standard.param.stokes, 
            widgets={
                'stokes': pn.widgets.TextInput,
        })

        self.outframe_widget = pn.Param(
            self.standard.param.deconvolver, 
            widgets={
                'outframe': {'widget_type': pn.widgets.Select, 'options': ['LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP', 'CMB']}
        })

        self.gridder_widget = pn.Param(
            self.standard.param.deconvolver, 
            widgets={
                'gridder': {'widget_type': pn.widgets.Select, 'options': ['standard', 'wproject', 'widefield', 'mosaic', 'awproject']}
        })

        self.weighting_widget = pn.Param(
            self.standard.param.weighting, 
            widgets={
                'weighting': {'widget_type': pn.widgets.Select, 'options': ['natural', 'uniform', 'briggs', 'superuniform', 'radial', 'briggsabs', 'briggsbwtaper']}
        })

        self.restoringbeam_widget = pn.Param(
            self.standard.param.restoringbeam, 
            widgets={
                'restoringbeam': {'widget_type': pn.widgets.LiteralInput, 'type':  list}
        })

        self.robust_widget = pn.Param(
            self.standard.param.robust, 
            widgets={
                'robust': pn.widgets.FloatInput,
        })

        self.threshold_widget = pn.Param(
            self.standard.param.threshold, 
            widgets={
                'threshold': pn.widgets.TextInput,
        })

        self.nsigma_widget = pn.Param(
            self.standard.param.nsigma, 
            widgets={
                'nsigma': pn.widgets.FloatInput,
        })

        self.usemask_widget = pn.Param(
            self.standard.param.usemask, 
            widgets={
                'usemask': pn.widgets.TextInput,
        })

        self.sidelobethreshold_widget = pn.Param(
            self.standard.param.sidelobethreshold, 
            widgets={
                'sidelodethreshold': pn.widgets.FloatInput,
        })

        self.noisethreshold_widget = pn.Param(
            self.standard.param.noisethreshold, 
            widgets={
                'noisethreshold': pn.widgets.FloatInput,
        })

        self.lownoisethreshold_widget = pn.Param(
            self.standard.param.lownoisethreshold, 
            widgets={
                'lownoisethreshold': pn.widgets.FloatInput,
        })

        self.negativethreshold_widget = pn.Param(
            self.standard.param.negativethreshold, 
            widgets={
                'negativethreshold': pn.widgets.FloatInput,
        })

        self.minbeamfrac_widget = pn.Param(
            self.standard.param.minbeamfrac, 
            widgets={
                'minbeamfrac': pn.widgets.FloatInput,
        })

        self.growiterations_widget = pn.Param(
            self.standard.param.growiterations, 
            widgets={
                'growiterations': pn.widgets.IntInput,
        })

        self.minpercentchange_widget = pn.Param(
            self.standard.param.minpercentchange, 
            widgets={
                'minpercentchange': pn.widgets.FloatInput,
        })

        self.savemodel_widget = pn.Param(
            self.standard.param.savemodel, 
            widgets={
                'savemodel': pn.widgets.TextInput,
        })

        # Iteration

        self.nchan_widget = pn.Param(
            self.standard.param.nchan, 
            widgets={
                'nchan': pn.widgets.IntSlider
            })
        
        self.niter_widget = pn.Param(
            self.standard.param.niter, 
            widgets={
                'niter': pn.widgets.IntInput,
        })

        self.terminal_widget = pn.widgets.Terminal(
            "CASA TClean Terminal\n\n",
            options={"cursorBlink": True},
            height=1100, width=900,
            name='Terminal'
        )

        # Boolean

        self.interactive_widget = pn.widgets.Toggle(
            name='Interactive', 
            button_type='primary'
        )
        
        self.interactive_widget.param.watch(self._set_interactive, 'value')

        self.perchanweightdensity_widget = pn.widgets.Toggle(
            name='perchanweightdensity', 
            button_type='primary'
        )
        self.perchanweightdensity_widget.param.watch(self._set_perchanweightdensity, 'value')

        self.mosweight_widget = pn.widgets.Toggle(
            name='mosweight', 
            button_type='primary'
        )
        self.mosweight_widget.param.watch(self._set_mosweight, 'value')

        self.usepointing_widget = pn.widgets.Toggle(
            name='usepointing', 
            button_type='primary'
        )
        self.usepointing_widget.param.watch(self._set_usepointing, 'value')

        self.restoration_widget = pn.widgets.Toggle(
            name='restoration', 
            button_type='primary'
        )
        self.restoration_widget.param.watch(self._set_restoration, 'value')

        self.dogrowprune_widget = pn.widgets.Toggle(
            name='dogrowprune', 
            button_type='primary'
        )
        self.dogrowprune_widget.param.watch(self._set_dogrowprune, 'value')

        self.fastnoise_widget = pn.widgets.Toggle(
            name='fastnoise', 
            button_type='primary'
        )
        self.fastnoise_widget.param.watch(self._set_fastnoise, 'value')

        self.parallel_widget = pn.widgets.Toggle(
            name='parallel', 
            button_type='primary'
        )
        self.parallel_widget.param.watch(self._set_parallel, 'value')

        self.verbose_widget = pn.widgets.Toggle(
            name='verbose', 
            button_type='primary'
        )
        self.verbose_widget.param.watch(self._set_verbose, 'value')

        # ------------------------------------ #
    
        self.play_button = pn.widgets.Button(
            name="Play", 
            button_type="success",  
            width=300
        )

        self.exit_button = pn.widgets.Button(
            name="Exit", 
            button_type="danger",  
            width=300
        )    
        
        self.play_button.on_click(self._clean)
        self.exit_button.on_click(self._exit)
    
    
        image_controls = pn.Column(
            self.imsize_widget,
            self.cell_widget,
            self.specmode_widget,
            self.npixels_widget,
            self.start_widget,
            self.width_widget,
            self.stokes_widget,
            self.nchan_widget,
            self.niter_widget,
        )

        analysis_controls = pn.Column(
            pn.Row(
                pn.Column(
                    self.interpolation_widget,
                    self.deconvolver_widget,
                    self.scales_widget,
                    self.field_widget,
                    self.spw_widget,
                    self.antenna_widget,
                    self.scan_widget,
                    self.intent_widget,
                    self.datacolumn_widget,
                    self.outframe_widget,
                    self.gridder_widget,
                ),
                pn.Column(
                    self.weighting_widget,
                    self.restoringbeam_widget,
                    self.robust_widget,
                    self.threshold_widget,
                    self.nsigma_widget,
                    self.usemask_widget,
                    self.sidelobethreshold_widget,
                    self.noisethreshold_widget,
                    self.lownoisethreshold_widget,
                    self.negativethreshold_widget,
                    self.minbeamfrac_widget,
                    self.growiterations_widget,
                    self.savemodel_widget
                )
            ), 
        )

        boolean_controls = pn.Column(
            self.perchanweightdensity_widget,
            self.mosweight_widget,
            self.usepointing_widget,
            self.dogrowprune_widget,
            self.fastnoise_widget,
            self.parallel_widget,
            self.verbose_widget,
            
        )

        if self.terminal is True:
            sys.stdout = self.terminal_widget
            sys.stderr = self.terminal_widget
            self.layout.sidebar.append(self.file_widget)
            self.layout.sidebar.append(pn.Tabs(
                    ('Image', image_controls), 
                    ('Analysis', analysis_controls),
                    ('Boolean', boolean_controls)
                )
            )
            self.layout.main.append(pn.Column(
                self.terminal_widget,
                    pn.Row(
                        self.play_button, 
                        self.exit_button
                    )
                )
            )
            
            self.layout.show()
        else:
            self.play_button.visible = False
            self.exit_button.visible = False
            
            self.layout = pn.Row(
            pn.Column(
                pn.Card(
                    self.file_widget,
                    width=900,
                    header_background = '#21618C',
                    header_color = 'white',
                    title='File Selector'),
                pn.Card(
                    pn.Tabs(
                        ('Image', image_controls),
                        ('Analysis', analysis_controls),
                        ('Boolean', boolean_controls)),
                    width=900,
                    header_background=' #21618C',
                    header_color = 'white',
                    title='TClean Controls'
                ),
                pn.Row(
                    self.play_button, self.exit_button,
                )
            )
        )
            

    # Private utility functions
    
    def _update_file(self, event):
            self.standard.vis = self.file_widget.value[0]
            print('Selected file: ' + self.standard.vis)

    def _clean(self, event):
        self.standard.test_standard_cube()

    def _exit(self, event):
        sys.exit('Exiting interactive stakeholder test.')
    
    def _set_interactive(self, event):
            if self.interactive_widget.value is True:
                self.standard.interactive = 1
            elif self.interactive_widget.value is False:
                self.standard.interactive = 0
            else:
                pass
    
    def _set_perchanweightdensity(self, event):
        self.standard.perchanweightdensity = self.perchanweightdensity_widget.value

    def _set_mosweight(self, event):
        self.standard.mosweight = self.mosweight_widget.value

    def _set_usepointing(self, event):
        self.standard.usepointing = self.usepointing_widget.value

    def _set_restoration(self, event):
        self.standard.restoration = self.restoration_widget.value

    def _set_dogrowprune(self, event):
        self.standard.dogrowprune = self.dogrowprune_widget.value

    def _set_fastnoise(self, event):
        self.standard.fastnoise = self.fastnoise_widget.value

    def _set_parallel(self, event):
        self.standard.parallel = self.parallel_widget.value

    def _set_verbose(self, event):
        self.standard.verbose = self.verbose_widget.value

    # Public utility functions

    def clean(self):
        self.standard.test_standard_cube()

    def show(self):
        if self.terminal is True:
            return self.layout.show()
        else:
            return self.layout
        
