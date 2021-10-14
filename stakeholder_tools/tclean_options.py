import yaml
import param
import numpy as np

class TCleanOptionsBaseClass(param.Parameterized):
    vis = param.String(default="E2E6.1.00034.S_tclean.ms", doc="Measurement file.") 
    imagename = param.String(default="standard_cube", doc="Name stub of output file.") 
    imsize = param.Array(np.array([80, 80]))
    cell = param.String(default="1.1arcsec", doc="Cell size") 
    specmode = param.String(default="cube", doc="Specmode") 
    interpolation = param.String(default="nearest", doc="") 
    nchan = param.Integer(508, bounds=(1, 580))
    start = param.String(default="220.2526743594GHz", doc="") 
    width = param.String(default="0.2441741MHz", doc="") 
    pblimit = param.Number(0.2)
    deconvolver = param.String(default="hogbom", doc="") 
    niter = param.Integer(1, bounds=(0, None))
    cyclefactor = param.Integer(2, bounds=(1, 50))
    scales = param.List([0, 3, 10])
    interactive = param.Integer(0, bounds=(0,1), doc="Interactive mode")
    field = param.String(default='1')
    spw = param.ListSelector(default=['0'])
    antenna = param.ListSelector(default=['0,1,2,3,4,5,6,7,8'])
    scan = param.ListSelector(default=['8,12,16'])
    intent = param.List(['OBSERVE_TARGET#ON_SOURCE'])
    datacolumn = param.String(default='data')
    phasecenter = param.String(default='ICRS 00:45:54.3836 -073.15.29.413')
    stokes= param.String(default='I')
    outframe = param.String(default='LSRK') 
    perchanweightdensity= param.Boolean(False)
    gridder = param.String(default='standard')  
    mosweight = param.Boolean(False)
    usepointing = param.Boolean(False) 
    restoration = param.Boolean(False)
    pbcor = param.Boolean(False) 
    weighting = param.String(default='briggs') 
    restoringbeam = param.List(['common'])
    robust = param.Number(default=0.5) 
    npixels = param.Integer(default=0) 
    threshold = param.String(default='0.0mJy')
    nsigma = param.Number(default=0.0)
    usemask = param.String(default='auto-multithresh')
    sidelobethreshold = param.Number(1.25)
    noisethreshold = param.Number(5.0)
    lownoisethreshold = param.Number(2.0)
    negativethreshold = param.Number(0.0) 
    minbeamfrac = param.Number(0.1)
    growiterations = param.Integer(75)
    dogrowprune = param.Boolean(True)
    minpercentchange = param.Number(1.0)
    fastnoise = param.Boolean(False)
    savemodel = param.String(default='none')
    parallel = param.Boolean(False)
    verbose = param.Boolean(True)
    restart = param.Boolean(False)
    calcres = param.Boolean(True)
    calcpsf = param.Boolean(True)

    config_file = param.String(default="config/tclean.yaml")
    
    def read_configuration(self, config:str)->None:
        with open(config) as file:
            config = yaml.full_load(file)
        
            self.niter = config['tclean-second-cycle']['niter']
            self.pbcor = config['tclean-second-cycle']['pbcor']
            self.restart = config['tclean-second-cycle']['restart']
            self.calcres = config['tclean-second-cycle']['calcres']
            self.calcpsf = config['tclean-second-cycle']['calcpsf']
            self.threshold = config['tclean-second-cycle']['threshold']
            self.restoration = config['tclean-second-cycle']['restoration']
