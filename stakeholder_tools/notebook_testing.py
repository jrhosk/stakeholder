##########################################################################
##########################################################################
# test_stk_alma_pipeline_imaging.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA.
#
# This script is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# [https://open-jira.nrao.edu/browse/CAS-12428]
#
#
##########################################################################

'''
Datasets (MOUS)
E2E6.1.00034.S (uid://A002/Xcff05c/X1ec)
2018.1.00879.S (uid://A001/X133d/X169f)
E2E6.1.00020.S (uid://A002/Xcff05c/Xe5)
2017.1.00750.T (uid://A001/X131b/X57)

Test list - 21 total
1a.  Single field(SF) cube with perchanweightdensity=False(pcwdF), weighting=briggs - E2E6.1.00034.S
1b.  SF cube with pcwdT, weighting=briggs - E2E6.1.00034.S
1c.  SF cube with pcwdT, weighting=briggsbwtaper - E2E6.1.00034.S
2.   SF MFS - E2E6.1.00020.S
3.   SF mtmfs - E2E6.1.00020.S
4a.  SF ephemeris cube (multi-EB) with pcwdF+briggs - 2017.1.00750.T
4b.  SF ephemeris cube (multi-EB) with pcwdT+briggs - 2017.1.00750.T
4c.  SF ephemeris cube (multi-EB) with pcwdT+briggsbwtaper - 2017.1.00750.T
5.   SF ephemeris MFS - 2018.1.00879.S
6.   SF ephemeris mtmfs - 2018.1.00879.S
7.   SF Calibrator - E2E6.1.00034.S
8.   SF ephemeris Calibrator - 2018.1.00879.S  
9a.  Mosaic cube with pcwdF, briggs- E2E6.1.00034.S
9b.  Mosaic cube with pcwdT+brigs- E2E6.1.00034.S
9c.  Mosaic cube with pcwdT+briggsbwtaper- E2E6.1.00034.S
10.  Mosaic MFS - E2E6.1.00020.S
11.  Mosaic mtmfs - E2E6.1.00020.S
12a. Mosaic ephemeris cube with pcwdF- 2018.1.00879.S
12b. Mosaic ephemeris cube with pcwdT+briggs - 2018.1.00879.S
12c. Mosaic ephemeris cube with pcwdT+briggsbwtaper - 2018.1.00879.S
13.  Mosaic ephemeris MFS - 2018.1.00879.S
14.  Mosaic ephemeris mtmfs - 2018.1.00879.S

Each test stores reference values in dictionaries for the metrics 
to be tested:
The following variable names are used:
the standard default sets are,
exp_im_stats, exp_mask_stats, exp_pb_stats, exp_psf_stats, 
exp_model_stats, exp_resid_stats, exp_sumwt_stats
exp_wt_stats (for mosaic)
Addtionally, for cube imaging (if self.parallel=True), 
exp_bmin_dict, exp_bmaj_dict, exp_pa_dict
And for mtmfs
exp_im1_stats, exp_model1_stats, exp_resid1_stats, exp_sumwt1_stats

'''

import os
import glob
import sys
import subprocess
import unittest
import numpy
import shutil
import inspect
import scipy
import matplotlib.pyplot as pyplot
import json
import casatools
import casatasks
import copy

from casatestutils.imagerhelpers import TestHelpers

th = TestHelpers()

from casatestutils import generate_weblog
from casatestutils import add_to_dict
from casatestutils import stats_dict

from .tclean_options import TCleanOptionsBaseClass

CASA6 = False
try:
    from casatools import ctsys, quanta, measures, image, vpmanager, calibrater
    from casatasks import casalog, delmod, imsubimage, tclean, uvsub, imhead, imsmooth, immath, widebandpbcor, immoments
    from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
    from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper
    from casaviewer import imview

    CASA6 = True
    _ia = image()
    ctsys_resolve = ctsys.resolve

except ImportError:
    from __main__ import default  # reset given task to its default values
    from tasks import *           # Imports all casa tasks
    from taskinit import *        # Imports all casa tools
    from parallel.parallel_task_helper import ParallelTaskHelper

    _ia = iatool()
    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0], 'casatestdata/')
        return os.path.join(dataPath,apath)

    #data_path = ctsys_resolve('./') This should be added back for unittest
data_path = './data/'

    # save the dictionaries of the metrics to files (per tests) 
    # mostly useful for the maintenace (updating the expected metric i
    # parameters based
    # on the current metrics)
savemetricdict=True

class test_tclean_base(unittest.TestCase, TCleanOptionsBaseClass):
    def setUp(self):
        self._myia = _ia
        self.epsilon = 0.01 # sets epsilon as a percentage (1%)
        self.msfile = ""
        self.img_subdir = 'testdir'
        self.parallel = False
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True

        print('Finished CasaTest initialization.')

    def tearDown(self):
        #print("TEST_DICT=",test_dict)
        generate_weblog("tclean_ALMA_pipeline",test_dict)
        print("Closing ia tool")
        self._myia.done()
        """ don't delete it all """
#        self.delData()

    # Separate functions here, for special-case tests that need their own MS.
    def prepData(self, msname=[""]):
        if msname != [""]:
            self.msfile=msname

    def delData(self, msname=[""]):
        del_files = [self.img_subdir]
        if msname != [""]:
            self.msfile=msname
    
        if (os.path.exists(self.msfile)):
            del_files.append(self.msfile)
      
        img_files = glob.glob(self.img+'*')
        del_files += img_files
      
        for f in del_files:
            shutil.rmtree(f)

    def prepInputmask(self, maskname=""):
        if maskname!="":
            self.maskname=maskname
        if (os.path.exists(self.maskname)):
            shutil.rmtree(self.maskname)
    
        shutil.copytree(refdatapath+self.maskname, self.maskname, symlinks=True)

    def checkfinal(self, pstr=""):
        pstr += "["+inspect.stack()[1][3]+"] : To re-run this test : " 
        "runUnitTest.main(['test_tclean["+ inspect.stack()[1][3] +"]'])"
        
        casalog.post(pstr,'INFO')
        if(pstr.count("(Fail") > 0 ):
            self.fail("\n"+pstr)

    def check_dict_vals_beam(self, exp_dict, act_dict, suffix, epsilon=0.01):
        """ Compares expected dictionary with actual dictionary. Useful for comparing the restoring beam.

            Parameters
            ----------
            exp_dict: dictionary
                Expected values, as key:value pairs.
                Keys must match between exp_dict and act_dict.
                Values are compared between exp_dict and act_dict. A summary
                line is returned with the first mismatched value, or the last
                successfully matched value.
            act_dict: dictionary
                Actual values to compare to exp_dict (and just the values).
            suffix: string
                For use with summary print statements.
        """

        report = ''
        eps = epsilon
        passed = True
        chans = 0
    
        for key in exp_dict:
            result = th.check_val(act_dict[key], exp_dict[key],
                                  valname=suffix+' chan'+str(chans), epsilon=eps)[1]
            chans += 1
            if 'Fail' in result:
                passed = False
                break
        
        report += th.check_val(passed, True, valname=suffix+' chan'+str(chans), exact=True)[1]

        return report

    def copy_products(self, old_pname, new_pname, ignore=None):
        """ function to copy iter0 images to iter1 images
            (taken from pipeline)
        """

        imlist = glob.glob('%s.*' % old_pname)
        imlist = [xx for xx in imlist if ignore is None or ignore not in xx]
        
        for image_name in imlist:
            newname = image_name.replace(old_pname, new_pname)
      
            if image_name == old_pname + '.workdirectory':
                mkcmd = 'mkdir '+ newname
                os.system(mkcmd)
                self.copy_products(os.path.join(image_name, old_pname), 
                                   os.path.join(newname, new_pname))
            else:
                shutil.copytree(image_name, newname, symlinks=True)

    def cube_beam_stats(self, image):
        """ function to return per-channel beam statistics
            will be deprecated and combined into image_stats 
            once CASA beam issue is fixed
        """

        self._myia.open(image)

        bmin_dict = {} 
        bmaj_dict = {} 
        pa_dict = {}

        beam_dict = self._myia.restoringbeam()['beams']
        for item in beam_dict.keys():
            bmin_dict[item] = beam_dict[item]['*0']['minor']['value']
            bmaj_dict[item] = beam_dict[item]['*0']['major']['value']
            pa_dict[item] = beam_dict[item]['*0']['positionangle']['value']

        self._myia.close()

        return bmin_dict, bmaj_dict, pa_dict

    def image_stats(self, image, fit_region=None, field_regions=None, masks=None):
        """ function that takes an image file and returns a statistics
            dictionary
        """

        self._myia.open(image)
        imagename=os.path.basename(image)
        stats_dict = {}

        statistics = self._myia.statistics()
            # Return data chunk; transpose to make channel selection easier
        chunk = numpy.transpose(self._myia.getchunk(dropdeg=True))

            # stats returned for all images
        im_size = self._myia.boundingbox()['imageShape'].tolist()
        stats_dict['npts'] = im_size[0]*im_size[1]*im_size[3]
        stats_dict['npts_unmasked'] = statistics['npts'][0]
        stats_dict['npts_real'] = numpy.count_nonzero(~numpy.isnan(chunk))
        stats_dict['freq_bin'] = self._myia.summary()['incr'][3]
        stats_dict['start'] = float( 
            statistics['blcf'].split(', ')[3].split('Hz')[0])
        stats_dict['end'] = float( 
            statistics['trcf'].split(', ')[3].split('Hz')[0])
        stats_dict['start_delta'] = stats_dict['start']
        stats_dict['end_delta'] = stats_dict['end']
        stats_dict['nchan'] = im_size[3]


        # stats returned for all images except .mask
        if not image.endswith('.mask'):
            stats_dict['max_val'] = statistics['max'][0]
            stats_dict['max_val_pos'] = statistics['maxpos'].tolist()
            max_loc = [stats_dict['max_val_pos'][0], 
                       stats_dict['max_val_pos'][1]]
            stats_dict['min_val'] = statistics['min'][0]
            stats_dict['min_val_pos'] = statistics['minpos'].tolist()
            stats_dict['im_rms'] = statistics['rms'][0]

            # stats returned if a region file is given
        if fit_region != None:
            if '_cube' in imagename:
                if '.pb' in imagename:
                    fit_region = fit_region + ', range=[%schan,%schan]' % (int(im_size[3]/2), int(im_size[3]/2))
                
            if '.psf' in imagename:
                    # using chan 1 as first because ia.fitcomponents fits
                    # every channel if chan=0
                fit_regions = [(fit_region + ', range=[%schan,%schan]' % (1, 1)), 
                               (fit_region + ', range=[%schan,%schan]' % (int(im_size[3]/2), 
                                                                          int(im_size[3]/2))), 
                                (fit_region + ', range=[%schan,%schan]' % ((im_size[3]-1), 
                                                                           (im_size[3]-1)))]
                i = 0
                for region in fit_regions:
                    try:
                        fit_dict = self._myia.fitcomponents( 
                        region=region)['results']['component0']
                
                        stats_dict['fit_'+str(i)] = [fit_dict['peak']['value'], 
                                             fit_dict['shape']['majoraxis']['value'], 
                                             fit_dict['shape']['minoraxis']['value']]
                        stats_dict['fit_loc_chan_'+str(i)] = fit_dict['spectrum']['channel']
                        stats_dict['fit_loc_freq_'+str(i)] = fit_dict['spectrum']['frequency']['m0']['value']
                        stats_dict['fit_pix_'+str(i)] = fit_dict['pixelcoords'].tolist()
                    except KeyError:
                        stats_dict['fit_'+str(i)] = [1.0, 1.0, 1.0]
                        stats_dict['fit_loc_chan_'+str(i)] = 1.0
                        stats_dict['fit_loc_freq_'+str(i)] = 1.0
                        stats_dict['fit_pix_'+str(i)] = [1.0, 1.0]
                        
                    i += 1
            
                if '.model' in imagename:
                    fit_region = fit_region
                
                if '.model' not in imagename and '.pb' not in imagename and '.psf' not in imagename:
                    # WARN: If max value channel is 0, tool fits all channels
                    fit_region = fit_region + ', range=[%schan,%schan]' % (stats_dict['max_val_pos'][3], 
                                                                           stats_dict['max_val_pos'][3])
            
            if '.psf' in imagename and '_cube' in imagename:
                stats_dict['regn_sum'] = self._myia.statistics(region=fit_regions[1])['sum'][0]
            
            else:
                stats_dict['regn_sum'] = self._myia.statistics(region=fit_region)['sum'][0]
            
            if ('image' in imagename and 'mosaic_cube_eph' not in imagename) or 'pb' in imagename or ('psf' in imagename and 'cube' not in imagename):
                try:
                    fit_dict = self._myia.fitcomponents( region=fit_region)['results']['component0']    
                    stats_dict['fit'] = [
                        fit_dict['peak']['value'], 
                        fit_dict['shape']['majoraxis']['value'], 
                        fit_dict['shape']['minoraxis']['value']
                    ]
                    stats_dict['fit_loc_chan'] = fit_dict['spectrum']['channel']
                    stats_dict['fit_loc_freq'] = fit_dict['spectrum']['frequency']['m0']['value']
                    stats_dict['fit_pix'] = fit_dict['pixelcoords'].tolist()
                
                except KeyError:
                    stats_dict['fit'] = [1.0, 1.0, 1.0]
                    stats_dict['fit_loc_chan'] = 1.0
                    stats_dict['fit_loc_freq'] = 1.0
                    stats_dict['fit_pix'] = [1.0, 1.0]

        # stats returned for .image(.tt0)
        if 'image' in imagename:
            commonbeam = self._myia.commonbeam()
            stats_dict['com_bmin'] = commonbeam['minor']['value']
            stats_dict['com_bmaj'] = commonbeam['major']['value']
            stats_dict['com_pa'] = commonbeam['pa']['value']
      
            if 'cube' in imagename:
                stats_dict['rms_per_chan'] = self._myia.statistics(axes=[0,1])['rms'].tolist()
                stats_dict['profile'] = self.cube_profile_fit( image, max_loc, stats_dict['nchan'])
            if 'mosaic' in imagename:
                stats_dict['rms_per_field'] = []
                for region in field_regions:
                    stats_dict['rms_per_field'].append(self._myia.statistics(region=region)['rms'][0])

            # stats returned if not .pb(.tt0), .sumwt(.tt0), or .mask
            # if 'pb' not in image and 'sumwt' not in image and not image.endswith('.mask'):
        stats_dict['im_sum'] = statistics['sum'][0]

        if image.endswith('.mask'):
            stats_dict['mask_pix'] = numpy.count_nonzero(chunk)
            stats_dict['mask_regns'] = scipy.ndimage.label(chunk)[1]
            stats_dict['mask'] = ~numpy.array(chunk, dtype=bool)

        if 'pb' in imagename:
            pb_mask_02 = chunk>0.2
            pb_mask_05 = chunk>0.5
      
            if 'cube' in image:
                pb_02_list = []
                pb_05_list = []
                i = 0
        
                for chan in chunk:
                    pb_02_list.append(numpy.count_nonzero(chan*pb_mask_02[i]))
                    pb_05_list.append(numpy.count_nonzero(chan*pb_mask_05[i]))
                    i += 1
                stats_dict['npts_0.2'] = pb_02_list
                stats_dict['npts_0.5'] = pb_05_list
            
            else:
                stats_dict['npts_0.2'] = numpy.count_nonzero(pb_mask_02)
                stats_dict['npts_0.5'] = numpy.count_nonzero(pb_mask_05)
      
            if 'mosaic' in imagename:
                stats_dict['pb_mask_0.2'] = pb_mask_02
                stats_dict['pb_mask_0.5'] = pb_mask_05

        if 'model' in imagename or image.endswith('.alpha'):
            stats_dict['mask_non0'] = numpy.count_nonzero(chunk*masks)

        if 'weight' in imagename:
            if 'cube' in imagename:
                wt_02_list = []
                wt_05_list = []
                i = 0
        
                for chan in chunk:
                    wt_02_list.append(numpy.count_nonzero(chan*masks[0][i]))
                    wt_05_list.append(numpy.count_nonzero(chan*masks[1][i]))
                    i += 1
                
                stats_dict['npts_0.2'] = wt_02_list
                stats_dict['npts_0.5'] = wt_05_list
      
            else:
                stats_dict['npts_0.2'] = numpy.count_nonzero(chunk*masks[0])
                stats_dict['npts_0.5'] = numpy.count_nonzero(chunk*masks[1])

        self._myia.close()

        return stats_dict

    def image_list(self, image, mode):
        """ function used to return expected imaging output files """

        standard = [image+'.psf', image+'.residual', image+'.image', 
                    image+'.image.pbcor', image+'.mask', image+'.pb', image+'.model', 
                    image+'.sumwt']
        
        mosaic = [image+'.weight']
        mtmfs = [image+'.alpha', image+'.alpha.error', image+'.alpha.pbcor', \
          image+'.psf.tt0', image+'.psf.tt1', image+'.psf.tt2', \
          image+'.residual.tt0', image+'.residual.tt1', image+'.image.tt0',\
          image+'.image.tt1', image+'.image.tt0.pbcor', image+'.image.tt1.pbcor', \
          image+'.mask', image+'.pb.tt0', image+'.model.tt0', image+'.model.tt1', \
          image+'.sumwt.tt0', image+'.sumwt.tt1', image+'.sumwt.tt2']
        mos_mtmfs = [image+'.weight.tt0', image+'.weight.tt1', image+'.weight.tt2']

        if mode == 'standard':
            img_list = standard
        if mode == 'mosaic':
            img_list = standard+mosaic
        if mode == 'mtmfs':
            img_list = mtmfs
        if mode == 'mos_mtmfs':
            img_list = mtmfs+mos_mtmfs

        return img_list

    def mom8_creator(self, image, range_list):
        """ function that takes and image and turns it into a .png for
            weblog
        """

        immoments(imagename = image, moments = 8, outfile = image+'.moment8')
        imview(raster={'file': image+'.moment8', 'range': range_list}, out = {'file': image+'.moment8.png'})
        subprocess.call('mogrify -trim '+image+'.moment8.png', shell=True)

    def cube_profile_fit(self, image, max_loc, nchan):
        """ function that will retrieve a profile for cubes at the max position
            and create a png showing the profile plot; must be called with
            image already opened
        """

        pyplot.clf()
        box = str(max_loc[0])+','+str(max_loc[1])+','+str(max_loc[0])+','+str(max_loc[1])
        profile = self._myia.fitprofile(box=box)['gs']['amp'][0][0][0][0][0]
        X = self._myia.getchunk(blc=max_loc, trc=max_loc, axes=[0,1])[0][0][0]
        
        pyplot.title('Frequency Profile at Max Value Position')
        pyplot.xlabel('Channel Number')
        pyplot.xlim(0,(nchan+1))
        pyplot.ylabel('Amplitude (Jy/Beam)')
        pyplot.plot(X)
        pyplot.savefig(image+'.profile.png')
        pyplot.clf()

        return profile

    def filter_report(self, report, showonlyfail=True):
        """ 
          function to filter the test report, the input report is expected to be a string with the newline code 
        """

        ret = ''
        if showonlyfail:
            filter='Fail'
        else:
            filter='Pass' 
        
        if report!='':
            testItems = report.split('\n')
            retitems=[]
            for testitem in testItems:
                if '[ check_ims ]' in testitem or '[ check_pixmask ]' in testitem or '[ check_val ]' in testitem:
                    if '( '+filter in testitem:
                        retitems.append(testitem)
      
            nfail = len(retitems)
            msg = str(nfail)+' individual test failure(s) '
            ret = '\n' + '\n'.join(retitems)
            ret += '\n' + msg
      
            return ret

    def save_dict_to_file(self, topkey, indict, outfilename, appendversion=True, outformat='JSON'):
        """ function that will save input Python dictionaries to a JSON file (default)
            or pickle file. topkey is will be added as a top key for output (nested) dictionary 
            and indict is stored under the key.  
            Create a separate file with outfilename if appendversion=True casa version (based on
            casatasks version) will be appended to the output file name.
        """
        try: 
            import casatasks as __casatasks
            casaversion = __casatasks.version_string()
            del __casatasks
        except:
            casaversion = ''
            
        if casaversion !='':
            casaversion = '_' + casaversion
        if type(indict) != dict:
            print("indict is not a dict. Saved file may not be in correct format") 
        nestedDict={}
        nestedDict[topkey]=indict 
        print("Saving %s dictionaries", len(indict))
        if outformat == 'pickle':
            # writing to pickle: note if writing this way (without protocol=2) 
            # in casa6 and read in casa5 it will fail 
            with open(outfilename+casaversion+'.pickle', 'wb') as outf:
                pickle.dump(nestedDict, outf)
        elif outformat== 'JSON':
            with open(outfilename+casaversion+'.json', 'w') as outf:
                json.dump(nestedDict, outf)
        else:
            print("no saving with format:", outformat)
    
    def modify_dict(self, output=None, testname=None, parallel=None):
        ''' Modified test_dict costructed by casatestutils add_to_dict to include only 
            the task commands executed and also add self.parallel value to the dictionary.
            The cube imaging cases usually have if-else conditional based on parallel mode is on or not
            to trigger different set of tclean commands.
            Assumption: self.parallel is used to trigger different tclean commands at iter1 step.
            For self.parallel=True, iter1 has two tclean commands (2nd and 3rd tclean commands within
            each relevante test(cube) and so in test_dict['taskcall'], 1st(iter0) and 2nd and 3rd commands
            are the ones acutually executed and should remove 4th (self.parallel=False) case.
        '''

        if testname in output:
            if 'taskcall' in output[testname] and len(output[testname]['taskcall'])==4: 
                if parallel:
                  # 0,1,2th in the list are used pop last one
                    output[testname]['taskcall'].pop()
            else:
                output[testname]['taskcall'].pop(1)
                output[testname]['taskcall'].pop(1)
        
            output[testname]['self.parallel']=parallel


    def save_dict_to_disk(self, indict, outfilename):
        """ 
          function that will save input Python dictionary to file (json) 
        """
    
        with open(outfilename+'.json', 'w') as outf:
            json.dump(indict, outf)

    def remove_prefix(self,string, prefix):
        ''' Remove a specified prefix string from string '''
        return string[string.startswith(prefix) and len(prefix):]

    def analysis_report(self):
        return print(self.report)

test_dict = {}
class Test_standard(test_tclean_base):

    @stats_dict(test_dict)
    def test_standard_cube(self):
        ''' Standard (single field) cube imaging - central field of SMIDGE_NWCloud (field 3), spw 22 '''

        self._testMethodName = 'test_standard_cube'
        test_name = self._testMethodName
        file_name = self.remove_prefix(test_name, 'test_')+'.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path + 'E2E6.1.00034.S_tclean.ms')

        

        print("\nSTARTING: iter0 routine")
        
        tclean(
            vis=self.vis, 
            imagename=file_name+'0', 
            field=self.field, 
            spw=self.spw, 
            imsize=self.imsize, 
            antenna=self.antenna, 
            scan=self.scan, 
            intent=self.intent, 
            datacolumn=self.datacolumn, 
            cell=self.cell, 
            phasecenter=self.phasecenter, 
            stokes=self.stokes, 
            specmode=self.specmode, 
            nchan=self.nchan, 
            start=self.start, 
            width=self.width, 
            outframe=self.outframe, 
            pblimit=self.pblimit, 
            perchanweightdensity=self.perchanweightdensity, 
            gridder=self.gridder,  
            mosweight=self.mosweight, 
            deconvolver=self.deconvolver, 
            usepointing=self.usepointing, 
            restoration=self.restoration, 
            pbcor=self.pbcor, 
            weighting=self.weighting, 
            restoringbeam=self.restoringbeam, 
            robust=self.robust, 
            npixels=self.npixels, 
            niter=0, 
            threshold=self.threshold, 
            nsigma=self.nsigma, 
            interactive=self.interactive, 
            usemask=self.usemask, 
            sidelobethreshold=self.sidelobethreshold, 
            noisethreshold=self.noisethreshold, 
            lownoisethreshold=self.lownoisethreshold, 
            negativethreshold=self.negativethreshold, 
            minbeamfrac=self.minbeamfrac, 
            growiterations=self.growiterations, 
            dogrowprune=self.dogrowprune, 
            minpercentchange=self.minpercentchange, 
            fastnoise=self.fastnoise, 
            savemodel=self.savemodel, 
            parallel=self.parallel,
            verbose=self.verbose)

        # iter0 routine
        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_products(file_name+'0', file_name+'1')

        print("STARTING: iter1 routine")
        self.read_configuration(self.config_file)

        # iter1 (restart)
        if self.parallel:
            tclean(
                vis=self.vis, 
                imagename=file_name+'1', 
                field=self.field,
                spw=self.spw, 
                imsize=self.imsize, 
                antenna=self.antenna, 
                scan=self.scan, 
                intent=self.intent, 
                datacolumn=self.datacolumn, 
                cell=self.cell, 
                phasecenter=self.phasecenter, 
                stokes=self.stokes, 
                specmode=self.specmode, 
                nchan=self.nchan, 
                start=self.start, 
                width=self.width,
                outframe=self.outframe, 
                perchanweightdensity=self.perchanweightdensity, 
                usepointing=self.usepointing,
                pblimit=self.pblimit, 
                nsigma=self.nsigma, 
                gridder=self.gridder, 
                mosweight=self.mosweight, 
                deconvolver=self.deconvolver, 
                restoration=self.restoration, 
                pbcor=self.pbcor, 
                weighting=self.weighting, 
                robust=self.robust, 
                npixels=self.npixels, 
                niter=self.niter, 
                threshold=self.threshold, 
                interactive=self.interactive, 
                usemask=self.usemask, 
                sidelobethreshold=self.sidelobethreshold, 
                noisethreshold=self.noisethreshold, 
                lownoisethreshold=self.lownoisethreshold, 
                negativethreshold=self.negativethreshold, 
                minbeamfrac=self.minbeamfrac, 
                growiterations=self.growiterations,
                dogrowprune=self.dogrowprune, 
                minpercentchange=self.minpercentchange, 
                fastnoise=self.fastnoise, 
                restart=self.restart, 
                calcres=self.calcres, 
                calcpsf=self.calcpsf, 
                savemodel=self.savemodel, 
                restoringbeam=self.restoringbeam, 
                parallel=self.parallel, 
                verbose=self.verbose)

            # retrieve per-channel beam statistics (only in parallel)
            bmin_dict, bmaj_dict, pa_dict = \
                self.cube_beam_stats(image=img+'.image')

        else:
            tclean(
                vis=self.vis, 
                imagename=file_name+'1', 
                field=self.field,
                spw=self.spw, 
                imsize=self.imsize, 
                antenna=self.antenna, 
                scan=self.scan, 
                intent=self.intent, 
                datacolumn=self.datacolumn, 
                cell=self.cell, 
                phasecenter=self.phasecenter, 
                stokes=self.stokes, 
                specmode=self.specmode, 
                nchan=self.nchan, 
                start=self.start, 
                width=self.width,
                outframe=self.outframe, 
                perchanweightdensity=self.perchanweightdensity, 
                usepointing=self.usepointing,
                pblimit=self.pblimit, 
                nsigma=self.nsigma, 
                gridder=self.gridder, 
                mosweight=self.mosweight, 
                deconvolver=self.deconvolver, 
                restoration=self.restoration, 
                pbcor=self.pbcor, 
                weighting=self.weighting, 
                robust=self.robust, 
                npixels=self.npixels, 
                niter=self.niter, 
                threshold=self.threshold, 
                interactive=self.interactive, 
                usemask=self.usemask, 
                sidelobethreshold=self.sidelobethreshold, 
                noisethreshold=self.noisethreshold, 
                lownoisethreshold=self.lownoisethreshold, 
                negativethreshold=self.negativethreshold, 
                minbeamfrac=self.minbeamfrac, 
                growiterations=self.growiterations,
                dogrowprune=self.dogrowprune, 
                minpercentchange=self.minpercentchange, 
                fastnoise=self.fastnoise, 
                restart=self.restart, 
                calcres=self.calcres, 
                calcpsf=self.calcpsf, 
                savemodel=self.savemodel, 
                restoringbeam=self.restoringbeam, 
                parallel=self.parallel, 
                verbose=self.verbose)

        # Check that images, history, and keywords exist th.checkall > th.check_imexist()
        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        # Calculate stats dictionary for cleaned image.
        im_stats_dict = self.image_stats(image=img+'.image', 
                                         fit_region = 'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        exp_im_stats = {
            'com_bmaj': [False, 8.509892605313942],
            'com_bmin': [False, 5.950050676606115],
            'com_pa': [False, 72.54607919421503],
            'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.94608676433563232],
            'max_val_pos': [True, [38, 36, 0, 254]],
            'min_val': [False, -0.70467984676361084],
            'min_val_pos': [True, [18, 57, 0, 374]],
            'im_rms': [False, 0.143986095161],
            'rms_per_chan': [False, 
                             [0.12634143789549535, 0.14184103443455054, 0.146980848491423, 0.11560714721775701, 
                              0.16661204089962214, 0.14476185721954618, 0.1316304465304373, 0.11433992978005528, 
                              0.1442336908847237, 0.16543686095641913, 0.11873031602382604, 0.15284339493903618, 
                              0.17195923304927352, 0.15046376352374186, 0.14070242684409134, 0.13744696370230045, 
                              0.11158401597402363, 0.12687194824894843, 0.14295610599459097, 0.16862585818094997,
                              0.13008023644212954, 0.14186757490541813, 0.1169917541756216, 0.13693402712045005, 
                              0.16194773183534902, 0.13634870573065122, 0.13860445389090703, 0.1309492701136035, 
                              0.14819650092974662, 0.15030436484700252, 0.14931368127692368, 0.11984058396768074, 
                              0.13119726238629406, 0.1278997483823513, 0.15680618961802364, 0.14782343113879803, 
                              0.1452811146145065, 0.14962350388870774, 0.12727392822661138, 0.12403611951801675, 
                              0.13971310504808565, 0.14443747442976043, 0.13947457857066817, 0.14603448593891352,
                              0.1487357653330162, 0.13728792717834695, 0.12754218448487697, 0.13472363429296094, 
                              0.17318897000268654, 0.15007875971445414, 0.1210452212469422, 0.15977553256440455, 
                              0.13077138200444186, 0.12679151267047647, 0.12091204082027505, 0.1338966333695089, 
                              0.13556991771575277, 0.15456345918134376, 0.12044465611280492, 0.14836982561861023,
                              0.1349896116866282, 0.15311214064438922, 0.11251497655504887, 0.134867796496014,
                              0.13574313457554038, 0.14582224580240324, 0.12753531221719416, 0.15335445312643003,
                              0.13482732612307324, 0.1622050903445585, 0.13260306174268546, 0.1345326608100535,
                              0.16404765102131583, 0.13449430188802702, 0.14543809289295098, 0.1606584196112734,
                              0.12484651484486906, 0.16251383851634701, 0.13756025624117688, 0.13165353467440083,
                              0.1308248320448295, 0.14752778635690292, 0.1274645256107852, 0.16421712463271607,
                              0.15255317243782812, 0.1497707840063393, 0.11911825364867326, 0.14541033702618353,
                              0.1659723426787793, 0.1554971226410762, 0.14703675741501698, 0.12325846980328654,
                              0.15070706791866434, 0.1243073669840061, 0.13646642844468243, 0.1301143392639293,
                              0.12734602178400867, 0.1553600823593344, 0.15035594210430997, 0.11530605847413075,
                              0.1611567346343003, 0.12221832384850957, 0.14207389319672978, 0.14522516033398006,
                              0.1345322788758837, 0.1486176245373929, 0.15765848896613346, 0.1308440759384876,
                              0.1466820831226493, 0.13598865468593319, 0.15187538855740168, 0.1478468013010444,
                              0.1383732515889412, 0.1276861625889527, 0.11697230161534232, 0.1373960738818452, 
                              0.11303259344169146, 0.1361001584583741, 0.12857356426667815, 0.1437570752313611,
                              0.13169397143643052, 0.15326431411050365, 0.12383180315967929, 0.1526310794015497, 
                              0.14746177769245866, 0.15194893390457265, 0.1421630320154613, 0.15662308690272084, 
                              0.12239198421329735, 0.12071542153915982, 0.14268554321174182, 0.13489697242976567,
                              0.15127855443293006, 0.1542443819690316, 0.15752918577920158, 0.11478434733366248,
                              0.17298964180575135, 0.13177526480150695, 0.12236732291938952, 0.15625856947990782,
                              0.13687165189461548, 0.1536631153928859, 0.14669563803395924, 0.1277170908624889, 
                              0.14966567842171496, 0.12823515897560267, 0.13577828413547297, 0.16140169123660877,
                              0.13133284404676335, 0.14223570583416104, 0.1603292311222728, 0.10759630495294702, 
                              0.15787039978749143, 0.1327200609847152, 0.14655899389809018, 0.14008820956915727,
                              0.1442107348583108, 0.1317943450568934, 0.12972989243424995, 0.1625036947147829, 
                              0.12241712383574781, 0.14998173521745944, 0.13352731228428555, 0.1741676258276787, 
                              0.15545996482656257, 0.13121844421079562, 0.1389256768353536, 0.1475992903718036, 
                              0.14205849908080379, 0.14975427804440275, 0.1532491403618113, 0.12531915969323904, 
                              0.14153689035122896, 0.16741877503811964, 0.1355536447212321, 0.12548585056941425, 
                              0.16334800417248366, 0.14220841606737944, 0.1376802362928535, 0.1394159389365598,
                              0.1533008119644231, 0.12568227593323275, 0.14138024496799945, 0.14688836279261966, 
                              0.12037367892758656, 0.12335138886587714, 0.16740640885840646, 0.11756235238942149, 
                              0.13221931449560975, 0.14605469946826174, 0.12287649136200192, 0.13900407591276098, 
                              0.1477935699475207, 0.14723640198504923, 0.12637771862286276, 0.14264989851200444, 
                              0.14188497863070984, 0.1517498748029243, 0.1745550071541481, 0.14693061119966988,
                              0.12180541963696558, 0.17178472812899895, 0.134842796032342, 0.1587769050427257,
                              0.16022475326023228, 0.12598385136025822, 0.12173065475536829, 0.1358700032273519, 
                              0.12249230371601251, 0.1320416693266833, 0.1380195667444624, 0.17036819494074398,
                              0.14449179298441997, 0.1363579047545357, 0.15814587607932587, 0.1387404461979901, 
                              0.13421674959986293, 0.1221729254232071, 0.15007074873391474, 0.1519841002224019, 
                              0.17405910974305452, 0.10810253208919626, 0.14404509620995673, 0.12925102011532486, 
                              0.13284702789448985, 0.16316517507291742, 0.18004246985230368, 0.12352109323053732, 
                              0.13438971701846103, 0.14110722423724795, 0.15240505247738928, 0.16299890474660164, 
                              0.13862726296963718, 0.13653417057201506, 0.15748574227626927, 0.13330507817933285, 
                              0.11630210517195279, 0.14310200319865532, 0.16947649357417122, 0.19276632]], 
            'im_sum': [False, 168.242297979], 
            'regn_sum': [False, 66.089371562004089],
            'npts_real': [True, 3251200],
            'profile': [False, 0.90354546054915941],
            'fit': [False, [0.9222914423989481, 14.289411729097703, 6.881365508161659]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [38.263402177385942, 37.306443753086633]]}

        report1 = th.checkall(
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 70, 0, 0]), 
                       (img+'.image', False, [40, 71, 0, 0]), 
                       (img+'.image', True, [10, 40, 0, 0]), 
                       (img+'.image', False, [9, 40, 0, 0])])

        # .image report
        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(image=img+'.mask')

        exp_mask_stats = {'npts': [True, 3251200],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'mask_pix': [False, 437],
            'mask_regns': [True, 1],
            'npts_real': [True, 3251200]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(image=img+'.pb', 
                                         fit_region = 'ellipse[[11.47659846deg, -73.25817055deg], [23.1086arcsec, 23.0957arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.200896695256],
            'im_rms': [False, 0.578238326026],
            'npts_0.2': [False, [2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997]],
            'npts_0.5': [False, [1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449]],
            'npts_real': [True, 3251200],
            'fit': [False, [1.0308127949041446, 46.61751391582679, 46.61253844001269]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [40.00032808200995, 40.00099739969875]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(image=img+'.psf', fit_region = 'ellipse[[11.47648725deg, -73.25812003deg], [8.0291arcsec, 6.8080arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 3251200.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.218764916062],
            'min_val_pos': [True, [1, 16, 0, 503]],
            'im_rms': [False, 0.136036099793],
            'im_sum': [False, 7472.57665916],
            'npts_real': [True, 3251200],
            'fit_0': [False, [1.0959520385253885, 7.675969776744627, 5.143545685744538]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 220.2529185335],
            'fit_pix_0': [False, [39.998799449215596, 39.99523672953224]],
            'fit_1': [False, [1.0959863390592945, 7.672871552789668, 5.141790170376213]],
            'fit_loc_chan_1': [True, 254],
            'fit_loc_freq_1': [1e-10, 220.31469458079383],
            'fit_pix_1': [False, [39.99880225653659, 39.99524870969922]],
            'fit_2': [False, [1.0960422882714267, 7.669928861314425, 5.140004109591353]],
            'fit_loc_chan_2': [True, 507],
            'fit_loc_freq_2': [1e-10, 220.37647062808767],
            'fit_pix_2': [False, [39.9988051116427, 39.995258207738075]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(image=img+'.residual', 
                                            fit_region = 'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        exp_resid_stats = {'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.785612404346],
            'max_val_pos': [True, [42, 43, 0, 256]],
            'min_val': [False, -0.704679846764],
            'min_val_pos': [True, [18, 57, 0, 374]],
            'im_rms': [False, 0.143918523224],
            'im_sum': [False, 124.317946204],
            'regn_sum': [False, 19.0307790947],
            'npts_real': [True, 3251200]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(image=img+'.model', fit_region = 'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 3251200.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.286023736],
            'max_val_pos': [True, [38, 36, 0, 254]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, 
                            [0, 0, 0, 0]],
            'im_rms': [False, 0.000249846621096],
            'im_sum': [False, 0.92636379227],
            'regn_sum': [False, 0.92636379227],
            'mask_non0': [True, 0],
            'npts_real': [True, 3251200]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(image=img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 508],
            'npts_unmasked': [True, 508.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 94.4766769409],
            'max_val_pos': [True, [0, 0, 0, 17]],
            'min_val': [False, 94.4766464233],
            'min_val_pos': [True, [0, 0, 0, 449]],
            'im_rms': [False, 94.476661707],
            'npts_real': [True, 508]}
            
        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + report6 + report7 + report8
        print('--> Reports finished.')

        if self.parallel:
            exp_bmin_dict = {'*0': 5.944103717803955, '*1': 5.94410514831543, '*10': 5.944015026092529, '*100': 5.94343900680542, '*101': 5.943426609039307, '*102': 5.943421840667725, '*103': 5.943414211273193, '*104': 5.943387985229492, '*105': 5.943387985229492, '*106': 5.943348407745361, '*107': 5.943348407745361, '*108': 5.94331693649292, '*109': 5.94331693649292, '*11': 5.943998336791992, '*110': 5.94331693649292, '*111': 5.94331693649292, '*112': 5.94331693649292, '*113': 5.94331693649292, '*114': 5.94331693649292, '*115': 5.94331693649292, '*116': 5.94331693649292, '*117': 5.943300247192383, '*118': 5.9432783126831055, '*119': 5.9432783126831055, '*12': 5.94395637512207, '*120': 5.943252086639404, '*121': 5.943252086639404, '*122': 5.943222999572754, '*123': 5.943180561065674, '*124': 5.943145275115967, '*125': 5.943161487579346, '*126': 5.943154811859131, '*127': 5.943154811859131, '*128': 5.943154811859131, '*129': 5.9431471824646, '*13': 5.94395637512207, '*130': 5.9431471824646, '*131': 5.9431471824646, '*132': 5.9431471824646, '*133': 5.943125247955322, '*134': 5.943125247955322, '*135': 5.94309663772583, '*136': 5.94309663772583, '*137': 5.943094253540039, '*138': 5.94307279586792, '*139': 5.94307279586792, '*14': 5.94395637512207, '*140': 5.943065166473389, '*141': 5.943065166473389, '*142': 5.943065166473389, '*143': 5.943065166473389, '*144': 5.943065166473389, '*145': 5.943065166473389, '*146': 5.943028926849365, '*147': 5.943028926849365, '*148': 5.943004131317139, '*149': 5.942957878112793, '*15': 5.94395637512207, '*150': 5.942957878112793, '*151': 5.942957878112793, '*152': 5.94295072555542, '*153': 5.94295072555542, '*154': 5.942929267883301, '*155': 5.942929267883301, '*156': 5.94288444519043, '*157': 5.94288444519043, '*158': 5.94288444519043, '*159': 5.94288444519043, '*16': 5.943938255310059, '*160': 5.94288444519043, '*161': 5.942834854125977, '*162': 5.942834854125977, '*163': 5.942787170410156, '*164': 5.942765712738037, '*165': 5.942748546600342, '*166': 5.942732810974121, '*167': 5.942732810974121, '*168': 5.94273042678833, '*169': 5.94273042678833, '*17': 5.94391393661499, '*170': 5.94273042678833, '*171': 5.94273042678833, '*172': 5.942715644836426, '*173': 5.942707538604736, '*174': 5.942705154418945, '*175': 5.942681789398193, '*176': 5.942651271820068, '*177': 5.9426374435424805, '*178': 5.942581653594971, '*179': 5.942581653594971, '*18': 5.943907260894775, '*180': 5.942537307739258, '*181': 5.942537307739258, '*182': 5.942537307739258, '*183': 5.942534923553467, '*184': 5.942537784576416, '*185': 5.942537784576416, '*186': 5.942508697509766, '*187': 5.942508697509766, '*188': 5.942508697509766, '*189': 5.942508697509766, '*19': 5.943899154663086, '*190': 5.942508697509766, '*191': 5.942508697509766, '*192': 5.942503929138184, '*193': 5.942503929138184, '*194': 5.942503929138184, '*195': 5.942503929138184, '*196': 5.942503929138184, '*197': 5.942503929138184, '*198': 5.942503452301025, '*199': 5.9424729347229, '*2': 5.9440813064575195, '*20': 5.943899154663086, '*200': 5.9424729347229, '*201': 5.9424729347229, '*202': 5.94244384765625, '*203': 5.942445278167725, '*204': 5.942437171936035, '*205': 5.942437171936035, '*206': 5.942437171936035, '*207': 5.942437171936035, '*208': 5.942437171936035, '*209': 5.942437171936035, '*21': 5.943899154663086, '*210': 5.9424333572387695, '*211': 5.9424333572387695, '*212': 5.9424333572387695, '*213': 5.9424333572387695, '*214': 5.9424333572387695, '*215': 5.9424333572387695, '*216': 5.942417144775391, '*217': 5.942417144775391, '*218': 5.942412376403809, '*219': 5.942412376403809, '*22': 5.943899154663086, '*220': 5.942407131195068, '*221': 5.942407131195068, '*222': 5.942407131195068, '*223': 5.942407131195068, '*224': 5.942375659942627, '*225': 5.94236946105957, '*226': 5.94236946105957, '*227': 5.942365646362305, '*228': 5.942312717437744, '*229': 5.942312717437744, '*23': 5.943899154663086, '*230': 5.942312717437744, '*231': 5.942311763763428, '*232': 5.942311763763428, '*233': 5.942311763763428, '*234': 5.942311763763428, '*235': 5.942282676696777, '*236': 5.942282676696777, '*237': 5.942237854003906, '*238': 5.942237854003906, '*239': 5.942210674285889, '*24': 5.9438910484313965, '*240': 5.942195415496826, '*241': 5.942195415496826, '*242': 5.942195415496826, '*243': 5.942178249359131, '*244': 5.942124366760254, '*245': 5.942124366760254, '*246': 5.942124366760254, '*247': 5.942124366760254, '*248': 5.942124366760254, '*249': 5.942124366760254, '*25': 5.943882465362549, '*250': 5.942102909088135, '*251': 5.942102909088135, '*252': 5.942080497741699, '*253': 5.942091464996338, '*254': 5.942091464996338, '*255': 5.94207763671875, '*256': 5.942057132720947, '*257': 5.942053318023682, '*258': 5.942053318023682, '*259': 5.942053318023682, '*26': 5.943882465362549, '*260': 5.942053318023682, '*261': 5.94204568862915, '*262': 5.94204568862915, '*263': 5.94204568862915, '*264': 5.94204568862915, '*265': 5.94204568862915, '*266': 5.941995143890381, '*267': 5.941966533660889, '*268': 5.941969394683838, '*269': 5.941969394683838, '*27': 5.943882465362549, '*270': 5.94195032119751, '*271': 5.941928386688232, '*272': 5.941928386688232, '*273': 5.941908359527588, '*274': 5.941908359527588, '*275': 5.941908359527588, '*276': 5.941908359527588, '*277': 5.941891670227051, '*278': 5.941891670227051, '*279': 5.941891193389893, '*28': 5.943875789642334, '*280': 5.941891193389893, '*281': 5.941891193389893, '*282': 5.941891193389893, '*283': 5.941833972930908, '*284': 5.941833972930908, '*285': 5.941833972930908, '*286': 5.941833972930908, '*287': 5.941833972930908, '*288': 5.941832542419434, '*289': 5.941824913024902, '*29': 5.943875789642334, '*290': 5.941824913024902, '*291': 5.941812038421631, '*292': 5.941808223724365, '*293': 5.941808223724365, '*294': 5.941808223724365, '*295': 5.941808223724365, '*296': 5.941801071166992, '*297': 5.941802501678467, '*298': 5.941776752471924, '*299': 5.941746711730957, '*3': 5.9440813064575195, '*30': 5.943875789642334, '*300': 5.9417266845703125, '*301': 5.9417266845703125, '*302': 5.941720008850098, '*303': 5.941673278808594, '*304': 5.941673278808594, '*305': 5.9416584968566895, '*306': 5.9416584968566895, '*307': 5.9416584968566895, '*308': 5.9416584968566895, '*309': 5.9416584968566895, '*31': 5.943875789642334, '*310': 5.941618919372559, '*311': 5.941612243652344, '*312': 5.941611289978027, '*313': 5.941611289978027, '*314': 5.941610336303711, '*315': 5.941610336303711, '*316': 5.941610336303711, '*317': 5.941610336303711, '*318': 5.941595077514648, '*319': 5.941595077514648, '*32': 5.943851947784424, '*320': 5.941573143005371, '*321': 5.941573143005371, '*322': 5.941573143005371, '*323': 5.94158935546875, '*324': 5.94158935546875, '*325': 5.941571235656738, '*326': 5.941502571105957, '*327': 5.941502571105957, '*328': 5.941502571105957, '*329': 5.941502571105957, '*33': 5.943851947784424, '*330': 5.941502571105957, '*331': 5.941502571105957, '*332': 5.941502571105957, '*333': 5.9414777755737305, '*334': 5.9414777755737305, '*335': 5.941437721252441, '*336': 5.941437721252441, '*337': 5.941437721252441, '*338': 5.94136381149292, '*339': 5.94136381149292, '*34': 5.943851947784424, '*340': 5.94134521484375, '*341': 5.941339015960693, '*342': 5.941339015960693, '*343': 5.941314220428467, '*344': 5.941311836242676, '*345': 5.941303253173828, '*346': 5.941303253173828, '*347': 5.941303253173828, '*348': 5.941303253173828, '*349': 5.941303253173828, '*35': 5.943843364715576, '*350': 5.941303253173828, '*351': 5.941265106201172, '*352': 5.941249847412109, '*353': 5.941249847412109, '*354': 5.941246509552002, '*355': 5.941246509552002, '*356': 5.941233158111572, '*357': 5.941233158111572, '*358': 5.941233158111572, '*359': 5.941233158111572, '*36': 5.943843364715576, '*360': 5.941233158111572, '*361': 5.9412522315979, '*362': 5.9412522315979, '*363': 5.9412522315979, '*364': 5.941198825836182, '*365': 5.941179275512695, '*366': 5.941154479980469, '*367': 5.941154479980469, '*368': 5.941154479980469, '*369': 5.9411540031433105, '*37': 5.943843364715576, '*370': 5.941145896911621, '*371': 5.941145896911621, '*372': 5.941145896911621, '*373': 5.941145896911621, '*374': 5.941145896911621, '*375': 5.941142559051514, '*376': 5.941142559051514, '*377': 5.941142559051514, '*378': 5.941142559051514, '*379': 5.941142559051514, '*38': 5.943840980529785, '*380': 5.941142559051514, '*381': 5.941142559051514, '*382': 5.941142559051514, '*383': 5.941102504730225, '*384': 5.941102504730225, '*385': 5.941101551055908, '*386': 5.941097259521484, '*387': 5.94110107421875, '*388': 5.94110107421875, '*389': 5.941096305847168, '*39': 5.943840980529785, '*390': 5.941096305847168, '*391': 5.941074848175049, '*392': 5.941074848175049, '*393': 5.941074848175049, '*394': 5.94106912612915, '*395': 5.941043376922607, '*396': 5.941043376922607, '*397': 5.941043376922607, '*398': 5.941043376922607, '*399': 5.941043376922607, '*4': 5.9440813064575195, '*40': 5.943840980529785, '*400': 5.941017150878906, '*401': 5.940994739532471, '*402': 5.940994739532471, '*403': 5.940974712371826, '*404': 5.940969467163086, '*405': 5.940962791442871, '*406': 5.940962791442871, '*407': 5.940912246704102, '*408': 5.940912246704102, '*409': 5.940912246704102, '*41': 5.943840980529785, '*410': 5.940889835357666, '*411': 5.940889835357666, '*412': 5.940889835357666, '*413': 5.940889835357666, '*414': 5.940862655639648, '*415': 5.940829277038574, '*416': 5.940829277038574, '*417': 5.940829277038574, '*418': 5.940817356109619, '*419': 5.9407854080200195, '*42': 5.943844795227051, '*420': 5.940733909606934, '*421': 5.940733909606934, '*422': 5.940733909606934, '*423': 5.940733909606934, '*424': 5.940709114074707, '*425': 5.940709114074707, '*426': 5.940701961517334, '*427': 5.940701961517334, '*428': 5.940701961517334, '*429': 5.940701961517334, '*43': 5.943844795227051, '*430': 5.940701961517334, '*431': 5.940697193145752, '*432': 5.940697193145752, '*433': 5.940697193145752, '*434': 5.940697193145752, '*435': 5.940667629241943, '*436': 5.940667629241943, '*437': 5.940667629241943, '*438': 5.940667629241943, '*439': 5.940667629241943, '*44': 5.943831920623779, '*440': 5.940667629241943, '*441': 5.9406633377075195, '*442': 5.9406633377075195, '*443': 5.9406633377075195, '*444': 5.9406633377075195, '*445': 5.9406633377075195, '*446': 5.9406633377075195, '*447': 5.9406633377075195, '*448': 5.9406633377075195, '*449': 5.940667629241943, '*45': 5.94381856918335, '*450': 5.94066858291626, '*451': 5.94066858291626, '*452': 5.940656661987305, '*453': 5.940609931945801, '*454': 5.940609931945801, '*455': 5.940609931945801, '*456': 5.940608024597168, '*457': 5.940608024597168, '*458': 5.940586566925049, '*459': 5.940587520599365, '*46': 5.943813800811768, '*460': 5.940587520599365, '*461': 5.940557479858398, '*462': 5.940558910369873, '*463': 5.940558910369873, '*464': 5.940558910369873, '*465': 5.940556049346924, '*466': 5.940556049346924, '*467': 5.940556049346924, '*468': 5.940556049346924, '*469': 5.940556049346924, '*47': 5.943813800811768, '*470': 5.940556049346924, '*471': 5.940539836883545, '*472': 5.940539836883545, '*473': 5.940534591674805, '*474': 5.940478324890137, '*475': 5.940452575683594, '*476': 5.940426826477051, '*477': 5.940426826477051, '*478': 5.9403886795043945, '*479': 5.940366744995117, '*48': 5.943813800811768, '*480': 5.940366744995117, '*481': 5.940366744995117, '*482': 5.940362453460693, '*483': 5.940362453460693, '*484': 5.940362453460693, '*485': 5.940317630767822, '*486': 5.940317630767822, '*487': 5.940317630767822, '*488': 5.940317630767822, '*489': 5.940317630767822, '*49': 5.943813800811768, '*490': 5.940311908721924, '*491': 5.940311908721924, '*492': 5.940311908721924, '*493': 5.940311908721924, '*494': 5.940311908721924, '*495': 5.940311908721924, '*496': 5.940311908721924, '*497': 5.940313816070557, '*498': 5.940313816070557, '*499': 5.940313816070557, '*5': 5.944080829620361, '*50': 5.943783283233643, '*500': 5.940237045288086, '*501': 5.940237045288086, '*502': 5.940237045288086, '*503': 5.9401655197143555, '*504': 5.9401655197143555, '*505': 5.940146446228027, '*506': 5.940146446228027, '*507': 5.940146446228027, '*51': 5.943783283233643, '*52': 5.943783283233643, '*53': 5.943783283233643, '*54': 5.943783283233643, '*55': 5.943783283233643, '*56': 5.943750381469727, '*57': 5.943746089935303, '*58': 5.943742275238037, '*59': 5.943742275238037, '*6': 5.944068908691406, '*60': 5.943739414215088, '*61': 5.943734169006348, '*62': 5.943727970123291, '*63': 5.943719863891602, '*64': 5.943719863891602, '*65': 5.943719863891602, '*66': 5.943719863891602, '*67': 5.943719863891602, '*68': 5.943719863891602, '*69': 5.943702697753906, '*7': 5.9440226554870605, '*70': 5.943700313568115, '*71': 5.943700313568115, '*72': 5.943695545196533, '*73': 5.943671703338623, '*74': 5.943672180175781, '*75': 5.943672180175781, '*76': 5.943665504455566, '*77': 5.943638801574707, '*78': 5.943638801574707, '*79': 5.943638801574707, '*8': 5.9440226554870605, '*80': 5.943630695343018, '*81': 5.943630695343018, '*82': 5.943585395812988, '*83': 5.943577289581299, '*84': 5.943577289581299, '*85': 5.943577289581299, '*86': 5.94357442855835, '*87': 5.94357442855835, '*88': 5.94357442855835, '*89': 5.943568229675293, '*9': 5.944015026092529, '*90': 5.943568229675293, '*91': 5.943568229675293, '*92': 5.943528175354004, '*93': 5.94349479675293, '*94': 5.94349479675293, '*95': 5.943480014801025, '*96': 5.943466663360596, '*97': 5.943450450897217, '*98': 5.943450450897217, '*99': 5.943450450897217}

            exp_bmaj_dict = {'*0': 8.501383781433105, '*1': 8.5013427734375, '*10': 8.501267433166504, '*100': 8.49976634979248, '*101': 8.499771118164062, '*102': 8.499725341796875, '*103': 8.499642372131348, '*104': 8.499631881713867, '*105': 8.499631881713867, '*106': 8.49958610534668, '*107': 8.49958610534668, '*108': 8.499568939208984, '*109': 8.499568939208984, '*11': 8.50125789642334, '*110': 8.499568939208984, '*111': 8.499568939208984, '*112': 8.499568939208984, '*113': 8.499568939208984, '*114': 8.499568939208984, '*115': 8.499568939208984, '*116': 8.499568939208984, '*117': 8.499425888061523, '*118': 8.499419212341309, '*119': 8.499419212341309, '*12': 8.501240730285645, '*120': 8.49941635131836, '*121': 8.49941635131836, '*122': 8.49940013885498, '*123': 8.499324798583984, '*124': 8.499231338500977, '*125': 8.499188423156738, '*126': 8.499133110046387, '*127': 8.499133110046387, '*128': 8.499133110046387, '*129': 8.499103546142578, '*13': 8.501240730285645, '*130': 8.499103546142578, '*131': 8.499103546142578, '*132': 8.499103546142578, '*133': 8.499098777770996, '*134': 8.499098777770996, '*135': 8.49912166595459, '*136': 8.49912166595459, '*137': 8.499106407165527, '*138': 8.498995780944824, '*139': 8.498995780944824, '*14': 8.501240730285645, '*140': 8.499002456665039, '*141': 8.499002456665039, '*142': 8.499002456665039, '*143': 8.499002456665039, '*144': 8.499002456665039, '*145': 8.499002456665039, '*146': 8.498955726623535, '*147': 8.498955726623535, '*148': 8.49891185760498, '*149': 8.498912811279297, '*15': 8.501240730285645, '*150': 8.498912811279297, '*151': 8.498912811279297, '*152': 8.498917579650879, '*153': 8.498917579650879, '*154': 8.498931884765625, '*155': 8.498931884765625, '*156': 8.49893569946289, '*157': 8.49893569946289, '*158': 8.49893569946289, '*159': 8.49893569946289, '*16': 8.501228332519531, '*160': 8.49893569946289, '*161': 8.498952865600586, '*162': 8.498952865600586, '*163': 8.498932838439941, '*164': 8.498906135559082, '*165': 8.498762130737305, '*166': 8.498750686645508, '*167': 8.498750686645508, '*168': 8.498736381530762, '*169': 8.498736381530762, '*17': 8.501107215881348, '*170': 8.498736381530762, '*171': 8.498736381530762, '*172': 8.498745918273926, '*173': 8.498677253723145, '*174': 8.498661994934082, '*175': 8.498669624328613, '*176': 8.498653411865234, '*177': 8.49864387512207, '*178': 8.498629570007324, '*179': 8.498629570007324, '*18': 8.50111198425293, '*180': 8.498604774475098, '*181': 8.498604774475098, '*182': 8.498604774475098, '*183': 8.498589515686035, '*184': 8.498554229736328, '*185': 8.498554229736328, '*186': 8.498551368713379, '*187': 8.498551368713379, '*188': 8.498551368713379, '*189': 8.498551368713379, '*19': 8.501042366027832, '*190': 8.498551368713379, '*191': 8.498551368713379, '*192': 8.498496055603027, '*193': 8.498496055603027, '*194': 8.498496055603027, '*195': 8.498496055603027, '*196': 8.498496055603027, '*197': 8.498496055603027, '*198': 8.498467445373535, '*199': 8.49845027923584, '*2': 8.501371383666992, '*20': 8.501042366027832, '*200': 8.49845027923584, '*201': 8.49845027923584, '*202': 8.498435020446777, '*203': 8.498387336730957, '*204': 8.498392105102539, '*205': 8.498392105102539, '*206': 8.498392105102539, '*207': 8.498392105102539, '*208': 8.498392105102539, '*209': 8.498392105102539, '*21': 8.501042366027832, '*210': 8.498343467712402, '*211': 8.498343467712402, '*212': 8.498343467712402, '*213': 8.498343467712402, '*214': 8.498343467712402, '*215': 8.498343467712402, '*216': 8.49832534790039, '*217': 8.49832534790039, '*218': 8.498236656188965, '*219': 8.498236656188965, '*22': 8.501042366027832, '*220': 8.498186111450195, '*221': 8.498186111450195, '*222': 8.498186111450195, '*223': 8.498186111450195, '*224': 8.49821949005127, '*225': 8.498163223266602, '*226': 8.498163223266602, '*227': 8.498127937316895, '*228': 8.498111724853516, '*229': 8.498111724853516, '*23': 8.501042366027832, '*230': 8.498111724853516, '*231': 8.498093605041504, '*232': 8.498093605041504, '*233': 8.49803638458252, '*234': 8.49803638458252, '*235': 8.498059272766113, '*236': 8.498059272766113, '*237': 8.498035430908203, '*238': 8.498035430908203, '*239': 8.498000144958496, '*24': 8.500977516174316, '*240': 8.497992515563965, '*241': 8.497992515563965, '*242': 8.497992515563965, '*243': 8.498001098632812, '*244': 8.497981071472168, '*245': 8.497981071472168, '*246': 8.497981071472168, '*247': 8.497981071472168, '*248': 8.497981071472168, '*249': 8.497981071472168, '*25': 8.500937461853027, '*250': 8.497970581054688, '*251': 8.497970581054688, '*252': 8.497958183288574, '*253': 8.497900009155273, '*254': 8.497900009155273, '*255': 8.497868537902832, '*256': 8.497882843017578, '*257': 8.497847557067871, '*258': 8.497847557067871, '*259': 8.497847557067871, '*26': 8.500937461853027, '*260': 8.497847557067871, '*261': 8.497831344604492, '*262': 8.497831344604492, '*263': 8.497831344604492, '*264': 8.497831344604492, '*265': 8.497831344604492, '*266': 8.49786376953125, '*267': 8.497773170471191, '*268': 8.497733116149902, '*269': 8.497733116149902, '*27': 8.500937461853027, '*270': 8.4977445602417, '*271': 8.497739791870117, '*272': 8.497739791870117, '*273': 8.497754096984863, '*274': 8.497754096984863, '*275': 8.497754096984863, '*276': 8.497712135314941, '*277': 8.497693061828613, '*278': 8.497693061828613, '*279': 8.497665405273438, '*28': 8.500882148742676, '*280': 8.497665405273438, '*281': 8.497665405273438, '*282': 8.497665405273438, '*283': 8.49767017364502, '*284': 8.49767017364502, '*285': 8.49767017364502, '*286': 8.49767017364502, '*287': 8.49767017364502, '*288': 8.497586250305176, '*289': 8.497503280639648, '*29': 8.500882148742676, '*290': 8.497503280639648, '*291': 8.4975004196167, '*292': 8.497406005859375, '*293': 8.497406005859375, '*294': 8.497406005859375, '*295': 8.497406005859375, '*296': 8.497410774230957, '*297': 8.497359275817871, '*298': 8.497355461120605, '*299': 8.497365951538086, '*3': 8.501371383666992, '*30': 8.500882148742676, '*300': 8.49736213684082, '*301': 8.49736213684082, '*302': 8.49730396270752, '*303': 8.497294425964355, '*304': 8.497294425964355, '*305': 8.497184753417969, '*306': 8.497184753417969, '*307': 8.497184753417969, '*308': 8.497184753417969, '*309': 8.497184753417969, '*31': 8.500882148742676, '*310': 8.497190475463867, '*311': 8.497152328491211, '*312': 8.49712085723877, '*313': 8.49712085723877, '*314': 8.497085571289062, '*315': 8.497085571289062, '*316': 8.497085571289062, '*317': 8.497085571289062, '*318': 8.497077941894531, '*319': 8.497077941894531, '*32': 8.500880241394043, '*320': 8.49706745147705, '*321': 8.49706745147705, '*322': 8.49706745147705, '*323': 8.496978759765625, '*324': 8.496978759765625, '*325': 8.496988296508789, '*326': 8.49695873260498, '*327': 8.49695873260498, '*328': 8.49695873260498, '*329': 8.49695873260498, '*33': 8.500880241394043, '*330': 8.49695873260498, '*331': 8.49695873260498, '*332': 8.49695873260498, '*333': 8.496955871582031, '*334': 8.496955871582031, '*335': 8.496960639953613, '*336': 8.496960639953613, '*337': 8.496960639953613, '*338': 8.496953964233398, '*339': 8.496953964233398, '*34': 8.500880241394043, '*340': 8.496964454650879, '*341': 8.496933937072754, '*342': 8.496933937072754, '*343': 8.496898651123047, '*344': 8.49687671661377, '*345': 8.496785163879395, '*346': 8.496785163879395, '*347': 8.496785163879395, '*348': 8.496785163879395, '*349': 8.496785163879395, '*35': 8.500884056091309, '*350': 8.496785163879395, '*351': 8.49678897857666, '*352': 8.496676445007324, '*353': 8.496676445007324, '*354': 8.496588706970215, '*355': 8.496588706970215, '*356': 8.496580123901367, '*357': 8.496580123901367, '*358': 8.496580123901367, '*359': 8.496580123901367, '*36': 8.500884056091309, '*360': 8.496580123901367, '*361': 8.49642276763916, '*362': 8.49642276763916, '*363': 8.49642276763916, '*364': 8.496406555175781, '*365': 8.496430397033691, '*366': 8.496448516845703, '*367': 8.496448516845703, '*368': 8.496448516845703, '*369': 8.496413230895996, '*37': 8.500884056091309, '*370': 8.49639892578125, '*371': 8.49639892578125, '*372': 8.49639892578125, '*373': 8.49639892578125, '*374': 8.49639892578125, '*375': 8.496352195739746, '*376': 8.496352195739746, '*377': 8.496352195739746, '*378': 8.496352195739746, '*379': 8.496352195739746, '*38': 8.500860214233398, '*380': 8.496352195739746, '*381': 8.496352195739746, '*382': 8.496352195739746, '*383': 8.496353149414062, '*384': 8.496353149414062, '*385': 8.496321678161621, '*386': 8.496297836303711, '*387': 8.496237754821777, '*388': 8.496237754821777, '*389': 8.49616813659668, '*39': 8.500860214233398, '*390': 8.49616813659668, '*391': 8.4961576461792, '*392': 8.4961576461792, '*393': 8.4961576461792, '*394': 8.49609375, '*395': 8.496075630187988, '*396': 8.496075630187988, '*397': 8.496075630187988, '*398': 8.496075630187988, '*399': 8.496075630187988, '*4': 8.501371383666992, '*40': 8.500860214233398, '*400': 8.496064186096191, '*401': 8.496079444885254, '*402': 8.496079444885254, '*403': 8.496061325073242, '*404': 8.49602222442627, '*405': 8.49592113494873, '*406': 8.49592113494873, '*407': 8.495946884155273, '*408': 8.495946884155273, '*409': 8.495946884155273, '*41': 8.500860214233398, '*410': 8.495942115783691, '*411': 8.495942115783691, '*412': 8.495942115783691, '*413': 8.495942115783691, '*414': 8.495928764343262, '*415': 8.495914459228516, '*416': 8.495914459228516, '*417': 8.495914459228516, '*418': 8.495862007141113, '*419': 8.495800971984863, '*42': 8.500804901123047, '*420': 8.495829582214355, '*421': 8.495829582214355, '*422': 8.495829582214355, '*423': 8.495829582214355, '*424': 8.495783805847168, '*425': 8.495783805847168, '*426': 8.495759010314941, '*427': 8.495759010314941, '*428': 8.495759010314941, '*429': 8.495759010314941, '*43': 8.500804901123047, '*430': 8.495759010314941, '*431': 8.495683670043945, '*432': 8.495683670043945, '*433': 8.495683670043945, '*434': 8.495683670043945, '*435': 8.49563980102539, '*436': 8.49563980102539, '*437': 8.49563980102539, '*438': 8.49563980102539, '*439': 8.49563980102539, '*44': 8.500791549682617, '*440': 8.49563980102539, '*441': 8.495579719543457, '*442': 8.495579719543457, '*443': 8.495579719543457, '*444': 8.495579719543457, '*445': 8.495579719543457, '*446': 8.495579719543457, '*447': 8.495579719543457, '*448': 8.495579719543457, '*449': 8.495530128479004, '*45': 8.500797271728516, '*450': 8.495490074157715, '*451': 8.495490074157715, '*452': 8.495403289794922, '*453': 8.495404243469238, '*454': 8.495404243469238, '*455': 8.495404243469238, '*456': 8.495386123657227, '*457': 8.495386123657227, '*458': 8.495345115661621, '*459': 8.495296478271484, '*46': 8.50074291229248, '*460': 8.495296478271484, '*461': 8.495307922363281, '*462': 8.495274543762207, '*463': 8.495274543762207, '*464': 8.495274543762207, '*465': 8.495259284973145, '*466': 8.495259284973145, '*467': 8.495259284973145, '*468': 8.495259284973145, '*469': 8.495259284973145, '*47': 8.50074291229248, '*470': 8.495259284973145, '*471': 8.495250701904297, '*472': 8.495250701904297, '*473': 8.495205879211426, '*474': 8.495189666748047, '*475': 8.495183944702148, '*476': 8.495180130004883, '*477': 8.495180130004883, '*478': 8.495171546936035, '*479': 8.495190620422363, '*48': 8.50074291229248, '*480': 8.495190620422363, '*481': 8.495190620422363, '*482': 8.495136260986328, '*483': 8.495136260986328, '*484': 8.495136260986328, '*485': 8.495139122009277, '*486': 8.495139122009277, '*487': 8.495139122009277, '*488': 8.495139122009277, '*489': 8.495139122009277, '*49': 8.50074291229248, '*490': 8.495108604431152, '*491': 8.495108604431152, '*492': 8.495108604431152, '*493': 8.495108604431152, '*494': 8.495108604431152, '*495': 8.495108604431152, '*496': 8.495108604431152, '*497': 8.495058059692383, '*498': 8.495058059692383, '*499': 8.495058059692383, '*5': 8.501351356506348, '*50': 8.500726699829102, '*500': 8.49502944946289, '*501': 8.49502944946289, '*502': 8.49502944946289, '*503': 8.495006561279297, '*504': 8.495006561279297, '*505': 8.494993209838867, '*506': 8.494993209838867, '*507': 8.494993209838867, '*51': 8.500726699829102, '*52': 8.50068473815918, '*53': 8.50068473815918, '*54': 8.50068473815918, '*55': 8.50068473815918, '*56': 8.500676155090332, '*57': 8.500635147094727, '*58': 8.500604629516602, '*59': 8.500604629516602, '*6': 8.501274108886719, '*60': 8.500577926635742, '*61': 8.500565528869629, '*62': 8.50053882598877, '*63': 8.500545501708984, '*64': 8.500545501708984, '*65': 8.500545501708984, '*66': 8.500545501708984, '*67': 8.500545501708984, '*68': 8.500545501708984, '*69': 8.500402450561523, '*7': 8.501262664794922, '*70': 8.500380516052246, '*71': 8.500380516052246, '*72': 8.500340461730957, '*73': 8.500334739685059, '*74': 8.500287055969238, '*75': 8.500287055969238, '*76': 8.500186920166016, '*77': 8.500128746032715, '*78': 8.500128746032715, '*79': 8.500128746032715, '*8': 8.501262664794922, '*80': 8.500129699707031, '*81': 8.500129699707031, '*82': 8.500100135803223, '*83': 8.500035285949707, '*84': 8.500035285949707, '*85': 8.500035285949707, '*86': 8.499926567077637, '*87': 8.499926567077637, '*88': 8.499926567077637, '*89': 8.499892234802246, '*9': 8.501267433166504, '*90': 8.499892234802246, '*91': 8.499892234802246, '*92': 8.499894142150879, '*93': 8.499887466430664, '*94': 8.499887466430664, '*95': 8.499871253967285, '*96': 8.499770164489746, '*97': 8.499777793884277, '*98': 8.499777793884277, '*99': 8.499777793884277}

            exp_pa_dict = {'*0': 72.54618072509766, '*1': 72.5458984375, '*10': 72.54621124267578, '*100': 72.54096984863281, '*101': 72.54096221923828, '*102': 72.54052734375, '*103': 72.54006958007812, '*104': 72.54045867919922, '*105': 72.54045867919922, '*106': 72.54037475585938, '*107': 72.54037475585938, '*108': 72.54067993164062, '*109': 72.54067993164062, '*11': 72.54644775390625, '*110': 72.54067993164062, '*111': 72.54067993164062, '*112': 72.54067993164062, '*113': 72.54067993164062, '*114': 72.54067993164062, '*115': 72.54067993164062, '*116': 72.54067993164062, '*117': 72.54016876220703, '*118': 72.54037475585938, '*119': 72.54037475585938, '*12': 72.547119140625, '*120': 72.54056549072266, '*121': 72.54056549072266, '*122': 72.54084777832031, '*123': 72.54094696044922, '*124': 72.54084014892578, '*125': 72.54124450683594, '*126': 72.54103088378906, '*127': 72.54103088378906, '*128': 72.54103088378906, '*129': 72.5400161743164, '*13': 72.547119140625, '*130': 72.5400161743164, '*131': 72.5400161743164, '*132': 72.5400161743164, '*133': 72.54019927978516, '*134': 72.54019927978516, '*135': 72.54015350341797, '*136': 72.54015350341797, '*137': 72.53997039794922, '*138': 72.54032897949219, '*139': 72.54032897949219, '*14': 72.547119140625, '*140': 72.54021453857422, '*141': 72.54021453857422, '*142': 72.54021453857422, '*143': 72.54021453857422, '*144': 72.54021453857422, '*145': 72.54021453857422, '*146': 72.53984069824219, '*147': 72.53984069824219, '*148': 72.53958129882812, '*149': 72.53987884521484, '*15': 72.547119140625, '*150': 72.53987884521484, '*151': 72.53987884521484, '*152': 72.5394287109375, '*153': 72.5394287109375, '*154': 72.53932189941406, '*155': 72.53932189941406, '*156': 72.53958892822266, '*157': 72.53958892822266, '*158': 72.53958892822266, '*159': 72.53958892822266, '*16': 72.54747009277344, '*160': 72.53958892822266, '*161': 72.53949737548828, '*162': 72.53949737548828, '*163': 72.54027557373047, '*164': 72.54119110107422, '*165': 72.54067993164062, '*166': 72.54100036621094, '*167': 72.54100036621094, '*168': 72.54080963134766, '*169': 72.54080963134766, '*17': 72.54749298095703, '*170': 72.54080963134766, '*171': 72.54080963134766, '*172': 72.53984832763672, '*173': 72.53959655761719, '*174': 72.53941345214844, '*175': 72.53883361816406, '*176': 72.53913116455078, '*177': 72.53937530517578, '*178': 72.53985595703125, '*179': 72.53985595703125, '*18': 72.54743957519531, '*180': 72.54026794433594, '*181': 72.54026794433594, '*182': 72.54026794433594, '*183': 72.54010009765625, '*184': 72.54019165039062, '*185': 72.54019165039062, '*186': 72.54039001464844, '*187': 72.54039001464844, '*188': 72.54039001464844, '*189': 72.54039001464844, '*19': 72.54718780517578, '*190': 72.54039001464844, '*191': 72.54039001464844, '*192': 72.54012298583984, '*193': 72.54012298583984, '*194': 72.54012298583984, '*195': 72.54012298583984, '*196': 72.54012298583984, '*197': 72.54012298583984, '*198': 72.54006958007812, '*199': 72.54035949707031, '*2': 72.54558563232422, '*20': 72.54718780517578, '*200': 72.54035949707031, '*201': 72.54035949707031, '*202': 72.54064178466797, '*203': 72.54032135009766, '*204': 72.54026794433594, '*205': 72.54026794433594, '*206': 72.54026794433594, '*207': 72.54026794433594, '*208': 72.54026794433594, '*209': 72.54026794433594, '*21': 72.54718780517578, '*210': 72.54014587402344, '*211': 72.54014587402344, '*212': 72.54014587402344, '*213': 72.54014587402344, '*214': 72.54014587402344, '*215': 72.54014587402344, '*216': 72.54064178466797, '*217': 72.54064178466797, '*218': 72.54141235351562, '*219': 72.54141235351562, '*22': 72.54718780517578, '*220': 72.54093933105469, '*221': 72.54093933105469, '*222': 72.54093933105469, '*223': 72.54093933105469, '*224': 72.54061889648438, '*225': 72.54041290283203, '*226': 72.54041290283203, '*227': 72.54000854492188, '*228': 72.54046630859375, '*229': 72.54046630859375, '*23': 72.54718780517578, '*230': 72.54046630859375, '*231': 72.54094696044922, '*232': 72.54094696044922, '*233': 72.54092407226562, '*234': 72.54092407226562, '*235': 72.5408935546875, '*236': 72.5408935546875, '*237': 72.54131317138672, '*238': 72.54131317138672, '*239': 72.54127502441406, '*24': 72.54695129394531, '*240': 72.54150390625, '*241': 72.54150390625, '*242': 72.54150390625, '*243': 72.54146575927734, '*244': 72.54193115234375, '*245': 72.54193115234375, '*246': 72.54193115234375, '*247': 72.54193115234375, '*248': 72.54193115234375, '*249': 72.54193115234375, '*25': 72.54627990722656, '*250': 72.54222869873047, '*251': 72.54222869873047, '*252': 72.54244995117188, '*253': 72.54251098632812, '*254': 72.54251098632812, '*255': 72.5427474975586, '*256': 72.5426254272461, '*257': 72.54222106933594, '*258': 72.54222106933594, '*259': 72.54222106933594, '*26': 72.54627990722656, '*260': 72.54222106933594, '*261': 72.54263305664062, '*262': 72.54263305664062, '*263': 72.54263305664062, '*264': 72.54263305664062, '*265': 72.54263305664062, '*266': 72.54269409179688, '*267': 72.54336547851562, '*268': 72.54347229003906, '*269': 72.54347229003906, '*27': 72.54627990722656, '*270': 72.54338836669922, '*271': 72.5435791015625, '*272': 72.5435791015625, '*273': 72.54358673095703, '*274': 72.54358673095703, '*275': 72.54358673095703, '*276': 72.54357147216797, '*277': 72.54421997070312, '*278': 72.54421997070312, '*279': 72.54415893554688, '*28': 72.54607391357422, '*280': 72.54415893554688, '*281': 72.54415893554688, '*282': 72.54415893554688, '*283': 72.54443359375, '*284': 72.54443359375, '*285': 72.54443359375, '*286': 72.54443359375, '*287': 72.54443359375, '*288': 72.54430389404297, '*289': 72.54385375976562, '*29': 72.54607391357422, '*290': 72.54385375976562, '*291': 72.54399108886719, '*292': 72.54380798339844, '*293': 72.54380798339844, '*294': 72.54380798339844, '*295': 72.54380798339844, '*296': 72.54374694824219, '*297': 72.54345703125, '*298': 72.54364776611328, '*299': 72.54367065429688, '*3': 72.54558563232422, '*30': 72.54607391357422, '*300': 72.54383087158203, '*301': 72.54383087158203, '*302': 72.54360961914062, '*303': 72.54415130615234, '*304': 72.54415130615234, '*305': 72.54395294189453, '*306': 72.54395294189453, '*307': 72.54395294189453, '*308': 72.54395294189453, '*309': 72.54395294189453, '*31': 72.54607391357422, '*310': 72.54415893554688, '*311': 72.54348754882812, '*312': 72.54342651367188, '*313': 72.54342651367188, '*314': 72.54335021972656, '*315': 72.54335021972656, '*316': 72.54335021972656, '*317': 72.54335021972656, '*318': 72.54357147216797, '*319': 72.54357147216797, '*32': 72.54622650146484, '*320': 72.54387664794922, '*321': 72.54387664794922, '*322': 72.54387664794922, '*323': 72.54393005371094, '*324': 72.54393005371094, '*325': 72.54389190673828, '*326': 72.54483795166016, '*327': 72.54483795166016, '*328': 72.54483795166016, '*329': 72.54483795166016, '*33': 72.54622650146484, '*330': 72.54483795166016, '*331': 72.54483795166016, '*332': 72.54483795166016, '*333': 72.54502868652344, '*334': 72.54502868652344, '*335': 72.54524993896484, '*336': 72.54524993896484, '*337': 72.54524993896484, '*338': 72.54592895507812, '*339': 72.54592895507812, '*34': 72.54622650146484, '*340': 72.5458755493164, '*341': 72.54499053955078, '*342': 72.54499053955078, '*343': 72.5450210571289, '*344': 72.54485321044922, '*345': 72.54437255859375, '*346': 72.54437255859375, '*347': 72.54437255859375, '*348': 72.54437255859375, '*349': 72.54437255859375, '*35': 72.54571533203125, '*350': 72.54437255859375, '*351': 72.54459381103516, '*352': 72.54454803466797, '*353': 72.54454803466797, '*354': 72.54415130615234, '*355': 72.54415130615234, '*356': 72.54439544677734, '*357': 72.54439544677734, '*358': 72.54439544677734, '*359': 72.54439544677734, '*36': 72.54571533203125, '*360': 72.54439544677734, '*361': 72.54432678222656, '*362': 72.54432678222656, '*363': 72.54432678222656, '*364': 72.54479217529297, '*365': 72.54451751708984, '*366': 72.54434967041016, '*367': 72.54434967041016, '*368': 72.54434967041016, '*369': 72.54428100585938, '*37': 72.54571533203125, '*370': 72.5446548461914, '*371': 72.5446548461914, '*372': 72.5446548461914, '*373': 72.5446548461914, '*374': 72.5446548461914, '*375': 72.5445327758789, '*376': 72.5445327758789, '*377': 72.5445327758789, '*378': 72.5445327758789, '*379': 72.5445327758789, '*38': 72.5455322265625, '*380': 72.5445327758789, '*381': 72.5445327758789, '*382': 72.5445327758789, '*383': 72.5447769165039, '*384': 72.5447769165039, '*385': 72.54471588134766, '*386': 72.54437255859375, '*387': 72.54450225830078, '*388': 72.54450225830078, '*389': 72.54415893554688, '*39': 72.5455322265625, '*390': 72.54415893554688, '*391': 72.54446411132812, '*392': 72.54446411132812, '*393': 72.54446411132812, '*394': 72.54413604736328, '*395': 72.54474639892578, '*396': 72.54474639892578, '*397': 72.54474639892578, '*398': 72.54474639892578, '*399': 72.54474639892578, '*4': 72.54558563232422, '*40': 72.5455322265625, '*400': 72.54512786865234, '*401': 72.54499053955078, '*402': 72.54499053955078, '*403': 72.54547119140625, '*404': 72.5450210571289, '*405': 72.54446411132812, '*406': 72.54446411132812, '*407': 72.54449462890625, '*408': 72.54449462890625, '*409': 72.54449462890625, '*41': 72.5455322265625, '*410': 72.54468536376953, '*411': 72.54468536376953, '*412': 72.54468536376953, '*413': 72.54468536376953, '*414': 72.54508972167969, '*415': 72.54559326171875, '*416': 72.54559326171875, '*417': 72.54559326171875, '*418': 72.54500579833984, '*419': 72.54496765136719, '*42': 72.545654296875, '*420': 72.54496002197266, '*421': 72.54496002197266, '*422': 72.54496002197266, '*423': 72.54496002197266, '*424': 72.54521942138672, '*425': 72.54521942138672, '*426': 72.54429626464844, '*427': 72.54429626464844, '*428': 72.54429626464844, '*429': 72.54429626464844, '*43': 72.545654296875, '*430': 72.54429626464844, '*431': 72.54389953613281, '*432': 72.54389953613281, '*433': 72.54389953613281, '*434': 72.54389953613281, '*435': 72.54237365722656, '*436': 72.54237365722656, '*437': 72.54237365722656, '*438': 72.54237365722656, '*439': 72.54237365722656, '*44': 72.54613494873047, '*440': 72.54237365722656, '*441': 72.54206848144531, '*442': 72.54206848144531, '*443': 72.54206848144531, '*444': 72.54206848144531, '*445': 72.54206848144531, '*446': 72.54206848144531, '*447': 72.54206848144531, '*448': 72.54206848144531, '*449': 72.54217529296875, '*45': 72.54609680175781, '*450': 72.5418930053711, '*451': 72.5418930053711, '*452': 72.54092407226562, '*453': 72.54122161865234, '*454': 72.54122161865234, '*455': 72.54122161865234, '*456': 72.5416259765625, '*457': 72.5416259765625, '*458': 72.54183959960938, '*459': 72.5418472290039, '*46': 72.54585266113281, '*460': 72.5418472290039, '*461': 72.54179382324219, '*462': 72.54154205322266, '*463': 72.54154205322266, '*464': 72.54154205322266, '*465': 72.54129028320312, '*466': 72.54129028320312, '*467': 72.54129028320312, '*468': 72.54129028320312, '*469': 72.54129028320312, '*47': 72.54585266113281, '*470': 72.54129028320312, '*471': 72.5415267944336, '*472': 72.5415267944336, '*473': 72.54109191894531, '*474': 72.54157257080078, '*475': 72.54180145263672, '*476': 72.5419921875, '*477': 72.5419921875, '*478': 72.54230499267578, '*479': 72.54096984863281, '*48': 72.54585266113281, '*480': 72.54096984863281, '*481': 72.54096984863281, '*482': 72.54082489013672, '*483': 72.54082489013672, '*484': 72.54082489013672, '*485': 72.54109191894531, '*486': 72.54109191894531, '*487': 72.54109191894531, '*488': 72.54109191894531, '*489': 72.54109191894531, '*49': 72.54585266113281, '*490': 72.54019927978516, '*491': 72.54019927978516, '*492': 72.54019927978516, '*493': 72.54019927978516, '*494': 72.54019927978516, '*495': 72.54019927978516, '*496': 72.54019927978516, '*497': 72.53990936279297, '*498': 72.53990936279297, '*499': 72.53990936279297, '*5': 72.54608154296875, '*50': 72.546142578125, '*500': 72.54069519042969, '*501': 72.54069519042969, '*502': 72.54069519042969, '*503': 72.5416488647461, '*504': 72.5416488647461, '*505': 72.54203033447266, '*506': 72.54203033447266, '*507': 72.54203033447266, '*51': 72.546142578125, '*52': 72.5461196899414, '*53': 72.5461196899414, '*54': 72.5461196899414, '*55': 72.5461196899414, '*56': 72.54659271240234, '*57': 72.54618835449219, '*58': 72.54581451416016, '*59': 72.54581451416016, '*6': 72.54618072509766, '*60': 72.54548645019531, '*61': 72.54580688476562, '*62': 72.5449447631836, '*63': 72.54444122314453, '*64': 72.54444122314453, '*65': 72.54444122314453, '*66': 72.54444122314453, '*67': 72.54444122314453, '*68': 72.54444122314453, '*69': 72.54393768310547, '*7': 72.54672241210938, '*70': 72.54375457763672, '*71': 72.54375457763672, '*72': 72.54331970214844, '*73': 72.54352569580078, '*74': 72.54351806640625, '*75': 72.54351806640625, '*76': 72.54296875, '*77': 72.5428695678711, '*78': 72.5428695678711, '*79': 72.5428695678711, '*8': 72.54672241210938, '*80': 72.54236602783203, '*81': 72.54236602783203, '*82': 72.54192352294922, '*83': 72.54167938232422, '*84': 72.54167938232422, '*85': 72.54167938232422, '*86': 72.54149627685547, '*87': 72.54149627685547, '*88': 72.54149627685547, '*89': 72.5408706665039, '*9': 72.54621124267578, '*90': 72.5408706665039, '*91': 72.5408706665039, '*92': 72.5411148071289, '*93': 72.5413818359375, '*94': 72.5413818359375, '*95': 72.54183197021484, '*96': 72.54070281982422, '*97': 72.54064178466797, '*98': 72.54064178466797, '*99': 72.54064178466797}

            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed=self.filter_report(report)      
        self.report = report

        add_to_dict(self, output = test_dict, dataset = "E2E6.1.00034.S_tclean.ms")
        
        self.modify_dict(test_dict, 'test_standard_cube', self.parallel)
        
        print('--> check_dict_vals finished.')
        
        test_dict['test_standard_cube']['report'] = report
        test_dict['test_standard_cube']['images'] = []
        
        print('--> dict_mods finished.')

        img = shutil._basename(img)
        self.mom8_creator(image=img+'.image', range_list=[0.3, 1.0])    
        self.mom8_creator(image=img+'.residual', range_list=[0.3, 1.0])

        print('--> mom8_creator finished.')

        test_dict['test_standard_cube']['images'].extend((img+'.image.moment8.png',img+'.residual.moment8.png'))
        test_dict['test_standard_cube']['images'].append(img+'.image.profile.png')

        if savemetricdict:
            ### serialize ndarray in mask_stats_dcit
            mask_stats_mod_dict = copy.deepcopy(mask_stats_dict)
            mask_stats_mod_dict['mask'] = mask_stats_dict['mask'].tolist()
            #create a nested dictionary containing exp dictionaries to save
            savedict = {}
            # list of stats to save
            # im_stats, mask_stats, pb_stats, psf_stats,\
            # model_stats, resid_stats, sumwt_stats]

            savedict['im_stats_dict']=im_stats_dict
            savedict['exp_im_stats_dict']=exp_im_stats
#            savedict['mask_stats_dict']=mask_stats_mod_dict
#            savedict['pb_stats_dict']=pb_stats_dict
#            savedict['psf_stats_dict']=psf_stats_dict
#            savedict['model_stats_dict']=model_stats_dict
#            savedict['resid_stats_dict']=resid_stats_dict
#            savedict['sumwt_stats_dict']=sumwt_stats_dict
            
            if self.parallel:
                savedict['bmin_dict']=bmin_dict
                savedict['bmaj_dict']=bmaj_dict
                savedict['pa_dict']=pa_dict

            self.save_dict_to_file(test_name,savedict, test_name+'_current_metrics')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

def suite():
     return [Test_standard, Test_mosaic]

def test():
    standard = Test_standard()
    standard.setUp()
    standard.test_standard_cube()


# Main #
if __name__ == '__main__':
    #unittest.main()
    test()
