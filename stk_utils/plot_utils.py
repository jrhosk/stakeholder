import casatasks
import casatools
import pylab as pl
import numpy as np
import subprocess

from astropy.io import fits
from astropy.wcs import WCS

def make_moment_plot(imname='',chan=0):
    casatasks.immoments(imagename = imname, moments = 8, outfile = imname+'.moment8')

    ia = casatools.image()

    ia.open(imname + '.moment8')
    pix = ia.getchunk()[:,:,0,chan]
    csys = ia.coordsys()
    ia.close()
    shp = pix.shape

    rad_to_deg =  180/np.pi
    w = WCS(naxis=2)
    w.wcs.crpix = csys.referencepixel()['numeric'][0:2]
    w.wcs.cdelt = csys.increment()['numeric'][0:2]*rad_to_deg
    w.wcs.crval = csys.referencevalue()['numeric'][0:2]*rad_to_deg
    w.wcs.ctype = ['RA---SIN', 'DEC--SIN']

    pl.subplot(projection=w)

    p1 = int(shp[0]*0.25)
    p2 = int(shp[0]*0.75)

    pl.imshow(pix[p1:p2,p1:p2].transpose(), origin='lower',  cmap=pl.cm.viridis)
    pl.xlabel('Right Ascension')
    pl.ylabel('Declination')
    
    pl.imsave(imname + '.moment8.png', pix[p1:p2,p1:p2].transpose())
    subprocess.call('mogrify -trim '+imname+'.moment8.png', shell=True)