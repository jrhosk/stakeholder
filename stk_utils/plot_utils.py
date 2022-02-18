import casatasks
import casatools
import pylab as pl
import numpy as np
import subprocess

from astropy.io import fits
from astropy.wcs import WCS

def plot_image(imname:str, type="", chan=0, trim=False):
    """ 
        Utility function to produce a .png plot from a .image file.

    Args:
        imname (str, required): Image file name.
        type (str, optional): Additional image type qualifiers. Defaults to "".
        chan (int, optional): Channel number. Defaults to 0.
        trim (bool, optional): Apply imagemagik mogrify -trim to plot. Defaults to False.
    """

    ia = casatools.image()

    ia.open(imname + type)
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
    
    pl.imsave(imname + type + '.png', pix[p1:p2,p1:p2].transpose())

    if trim is True:
        subprocess.call('mogrify -trim ' + imname + type + '.png', shell=True)