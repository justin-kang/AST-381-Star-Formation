import warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from pvextractor import Path
from pvextractor import extract_pv_slice

'''
# PROBLEM 1
# (a) sketch how the disk looks
plt.figure()
x = np.linspace(0,0.01,100)
plt.xlim((0,0.01))
plt.plot(x,x**(3/2),c='C0')
plt.plot(x,-x**(3/2),c='C0')
plt.xlabel(r'$h/r$')s
plt.ylabel('Relative disk height')
plt.tight_layout()
plt.savefig('data/disk.png')
'''


# PROBLEM 2
# (a) produce an integrated intensity map of the outflow
# convert the data from Jy/beam to K (brightness temperature)
def _brightness_temperature():
    # ignore errors from FITS file being "incomplete"
    warnings.simplefilter('ignore')
    with fits.open('data/hh114_12co.fits') as file:
        warnings.simplefilter('default')
        # try to make no overflow errors happen
        data = file[0].data.astype(np.float64)
        hdr = file[0].header
    # github.com/astropy/astropy/blob/master/docs/units/equivalencies.rst
    freq = hdr['RESTFREQ']*u.Hz
    beam_area = (4*u.arcsec)**2
    return ((data*u.Jy/u.beam).to(u.K,
        u.brightness_temperature(freq,beam_area))).value
def intensity():
    data = _brightness_temperature()
    # ignore errors from FITS file being "incomplete"
    warnings.simplefilter('ignore')
    with fits.open('data/hh114_12co.fits') as file:
        warnings.simplefilter('default')
        hdr = file[0].header
    stp = hdr['CDELT3'] # [m/s]
    # get the intensity map
    imap = np.zeros((data.shape[-2],data.shape[-1]))
    for i in range(imap.shape[0]):
        for j in range(imap.shape[1]):
            imap[i,j] = np.trapz(data[0,:,i,j],dx=stp/1000)
    return imap
imap = intensity()
plt.figure()
plt.imshow(imap,origin='lower')
plt.set_cmap('gray')
plt.tight_layout()
plt.savefig('data/intensity.png')

# (b) make a "Hubble diagram" (p-v) of the outflow
# ignore the many warnings that come with pvextractor
warnings.simplefilter('ignore')
with fits.open('data/hh114_12co.fits') as file:
    data = file[0].data.astype(np.float64)
path = Path([(103,139),(133,127),(173,98)],width=10.)
pvslice = extract_pv_slice(data[0],path)
warnings.simplefilter('default')
plt.figure()
plt.imshow(pvslice.data,origin='lower')
plt.set_cmap('gray')
plt.tight_layout()
plt.savefig('data/hubble.png')

# (c) use the CO emission to calculate the outflow mass
def temperature():
    def _Tex(Tb):
        return 5.5/np.log(1+5.5/(Tb+0.82))
    data = _brightness_temperature()
    temp = np.zeros((data.shape[-2],data.shape[-1]))
    for i in range(temp.shape[0]):
        for j in range(temp.shape[1]):
            Tb = max(data[0,:,i,j])
            temp[i,j] = _Tex(Tb)
    return temp
def mass():
    mH2 = (const.m_p.value + const.m_e.value)*2
    Msun = const.M_sun.value
    pc = const.pc.cgs.value
    Nco = 3.0e14*intensity()/(1-np.exp(-5.3/temperature()))
    ratio = 1/(60*2e-6)
    num = np.sum(ratio*Nco*(4*460*pc/(180/np.pi*60*60))**2)
    return num*mH2/Msun
print('Outflow mass:',mass(),'solar masses')