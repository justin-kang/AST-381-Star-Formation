import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy.io import fits
import warnings

warnings.filterwarnings("ignore")
SHOW = [0,0,1,1,1]
molecules = ['12co','13co','c18o']

def pxrange(hdul):
    stp = hdul[0].header['CDELT3']
    pix = hdul[0].header['CRPIX3']
    val = hdul[0].header['CRVAL3']
    vmin = (val-3500)/stp
    vmax = (12000-val)/stp
    pxmin = int(np.ceil(pix-vmin))
    pxmax = int(np.ceil(pix+vmax))
    return pxmin,pxmax

# Make integrated intensity maps from the different line maps for the velocity 
# range 3.5 < v_LSR < 12.0 km/s.

def intensity(molecule):
    with fits.open('ngc1333_'+molecule+'.fits') as file:
        data = file[0].data
        stp = file[0].header['CDELT3']
        pxmin,pxmax = pxrange(file)
        imap = np.zeros((file[0].header['NAXIS2'],file[0].header['NAXIS1']))
        for i in range(imap.shape[0]):
            for j in range(imap.shape[1]):
                imap[i,j] = np.trapz(2*data[pxmin:pxmax,i,j],dx=stp/1000)
        return imap

if SHOW[0]:
    for molecule in molecules:
        I = intensity(molecule)
        plt.figure()
        plt.imshow(I)
        plt.colorbar()
        plt.tight_layout()
        plt.savefig('img/'+molecule+'_intensity.png')

# Obtain an excitation temperature map of the cloud using the peak emission of 
# the 12CO in each pixel.

def temperature():
    def _Tex(Tb):
        return 5.5/np.log(1+5.5/(Tb+0.82)) # (15.30)
    with fits.open('ngc1333_12co.fits') as file:
        data = file[0].data
        pxmin,pxmax = pxrange(file)
        temp = np.zeros((file[0].header['NAXIS2'],file[0].header['NAXIS1']))
        for i in range(temp.shape[0]):
            for j in range(temp.shape[1]):
                Tb = max(2*data[pxmin:pxmax,i,j])
                temp[i,j] = _Tex(Tb)
        return temp

if SHOW[1]:
    T = temperature()
    plt.figure()
    plt.imshow(T)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig('img/temp.png')

# Calculate the mass of the cloud for the same velocity range as above, using 
# each of the three maps and assuming all lines are optically thin.

def mass(molecule):
    mH2 = (const.m_p.value + const.m_e.value)*2
    Msun = const.M_sun.value
    pc = const.pc.cgs.value
    Nco = 3.0e14*intensity(molecule)/(1-np.exp(-5.3/temperature())) # (15.37)
    if molecule == '12co':
        ratio = 1/(60*2e-6)
    elif molecule == '13co':
        ratio = 1/(2e-6)
    elif molecule == 'c18o':
        ratio = 12.3/(2e-6)
    num = np.sum(ratio*Nco*(25*280*pc/(180/np.pi*60*60))**2)
    return num*mH2/Msun

if SHOW[2]:
    for molecule in molecules:
        print('Mass of',molecule,':',int(round(mass(molecule))),'Msun')

# Now obtain the mass using the 13CO (1-0) map, taking into consideration the 
# line optical depth, and compare any differences or similarities with the 
# mass estimate from the C18O (1-0) map (above).

def corrected_mass():
    X = 12.3
    F = 0
    with fits.open('ngc1333_13co.fits') as co13, \
        fits.open('ngc1333_c18o.fits') as c18o:
        data13 = co13[0].data
        data18 = c18o[0].data
        stp13 = co13[0].header['CDELT3']
        stp18 = c18o[0].header['CDELT3']
        pxmin13,pxmax13 = pxrange(co13)
        pxmin18,pxmax18 = pxrange(c18o)
        F = np.zeros((co13[0].header['NAXIS2'],co13[0].header['NAXIS1']))
        for i in range(F.shape[0]):
            for j in range(F.shape[1]):
                T13 = np.trapz(2*data13[pxmin13:pxmax13,i,j],dx=stp13/1000)
                T18 = np.trapz(2*data18[pxmin18:pxmax18,i,j],dx=stp18/1000)
                f = X*T18/T13
                if T13 == 0:
                    f = 0
                if f < 0:
                    f = 0
                F[i,j] = f
        F = np.mean(F)
    return F*mass('13co')

if SHOW[3]:
    print('Mass accounting for τ:',int(round(corrected_mass())),'Msun')

# Use your mass estimates and the velocity dispersion derived in HW 1 to 
# compute the virial parameter under the optically thick and thin assumptions. 

def alpha(M,R):
    G = const.G.value
    M *= const.M_sun.value
    sigma = 1700
    return 5*(sigma**2)*R/(G*M)

if SHOW[4]:
    R = 30*25*280*const.pc.value/(180/np.pi*60*60)
    Ml = mass('12co')
    Mu = corrected_mass()
    print('Optically thin virial parameter:',alpha(Ml,R))
    print('Optically thick virial parameter:',alpha(Mu,R))