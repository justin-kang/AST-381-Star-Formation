import os
import re # regex
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astrodendro import Dendrogram
from astropy import constants as const
from astropy.io import fits
from scipy.optimize import curve_fit
# ignore 'BLANK' from FITS reading
import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter('ignore', category=VerifyWarning)



# PROBLEM 1a
# get the cores from the data
data = fits.getdata('data/ngc1333_c18o.fits')
# http://iopscience.iop.org/article/10.1086/375455/pdf
d = Dendrogram.compute(data,min_value=0.14*3,min_npix=66/2)
cores = []
descendants = []
for i in d.trunk:
    for j in i.descendants:
        descendants.append(j)
for child in descendants:
    if child.is_leaf:
        cores.append(child)
# plot the cores on top of an integrated intensity map
p = d.plotter()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
def _intensity():
    with fits.open('data/ngc1333_c18o.fits') as file:
        stp = file[0].header['CDELT3']/1000
    imap = np.zeros((data.shape[1],data.shape[2]))
    for i in range(imap.shape[0]):
        for j in range(imap.shape[1]):
            imap[i,j] = np.trapz(2*data[:,i,j],dx=stp)
    return imap
ax.imshow(_intensity(),cmap=plt.cm.Blues)
colors = ['r','g','b','c','m','y','k','w']
for i,core in enumerate(cores):
    p.plot_contour(ax,structure=core.idx,color=colors[i%len(colors)])
plt.tight_layout()
plt.savefig('data/cores.png')


# PROBLEM 1b
def get_core(index):
    x,y,z = cores[index].indices()
    cube = np.zeros((len(x),len(y),len(z)))
    for i,xi in enumerate(x):
        for j,yj in enumerate(y):
            for k,zk in enumerate(z):
                cube[i,j,k] = data[xi,yj,zk]
    return cube
# calculate the mass of the core based on the intensity
def mass(index):
    with fits.open('data/ngc1333_c18o.fits') as file:
        stp = file[0].header['CDELT3']/1000
    data = get_core(index)
    def _intensity():
        imap = np.zeros((data.shape[1],data.shape[2]))
        for i in range(imap.shape[0]):
            for j in range(imap.shape[1]):
                imap[i,j] = np.trapz(2*data[:,i,j],dx=stp)
        return imap
    def _temperature():
        def _Tex(Tb):
            return 5.5/np.log(1+5.5/(Tb+0.82))
        temp = np.zeros((data.shape[1],data.shape[2]))
        for i in range(temp.shape[0]):
            for j in range(temp.shape[1]):
                Tb = max(2*data[:,i,j])
                temp[i,j] = _Tex(Tb)
        return temp
    mH2 = (const.m_p.value + const.m_e.value)*2
    Msun = const.M_sun.value
    pc = const.pc.cgs.value
    Nco = 3.0e14*_intensity()/(1-np.exp(-5.3/_temperature())) * 12.3/(2e-6)
    num = np.sum(Nco*(47*350*pc/(180/np.pi*60*60))**2)
    return num*mH2/Msun/1000000
masses = []
for i in range(len(cores)):
    masses.append(mass(i))
hist,bins = np.histogram(np.log10(masses))
x = np.array(bins[:-1]+bins[1:])/2
# get the slope of the CMF
def f(x,m,b):
    return m*x+b
popt,pcov = curve_fit(f,x,hist)
print(popt[0],np.sqrt(np.diag(pcov))[0])
# graph out data + fit
plt.figure()
plt.plot(x,hist,'.')
plt.plot(x,popt[0]*x+popt[1])
plt.xlabel(r'Core Mass ($\log(M_{\odot})$)')
plt.ylabel('Number of cores')
plt.tight_layout()
plt.savefig('data/cmf.png')



# PROBLEM 2b
time = []
lumm = []
with open('data/const_imf.quanta1') as file:
    counter = 0
    for line in file:
        # ignore non numbers in the file
        if counter < 7:
            counter += 1
            continue
        # ignore spacing and newline, just get the data
        line = re.sub(' +',' ',line)[1:]
        line = re.sub('\n','',line)
        line = line.split(' ')
        time.append(float(line[0]))
        lumm.append(float(line[-1]))
plt.figure()
plt.plot(np.log10(time),lumm,'.')
plt.plot(np.log10(time),[43.41]*len(time))
plt.xlabel(r'$\log($Time$)$ (years)')
plt.ylabel('Exponent of luminosity (ergs/s)')
plt.tight_layout()
plt.savefig('data/const_lum.png')
lum1 = np.array(lumm)

# Problem 2c
time = []
lumm = []
with open('data/heavy_imf.quanta1') as file:
    counter = 0
    for line in file:
        # ignore non numbers in the file
        if counter < 7:
            counter += 1
            continue
        # ignore spacing and newline, just get the data
        line = re.sub(' +',' ',line)[1:]
        line = re.sub('\n','',line)
        line = line.split(' ')
        time.append(float(line[0]))
        lumm.append(float(line[-1]))
lum2 = np.array(lumm)
# fractional difference in exponent at 1 Gyr
print(np.mean(lum2/lum1))
# fractional difference in luminosity at 1 Gyr
print(10**(lum2[-1]-lum1[-1]))