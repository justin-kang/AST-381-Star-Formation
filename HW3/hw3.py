# http://docs.hyperion-rt.org/en/stable/tutorials/example_class1.html
import numpy as np
import matplotlib
import sys
matplotlib.use('agg')
import matplotlib.pyplot as plt
import astropy.constants as const
from scipy.integrate import simps
from hyperion.model import AnalyticalYSOModel, ModelOutput
from hyperion.util.constants import rsun, lsun, au, pc, msun, yr, c

RUN = [1,1]
c = np.longdouble(const.c.cgs.value)
#sys.setrecursionlimit(1500)

def _bin_search(func,l,r,value):
    mid = (l + r)/2
    val = func(mid)
    dif = abs(val - value)
    avg = abs(val + value)/2
    print('val',val)
    print('value',value)
    if (dif/avg) < 1/100:
        return mid
    elif val > value:
        return _bin_search(func,l,mid,value)
    else:
        return _bin_search(func,mid,r,value)

def _planck(v,T):
    h = np.longdouble(const.h.cgs.value)
    k = np.longdouble(const.k_B.cgs.value)
    coeff = 2*h*np.longdouble(np.power(v,3))/(c**2)
    exp = np.longdouble(1/np.longdouble(np.exp(h*v/(k*T))-1))
    return coeff*exp

def _temperature(v,T):
    return simps(v*_planck(v,T),v)/simps(_planck(v,T),v)

# Initalize the model
m = AnalyticalYSOModel()
# Set the stellar parameters
m.star.radius = 2*rsun
m.star.luminosity = 5*lsun
m.star.temperature = 6200.
# Add a flared disk
disk = m.add_flared_disk()
disk.mass = 0.01*msun
disk.rmin = 10.*m.star.radius
disk.rmax = 200.*au
disk.r_0 = 100.*au#m.star.radius
disk.h_0 = 5.*au
disk.p = -1.0
disk.beta = 1.25
disk.dust = 'kmh94_3.1_full.hdf5'
# Add an envelope
envelope = m.add_power_law_envelope()
envelope.mass = 0.4*msun
envelope.rmin = 200.*au
envelope.rmax = 1e4*au
envelope.power = -2
envelope.r_0 = disk.rmax
envelope.dust = 'kmh94_3.1_full.hdf5'
# Use raytracing to improve s/n of thermal/source emission
m.set_raytracing(True)
# Use the modified random walk
m.set_mrw(True, gamma=2.)
# Set up grid
m.set_spherical_polar_grid_auto(100,5,5)
# Set up SED
sed = m.add_peeled_images(sed=True, image=False)
sed.set_viewing_angles(np.linspace(0.,90.,10), np.repeat(45.,10))
sed.set_wavelength_range(150,0.02,2000.)
# Set number of photons
m.set_n_photons(initial=1e4, imaging=1e4,
    raytracing_sources=1e4, raytracing_dust=1e4)
# Set number of temperature iterations and convergence criterion
m.set_n_initial_iterations(10)
m.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)
# Write out file
file = open('temp.txt','w')
for mass in [0.001,0.01,0.1,1.]:
    disk.mass = mass*msun
    if RUN[0]:
        m.write('model_'+str(mass)+'.rtin')
        m.run('model_'+str(mass)+'.rtout', mpi=True)

    mo = ModelOutput('model_'+str(mass)+'.rtout')
    sed = mo.get_sed(inclination='all',aperture=-1,distance=300.*pc)
    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot(1,1,1)
    for i in range(sed.val.shape[0]):
        #ax.loglog(sed.wav,sed.val.transpose(),color='black')
        ax.loglog(sed.wav,sed.val[i,:],color='black')
    #ax.set_xlim(0.03,2000.)
    ax.set_xlim(0.1,1500.)
    #ax.set_ylim(2.e-15,1e-8)
    ax.set_ylim(2e-16,2e-9)
    ax.set_title(r'Model SED with '+str(mass)+' $M_{\odot}$ envlope')
    ax.set_xlabel(r'$\lambda$ [$\mu$m]')
    ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/cm$^2/s$]')
    if RUN[1]:
        fig.savefig('model_'+str(mass)+'_sed.png',bbox_inches='tight')
    
    v = c / np.longdouble(sed.wav*1e-4)
    max_index = sed.val.sum(axis=1).argmax()
    def _temp(T):
        return _temperature(v,T)
    ratio = simps(sed.val[max_index,:],v) \
        / simps(sed.val[max_index,:]/v,v)
    t = _bin_search(_temp,0.,m.star.temperature,ratio)
    file.write(str(t)+'\n')

file.close()