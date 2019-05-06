import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy.io import fits
from scipy import linalg, integrate, optimize
import warnings

warnings.filterwarnings("ignore")
mpl.rcParams['agg.path.chunksize'] = 10000
plt.rc('text',usetex=True)
OUTPUT = [0,0,0,0,1]


# Problem 1b
def populations(dim,T,N,e):
    def _make_matrix(T,n,e):
        def _A(u,l):
            A10 = 7.67*10**(-8)
            return 3*(u**4)/(2*u+1)*A10
        gamma = 10**(-11)*np.sqrt(T)
        matrix = np.zeros((dim+1,dim+1))
        for i in range(dim+1):
            matrix[dim,i] = 1
        for r in range(dim):
            for c in range(dim+1):
                if r == c:
                    matrix[r,c] = n*gamma
                elif r == c-1:
                    matrix[r,c] = -(n*gamma+e*_A(c,r))
        return matrix
    b = [0 for i in range(dim)]+[1]
    curves = [None]*(dim+1)
    for i in range(len(curves)):
        if isinstance(T,int):
            curves[i] = np.ones((1,1))
        else:
            curves[i] = np.ones((len(T),len(N)))
    if isinstance(T,int):
        a = _make_matrix(T,N,e)
        x = linalg.lstsq(a,b)[0]
        for i in range(dim+1):
            curves[i] = x[i]
    else:
        for t in range(len(T)):
            for n in range(len(N)):
                a = _make_matrix(T[t],N[n],e)
                x = linalg.lstsq(a,b)[0]
                for i in range(dim+1):
                    curves[i][t,n] = x[i]
    return curves
if OUTPUT[0]:
    T = np.linspace(0.1,200,1000,endpoint=False)
    N = np.linspace(10**3,10**5,1000,endpoint=False)
    e = [0.01,1]
    for i in e:
        curves = populations(5,T,N,i)
        for j,pop in enumerate(curves):
            plt.figure()
            contour = plt.contour(T,N,pop)
            cbar = plt.colorbar(contour)
            cbar.ax.set_ylabel('Fractional level population')
            cbar.add_lines(contour)
            plt.title(r'Fractional level populations of the $J =\ $'+str(j)
                +r' state for $\epsilon =\ $'+str(i))
            plt.xlabel('Temperature (K)')
            plt.ylabel(r'H$_{2}$ Density (cm$^{-3}$)')
            plt.xscale('log')
            plt.yscale('log')
            plt.tight_layout()
            temp = '0'
            if i == 1:
                temp = '1'
            plt.savefig('img/1/'+str(j)+'_'+temp+'.png')


# Problem 1d
def optical_depth(u,l,e):
    def _A(u,l):
        A10 = 7.67*10**(-8)
        return 3*(u**4)/(2*u+1)*A10 
    def _E(J):
        hbar = const.hbar.value
        amu = const.u.value
        return (hbar)**2/(2*(6.86*amu)*(112.8e-12)**2)*J*(J+1)
    def _nu(u,l):
        h = const.h.value
        return (_E(u)-_E(l))/h
    c = const.c.value*100
    N = 10**18
    T = 10
    nH = 300
    n = populations(215,T,nH,e)
    nl = n[l]
    nu = n[u]
    b = 300000
    def _g(i):
        return 2*i+1
    return _A(u,l)*c**3/(8*np.pi*_nu(u,l)**3)*N*nu/b*(nl*_g(u)/(nu*_g(l))-1)
if OUTPUT[1]:
    print(optical_depth(1,0,1))
    print(optical_depth(2,1,1))


PATH = 'Turbulence_Archive/'
# Problem 2a
if OUTPUT[2]:
    # make the gas density PDF
    rho = fits.getdata(PATH+'data.1998.3d.hdf5_flatrho_256.fits').flatten()
    y,bins = np.histogram(np.log10(rho),bins='auto',density=True)
    x = np.array(bins[:-1]+bins[1:])/2
    plt.figure()
    plt.plot(x,np.log10(y),'.')
    plt.title('Gas Density Probability Distribution Function')
    plt.xlabel(r'$\log_{10}\rho$ (g/cm$^{-3}$)')
    plt.ylabel(r'$\log_{10}PDF$')
    plt.tight_layout()
    plt.savefig('img/2/pdf.png')
    # calculate Mach numbers and dispersion
    T = 10
    kb = const.k_B.value
    mH = const.m_e.value + const.m_p.value
    cs = np.sqrt(kb*T/(2.8*mH))
    xmax = x[np.where(y>0.5*max(y))[0][-1]]
    xdiff = xmax - x[np.where(y==max(y))[0][0]]
    xmin = xmax - 2*xdiff
    fwhm = xmax - xmin
    std = fwhm/(2*np.sqrt(2*np.log(2)))
    b = 0.2
    mach = np.sqrt(np.exp(std**2)-1)/b
    print('Mach number for b=0.2:',mach)
    print('Velocity dispersion for b=0.2:',mach*cs,'m/s')
    b = 0.5
    mach = np.sqrt(np.exp(std**2)-1)/b
    print('Mach number for b=0.5:',mach)
    print('Velocity dispersion for b=0.2:',mach*cs,'m/s')


# Problem 2b
if OUTPUT[3]:
    # generate the velocity power spectrum
    vx = fits.getdata(PATH+'data.1998.3d.hdf5_flatvx_256.fits')
    vy = fits.getdata(PATH+'data.1998.3d.hdf5_flatvy_256.fits')
    vz = fits.getdata(PATH+'data.1998.3d.hdf5_flatvz_256.fits')
    pc = const.pc.cgs.value
    L = 5*pc
    n = vx.shape[0]
    vkx = L*np.fft.fftn(vx)[0:n//2+1,0:n//2+1,0:n//2+1]
    vkx = np.abs(vkx/(n**3))**2
    vky = L*np.fft.fftn(vy)[0:n//2+1,0:n//2+1,0:n//2+1]
    vky = np.abs(vky/(n**3))**2
    vkz = L*np.fft.fftn(vz)[0:n//2+1,0:n//2+1,0:n//2+1]
    vkz = np.abs(vkz/(n**3))**2
    kx = np.fft.rfftfreq(n)*n/L*2
    ky = np.fft.rfftfreq(n)*n/L*2
    kz = np.fft.rfftfreq(n)*n/L*2
    kmin = 1.0/L
    kmax = 0.5*n/L
    bins = np.arange(kmin,kmax,kmin)
    kx,ky,kz = np.meshgrid(kx,ky,kz,indexing="ij")
    k = np.sqrt(kx**2+ky**2+kz**2)
    whichbin = np.digitize(k.flat, bins)
    ncount = np.bincount(whichbin)
    ps = np.zeros(len(ncount)-1)
    psk = 0.5/np.pi*k**2*(vkx+vky+vkz)
    for n in range(1,len(ncount)):
        ps[n-1] = np.mean(psk.flat[whichbin==n])
    k = np.array(bins[:-1]+bins[1:])/2
    ps = ps[1:len(bins)]
    x = np.log10(k/kmin)
    y = np.log10(ps)
    plt.figure()
    plt.plot(x,y)
    # slope of linear region
    imin = np.where(x>0.4)[0][0]
    imax = np.where(x>1.95)[0][0]
    def linear(x,m,b):
        return m*x+b
    popt, pcov = optimize.curve_fit(linear,x[imin:imax+1],y[imin:imax+1])
    plt.plot(x[np.where(x>0.25)[0][0]:np.where(x>2)[0][0]],
        popt[0]*x[np.where(x>0.25)[0][0]:np.where(x>2)[0][0]]+popt[1])
    plt.title('Velocity Power Spectrum')
    plt.xlabel(r'$\log_{10}(k/k_{min})$')
    plt.ylabel(r'$\log_{10}\Psi(k)$')
    plt.tight_layout()
    plt.savefig('img/2/ps.png')
    print('Slope of velocity power spectrum:',popt[0],
        '±',np.sqrt(np.diag(pcov))[0])
    # injection, dissipation scale
    print('k scales:',kmin*10**x[imin],kmin*10**x[imax])
    print('Energy injection scale:',2*np.pi/(kmin*10**x[imin])/pc,'pc')
    print('Energy dissipation scale:',2*np.pi/(kmin*10**x[imax])/pc,'pc')


# Problem 3
if OUTPUT[4]:
    cube = fits.getdata('mystery_12co.fits')
    disp = np.std(np.linspace(-9.99541,29.9948,cube.shape[0]))
    print('Cloud velocity dispersion:',disp,'km/s')
    T = 10
    kb = const.k_B.value
    mH = const.m_e.value + const.m_p.value
    cs = np.sqrt(kb*T/(2.8*mH))
    M = disp/cs
    print('Cloud Mach number:',M)
    def _virial(R,v):
        R = R*const.pc.value
        v *= 1000
        return int(5*R*v**2/(6*const.G.value)/const.M_sun.value)
    print('Virial mass: %.2E Msun'%_virial(10,disp))