import numpy as np
import matplotlib.pyplot as plt

plt.rc('text',usetex=True)
plt.rc('text.latex',preamble=r'\usepackage[version=4]{mhchem}'+
    r'\usepackage{siunitx}'+r'\usepackage{amsmath}')

with open('data/timeres.dat') as file:
    t = [float(i) for i in file.read().strip().split()]
t = np.log10(t[1:])

# PROBLEM 2a
with open('data/a/HCN.txt') as file:
    HCN = [float(i) for i in file.read().split('  ')][1:]
with open('data/a/HC3N.txt') as file:
    HC3N = [float(i) for i in file.read().split('  ')][1:]
with open('data/a/CH3CN.txt') as file:
    CH3CN = [float(i) for i in file.read().split('  ')][1:]

plt.figure()
ax1 = plt.subplot(211)
plt.plot(t,np.log10(np.divide(HC3N,HCN)))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log(\ce{HC3N}/\ce{HCN})$')
plt.setp(ax1.get_xticklabels(),visible=False)
ax2 = plt.subplot(212)
plt.plot(t,np.log10(np.divide(CH3CN,HCN)))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log(\ce{CH3CN}/\ce{HCN})$')
plt.xlabel(r'$\log{t}\ (\text{yr})$')
plt.tight_layout()
plt.savefig('data/a/ratios.png')



# PROBLEM 2b
HCN = []
HC3N = []
CH3CN = []
for i in range(-1,4):
    with open('data/b/temp/'+str(i)+'/HCN.txt') as file:
        HCN.append([float(i) for i in file.read().split('  ')][1:])
    with open('data/b/temp/'+str(i)+'/HC3N.txt') as file:
        HC3N.append([float(i) for i in file.read().split('  ')][1:])
    with open('data/b/temp/'+str(i)+'/CH3CN.txt') as file:
        CH3CN.append([float(i) for i in file.read().split('  ')][1:])

plt.figure()
ax1 = plt.subplot(311)
for hcn in HCN:
    plt.plot(t,np.log10(hcn))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{HCN}}$')
plt.legend((r'$10^{-1}$ K',r'$10^0$ K',r'$10^1$ K',r'$10^2$ K',r'$10^3$ K'),
    loc='lower left',ncol=2,frameon=False,fontsize='xx-small')
plt.setp(ax1.get_xticklabels(),visible=False)
ax2 = plt.subplot(312)
for hc3n in HC3N:
    plt.plot(t,np.log10(hc3n))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{HC3N}}$')
plt.legend((r'$10^{-1}$ K',r'$10^0$ K',r'$10^1$ K',r'$10^2$ K',r'$10^3$ K'),
    loc='lower left',ncol=2,frameon=False,fontsize='xx-small')
plt.setp(ax2.get_xticklabels(),visible=False)
ax3 = plt.subplot(313)
for ch3cn in CH3CN:
    plt.plot(t,np.log10(ch3cn))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{CH3CN}}$')
plt.legend((r'$10^{-1}$ K',r'$10^0$ K',r'$10^1$ K',r'$10^2$ K',r'$10^3$ K'),
    loc='lower left',ncol=2,frameon=False,fontsize='xx-small')
plt.xlabel(r'$\log{t}\ (\text{yr})$')
plt.tight_layout()
plt.savefig('data/b/t_abundances.png')

plt.figure()
ax1 = plt.subplot(211)
for i in range(len(HCN)):
    plt.plot(t,np.log10(np.divide(HC3N[i],HCN[i])))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{(\ce{HC3N}/\ce{HCN})}$')
plt.legend((r'$10^{-1}$ K',r'$10^0$ K',r'$10^1$ K',r'$10^2$ K',r'$10^3$ K'),
    loc='lower left',ncol=2,frameon=False,fontsize='xx-small')
plt.setp(ax1.get_xticklabels(),visible=False)
ax2 = plt.subplot(212)
for i in range(len(HCN)):
    plt.plot(t,np.log10(np.divide(CH3CN[i],HCN[i])))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{(\ce{CH3CN}/\ce{HCN})}$')
plt.legend((r'$10^{-1}$ K',r'$10^0$ K',r'$10^1$ K',r'$10^2$ K',r'$10^3$ K'),
    loc='lower left',ncol=2,frameon=False,fontsize='xx-small')
plt.xlabel(r'$\log{t}\ (\text{yr})$')
plt.tight_layout()
plt.savefig('data/b/t_ratios.png')

HCN = []
HC3N = []
CH3CN = []
for i in range(3,8):
    with open('data/b/size/'+str(i)+'/HCN.txt') as file:
        HCN.append([float(i) for i in file.read().split('  ')][1:])
    with open('data/b/size/'+str(i)+'/HC3N.txt') as file:
        HC3N.append([float(i) for i in file.read().split('  ')][1:])
    with open('data/b/size/'+str(i)+'/CH3CN.txt') as file:
        CH3CN.append([float(i) for i in file.read().split('  ')][1:])

plt.figure()
ax1 = plt.subplot(311)
for hcn in HCN:
    plt.plot(t,np.log10(hcn))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{HCN}}$')
plt.legend((r'$10^{-3}$ cm',r'$10^{-4}$ cm',r'$10^{-5}$ cm',r'$10^{-6}$ cm',
    r'$10^{-7}$ cm'),loc='lower left',ncol=2,frameon=False,fontsize='xx-small')
plt.setp(ax1.get_xticklabels(),visible=False)
ax2 = plt.subplot(312)
for hc3n in HC3N:
    plt.plot(t,np.log10(hc3n))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{HC3N}}$')
plt.legend((r'$10^{-3}$ cm',r'$10^{-4}$ cm',r'$10^{-5}$ cm',r'$10^{-6}$ cm',
    r'$10^{-7}$ cm'),loc='lower left',ncol=2,frameon=False,fontsize='xx-small')
plt.setp(ax2.get_xticklabels(),visible=False)
ax3 = plt.subplot(313)
for ch3cn in CH3CN:
    plt.plot(t,np.log10(ch3cn))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{CH3CN}}$')
plt.legend((r'$10^{-3}$ cm',r'$10^{-4}$ cm',r'$10^{-5}$ cm',r'$10^{-6}$ cm',
    r'$10^{-7}$ cm'),loc='lower left',ncol=2,frameon=False,fontsize='xx-small')
plt.xlabel(r'$\log{t}\ (\text{yr})$')
plt.tight_layout()
plt.savefig('data/b/s_abundances.png')

plt.figure()
ax1 = plt.subplot(211)
for i in range(len(HCN)):
    plt.plot(t,np.log10(np.divide(HC3N[i],HCN[i])))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{(\ce{HC3N}/\ce{HCN})}$')
plt.legend((r'$10^{-3}$ cm',r'$10^{-4}$ cm',r'$10^{-5}$ cm',r'$10^{-6}$ cm',
    r'$10^{-7}$ cm'),loc='lower left',ncol=2,frameon=False,fontsize='xx-small')
plt.setp(ax1.get_xticklabels(),visible=False)
ax2 = plt.subplot(212)
for i in range(len(HCN)):
    plt.plot(t,np.log10(np.divide(CH3CN[i],HCN[i])))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{(\ce{CH3CN}/\ce{HCN})}$')
plt.legend((r'$10^{-3}$ cm',r'$10^{-4}$ cm',r'$10^{-5}$ cm',r'$10^{-6}$ cm',
    r'$10^{-7}$ cm'),loc='lower left',ncol=2,frameon=False,fontsize='xx-small')
plt.xlabel(r'$\log{t}\ (\text{yr})$')
plt.tight_layout()
plt.savefig('data/b/s_ratios.png')



# PROBLEM 2c
HCN = []
HC3N = []
CH3CN = []
for i in range(15,20):
    with open('data/c/'+str(i)+'/HCN.txt') as file:
        HCN.append([float(i) for i in file.read().split('  ')][1:])
    with open('data/c/'+str(i)+'/HC3N.txt') as file:
        HC3N.append([float(i) for i in file.read().split('  ')][1:])
    with open('data/c/'+str(i)+'/CH3CN.txt') as file:
        CH3CN.append([float(i) for i in file.read().split('  ')][1:])

plt.figure()
ax1 = plt.subplot(311)
for hcn in HCN:
    plt.plot(t,np.log10(hcn))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{HCN}}$')
plt.legend((r'$10^{-15}$ \si{\per\second}',r'$10^{-16}$ \si{\per\second}',
    r'$10^{-17}$ \si{\per\second}',r'$10^{-18}$ \si{\per\second}',
    r'$10^{-19}$ \si{\per\second}'),loc='lower left',ncol=3,frameon=False,
    fontsize='xx-small')
plt.setp(ax1.get_xticklabels(),visible=False)
ax2 = plt.subplot(312)
for hc3n in HC3N:
    plt.plot(t,np.log10(hc3n))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{HC3N}}$')
plt.legend((r'$10^{-15}$ \si{\per\second}',r'$10^{-16}$ \si{\per\second}',
    r'$10^{-17}$ \si{\per\second}',r'$10^{-18}$ \si{\per\second}',
    r'$10^{-19}$ \si{\per\second}'),loc='lower left',ncol=3,frameon=False,
    fontsize='xx-small')
plt.setp(ax2.get_xticklabels(),visible=False)
ax3 = plt.subplot(313)
for ch3cn in CH3CN:
    plt.plot(t,np.log10(ch3cn))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{CH3CN}}$')
plt.legend((r'$10^{-15}$ \si{\per\second}',r'$10^{-16}$ \si{\per\second}',
    r'$10^{-17}$ \si{\per\second}',r'$10^{-18}$ \si{\per\second}',
    r'$10^{-19}$ \si{\per\second}'),loc='lower left',ncol=3,frameon=False,
    fontsize='xx-small')
plt.xlabel(r'$\log{t}\ (\text{yr})$')
plt.tight_layout()
plt.savefig('data/c/abundances.png')

plt.figure()
ax1 = plt.subplot(211)
for i in range(len(HCN)):
    plt.plot(t,np.log10(np.divide(HC3N[i],HCN[i])))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{(\ce{HC3N}/\ce{HCN})}$')
plt.legend((r'$10^{-15}$ \si{\per\second}',r'$10^{-16}$ \si{\per\second}',
    r'$10^{-17}$ \si{\per\second}',r'$10^{-18}$ \si{\per\second}',
    r'$10^{-19}$ \si{\per\second}'),loc='lower left',ncol=3,frameon=False,
    fontsize='xx-small')
plt.setp(ax1.get_xticklabels(),visible=False)
ax2 = plt.subplot(212)
for i in range(len(HCN)):
    plt.plot(t,np.log10(np.divide(CH3CN[i],HCN[i])))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{(\ce{CH3CN}/\ce{HCN})}$')
plt.legend((r'$10^{-15}$ \si{\per\second}',r'$10^{-16}$ \si{\per\second}',
    r'$10^{-17}$ \si{\per\second}',r'$10^{-18}$ \si{\per\second}',
    r'$10^{-19}$ \si{\per\second}'),loc='lower left',ncol=3,frameon=False,
    fontsize='xx-small')
plt.xlabel(r'$\log{t}\ (\text{yr})$')
plt.tight_layout()
plt.savefig('data/c/ratios.png')



# PROBLEM 2d
H2O = []
for i in range(-1,4):
    with open('data/d/temp/'+str(i)+'/H2O.txt') as file:
        H2O.append([float(i) for i in file.read().split('  ')][1:])

plt.figure()
for h2o in H2O:
    plt.plot(t,np.log10(h2o))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{H2O}}$')
plt.legend((r'$10^{-1}$ K',r'$10^0$ K',r'$10^1$ K',r'$10^2$ K',r'$10^3$ K'),
    loc='lower left',ncol=3,frameon=False,fontsize='small')
plt.xlabel(r'$\log{t}\ (\text{yr})$')
plt.tight_layout()
plt.savefig('data/d/t_abundances.png')

H2O = []
for i in range(3,8):
    with open('data/d/size/'+str(i)+'/H2O.txt') as file:
        H2O.append([float(i) for i in file.read().split('  ')][1:])

plt.figure()
for h2o in H2O:
    plt.plot(t,np.log10(h2o))
plt.xlim(t[0],t[-1])
plt.ylabel(r'$\log{\ce{H2O}}$')
plt.legend((r'$10^{-3}$ cm',r'$10^{-4}$ cm',r'$10^{-5}$ cm',r'$10^{-6}$ cm',
    r'$10^{-7}$ cm'),loc='lower left',ncol=3,frameon=False,fontsize='small')
plt.xlabel(r'$\log{t}\ (\text{yr})$')
plt.tight_layout()
plt.savefig('data/d/s_abundances.png')