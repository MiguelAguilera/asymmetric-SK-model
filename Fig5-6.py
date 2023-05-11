import numpy as np
from matplotlib import pyplot as plt
from SK import *

plt.rc('text', usetex=True)
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 20})
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')



size = 1000     # Network size

J0 = 1.0        # Average value of couplings
Js = 0.2    # standard deviation of couplings

interval = 1000  # interval of beta
beta_c = np.zeros(interval) # critical points
H0_list = np.linspace(0, 2, interval)

grid_num = 100  # number of grids

M = np.zeros((grid_num,grid_num))     # magnetization
Q = np.zeros((grid_num,grid_num))     # self-correlation
sig = np.zeros((grid_num,grid_num))   # entropy production
S = np.zeros((grid_num,grid_num))   # entropy rate
S_r = np.zeros((grid_num,grid_num))   # reverse entropy rate
sigb = np.zeros((grid_num,grid_num))

for ii in range(interval):

    error=1
    beta_max= 5
    beta_min = 0

    H0=0
    H0=H0_list[ii]
    
    while error >1E-10:
        beta = 0.5 * (beta_max + beta_min)
        if H0 == 0:
            f = Gaussian_integral(Gtanh2,beta*H0,beta*Js,Nint=50000,xmax=15) * beta * J0
        elif Js == 0:
            f = np.tanh(beta*H0) * (J0 / H0)
        else:
            f = Gaussian_integral(Gtanh,beta*H0,beta*Js,Nint=50000,xmax=15) * (J0 / H0)

        if f > 1:
            beta_max = beta
        else:
            beta_min = beta
            
        error = beta_max-beta_min
        
    print(ii)
    beta_c[ii] = beta


M, Q, sig, S, S_r, sigb= np.load('data/data-heatmap-Js-'+str(Js)+'.npy')

plt.rc('text', usetex=True)
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 20})

plt.figure()
plt.axis([0,4,0,1.5])
plt.imshow(M,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,2],vmin=0,vmax=1)
plt.colorbar()
plt.plot(beta_c, H0_list,'--',color='0.5',lw=2.5)
plt.plot([0,4], np.ones(2),':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta H$', rotation=0, labelpad=20)
plt.title(r'$m$', pad=8)
plt.savefig('img/m2.pdf', bbox_inches='tight')

plt.figure()
plt.axis([0,4,0,1.5])
plt.imshow(Q,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,2],vmin=0,vmax=1)
plt.colorbar()
plt.plot(beta_c, H0_list,'--',color='0.5',lw=2.5)
plt.plot([0,4], np.ones(2),':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta H$', rotation=0, labelpad=20)
plt.title(r'$q$', pad=8)
plt.savefig('img/q2.pdf', bbox_inches='tight')


plt.figure()
plt.axis([0,4,0,1.5])
plt.imshow(S,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,2],vmin=0)
plt.colorbar()
#plt.plot(beta_c, H0_list,'--',color='0.5',lw=2.5)
#plt.plot([0,4], np.ones(2),':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta H$', rotation=0, labelpad=20)
plt.title(r'$S_{u|u-1}/N$', pad=8)
plt.savefig('img/Scond2.pdf', bbox_inches='tight')


plt.figure()
plt.axis([0,4,0,1.5])
plt.imshow(S_r,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,2],vmin=0)
plt.colorbar()
#plt.plot(beta_c, H0_list,'--',color='0.5',lw=2.5)
#plt.plot([0,4], np.ones(2),':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta H$', rotation=0, labelpad=20)
plt.title(r'$S_{u|u-1}^r/N$', pad=8)
plt.savefig('img/Scond_r2.pdf', bbox_inches='tight')

plt.figure()
plt.axis([0,4,0,1.5])
plt.imshow(sig,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,2],vmin=0)
plt.colorbar()
#plt.plot(beta_c, H0_list,'--',color='0.5',lw=2.5)
#plt.plot([0,4], np.ones(2),':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta H$', rotation=0, labelpad=20)
plt.title(r'$\dfrac{1}{N}\left[\sigma_{u}\right]_{\mathbf{J}}$', pad=8)
plt.savefig('img/sigma_bath2.pdf', bbox_inches='tight')
#plt.figure()
#plt.axis([0,4,0,2])
#plt.imshow(sigb,cmap='inferno',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,2],vmin=0)
#plt.colorbar()
#plt.plot(beta_c, H0_list)
#plt.xlabel(r'$\beta$')
#plt.ylabel(r'$H_0$', rotation=0)

plt.show()
