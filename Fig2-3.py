import numpy as np
from matplotlib import pyplot as plt
from SK import *
import matplotlib.colors as colors

plt.rc('text', usetex=True)
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 20})



size = 1000     # Network size
H0 = 0       # Uniform distribution of fields parameter
J0 = 1.0        # Average value of couplings

interval = 100
Js_list = np.linspace(0, 1, interval)      # Standard deviation of couplings
beta_c = np.zeros(interval)

grid_num = 100

M = np.zeros((grid_num,grid_num))     # magnetization
Q = np.zeros((grid_num,grid_num))     # self-correlation
sig = np.zeros((grid_num,grid_num))   # entropy production
S = np.zeros((grid_num,grid_num))   # entropy rate
S_r = np.zeros((grid_num,grid_num))   # reverse entropy rate
sigb = np.zeros((grid_num,grid_num))

for ii in range(interval): # Compute critical points
    error=1
    beta_max= 10
    beta_min = 0

    Js = 0
    Js = Js_list[ii]

    while error >1E-10:
        beta = 0.5*(beta_max + beta_min)
        if Js==0:
            f = np.tanh(beta*H0)
        else:
            f = Gaussian_integral(Gtanh2,0,beta*Js,Nint=50000,xmax=15)
        if f > 1/beta/J0:
            beta_max = beta
        else:
            beta_min = beta
            
        error = beta_max-beta_min
    
    beta_c[ii] = beta

M, Q, sig, S, S_r,_= np.load('data/data-heatmap-H0-'+str(H0)+'.npy')

plt.rc('text', usetex=True)
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 20})

plt.figure()
plt.axis([0,4,0,1])
plt.imshow(M,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,1],vmin=0,vmax=1)
plt.colorbar()
plt.plot(beta_c[beta_c<beta_max/2], Js_list[beta_c<beta_max/2],'--',color='0.5',lw=2.5)
plt.plot([0,4], np.ones(2)*0.795019386059721,':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta J$', rotation=0, labelpad=20)
plt.title(r'$m$', pad=8)
plt.savefig('img/Fig2a.pdf', bbox_inches='tight')

plt.figure()
plt.axis([0,4,0,1])
plt.imshow(Q,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,1],vmin=0,vmax=1)
plt.colorbar()
plt.plot(beta_c[beta_c<beta_max/2], Js_list[beta_c<beta_max/2],'--',color='0.5',lw=2.5)
plt.plot([0,4], np.ones(2)*0.795019386059721,':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta J$', rotation=0, labelpad=20)
plt.title(r'$q$', pad=8)
plt.savefig('img/Fig2b.pdf', bbox_inches='tight')

plt.figure()
plt.axis([0,4,0,1])
plt.imshow(Q-M*M,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,1], norm=colors.PowerNorm(gamma=0.5, vmin=0,vmax=1))
plt.colorbar()
plt.plot(beta_c[beta_c<beta_max/2], Js_list[beta_c<beta_max/2],'--',color='0.5',lw=2.5)
plt.plot([0,4], np.ones(2)*0.795019386059721,':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta J$', rotation=0, labelpad=20)
plt.title(r'$q-m^2$', pad=8)
plt.savefig('img/FigS5.pdf', bbox_inches='tight')

plt.figure()
plt.axis([0,4,0,1])
plt.imshow(S,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,1],vmin=0)
plt.colorbar()
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta J$', rotation=0, labelpad=20)
plt.title(r'$S_{u|u-1}/N$', pad=8)
plt.savefig('img/Fig3a.pdf', bbox_inches='tight')


plt.figure()
plt.axis([0,4,0,1])
plt.imshow(S_r,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,1],vmin=0)
plt.colorbar()
#plt.plot(beta_c[beta_c<beta_max/2], Js_list[beta_c<beta_max/2],'--',color='0.5',lw=2.5)
#plt.plot([0,4], np.ones(2)*0.795019386059721,':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta J$', rotation=0, labelpad=20)
plt.title(r'$S_{u|u-1}^r/N$', pad=8)
plt.savefig('img/Fig3b.pdf', bbox_inches='tight')

plt.figure()
plt.axis([0,4,0,1])
plt.imshow(sig,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,1], norm=colors.PowerNorm(gamma=0.5))
plt.colorbar()
plt.plot(beta_c[beta_c<5], Js_list[beta_c<5],linestyle='--', dashes=(6, 6) ,color='1.',lw=1)
#plt.plot(beta_c[beta_c<beta_max/2], Js_list[beta_c<beta_max/2],linestyle='--', dashes=(4, 10) ,color='0.5',lw=1)
#plt.plot([0,4], np.ones(2)*0.795019386059721,':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta J$', rotation=0, labelpad=20)
plt.title(r'$\sigma_{u}/N$', pad=8)
plt.savefig('img/Fig3c.pdf', bbox_inches='tight')


plt.show()
