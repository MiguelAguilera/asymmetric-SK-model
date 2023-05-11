import numpy as np
from matplotlib import pyplot as plt
from SK import *

plt.rc('text', usetex=True)
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 20})
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


size = 1000     # Network size
H0 = 0.5       # Uniform distribution of fields parameter
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
            f = Gaussian_integral(Gtanh,beta*H0,beta*Js,Nint=50000,xmax=15)
        if f > H0/J0:
            beta_max = beta
        else:
            beta_min = beta
            
        error = beta_max-beta_min
        
        print(ii,beta, f, error)
    
    beta_c[ii] = beta


M, Q, sig, S, S_r= np.load('data-heatmap-H0-'+str(H0)+'.npy')


plt.rc('text', usetex=True)
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 20})



plt.figure()
plt.axis([0,4,0,1])
plt.imshow(M,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,1],vmin=0,vmax=1)
plt.colorbar()
plt.plot(beta_c[beta_c<beta_max/2], Js_list[beta_c<beta_max/2],'--',color='0.5',lw=2.5)
plt.plot([0,4], np.ones(2)*0.741056711631245,':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta J$', rotation=0, labelpad=20)
plt.title(r'$m$', pad=8)
plt.savefig('img/m1.pdf', bbox_inches='tight')

plt.figure()
plt.axis([0,4,0,1])
plt.imshow(Q,cmap='inferno_r',interpolation='None',aspect='auto',origin='lower',extent=[0,4,0,1],vmin=0,vmax=1)
plt.colorbar()
plt.plot(beta_c[beta_c<beta_max/2], Js_list[beta_c<beta_max/2],'--',color='0.5',lw=2.5)
plt.plot([0,4], np.ones(2)*0.741056711631245,':',color='0',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\Delta J$', rotation=0, labelpad=20)
plt.title(r'$q$', pad=8)
plt.savefig('img/q1.pdf', bbox_inches='tight')



plt.show()
