#!/usr/bin/env python3
"""
GPLv3 2020 Miguel Aguilera

This code computes the solution of the asymmetric SK model in the thermodynamic limit
"""

import numpy as np
import time
from matplotlib import pyplot as plt
from SK import *

plt.rc('text', usetex=True)
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 20})
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')



def nsf(num, n=4):
    """n-Significant Figures"""
    numstr = ("{0:.%ie}" % (n - 1)).format(num)
    return float(numstr)

size = 1000     # Network size
H0 = 1       # Uniform distribution of fields parameter
J0 = 1.0        # Average value of couplings
Js = 0.2      # Standard deviation of couplings
#Js=0.7410552838700823 
#Js=0.5

B = 201                    # Number of values of beta
T = 100000                    # Number of simulation time steps
alphas = (1 + np.linspace(-1, 1, B)*1) * 1.5/2     # Inverse temperature

M = np.zeros(B)     # magnetization
Q = np.zeros(B)     # self-correlation
sig = np.zeros(B)   # entropy production

steps = np.arange(T)
for ib in range(B):
    alpha=round(alphas[ib], 10)
    if True:
        q=0
        m=0.0001
        m_prev=m
        for t in steps:
            m_prev2=m_prev
            m_prev=m
            q_prev=q

            if alpha==0:
                m = np.sign(J0*m_prev)
                q = np.sign(J0*m_prev)**2
                k =0
            else:
                m = (Gaussian_integral(Gm0,H0*alpha + J0*m_prev,Js)-Gaussian_integral(Gm0,-H0*alpha + J0*m_prev,Js))/(2*H0*alpha)
                if t%10==9 or (t%3==2 and t<100):
                    q = (Gaussian_integral2D(Gmm0, H0*alpha, J0*m_prev, J0*m_prev2,Js, q_prev) - Gaussian_integral2D(Gmm0, -H0*alpha, J0*m_prev,J0*m_prev2,Js, q_prev) )/(2*H0*alpha)
                k = (Gaussian_integral(Gtanh0,H0*alpha + J0*m_prev,Js)-Gaussian_integral(Gtanh0,-H0*alpha + J0*m_prev,Js))/(2*H0*alpha)
            update_error = np.abs(m-m_prev)
            if update_error < 1E-14 and t>100:
                break
#            if t == steps[-1]:
#                print('end')
        M[ib] =m
        Q[ib] =q
        sig[ib] = (Js)**2*(1-q)*k
    
        print(alpha, M[ib], Q[ib], sig[ib],t )

dM=np.gradient(M,alphas)
dQ=np.gradient(Q,alphas)

plt.figure()
plt.plot(alphas,M)
plt.ylabel(r'$m$',fontsize=18, rotation=0, labelpad=10)
plt.xlabel(r'$H_0$',fontsize=18)
plt.axis([alphas[0],alphas[-1],0,1])
plt.figure()
plt.plot(alphas,dM)

plt.figure()
plt.plot(alphas,Q)
plt.ylabel(r'$q$',fontsize=18, rotation=0, labelpad=10)
plt.xlabel(r'$\Delta H$',fontsize=18)
plt.figure()
plt.plot(alphas,dQ)

plt.figure()
plt.plot(alphas,sig,'k')
plt.ylabel(r'$\dfrac{1}{N}\left[\sigma_{u}\right]_{\mathbf{J}}$',labelpad=20,rotation=0)
plt.xlabel(r'$\Delta H$',fontsize=18)
plt.axis([0,np.max(alphas),0,np.max(sig)*1.05])
plt.savefig('img/Fig6e.pdf',bbox_inches='tight')

plt.show()
        
