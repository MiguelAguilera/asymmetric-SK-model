#!/usr/bin/env python3
"""
GPLv3 2020 Miguel Aguilera

This code computes the solution of the asymmetric SK model in the thermodynamic limit
"""

import numpy as np
import time
from matplotlib import pyplot as plt
from matplotlib import cm
from SK import *

def nsf(num, n=4):
    """n-Significant Figures"""
    numstr = ("{0:.%ie}" % (n - 1)).format(num)
    return float(numstr)

plt.rc('text', usetex=True)
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 14})



size = 1000     # Network size
H0 = 0.0       # Uniform distribution of fields parameter
J0 = 1.0        # Average value of couplings
#Js = 0.75      # Standard deviation of couplings
#Js=0.7410552838700823
Js=0.2

B = 201                    # Number of values of beta
T = 10000                    # Number of simulation time steps
#T = 2                    # Number of simulation time steps
betas = (1 + np.linspace(-1, 1, B)) *2     # Inverse temperature

M = np.zeros(B)     # magnetization
Q = np.zeros(B)     # self-correlation
sig = np.zeros(B)   # entropy production
S = np.zeros(B)   # entropy rate
S_r = np.zeros(B)   # reverse entropy rate

steps = np.arange(T)

plt.figure()

vS=[0.4,0.5,0.6,0.7,0.795019386059721]

cmap = cm.get_cmap('viridis_r')
colors=[]
Nv=len(vS)
for i in range(Nv):
	colors+=[cmap(float(i+1)/(Nv+1-1))]

maxsig=0
for ind,v in enumerate(vS):
    H0=v
    for ib in range(B):
        beta=round(betas[ib], 3)
        if beta==0:
            M[ib] = 0
            Q[ib] = 0
            S[ib] = np.log(2)
            S_r[ib] = np.log(2)
        else:
            q=0
            m=0.0001
            m_prev=m
            for t in steps:
                m_prev2=m_prev
                m_prev=m
                q_prev=q

                if H0==0:
                    m = Gaussian_integral(Gtanh,beta*J0*m_prev,beta*Js)
                    q = Gaussian_integral2D(Gtanhtanh,0, beta*J0*m_prev, beta*J0*m_prev2,beta*Js, q_prev)
                elif Js==0:
                    m = (np.log(np.cosh(beta*H0 + beta*J0*m_prev))-np.log(np.cosh(-beta*H0 + beta*J0*m_prev)))/(2*beta*H0)
                    q = (beta*H0*2 -np.tanh(beta*H0+ beta*J0*m_prev) + np.tanh(-beta*H0+ beta*J0*m_prev))/(2*beta*H0)
                else:
                    m = (Gaussian_integral(Gm, beta*H0 + beta*J0*m_prev, beta*Js)-Gaussian_integral(Gm,-beta*H0 + beta*J0*m_prev, beta*Js))/(2*beta*H0)
                    q = (Gaussian_integral2D(Gmm, beta*H0, beta*J0*m_prev, beta*J0*m_prev2, beta*Js, q_prev) - Gaussian_integral2D(Gmm, -beta*H0, beta*J0*m_prev, beta*J0*m_prev2, beta*Js, q_prev) )/(2*beta*H0)

                update_error = np.abs(m-m_prev)
                if update_error < 1E-8:
#                if update_error < 1E-6:
                    break
                    
            if H0==0:
                k = Gaussian_integral(Gtanh2,beta*J0*m,beta*Js)
                p = Gaussian_integral(Phi1,beta*H0+beta*J0*m,beta*Js)
                S[ib] = p
                sig[ib] = (beta*Js)**2*(1-q)*k
            elif Js==0:
                k = (1-np.tanh(beta*H0 + beta*J0*m)**2)
                p = ( -((beta*H0 + beta*J0*m) * np.log(1+np.exp(2*(beta*H0 + beta*J0*m))) -fermi_poly2(2*(beta*H0 + beta*J0*m))) + ((-beta*H0 + beta*J0*m)* np.log(1+np.exp(2*(-beta*H0 + beta*J0*m))) -fermi_poly2(2*(-beta*H0 + beta*J0*m))))/(2*beta*H0)
                S[ib] = p
                sig[ib] = 0
            else:
                k= (Gaussian_integral(Gtanh,beta*H0 + beta*J0*m,beta*Js)-Gaussian_integral(Gtanh,-beta*H0 + beta*J0*m,beta*Js))/(2*beta*H0)
                p = (Gaussian_integral(Phi,beta*H0+beta*J0*m,beta*Js)-Gaussian_integral(Phi, -beta*H0+beta*J0*m,beta*Js))/(2*beta*H0)
                S[ib] = p
                sig[ib] = (beta*Js)**2*(1-q)*k
                
                
            M[ib] =m
            Q[ib] =q
            S_r[ib] = S[ib] + sig[ib]
            
       
        
        print(ind,v,beta, M[ib], Q[ib], sig[ib] , S[ib],S_r[ib])

    maxsig=max(maxsig,np.max(sig))
    print(colors[ind])
    plt.plot(betas,sig,color=colors[ind],label=r'$\Delta J = '+str(v)+'$')
plt.ylabel(r'$\sigma_{u}/N$',labelpad=20,rotation=0)
plt.xlabel(r'$\beta$')
plt.axis([np.min(betas),np.max(betas),0,maxsig])
plt.legend()
plt.savefig('img/Fig3d.pdf', bbox_inches='tight')
plt.show()
        
