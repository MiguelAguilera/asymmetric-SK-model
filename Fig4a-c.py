#!/usr/bin/env python3
"""
GPLv3 2020 Miguel Aguilera

This code computes the solution of the asymmetric SK model in the thermodynamic limit
"""

import numpy as np
import time
from matplotlib import pyplot as plt
from matplotlib import colormaps as cm
from SK import *


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 15})





H0=0.0
Js=1.0
J0=1
R=400000
Tr=1

#R=10000
#Tr=6

B=101
#betas = 1 + np.linspace(-1, 1, B) * 0.3
betas = np.linspace(0, 4, B)

#T=100000

run = False


if run:
    T=10000
    M = np.zeros(B)
    sig = np.zeros(B)
    K = np.zeros(B)
    Q = np.zeros(B)

    for ib in range(len(betas)):
        print(ib)
        steps = np.arange(T)
        beta=betas[ib]

        m=0.5
        m_prev=0.5
        q=m**2
        k=0
        for t in steps:
            m_prev2=m_prev
            m_prev=m
            q_prev=q

            m = Gaussian_integral(Gtanh,beta*J0*m_prev,beta*Js)
            q = Gaussian_integral2D(Gtanhtanh,0, beta*J0*m_prev, beta*J0*m_prev2, beta*Js, q_prev)

        k = Gaussian_integral(Gtanh2,beta*J0*m,beta*Js)
        print(m,q,(beta*Js)**2*(1-q)*k)
        sig[ib] = (beta*Js)**2*(1-q)*k
        M[ib] =m
        Q[ib] =q

    np.savez('data/analytical_DeltaJ'+str(Js)+'.npz',M=M,Q=Q, sig=sig)

else:
    data = np.load('data/analytical_DeltaJ'+str(Js)+'.npz')
    M=data['M']
    Q=data['Q']
    sig=data['sig']

sizes=[32,64,128,256,512,1024]
S=len(sizes)
cmap = cm.get_cmap('plasma_r')
cm.get_cmap('plasma_r')
colors=[]
for i in range(S):
	colors+=[cmap((i)/(S-1))]
	
ms=[]
Cs=[]
Ds=[]
sigmas=[]

for size in sizes:

    filename = 'data/data-gamma1-'+str(H0)+'-gamma2-'+str(Js)+'-s-'+str(size)+'-R-'+str(R)+'-Tr-'+str(Tr)+'.npz'
    data=np.load(filename)
    ms+=[data['m']]
    Cs+=[data['C']]
    Ds+=[data['D']]
    sigmas+=[data['sigma']/size]


plt.figure(figsize=(6.4,4.8*2/3))
if Js==0.5:
    plt.axis([0,4,0,1])
else:
    plt.axis([0,4,0,np.max(ms)*1.1])
for i in range(len(sizes)):
    plt.plot(betas,ms[i],lw=1.5,color=colors[i],label=r'$N='+str(sizes[i])+'$')
plt.plot(betas,M,'k',lw=2,label=r'$N\to\infty$')
plt.xlabel(r'$\beta$')
plt.legend(handlelength=1.5)
plt.ylabel(r'$m$', rotation=0, labelpad=15)
plt.savefig('img/simulation_m_DeltaJ'+str(Js)+'.pdf', bbox_inches='tight')


plt.figure()
for i in range(len(sizes)):
    plt.plot(betas,Cs[i]-ms[i]**2,lw=1.5,color=colors[i])
    
plt.figure(figsize=(6.4,4.8*2/3))
if Js==0.5:
    plt.axis([0,4,0,1])
else:
    plt.axis([0,4,0,np.max(Ds)*1.1])
for i in range(len(sizes)):
    plt.plot(betas,Ds[i],lw=1.5,color=colors[i])
plt.plot(betas,Q,'k',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$q$', rotation=0, labelpad=20)
plt.savefig('img/simulation_q_DeltaJ'+str(Js)+'.pdf', bbox_inches='tight')


plt.figure(figsize=(6.4,4.8*2/3))
plt.axis([0,4,0,np.max(Ds[0]-ms[0]**2)*1.05])
for i in range(len(sizes)):
    plt.plot(betas,Ds[i]-ms[i]**2,lw=1.5,color=colors[i])
plt.plot(betas,Q-M**2,'k',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$q-m^2$', rotation=0, labelpad=30)
plt.savefig('img/simulation_q-m2_DeltaJ'+str(Js)+'.pdf', bbox_inches='tight')

plt.figure(figsize=(6.4,4.8*2/3))
plt.axis([0,4,0,np.max(sig)*1.05])
for i in range(len(sizes)):
    plt.plot(betas,sigmas[i],lw=1.5,color=colors[i])
plt.plot(betas,sig,'k',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\dfrac{1}{N}\left[\sigma\right]_{\mathbf{J}}$', rotation=0, labelpad=35)
plt.savefig('img/simulation_sigma_DeltaJ'+str(Js)+'.pdf', bbox_inches='tight')


plt.show()
