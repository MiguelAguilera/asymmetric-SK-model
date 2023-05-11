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

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 15})
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


H0=0.0
Js=0.5
Js=1.0
J0=1
R=400000
Tr=1

#R=10000
#Tr=6

B=101
#betas = 1 + np.linspace(-1, 1, B) * 0.3
betas = np.linspace(0, 4, B)

#T=1000
#M = np.zeros(B)
#sig = np.zeros(B)
#K = np.zeros(B)
#Q = np.zeros(B)
#DF = np.zeros(B)


#cmap = cm.get_cmap('plasma_r')
#colors=[]
#for i in range(4):
#	colors+=[cmap((i)/(4-1))]

#alpha=0.01
#1
#L=20
#T1=20
##L=100
##T1=100
##L=400
##T1=400
#L=1000
#T1=1000
#for ib in range(30,len(betas)):

#    print(ib)
#    steps = np.arange(T)
#    beta=betas[ib]
#    beta=2

#    m=0.5
#    m_prev=0.5
#    q=m**2
#    k=0
#    for t in steps:
#        m_prev2=m_prev
#        m_prev=m
#        q_prev=q

#        m = Gaussian_integral(Gtanh,beta*J0*m_prev,beta*Js)
#        
#        
#    q0=m**2
#    for t in range(T1):
#        q0 =  Gaussian_integral2D(Gtanhtanh, 0, beta*J0*m, beta*J0*m, beta*Js,q0)
#    
#    q = np.ones(L+1)*q0
#    q[0]=1
#    for rep in range(T1):
#        print(rep,T1)
#        q_prev=q.copy()
#        q_ = np.zeros(L)
#        q = np.ones(L+1)*q0
#        q[0]=1
#        q_[0]=1
#        for l in range(1,L):
#            q_[l] = (1-alpha) * q_[l-1] + alpha * Gaussian_integral2D(Gtanhtanh, 0, beta*J0*m, beta*J0*m, beta*Js, q_prev[l])
#            q[l] = (1-alpha) * q_prev[l+1] + alpha * q_[l]
#    print(beta, m)

##    q2 = m**2
##    q2_ = m**2
##    for r in range(100):
##        q2_ = 1 - alpha *2  + 2*alpha * Gaussian_integral2D(Gtanhtanh, 0, beta*J0*m, beta*J0*m, beta*Js, q2)
##    q2 = np.ones(L)
##    for r in range(100):
##        q2[0]=1
##        q2 = (1-alpha)*q2 + alpha*q2_
##        
##    print(q)
##    print(q_)
##    print(q2,q2_)
##    plt.show()

#    k = Gaussian_integral(Gtanh2,beta*J0*m,beta*Js)
##    S[ib] = Gaussian_integral(Phi1,beta*H0+beta*J0*m,beta*Js)
#    sig[ib] = (beta*Js)**2*(1-q[1])*k
#    M[ib] =m
#    Q[ib] =q[-1]
#    DF[ib] = (beta*Js)**2*(q[1]-q[2])/alpha*k
#    print(beta,m,k,q[1],q0,sig[ib],DF[ib])
#    
#    fig,ax = plt.subplots(figsize=(6.4*1.3,4.8))
#    ax.set_xlabel(r"$d$", fontsize = 18)
#    # set y-axis label
#    ax.set_ylabel(r"$q$",rotation=0, fontsize = 18, labelpad=15)
##    ax2=ax.twinx()
#    ax.plot(np.arange(L)*alpha,q_,color=colors[-2],label=r'$q_{t+d,t}^{*,1}$')
#    ax.plot(np.arange(L+1)*alpha,q,color=colors[-1],label=r'$q_{t+d,t}$')
#    ax.plot(np.arange(L+1)*alpha,q*0+q0,'--',color=colors[-3],label=r'$q_{\infty,t}$')
#    ax.plot(np.arange(L+1)*alpha,q,color=colors[-1])
##    ax2.set_ylabel(r"$q^{1,*}$",color=colors[-2],rotation=0, labelpad=20)
#    plt.legend(fontsize = 18)
#    plt.axis([0,L*alpha,0.55,1.05])
#    plt.savefig('img/q_decay.pdf', bbox_inches='tight')
#    plt.show()
#    


#np.savez('analytical_a_DeltaJ'+str(Js)+'.npz',M=M,Q=Q, sig=sig,DF=DF)

data = np.load('data/analytical_a_DeltaJ'+str(Js)+'.npz')
M=data['M']
Q=data['Q']
sig=data['sig']


sizes=[32,64,128,256,512,1024]
#sizes=[32,64,128,256,512]
#sizes=[32,64,128,256]

S=len(sizes)
cmap = cm.get_cmap('plasma_r')
colors=[]
for i in range(S):
	colors+=[cmap((i)/(S-1))]
	
	
ms=[]
Cs=[]
Ds=[]
sigmas=[]

for size in sizes:

    filename = 'data/data_a_-gamma1-'+str(H0)+'-gamma2-'+str(Js)+'-s-'+str(size)+'-R-'+str(R)+'-Tr-'+str(Tr)+'.npz'
#    filename = 'data-gamma1-'+str(H0)+'-gamma2-'+str(Js)+'-s-'+str(size)+'-R-'+str(R)+'.npz'
    data=np.load(filename)
    ms+=[data['m']]
    Cs+=[data['C']]
    Ds+=[data['D']]
    sigmas+=[data['sigma']]

plt.figure(figsize=(6.4,4.8*2/3))
if Js==0.5:
    plt.axis([0,4,0,1])
else:
    plt.axis([0,4,0,np.max(ms)*1.1])
#plt.axis([0,4,0,np.max(ms)*2])
for i in range(len(sizes)):
    plt.plot(betas,ms[i],lw=1.5,color=colors[i],label=r'$N='+str(sizes[i])+'$')
plt.plot(betas,M,'k',lw=2,label=r'$N\to\infty$')
plt.xlabel(r'$\beta$')
plt.legend(handlelength=1.5)
plt.ylabel(r'$m$', rotation=0, labelpad=15)
plt.savefig('img/simulation_a_m_DeltaJ'+str(Js)+'.pdf', bbox_inches='tight')


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
plt.savefig('img/simulation_a_q_DeltaJ'+str(Js)+'.pdf', bbox_inches='tight')


plt.figure(figsize=(6.4,4.8*2/3))
plt.axis([0,4,0,np.max(Ds[0]-ms[0]**2)*1.05])
for i in range(len(sizes)):
    plt.plot(betas,Ds[i]-ms[i]**2,lw=1.5,color=colors[i])
plt.plot(betas,Q-M**2,'k',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$q-m^2$', rotation=0, labelpad=30)
plt.savefig('img/simulation_a_q-m2_DeltaJ'+str(Js)+'.pdf', bbox_inches='tight')

plt.figure(figsize=(6.4,4.8*2/3))
plt.axis([0,4,0,np.max(sig)*1.05])
for i in range(len(sizes)):
    plt.plot(betas,sigmas[i],lw=1.5,color=colors[i])
plt.plot(betas,sig,'k',lw=2)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\dfrac{1}{N}\left[\dfrac{d\sigma}{dt}\right]_{\mathbf{J},\boldsymbol{\tau}}$', rotation=0, labelpad=40)
plt.savefig('img/simulation_a_sigma_DeltaJ'+str(Js)+'.pdf', bbox_inches='tight')


#plt.figure(figsize=(6.4,4.8*2/3))
#plt.axis([0,4,0,np.max(sig)*1.05])
##for i in range(len(sizes)):
##    plt.plot(betas,sigmas[i]*sizes[i],lw=1.5,color=colors[i])
#plt.plot(betas,DF,'k',lw=2)
#plt.xlabel(r'$\beta$')
#plt.ylabel(r'$\beta \Delta F,\frac{d\sigma}{dt}$', rotation=0, labelpad=30)
#plt.savefig('img/simulation_a_sigma.pdf', bbox_inches='tight')

#plt.figure(figsize=(6.4,4.8*2/3))
#plt.axis([0,4,0,np.max(sig)*1.05])
##for i in range(len(sizes)):
##    plt.plot(betas,sigmas[i]*sizes[i],lw=1.5,color=colors[i])
#plt.plot(betas,DF-sig,'k',lw=2)
#plt.xlabel(r'$\beta$')
#plt.ylabel(r'$\frac{d\sigma}{dt}$', rotation=0, labelpad=20)
#plt.savefig('img/simulation_a_sigma.pdf', bbox_inches='tight')
plt.show()
