import numpy as np
from mpmath import polylog
from SK import *

size = 1000     # Network size
J0 = 1.0        # Average value of couplings
Js = 0.8

grid_num = 100
T = 10000                    # Number of simulation time steps
beta_list = np.linspace(0, 4, grid_num)     # Inverse temperature
H0_list = np.linspace(0, 2, grid_num)

M = np.zeros((grid_num,grid_num))     # magnetization
Q = np.zeros((grid_num,grid_num))     # self-correlation
sig = np.zeros((grid_num,grid_num))   # entropy production
sigb = np.zeros((grid_num,grid_num))    # EP over beta
S = np.zeros((grid_num,grid_num))   # entropy rate
S_r = np.zeros((grid_num,grid_num))   # reverse entropy rate

steps = np.arange(T)

def function(H0, beta):
    if beta==0:
        m = 0
        q = 0
        S = np.log(2)
        S_r = np.log(2)
        sig = 0
        sigb = 0
    else:
        q=0
        m=0.1
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
            if update_error < 1E-10:
                break
                
        if H0==0:
            k = Gaussian_integral(Gtanh2,beta*J0*m,beta*Js)
            p = Gaussian_integral(Phi1,beta*H0+beta*J0*m,beta*Js)
            S = p
            sig = (beta*Js)**2*(1-q)*k
        elif Js==0:
            k = (1-np.tanh(beta*H0 + beta*J0*m)**2) 
            p = ( -((beta*H0 + beta*J0*m) * np.log(1+np.exp(2*beta*H0 + 2*beta*J0*m)) + polylog(2,-np.exp(2*beta*H0 + 2*beta*J0*m))) + ((-beta*H0 + beta*J0*m)* np.log(1+np.exp(-2*beta*H0 + 2*beta*J0*m)) + polylog(2,-np.exp(-2*beta*H0 + 2*beta*J0*m))))/(2*beta*H0)
            S = p
            sig = 0
        else:
            k= (Gaussian_integral(Gtanh,beta*H0 + beta*J0*m,beta*Js)-Gaussian_integral(Gtanh,-beta*H0 + beta*J0*m,beta*Js))/(2*beta*H0)
            p = (Gaussian_integral(Phi,beta*H0+beta*J0*m,beta*Js)-Gaussian_integral(Phi, -beta*H0+beta*J0*m,beta*Js))/(2*beta*H0)
            S = p
            sig = (beta*Js)**2*(1-q)*k

        S_r = S + sig

        sigb = sig / beta

    print(m, q, sig, S, S_r, sigb)

    return m, q, sig, S, S_r, sigb

if __name__ == '__main__':

    for ii in range(grid_num):
        for jj in range(grid_num):
            M[ii][jj],Q[ii][jj],sig[ii][jj],S[ii][jj],S_r[ii][jj],sigb[ii][jj] = function(H0_list[ii], beta_list[jj])

    np.save('data-heatmap-Js-'+str(Js)+'.npy', (M, Q, sig, S, S_r, sigb))
