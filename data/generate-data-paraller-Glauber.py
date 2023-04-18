"""
GPLv3 2020 Miguel Aguilera

This code allows to run simulations of the kinetic Ising model, with asymmetric weights and parallel updates.
"""
from kinetic_ising import ising
import numpy as np
#from matplotlib import pyplot as plt
from sys import argv

if len(argv) < 4:
    print("Usage: " + argv[0] + " <network size> <repetitions> <gamma1> <gamma2> <b> <trial>")
    exit(1)

size=int(argv[1])
R=int(argv[2])
gamma1=float(argv[3])
gamma2=float(argv[4])
ind=int(argv[5])
ib = ind % 101
trial = ind // 101

print('Starting ib '+str(ib)+', trial '+str(trial))
random_s0 = False
stationary = True
I = ising(size)

B = 101
T = 2**7
#T = 2**8
beta0=1
betas =  np.linspace(0,4,B) 


print(ib)
beta_ref = round(betas[ib],3)
beta = beta_ref * beta0
print(beta_ref,str(ib)+'/'+str(len(betas)),gamma1,gamma2,size)

I.H = beta * gamma1 * (np.random.rand(size)*2-1)
if np.mean(I.H)<0:
    I.H*=-1
I.J = beta * (1/size + gamma2 * np.random.randn(size, size) / np.sqrt(size) )

m_exp = np.zeros((size))
m_exp_prev = np.zeros((size))
C_exp = np.zeros((size, size))
D_exp = np.zeros((size, size))
sigma=0


s0 = np.ones(size)

print('generate data')
for rep in range(R):
    I.H = beta * gamma1 * (np.random.rand(size)*2-1)
    if np.mean(I.H)<0:
        I.H*=-1
    I.J = beta * (1/size + gamma2 * np.random.randn(size, size) / np.sqrt(size) )

    I.s = s0.copy()
    for t in range(T-1):
        I.ParallelUpdate()
    h = I.H + np.dot(I.J, I.s)
    m_exp += np.tanh(h) / R
    C_exp += np.einsum('i,j->ij', np.tanh(h), np.tanh(h),  optimize=True) / R
    C_exp[range(size),range(size)] += (1 - np.tanh(h)**2) / R
    d = np.einsum('i,j->ij', np.tanh(h), I.s, optimize=True)
    D_exp += d / R
    sigma += np.sum(d*(I.J-I.J.T)) / R

filename = 'parts/data-gamma1-' + str(gamma1) +'-gamma2-' + str(gamma2) + '-s-' +  str(size) + '-R-' + str(R) + '-beta-' + str(beta_ref) + '-trial-' + str(trial) + '.npz'
np.savez_compressed(filename, C=np.mean(C_exp), Cdiag=np.mean(np.diag(C_exp)), m=np.mean(m_exp), D=np.mean(D_exp),Ddiag=np.mean(np.diag(D_exp)),sigma=sigma)
