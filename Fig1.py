import numpy as np
from matplotlib import pyplot as plt
from kinetic_ising import ising

plt.rc('text', usetex=True)
font = {'size': 18, 'family':'serif', 'serif': ['latin modern roman']}
plt.rc('font', **font)
plt.rc('legend', **{'fontsize': 20})
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


def bool2int(x):  # Transform bool array into positive integer
    y = 0
    for i, j in enumerate(np.array(x)[::-1]):
        y += j * 2**i
    return y


H0=0.
Js=0.6
beta=1.6281474038842134


N=64
I=ising(N)
J= 1/N + Js*np.random.randn(N,N) / np.sqrt(N)
H=np.zeros(N)
#D=J


I.H=H.copy()
I.J = J.copy()
I.beta=beta 
I.s=np.zeros(N)
T=4000
T=400*N
T=20001
T0=100000
S=np.zeros((N,T))
A=np.zeros(T)
for t in range(T0):
    I.GlauberStep()
for t in range(T):
    I.GlauberStep()
    S[:,t]=I.s
    A[t] = bool2int((I.s+1)//2)
    

plt.figure(figsize=(6,3))
plt.title(r'$\Delta J = ' +str(Js)+r',\beta='+str(beta)[0:5]+'$')
plt.imshow(S,cmap='inferno_r',aspect='auto',interpolation='None',origin='lower')

plt.ylabel(r'Neuron', labelpad=5)
plt.xlabel(r'$t$')
fig1 = plt.gcf()
plt.draw()
fig1.savefig('img/rasterplot_J_' +str(Js)+r'_beta_'+str(beta)[0:5]+'.pdf', bbox_inches='tight',dpi=200)

plt.show()

