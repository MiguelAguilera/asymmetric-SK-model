import numpy as np
from mpmath import polylog
from polylog import fermi_poly2
    
    
# Integral of tanh with Gaussian input
def Gm(x, g, sqrtD): return 1 / np.sqrt(2 * np.pi) * \
    np.exp(-0.5 * x**2) * np.log(np.cosh(g + x * sqrtD))
    
def Gm0(x, g, sqrtD): return 1 / np.sqrt(2 * np.pi) * \
    np.exp(-0.5 * x**2) * (np.abs(g + x * sqrtD))

# Integral of tanh with Gaussian input
def zGm(x, g, sqrtD): return 1 / np.sqrt(2 * np.pi) * \
    np.exp(-0.5 * x**2) * x*np.log(np.cosh(g + x * sqrtD))
    
# tanh (or integral of 1-tanh**2) with Gaussian input
def Gtanh(x, g, sqrtD): return 1 / np.sqrt(2 * np.pi) * \
    np.exp(-0.5 * x**2) * np.tanh(g + x * sqrtD)
    
# tanh (or integral of 1-tanh**2) with Gaussian input
def Gtanh0(x, g, sqrtD): return 1 / np.sqrt(2 * np.pi) * \
    np.exp(-0.5 * x**2) * np.sign(g + x * sqrtD)
    
    # tanh (or integral of 1-tanh**2) with Gaussian input
def Gsign(x, g, sqrtD): return 1 / np.sqrt(2 * np.pi) * \
    np.exp(-0.5 * x**2) * np.sign(g + x * sqrtD)

# 1-tanh**2 with Gaussian input
def Gtanh2(x, g, sqrtD): return 1 / np.sqrt(2 * np.pi) * \
    np.exp(-0.5 * x**2) * (1-np.tanh(g + x * sqrtD)**2)
    
# 1-tanh**2 with Gaussian input
def Gtanh20(x, g, sqrtD): return 1 / np.sqrt(2 * np.pi) * \
    np.exp(-0.5 * x**2) * (1-np.sign(g + x * sqrtD)**2)
#return 1 / np.sqrt(2 * np.pi) * \
#    np.exp(-0.5 * x**2) * (1-np.sign(g + x * sqrtD)**2)
    
# 1-tanh**2 with Gaussian input
def Phi1(x, g, sqrtD): return 1 / np.sqrt(2 * np.pi) * \
    np.exp(-0.5 * x**2) * (-np.tanh(g + x * sqrtD)*(g + x * sqrtD)+np.log(2*np.cosh(g + x * sqrtD)))
    
# phi with Gaussian input
def Phi(x, g, sqrtD): 
#    L=np.zeros_like(x)
#    for i in range(len(x)):
#        L[i] = polylog(2,-np.exp(2*(g + x[i] * sqrtD)))
    L= -fermi_poly2(2*(g + x * sqrtD))
    return 1 / np.sqrt(2 * np.pi) * \
    np.exp(-0.5 * x**2) *  (-1) * ((g + x * sqrtD) * np.log(1+np.exp(2*(g + x * sqrtD))) + L)

# 1D Integrator
def Gaussian_integral(F,g,sqrtD,Nint=10000,xmax=5,x0=0):
    if sqrtD==0:
        return np.sqrt(2 * np.pi) * F(0, g, 0)
    else:
        x1 = np.linspace(-1, 1, Nint)*xmax + x0
        f = F(x1, g, sqrtD)
        return np.sum(f)*(x1[1]-x1[0])


# Integral of tanh(x)*tanh(y) with Gaussian input
def Gmm(x,y,g,h1,h2,sqrtD,q):
    if 1-q**2>0:
        X, Y = np.meshgrid(x, y)
        return Gmm2(X,Y,g,h1,h2,sqrtD,q)
    else:
        return Gmm1(x,g,h1,h2,sqrtD)
        
# Integral of tanh(x)*tanh(y) with 2D inputs
def Gmm2(x,y,g,h1,h2,sqrtD,q):
    a=h1+sqrtD*x
    b=h2+sqrtD*(x*q + np.sqrt(1-q**2)*y)
    R=np.zeros_like(x)
    R[a==b] = 1 / (2 * np.pi) * np.exp(-0.5 * x[a==b]**2 -0.5 * y[a==b]**2) * (g-np.tanh(g+a[a==b]))
    R[a!=b] =  1 / (2 * np.pi) * np.exp(-0.5 * x[a!=b]**2 -0.5 * y[a!=b]**2) * \
    (g+((np.exp(2*b[a!=b])+np.exp(2*a[a!=b]))*(np.log(1+np.exp(2*g+2*a[a!=b]))-np.log(1+np.exp(2*g+2*b[a!=b]))))/(np.exp(2*b[a!=b])-np.exp(2*a[a!=b])))
    return R

# Integral of tanh(x) tanh(y) with 1D inputs
def Gmm1(x,g,h1,h2,sqrtD):
    a=h1+sqrtD*x
    b=h2+sqrtD*x
    if h1==h2:
        R= 1 / np.sqrt(2 * np.pi) * np.exp(-0.5 * x**2) * (g-np.tanh(g+a))
    else:
        R=  1 / np.sqrt(2 * np.pi) * np.exp(-0.5 * x**2) * \
    (g+(np.exp(2*b)+np.exp(2*a))*((np.exp(2*b)+np.exp(2*a))*(np.log(1+np.exp(2*g+2*a))-np.log(1+np.exp(2*g+2*b))))/(np.exp(4*b)+np.exp(4*a)))
    return R
    
## Integral of tanh(x)*tanh(y) with Gaussian input
# Integral of tanh(x)*tanh(y) with Gaussian input
def Gmm0(x,y,g,h1,h2,sqrtD,q):
    if 1-q**2>0:
        X, Y = np.meshgrid(x, y)
        return Gmm02(X,Y,g,h1,h2,sqrtD,q)
    else:
        return Gmm01(x,g,h1,h2,sqrtD)
        
# Integral of tanh(x)*tanh(y) with 2D inputs
def Gmm02(x,y,g,h1,h2,sqrtD,q):
    a=h1+sqrtD*x
    b=h2+sqrtD*(x*q + np.sqrt(1-q**2)*y)
    R=np.zeros_like(x)
    R =  1 / (2 * np.pi) * np.exp(-0.5 * x**2 -0.5 * y**2) * \
    (g+np.sign(b-a)*(np.abs(g+a) - np.abs(g+b)))
    return R

# Integral of tanh(x) tanh(y) with 1D inputs
def Gmm01(x,g,h1,h2,sqrtD):
    a=h1+sqrtD*x
    b=h2+sqrtD*x
    R =  1 / np.sqrt(2 * np.pi) * np.exp(-0.5 * x**2) * \
    (g+np.sign(b-a)*(np.abs(g+a) - np.abs(g+b)))
    return R
    
        
    
# Integral of tanh(x)*tanh(y) with Gaussian input
def Gtanhtanh(x,y,g,h1,h2,sqrtD,q):
    if 1-q**2>0:
        X, Y = np.meshgrid(x, y)
        return Gtanhtanh2(X,Y,g,h1,h2,sqrtD,q)
    else:
        return Gtanhtanh1(x,g,h1,h2,sqrtD)
        
# Integral of tanh(x)*tanh(y) with 2D inputs
def Gtanhtanh2(x,y,g,h1,h2,sqrtD,q):
    a=h1+sqrtD*x
    b=h2+sqrtD*(x*q + np.sqrt(1-q**2)*y)
    R=np.zeros_like(x)
    R[a==b] = 1 / (2 * np.pi) * np.exp(-0.5 * x[a==b]**2 -0.5 * y[a==b]**2) * np.tanh(g+a[a==b])*np.tanh(g+b[a==b])
    R[a!=b] =  1 / (2 * np.pi) * np.exp(-0.5 * x[a!=b]**2 -0.5 * y[a!=b]**2) * \
    np.tanh(g+a[a!=b])*np.tanh(g+b[a!=b])
    return R

# Integral of tanh(x) tanh(y) with 1D inputs
def Gtanhtanh1(x,g,h1,h2,sqrtD):
    a=h1+sqrtD*x
    b=h2+sqrtD*x
    R= 1 / np.sqrt(2 * np.pi) * np.exp(-0.5 * x**2) * np.tanh(g+a)*np.tanh(g+b)
    return R
     
# 2D Integrator
def Gaussian_integral2D(F,g, h1,h2, sqrtD, q, N=200):
    x = np.linspace(-1, 1, N) * 5
    y = np.linspace(-1, 1, N) * 5
    if 1-q**2>0:    # If correlation is not 1
        return np.sum(F(x, y, g, h1,h2, sqrtD, q)) * (x[1] - x[0]) * (y[1] - y[0])
    else:            # If correlation is 1
        return np.sum(F(x, y, g, h1,h2, sqrtD, q)) * (x[1] - x[0]) 
