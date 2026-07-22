import numpy as np
import matplotlib.pyplot as plt


def f(z,a,H0,T0,T1):
    return (T1+T0)/2 + 0.5*(T0-T1)*np.tanh(a*(z+0.5*H0)/H0)


H0 = 100
T0 = 13.25
T1 = 8.
a = 20.
z = np.linspace(0,-H0, 100)
Tini = f(z,a,H0,T0,T1)

fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
ax.plot(Tini,z)
ax.vlines(T1,z[-1],z[0])
ax.vlines(T0,z[-1],z[0])


def gaus(x,a,b):
    return a*np.exp(-0.5*x**2/b**2)

X = np.linspace(-5000,5000,200)
ZB = 50*gaus(X/5000,1,0.1)
fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
ax.plot(X,ZB)
ax.set_ylim([0,100])
plt.show()
