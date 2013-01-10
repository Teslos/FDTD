from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import matplotlib.pyplot as plt

NN = 50             # Number of points in the problem
MM = 50 
hbar = 1.054e-34    # Planck's constant
melec = 9.11e-31    # Mass of an electron
eV2J = 1.6e-19      # Energy conversion factors
J2eV = 1/eV2J

del_x = 4e-10       # The cells size
dt    = .25e-15     # Time steps
ra    = (0.5*hbar/melec)*(dt/del_x**2) # ra must be < .1
DX    = del_x * 1e9
XX    = np.linspace(0,DX*NN,NN)  # Length in nm for plotting
YY    = np.linspace(0,DX*MM,MM)  # Length in nm for plotting

NC = NN/2           # Starting position of the pulse
MC = MM/2  

V = np.zeros((NN,MM))

# ------ Input -------
nmode = input('X mode --->')
mmode = input('Y mode --->')

prl = np.zeros((NN,MM))
pim = np.zeros((NN,MM))
ptot = 0.

xx,yy = np.meshgrid(range(0,MM),range(0,NN))
print xx.shape, yy.shape
prl = np.sin(nmode*np.pi*xx/NN)*np.sin(mmode*np.pi*yy/MM)
print prl.shape
print pim.shape

ptot = np.sum(prl**2+pim**2)
pnorm = np.sqrt(np.sum(ptot))
print ptot

ptot1 = 0.
# Normalize wavefunction
prl = prl/pnorm
pim = pim/pnorm
ptot1 = np.sum(prl**2+pim**2)
print ptot1

#plt.subplot(2,2,1)
fig = plt.figure()
ax  = fig.add_subplot(2,2,1,projection='3d')
xx,yy = np.meshgrid(XX,YY)
surf = ax.plot_surface(yy,xx,prl, rstride=1, cstride = 1, cmap=cm.coolwarm,
    linewidth=0, antialiased = False)
ax.set_zlim(-.03,.03)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=.5,aspect=5)
ax1 = fig.add_subplot(2,2,2,projection='3d')
surf1 = ax1.plot_surface(xx,yy,pim, rstride=1, cstride=1, cmap=cm.coolwarm,
    linewidth=0, antialiased = False)

T = 0
n_step = 1
if n_step > 0:
    n_step = input('How many time steps --->')

# ------------- This is the core FDTD program -------------------
for T in range(1,n_step):
    for m in range(1,MM-1):
        for n in range(1,NN-1):
            prl[n,m] = prl[n,m] - ra * (-4*pim[n,m] + pim[n,m-1] + pim[n,m+1]
            + pim[n-1,m] + pim[n+1,m]) + (dt/hbar)*V[n,m]*pim[n,m]
        
    for m in range(1,MM-1):
        for n in range(1,NN-1):
            pim[n,m] = pim[n,m] + ra * (-4*prl[n,m] + prl[n,m-1] + prl[n,m+1]
            + prl[n-1,m] + prl[n+1,m]) - (dt/hbar)*V[n,m]*prl[n,m]
     
# Check normalization
ptot = 0.
psi = prl + 1j*pim
ptot = np.sum(prl**2 + pim**2)
print ptot

# Calculate the expected values
ke = 0. + 1j*0.
for m in range(2, MM-2):
    for n in range(2, NN-2):
        lap_p = psi[n+1,m] - 4 * psi[n,m] + psi[n-1,m] 
        + psi[n,m-1] + psi[n,m+1]

        ke = ke + lap_p * psi[n,m].conjugate()

KE = -J2eV*((hbar/del_x)**2/(2*melec))*ke.real
print "Kinetic energy: %2e" %  KE

surf = ax.plot_surface(yy,xx,prl, rstride=1, cstride = 1, cmap=cm.coolwarm,
    linewidth=0, antialiased = False)
ax.set_zlim(-.03,.03)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=.5,aspect=5)
surf1 = ax1.plot_surface(xx,yy,pim, rstride=1, cstride=1, cmap=cm.coolwarm,
    linewidth=0, antialiased = False)
ax1.set_xlabel('Y (nm)')
ax1.set_ylabel('X (nm)')
ax1.text(DX*MM,DX*NN/2,.035, " KE = {0} meV".format(1e3*KE))
plt.show()

