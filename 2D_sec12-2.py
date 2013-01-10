from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import matplotlib.pyplot as plt

NN = 50              # Number of points in the space
MM = 50
hbar = 1.054e-34     # Planck's constant
melec = 9.1e-31      # Mass of an electron
eV2J  = 1.6e-19      # Energy conversion factors
J2eV  = 1/eV2J

del_x = 2.e-10       # The cells size
dt    = .1e-15       # Time steps
ra    = (0.5*hbar/melec)*(dt/del_x**2) # ra must be < .1
print "Stability criteria ra: %2e" % ra
DX    = del_x * 1e9
XX    = np.linspace(0,DX*NN,NN)

# ---- Parameters for the FFT
Ntime = 2**16         # Size of the FFT buffer
Ptime = np.zeros(Ntime)
Pwin  = np.zeros(Ntime)
win   = np.zeros(Ntime)
PF    = np.zeros(Ntime)
FF    = np.zeros(Ntime)
del_F = 1./(Ntime*dt)
del_E = J2eV*2*np.pi*hbar*del_F
FF    = np.arange(0,del_E*Ntime,del_E)
delT  = 1e12*dt

PS    = np.zeros(Ntime)
PS    = np.arange(0,delT*Ntime,delT)

NC = NN/2                  # Starting position of the pulse
MC = NN/2 

# ---- Specify the potential -------------
V = np.zeros((NN,NN))

# Slanted potential
#for m in range(1,NN):
#    for n in range(10,NN):
    #    V(n,m) = k0 *( (NC-n)**2 + 0.75*(MC-m)**2 )*del_x**2
    #    V(n,m) = eV2J *(.1/NN)*n

#for n in range(1,9):
    # V(n,m) = eV2J *.1

# HO potential
# k0 = 0.005
# for m in range(1,NN):
#     for n in range(1,NN):
#         V[n,m] = k0*((NC-n)**2 + (MC-m)**2)*del_x**2

fig = plt.figure()
# print XX.shape
xx,yy = np.meshgrid(XX,XX)
# ax = fig.add_subplot(2,2,1,projection='3d')
# print xx.shape
# print yy.shape
# surf = ax.plot_surface(xx,yy,J2eV*V)    
# ax.set_xlabel('nm')
# ax.set_ylabel('nm')
# ax.set_zlabel('V (eV)')

# ------ Test function ---------------------------


win2D = np.zeros((NN,NN))
prl   = np.zeros((NN,NN))
pim = np.zeros((NN,NN))
sigma = 6.
ptot = 0.

 
nn,mm = np.meshgrid(range(0,NN),range(0,NN))
        # The test function is smoothed with a 2D Hanning window
win2D = 0.25 * (1. - np.cos(2*np.pi*mm/NN)) * (1. - np.cos(2*np.pi*nn/NN))

        # Single point
dist = np.sqrt((NC-nn)**2 + (MC-mm)**2)
prl  = win2D * np.exp(-(dist/sigma)**2)

        # Double point
dist1 = np.sqrt( (MC-mm)**2 + (NC+10-nn)**2 )
dist2 = np.sqrt( (MC-mm)**2 + (NC-10-nn)**2 )

# DEBUG
# prl = np.sin(4*np.pi*nn/NN)*np.sin(1*np.pi*mm/MM)
ptot  = np.sum(prl**2 + pim**2)
         
pnorm = np.sqrt( ptot )
msource = MC             # Position were the time-domain 
nsource = NC - 10        # data is monitored

# Normalize and check
ptot = 0.
prl = prl / pnorm
pim = pim / pnorm
ptot = np.sum( prl**2 + pim**2 )
print ptot

ax1 = fig.add_subplot(3,2,3, projection='3d')
surf1 = ax1.plot_surface(xx,yy,prl)
ax1.set_xlim3d(0, DX*NN)
ax1.set_ylim3d(0, DX*NN)
ax1.set_zlim3d(-0.03, 0.03)

ax2 = fig.add_subplot(3,2,4, projection='3d')
surf2 = ax2.plot_surface(xx,yy,pim)
ax2.set_xlim3d(0, DX*NN)
ax2.set_ylim3d(0, DX*NN)
ax2.set_zlim3d(-0.03, 0.03)

T = 0
n_step = 1
if n_step > 0:
    n_step = input('How many time steps --->')

# ------ This is the core FDTD program ---------------
for T in range(1, n_step):
    for m in range(1,MM-1):
        for n in range(1,NN-1):
            prl[n,m] = prl[n,m] - ra * (-4*pim[n,m] + pim[n,m-1] + pim[n,m+1]
            + pim[n-1,m] + pim[n+1,m]) + (dt/hbar)*V[n,m]*pim[n,m]
        
    for m in range(1,MM-1):
        for n in range(1,NN-1):
            pim[n,m] = pim[n,m] + ra * (-4*prl[n,m] + prl[n,m-1] + prl[n,m+1]
            + prl[n-1,m] + prl[n+1,m]) - (dt/hbar)*V[n,m]*prl[n,m]
        
    Ptime[T] = prl[msource,nsource]-1j*pim[msource,nsource]

# Check normalization
ptot = 0.
psi  = prl + 1j * pim
ptot = np.sum(prl**2 + pim**2)
print ptot

ax3 = fig.add_subplot(3,2,1, projection='3d')
surf3 = ax3.plot_surface(xx,yy,prl,rstride=1, cstride = 1, cmap=cm.coolwarm,
    linewidth=0, antialiased = False)
ax3.zaxis.set_major_locator(LinearLocator(5))
ax3.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax3.set_xlim3d(0, DX*NN)
ax3.set_ylim3d(0, DX*NN)
ax3.set_zlim3d(-.3, .3)
ax3.set_xlabel('nm')
ax3.set_xlabel('nm')
ax3.set_zlabel('Real (Psi)')
ax3.set_title('Se2d-find-E')

ax4 = fig.add_subplot(3,2,2, projection='3d')
surf4 = ax4.plot_surface(xx,yy,pim,rstride=1, cstride = 1, cmap=cm.coolwarm,
    linewidth=0, antialiased = False)
ax4.zaxis.set_major_locator(LinearLocator(5))
ax4.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax4.set_xlim3d(0, DX*NN)
ax4.set_ylim3d(0, DX*NN)
ax4.set_zlim3d(-.3, .3)
ax4.set_xlabel('nm')
ax4.set_xlabel('nm')
ax4.set_zlabel('Imag (Psi)')
ax4.set_title('Se2d-find-E')

# Plot the time domain data and FFT.

# ---- Create the Hanning window for the time-domain data
win = np.zeros(T)
n = np.arange(1,T)
win = 0.5*(1-np.cos(2*np.pi*n/T))
m = np.arange(1,T-1)
Pwin[m] = win[m] * Ptime[m]

ax5 = fig.add_subplot(3,2,5)
ax5.plot(PS,Pwin.real,'k')
ax5.set_xlabel('Time (ps)')
ax5.set_ylabel('Real Psi')

# ----- Take the FFT of the windowed time-domain data --------
PF = (1/np.sqrt(Ntime))*np.abs(np.fft.fft(Pwin))
Pmax = np.max(PF)

ax6 = fig.add_subplot(3,2,6)
ax6.plot(1e3*FF,PF,'k')
ax6.set_xlabel('E (meV)')
ax6.set_ylabel('FFT (Psi)')
ax6.set_xlim(0,50)
ax6.set_ylim(0,1.1*Pmax)
plt.savefig('Se2d-find-E-Inf-well.png')
plt.show()
