import numpy as np
import matplotlib.pyplot as plt

NN = 50             # Number of points in the problem
hbar = 1.054e-34    # Planck's constant
m0 = 9.1e-31        # Free space mass of an electron
meff = 1.0          # Effective mass: Si is 1.08, Ge is 0.067, GaAs is 0.55
melec = meff*m0     # Mass of an electron

eV2J = 1.6e-19      # Energy conversion factors
J2eV = 1/eV2J

del_x = .2e-9       # The cells size
dt = 1e-16          # Time steps
ra = (0.5*hbar/melec)*(dt/del_x**2) # ra must be < .1
DX = del_x * 1e9    # Cell size in nm
XX = np.linspace(0,DX*NN,NN)  # Length in nm for plotting

#-------Parameters for the FFT------------------------------------
Ntime = 2**16       # Number of points in the FFT buffer
Ptime = np.zeros(Ntime)
Pwin  = np.zeros(Ntime)
PF    = np.zeros(Ntime)
FF    = np.zeros(Ntime)
del_F = 1./(Ntime * dt)
del_E = J2eV*2*np.pi*hbar*del_F
FF    = np.arange(0,del_E*(Ntime),del_E)
delT  = 1e12*dt
PS    = np.zeros(Ntime)
PS    = np.arange(0,delT*Ntime,delT)

#------Specify the potential--------------------------------------

# Slanted potential
V = np.zeros(NN)
#V = eV2J*(.005)
n = np.arange(20,40)
#V[n] = eV2J*(.0005)*(np.abs(30-n))

# Plot V
plt.subplot(3,2,2)
plt.plot(XX, J2eV * V)
plt.ylabel('V (eV)')
plt.xlabel('nm')
Umax = max(J2eV * V)
Umax = 1
plt.title('Sel12-1')
plt.axis( [0, DX*NN, 0, 1.1 * Umax])

#-----Initialize the test function-------------------------------
sigma = 3.        # pulse width
nc = 25.          # Time domain data monitored at this point
prl = np.zeros(NN)
pim = np.zeros(NN)
ptot = 0.
for n in range(1,NN-1):
    prl[n] = np.exp(-1.*((nc-n)/sigma)**2)     # First test function
#    prl[n] = np.exp(-1.*((12.5-n)/sigma)**2)-np.exp(-1.*((37.5-n)/sigma)**2)
ptot = np.sum(prl**2+pim**2)
pnorm = np.sqrt(ptot)
# Normalize and check
ptot = 0.
prl = prl/pnorm
pim = pim/pnorm

ptot = np.sum(prl**2+pim**2)
print ptot

# Plot psi
#print XX.size
#print prl.size
plt.subplot(3,2,1)
plt.plot(XX,prl,'k')
plt.hold(True)
#plt.plot(XX,pim,'-.r')
plt.plot(XX,J2eV*V,'--k')
plt.hold(False)
plt.axis([0, DX*NN, -.4,.4])
plt.xlabel('nm')
plt.title('Sel12-1')
plt.savefig('se.png')

T = 0
n_step = 1
if n_step > 0:
    n_step = input("How many time steps--->")

#-----------This is the main FDTD loop------------------
for m in range(1,n_step):
    T = T + 1
    for n in range(1,NN-1):
        prl[n] = prl[n] - ra*(pim[n-1] - 2*pim[n] + pim[n+1]) 
        + (dt/hbar)*V[n]*pim[n]
      
    for n in range(1,NN-1):
        pim[n] = pim[n] + ra*(prl[n-1] - 2*prl[n] + prl[n+1])
        - (dt/hbar)*V[n]*prl[n]

    Ptime[T] = prl[nc]-1j*pim[nc]

# Check normalization
ptot1 = 0.
ptot1 = np.sum(prl**2 + pim**2)
ptot1

# Plot psi
plt.subplot(3,2,1)
plt.plot(XX,prl,'k')
plt.hold(True)
plt.plot(XX,pim,'--k')
#plt.plot(XX,J2eV*V,'--k')
plt.hold(False)
plt.axis([0, DX*NN, -.25,.5])
plt.text(7,.3,"{0} ps".format(T*dt*1e12))
plt.xlabel('nm')
plt.title('Sel2-1')

# FFT with windowing
# --------- Create the Hanning window for T points-----------------------
win = np.zeros(T)
#win = 1.
win = 0.5*(1-np.cos(2*np.pi*n/T))

# -------- Window the saved time-domain data ----------------------------
Pwin = win * Ptime

iin = np.zeros( T )
# --------- Plot the time-domain psi at one point -----------------------
#print PS.size
#print Pwin.size
plt.subplot(3,2,3)
plt.plot(PS, Pwin.real, 'k')
plt.axis([0, 1e12*dt*T, -.3, .3])
plt.xlabel('ps')
plt.ylabel('Windowed Psi')


# Take the FFT of the windowed time-domain data to display the eigenenergies
PF = (1/np.sqrt(Ntime))*abs(np.fft.fft(Pwin))
Pmax = np.max(PF)
#print FF.size
#print PF.size
plt.subplot(3,2,4)
plt.plot(1e3*FF, PF, 'k')
plt.axis([0,40,0, 1.1*Pmax])
plt.xlabel('E (meV)')
plt.ylabel('| fft (Psi) |')
plt.grid(True)
plt.show()


