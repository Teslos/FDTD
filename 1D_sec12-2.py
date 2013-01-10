import numpy as np
import matplotlib.pyplot as plt

NN = 50            # Number of points in the problem
hbar = 1.054e-34   # Plank's constant
melec = 9.1e-31    # Mass of electron
eV2J = 1.6e-19     # Energy conversion factors
J2eV = 1/eV2J

del_x = .2e-9      # The cells size
dt    = 1e-16      # Time steps
ra    = (0.5*hbar/melec)*(dt/del_x**2) # ra must be < .1
DX    = del_x * 1e9
XX    = np.linspace(0,DX*NN,NN)  # Length in nm for plotting

# -------------- Specify the potential ----------------------
V = np.zeros(NN)
# V = eV2J *(0.005)
n = np.arange(20,40)
# V = eV2J * (.0005)*(np.abs(30-n))

plt.subplot(3,2,2)
plt.plot(XX,J2eV*V,'k')
plt.ylabel('V (eV)')
plt.xlabel('nm')
Umax = np.max(J2eV*V)
plt.axis([0, DX*NN, 0, Umax])

# ----- Initialize the test function ------------------------
sigma = 3.                    # Pulse width
nc = NN/2                     # Starting position of the pulse
prl = np.zeros(NN)
pim = np.zeros(NN)
ptot = 0.
for n in range(1,NN-1):
    prl[n] = np.exp(-1.*((n-nc)/sigma)**2)
    #prl[n] = np.exp(-1.*((n-nc)/sigma)**2) - np.exp(-1.*((n-nc-25)/sigma)**2)

ptot  = np.sum(prl**2 + pim**2)
pnorm = np.sqrt(ptot)

# Normalize and check
ptot = 0.
prl = prl / pnorm
pim = pim / pnorm
ptot = np.sum(prl**2 + pim**2)
print ptot


# ---------- Parameters for the eigenfunction -------------
Ein = input('Eigenenergy (eV) --->')
freq = Ein / (J2eV*2*np.pi*hbar)
omega = 2 * np.pi *freq
arg = omega * dt
np.exp(-1j*arg)
Time_period = 1/(freq*dt)
# ------- Specify the Hanning window ----------------------
Than = input('Length of Hanning window --->')
win = np.zeros(2^16)
n = np.arange(0,Than)
win = 0.5*(1-np.cos(2*np.pi*n/Than))

plt.subplot(3,2,2)
plt.plot(win,'k')
plt.axis([0, Than, 0, 1])
plt.ylabel('Hanning')
plt.xlabel('Time steps')

phi = np.zeros(NN)
phi0_rl = np.zeros(NN)
phi0_im = np.zeros(NN)
psi     = np.zeros(NN)
phi_m   = np.zeros(NN)

T = 0
n_step = 1
if n_step > 0:
    n_step = input('How many time steps --->')

# ------------ This is the main FDTD loop --------------------
#print win
for m in range(1,n_step):
   T = T + 1
    
   for n in range(1,NN-1):
       prl[n] = prl[n] - ra*(pim[n-1] - 2*pim[n] + pim[n+1]) 
       + (dt/hbar)*V[n]*pim[n]
     
   for n in range(1,NN-1):
       pim[n] = pim[n] + ra*(prl[n-1] - 2*prl[n] + prl[n+1])
       - (dt/hbar)*V[n]*prl[n]
   
   for n in range(1,NN-1):
       psi[n] = prl[n] + 1j*pim[n]
       phi[n] = phi[n] + win[T] * np.exp(-1j*arg*T)*psi[n]

print np.dot(psi, psi.conjugate())

plt.subplot(3,2,1)
plt.plot(XX,prl,'k')
plt.hold(True)

plt.plot(XX,prl,'--k')
plt.hold(False)

plt.axis( [0, DX*NN, -.25, .5])
plt.text(6,.3, "{0} ps".format(T*dt*1e12))
plt.xlabel('nm')
plt.title('Se12-2')

# Normalize phi
ptot1 = np.dot(phi,phi.conjugate())
phi_m = phi / ptot1
angle = np.arctan2(phi.imag, phi.real)
print ptot1

# --- Plot the complex reconstructed eigenfunction ---------
pmax = np.max(np.abs(phi_m))
plt.subplot(3,2,3)
plt.plot(XX, phi_m.real, 'k')
plt.hold(True)
plt.plot(XX,phi_m.imag, 'k--')
plt.hold(False)

plt.text(2.5, .5*pmax, " E_i_n = {0} meV".format(Ein*1e3))

plt.axis([0, DX*NN, -1.2*pmax, 1.2*pmax])
plt.xlabel('nm')
plt.ylabel('phi')

# --- Set the angle at the source to zero. This gives and 
#     eignefunction that is all real.
ang0 = angle[25]
angle = angle - ang0
phi0_rl = np.abs(phi_m)*np.cos(angle)
phi0_im = np.abs(phi_m)*np.sin(angle)

ptot = np.dot(phi0_rl, phi0_rl.conjugate())

phi0_rl = phi0_rl/np.sqrt(ptot)

pmax = np.max(phi0_rl) +.01
pmin = np.min(phi0_rl)

plt.subplot(3,2,4)
plt.plot(XX, phi0_rl, 'k')
plt.hold(True)
plt.plot(XX, phi0_im, 'k--')
plt.hold(False)

plt.axis([0, DX*NN, 1.2*pmin, 1.2*pmax])
plt.text(2.5,.5*pmax, "E_i_n = {0} meV".format(Ein*1e3))
plt.grid(True)

plt.grid(True)
plt.xlabel('nm')
print np.dot(phi0_rl, phi0_rl.conjugate())

plt.show() 

