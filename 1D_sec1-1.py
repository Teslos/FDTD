import numpy as np
import matplotlib.pyplot as plt

NN = 400             # Number of points in the problem
hbar = 1.054e-34     # Planck's constant
m0 = 9.1e-31         # Free space mass of an electron
meff = 1.0           # Effective mass: Si is 1.08, Ge is 0.067, GaAs is 0.55
melec = meff*m0      # Mass of an electron
ecoul = 1.6e-19      # Charge of an electron
epsz  = 8.85e-9      # Dielectric of free space
eV2J  = 1.6e-19      # Energy conversion factor
J2eV  = 1/eV2J

del_x = .1e-9        # The cell size
dt = 2e-17           # Time step
ra = (0.5*hbar/melec)*(dt/del_x**2)   # ra must be <.15
DX = del_x * 1e9     # Cell size in nm
XX = np.linspace(0,DX*NN,NN) # Length in nm for plotting
# --- Specify the potential ------------------------------
V = np.zeros(NN)

# Barrier
n = range(NN/2,NN/2+50)
#V[n] = .15 * eV2J

# Semiconductor conduction band
n = range(1,NN/2)
#V[n] = .1*eV2J

n = range(NN/2+1,NN)
#V[n] = .2*eV2J

# Electric field
# n = range(1,NN)
# V[n] = -(0.2/400)*(n-400)*eV2J

# Initialize a sine wave in a gaussian envelope

mlambda = 50.0          # Pulse wavelength
#mlambda = 25.0 
sigma = 50.0           # Pulse width
nc = 150               # Starting position
prl = np.zeros(NN)     # Real part of the state variable
pim = np.zeros(NN)     # Imaginary part of the state variable
ptot = 0.
for n in range(1,NN-1):
    prl[n] = np.exp(-1.*((n-nc)/sigma)**2)*np.cos(2*np.pi*(n-nc)/mlambda)
    pim[n] = np.exp(-1.*((n-nc)/sigma)**2)*np.sin(2*np.pi*(n-nc)/mlambda)
ptot = np.sum(prl**2 + pim**2)
pnorm = np.sqrt(ptot)         # Normalization constant
print('Normalization constant is: %2e' % pnorm)
# Normalize and check
ptot = 0.
n = range(0,NN)
prl[n] = prl[n]/pnorm
pim[n] = pim[n]/pnorm
ptot = np.sum(prl**2 + pim**2)
print ('Total probability:%2e ' % ptot)

T = 0
n_step = 1
if n_step > 0:
    n_step = input("How many time steps -->")

# This is the core FDTD program
for m in range(1,n_step):
    T = T + 1
    
    for n in range(1,NN-1):
        prl[n] = prl[n] - ra*(pim[n-1] - 2*pim[n] + pim[n+1]) 
        + (dt/hbar)*V[n]*pim[n]
    
    for n in range(1,NN-1):
        pim[n] = pim[n] + ra*(prl[n-1] - 2*prl[n] + prl[n+1])
        - (dt/hbar)*V[n]*prl[n]
#    print "Finished iteration {0}".format(T)

PE = 0.
psi = prl + 1j*pim                     # Write as a complex function
#print psi
PE = np.dot(psi,psi.conjugate()*V).real
print PE
print('Normalization after FDTD: %2e' % np.dot(psi,psi.conjugate()).real)                    # This checks normalization
#print('Potential energy: %2e' %  PE * J2eV)


ke = 0. + 1j * 0.
for n in range(1,NN-1):
    lap_p = psi[n+1] - 2 * psi[n] + psi[n-1]
    ke = ke + lap_p * psi[n].conjugate()
print ke.real
KE = -J2eV *(( hbar/del_x)**2/(2*melec)) * ke.real  # Kinetic energy

#plt.subplot(2,1,1)
#fig = plt.figure()

plt.plot(XX,prl,'k')
plt.hold(True)
plt.plot(XX,pim,'-.k')
plt.plot(XX,J2eV*V,'--k')
plt.hold(False)
plt.axis([1,DX*NN, -.2, .3])
plt.text(5,.15, "{0} fs".format(T*dt*1e15))
plt.text(5,-.15, "KE: {0}".format(KE))
plt.text(350,-.15,"PE: {0}" .format(PE))
plt.xlabel('nm')
plt.title('Se1-1')
plt.draw()
plt.show()
     


