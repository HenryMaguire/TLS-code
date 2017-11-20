"""
Here I will code up the weak-coupling version of the electron-phonon interaction
"""

import numpy as np
import scipy as sp
import qutip as qt
from qutip import destroy, tensor, qeye, spre, spost, sprepost, ket, basis
import driving_liouv as EM
from sympy.functions import coth
import matplotlib.pyplot as plt
from utils import *

def integral_converge(f, a, omega):
    x = 1
    I = 0
    print abs(f(x)), 'f(x)'
    while abs(f(x))>0.0001:
        print a, x
        I += integrate.quad(f, a, x, weight='cauchy', wvar=omega)[0]
        a+=1
        x+=1
    return I # Converged integral

def K(omega, beta, J, alpha, wc, imag_part=True):
    G = 0
    # Here I define the functions which "dress" the integrands so they
    # have only 1 free parameter for Quad.
    F_0 = lambda x: J(omega, alpha, wc)
    w='cauchy'
    G = (np.pi/2)*(2*alpha/(wc*beta))
    # The limit as omega tends to zero is zero for superohmic case?
    if imag_part:
        G += -(1j)*integral_converge(F_0, -1e-12,0)
    #print (integrate.quad(F_0, -1e-12, 20, weight='cauchy', wvar=0)[0])
    return G

def L_phonon(J, alpha, beta, wc):
    XX = basis(2,1)*basis(2,1).dag()
    gamma_0 = np.pi*alpha/(beta*wc)
    #S = K(0, beta, J, alpha, wc).imag
    return -0.5*gamma_0*((spre(XX) + spost(XX)) - 2*sprepost(XX,XX)) + 0.5*1j*np.pi*alpha*(spre(XX) - spost(XX))
"""
def L_ph(XX, eps, beta, wc, w0, alpha_ph, kbT):
    dephasing_rate = 2*np.pi*alpha_ph*kbT
    return dephasing_rate*(sprepost(XX,XX)- 0.5*(spre(XX) + spost(XX)))
"""

if __name__ == '__main__':
    """
    Define all system and environment  parameters
    """
    G = ket([0])
    E = ket([1])
    XX = E*E.dag() # Definition of a sigma_- operator.
    time_units = 'cm'
    eps = 1.1*8065.5#*8065.5 # TLS splitting
    #eps = 2.*1519.3 # ps
    T_EM = 6000. # Optical bath temperature
    #alpha_EM = 0.3 # System-bath strength (optical)
    #Gamma_EM = 6.582E-4*8065.5 #bare decay of electronic transition in inv. cm
    #Gamma_EM = 6.582E-7*1519.3
    T_ph = 300. # Phonon bath temperature
    wc = 53. # Ind.-Boson frame phonon cutoff freq
    #wc = 53.*0.188
    w0 = 300. # underdamped SD parameter omega_0
    #w0 = 200.*0.188
    alpha_ph = 4/np.pi#(10./np.pi)# Ind.-Boson frame coupling
    kbT = (0.695*T_ph)
    mu = 1.
    H_S = eps*E*E.dag()
    L_wc = L_ph(XX, eps, T_ph, wc, w0, alpha_ph, kbT)
    #J = EM.J_multipolar
    #rho_init = 0.5*(G+E)*((E+G).dag())
    rho_init = E*G.dag()
    expects = [G*E.dag(), G*E.dag()+E*G.dag()]
    timelist = np.linspace(0,0.01,10000)#*0.188
    #print H_S, rho_init
    DATA = qt.mesolve(H_S, rho_init, timelist, [L_wc],  expects)
    #mesolve(H, rho_0, timelist, [L_RC+L_ns], expects, progress_bar=True)
    #plt.plot(timelist, DATA.expect[1].real)
    #plt.plot(timelist, DATA.expect[1].imag)

    C = mu**2*(DATA.expect[1])
    plt.plot(timelist, C)
    plt.figure()
    m = np.fft.fft(C)
    n = C.size
    freq = np.fft.fftfreq(n, 0.1)
    plt.scatter(freq,m)
    plt.show()
