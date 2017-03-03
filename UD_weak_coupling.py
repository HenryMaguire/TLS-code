"""
Here I will code up the weak-coupling version of the electron-phonon interaction
"""

import numpy as np
import scipy as sp
import qutip as qt
from qutip import destroy, tensor, qeye, spre, spost, sprepost, ket
import driving_liouv as EM
from sympy.functions import coth
import matplotlib.pyplot as plt

def J_UD_SB(omega, alpha, omega_0, Gamma):
    return (alpha*Gamma*(omega_0**2))*omega/((omega_0**2 - omega**2)**2+((Gamma**2)*(omega**2)))

def L_ph(sigma, eps, T_ph, wc, w0, alpha_ph, time_units='cm'):
    A = sigma.dag()*sigma
    Gamma = (w0**2)/wc
    beta = 1/(0.695*T_ph)
    dephasing_rate = np.pi*J_UD_SB(eps, alpha_ph, w0, Gamma)*(coth(0.5*beta*eps)+1)
    return float(dephasing_rate)*(sprepost(A,A)- 0.5*(spre(A) + spost(A)))


if __name__ == '__main__':
    """
    Define all system and environment  parameters
    """
    N = 6
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
    time_units = 'cm'
    eps = 1.*8065.5 # TLS splitting
    #eps = 2.*1519.3 # ps
    T_EM = 6000. # Optical bath temperature
    #alpha_EM = 0.3 # System-bath strength (optical)
    Gamma_EM = 6.582E-4*8065.5 #bare decay of electronic transition in inv. cm
    #Gamma_EM = 6.582E-7*1519.3
    T_ph = 300. # Phonon bath temperature
    wc = 53. # Ind.-Boson frame phonon cutoff freq
    #wc = 53.*0.188
    w0 = 300. # underdamped SD parameter omega_0
    #w0 = 200.*0.188
    alpha_ph = (10./np.pi)# Ind.-Boson frame coupling
    H_S = eps*E*E.dag()
    L_wc = L_ph(sigma, eps, T_ph, wc, w0, alpha_ph, time_units='cm')
    #J = EM.J_multipolar
    rho_init = 0.5*(G+E)*((E+G).dag())
    expects = [E*E.dag(), G*E.dag()]
    timelist = np.linspace(0,100.,15000)#*0.188
    #print H_S, rho_init
    DATA = qt.mesolve(H_S, rho_init, timelist, [L_wc],  expects)
    #mesolve(H, rho_0, timelist, [L_RC+L_ns], expects, progress_bar=True)
    plt.plot(timelist, DATA.expect[1])
    plt.plot(timelist, DATA.expect[0])
    plt.show()
