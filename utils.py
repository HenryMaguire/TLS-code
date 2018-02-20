import numpy as np
from numpy import pi
import scipy as sp
from qutip import spre, spost, sprepost
import qutip as qt
import pickle


ev_to_inv_cm = 8065.5
inv_ps_to_inv_cm = 5.309

def load_obj(name ):
    with open(name + '.pickle', 'rb') as f:
        return pickle.load(f)

def save_obj(obj, name ):
    with open(name + '.pickle', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def Coth(x):
    return (np.exp(2*x)+1)/(np.exp(2*x)-1)

def beta_f(T):
    conversion = 0.695
    beta = 0
    if T ==0.: # First calculate beta
        beta = np.infty
    else:
        # no occupation yet, make sure it converges
        beta = 1. / (conversion*T)
    return beta

def J_multipolar(omega, Gamma, omega_0):
    return Gamma*(omega**3)/(2*np.pi*(omega_0**3))

def J_minimal(omega, Gamma, omega_0):
    return Gamma*omega/(2*np.pi*omega_0)

def J_minimal_hard(omega, Gamma, omega_0, cutoff=10*omega_0):
    if omega <cutoff:
        return Gamma*omega/(omega_0) #2*np.pi*
    else:
        return 0.

def J_flat(omega, Gamma, omega_0):
    return Gamma

def J_overdamped(omega, alpha, wc):
    return alpha*wc*omega/(omega**2+wc**2)

def J_underdamped(omega, alpha, Gamma, omega_0):
    return alpha*Gamma*pow(omega_0,2)*omega/(pow(pow(omega_0,2)-pow(omega,2),2)+(Gamma**2 *omega**2))

def rate_up(w, n, gamma, J, w_0):
    rate = 0.5 * pi * n * J(w, gamma, w_0)
    return rate

def rate_down(w, n, gamma, J, w_0):
    rate = 0.5 * pi * (n + 1. ) * J(w, gamma, w_0)
    return rate

def lin_construct(O):
    Od = O.dag()
    L = 2. * spre(O) * spost(Od) - spre(Od * O) - spost(Od * O)
    return L
