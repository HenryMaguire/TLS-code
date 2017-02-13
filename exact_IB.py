import scipy as sp
from qutip import ket, qeye, tensor, destroy, Qobj
import numpy as np
import matplotlib.pyplot as plt
import sympy.functions.elementary.hyperbolic as hyp

def integrand(omega, t, alpha, KT, wc):
    return -4.*alpha*(omega*wc/(omega**2 + wc**2))*hyp.coth(omega/2*KT)*(1-np.cos(omega*t))

def exact_decay(t, alpha, KT, wc):
    I = np.array(sp.integrate.quad(integrand, 0, np.infty, args=(t, alpha, KT, wc)))
    return I

def exact_solution(t, eps, wc, KT, rho_init):
    Gamma = exact_decay(t, alpha, KT, wc)[0]
    print Gamma
    rho_00 = rho_init.matrix_element(G.dag(), G)
    rho_11 = rho_init.matrix_element(E.dag(), E)
    rho_01 = rho_init.matrix_element(G.dag(), E)*np.exp(Gamma)*np.exp(-1j*eps*t)
    rho_10 = rho_init.matrix_element(E.dag(), G)*np.exp(Gamma)*np.exp(1j*eps*t)
    rho = Qobj([[rho_00, rho_01],[rho_10, rho_11]])
    return rho


print integrand(1., 2., 0.1, 1., 10.)

G = ket([0])
E = ket([1])

x = np.arange(1.,5,0.5)
alpha = 1.
KT = 300.
wc = 53.
eps = 10.
t = 0.001

rho_1 = G*G.dag()
rho_2 = E*E.dag()
rho_3 = 0.5*(G+E)*(E.dag()+G.dag())

print exact_decay(0.1, alpha, KT, wc)
print exact_solution(t, eps, wc, KT, rho_1)
print exact_solution(t, eps, wc, KT, rho_2)
print exact_solution(t, eps, wc, KT, rho_3)
#plt.plot(x, integrand(x, 0, alpha, KT, wc))
