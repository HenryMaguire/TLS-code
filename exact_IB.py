from scipy.integrate import quad
from qutip import ket, qeye, tensor, destroy, Qobj
import numpy as np
import matplotlib.pyplot as plt
from sympy.functions import coth
from utils import *

def get_data():

    rho = [rho_01, rho_10]
    return


def integrand(omega, t, alpha, beta, Gamma, omega_0):
    if omega ==0:
        return 0.
    else:
        return (1/omega**2)*J_underdamped(omega, alpha, Gamma, omega_0)*coth(beta*omega/2)*(1-np.cos(omega*t))

def plot_integrand():
    omega = np.linspace(0.0001,20000,100000)
    y = []
    for i in omega:
        y.append(integrand(i, 1, alpha, beta, Gamma, omega_0))
    plt.plot(omega, y)
    plt.show()

def exact_decay(t, alpha, beta, Gamma, omega_0):
    I = -np.array(quad(integrand, 0.00001, 20000, args=(t, alpha, beta, Gamma, omega_0),limit=150))
    return I

def exact_solution_at_t(t, alpha, beta, Gamma, omega_0, rho_init):
    Gamma = exact_decay(t, alpha, beta, Gamma, omega_0)[0]
    rho_01 = rho_init.matrix_element(G.dag(), E)*np.exp(Gamma)*np.exp(-1j*eps*t)
    return rho_01


def exact_solution(eps, wc, beta, rho_init, a, b, points):
    T = np.linspace(a, b, points)
    rho_t = []
    for t in T:
        rho_t.append(exact_solution_at_t(t, alpha, beta, Gamma, omega_0, rho_init))
    return T, rho_t


if __name__ == "__main__":

    G = ket([0])
    E = ket([1])
    rho_init = 0.5*(G+E)*(E.dag()+G.dag())
    #x = np.arange(1.,5,0.5)
    alpha = 400/np.pi
    beta = 1/(0.695*300)
    eps = 8000.
    wc = 153.
    omega_0 = 2000.
    Gamma = omega_0**2/wc
    T, rho_t = exact_solution(eps, wc, beta, rho_init, 0, 0.02, 1000)
    reData = np.loadtxt("DATA/RealExactdata.dat")
    imData = np.loadtxt("DATA/ImagExactdata.dat")
    plt.plot(list(zip(*reData)[0]), list(zip(*reData)[1]), label='math real')
    #plt.plot(list(zip(*imData)[0]), list(zip(*imData)[1]), 'math imag')
    plt.plot(T, np.array(rho_t).real, label='real')
    #plt.plot(T, np.array(rho_t).imag, label='imag')
    plt.legend()
    plt.show()
    #t = 1
    print exact_decay(1., alpha, beta, Gamma, omega_0)
