from scipy.integrate import quad
from qutip import ket, qeye, tensor, destroy, Qobj
import numpy as np
import matplotlib.pyplot as plt
from sympy.functions import coth
from utils import *

G = ket([0])
E = ket([1])

def get_data():

    rho = [rho_01, rho_10]
    return

def integrand_OD(omega, t, alpha, beta, wc):
    if omega ==0:
        return 0.
    else:
        return (1/omega**2)*J_overdamped(omega, alpha,wc)*coth(beta*omega/2)*(1-np.cos(omega*t))



def integrand_UD(omega, t, alpha, beta, Gamma, w0):
    if omega ==0:
        return 0.
    else:
        return (1/omega**2)*J_underdamped(omega, alpha, Gamma, w0)*coth(beta*omega/2.)*(1-np.cos(omega*t))

def integral_converge(f, a):
    step = 110.
    x = step
    I = 0
    while abs(f(x))>5E-9:
        I += quad(f, a, x)[0]
        a+=step
        x+=step
    return I # Converged integral

def plot_integrand(beta, wc):
    omega = np.linspace(1E-11,20000,10000)
    y1 = []
    y2 = []
    for i in omega:
        y1.append(integrand_OD(i, 0.001, alpha, beta, wc))
        y2.append(integrand_OD(i, 0.002, alpha, beta, wc))
    plt.plot(omega, y1)
    plt.plot(omega, y2)
    plt.show()

def exact_decay(t, alpha, beta, Gamma, w0):
    f = lambda omega : integrand_UD(omega, t, alpha, beta, Gamma, w0)
    I = -integral_converge(f, 0.)
    return I

def exact_solution_at_t(t, eps, alpha, beta, Gamma, w0, rho_init):
    Decay = exact_decay(t, alpha, beta, Gamma, w0)
    rho_01 = rho_init.matrix_element(E.dag(), G)*np.exp(Decay)*np.exp(-1j*eps*t)
    return rho_01


def exact_dynamics(eps, alpha, wc, w0, Gamma, beta, rho_init, time_points, overdamped=False):
    eps = eps + np.pi*alpha/2.
    rho_t = []
    if overdamped:
        Gamma = w0**2/wc
    for t in time_points:
        rho_t.append(exact_solution_at_t(t, eps, alpha, beta, Gamma, w0, rho_init))
        if len(rho_t)%40 == 0:
            print "Exact solution {:0.3f} percent finished".format((float(len(rho_t))/len(time_points))*100)
    return rho_t

def absorption_integrand(t, omega, eps, shift, alpha, beta, Gamma, omega_0):
    #Similar to Integrand function but slightly different since was calculated in Heisenberg picture
    I = (1/omega**2)*J_underdamped(omega, alpha, Gamma, omega_0)
    I*= complex(coth(beta*omega/2)*(1-np.cos(omega*t))+1j*np.sin(omega*t))
    I = np.exp(-I)
    I*= np.exp(1j*(omega-eps)*t)
    return I

def absorption(omega, eps, shift, alpha, beta, Gamma, omega_0, mu=1.):
    I = np.array(quad(absorption_integrand, 0, 0.2, args=(omega, eps, shift, alpha, beta, Gamma, omega_0),limit=150))
    return mu**2*(I.real)


if __name__ == "__main__":

    G = ket([0])
    E = ket([1])
    rho_init = 0.5*(G+E)*(E.dag()+G.dag())
    #x = np.arange(1.,5,0.5)
    alpha = 1000/np.pi
    beta = 1/(0.695*300)
    eps = 8000.
    wc = 53.
    omega_0 = 500.
    Gamma = omega_0**2/wc
    plot_integrand(beta, wc)

    """
    T, rho_t = exact_dynamics(eps, wc, omega_0, Gamma, beta, rho_init, 0, 0.02, 1000)
    reData = np.loadtxt("DATA/RealExactdata.dat")
    imData = np.loadtxt("DATA/ImagExactdata.dat")
    plt.plot(list(zip(*reData)[0]), list(zip(*reData)[1]), label='math real')
    #plt.plot(list(zip(*imData)[0]), list(zip(*imData)[1]), 'math imag')
    plt.plot(T, np.array(rho_t).real, label='real')
    #plt.plot(T, np.array(rho_t).imag, label='imag')
    plt.legend()
    plt.show()
    """
    """
    m = []
    frequencies = np.linspace(-1000,1000,100)
    for w in frequencies:
        m.append(absorption(w, eps, 200, alpha, beta, Gamma, omega_0, mu=1.))
    plt.plot(frequencies, m)
    #print exact_decay(1., alpha, beta, Gamma, omega_0)
    """
