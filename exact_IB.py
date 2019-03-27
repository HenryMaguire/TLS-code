import time

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



def integrand_Re(omega, t, alpha, beta, Gamma, w0):
    if omega ==0:
        return 0.
    else:
        return (1/omega**2)*J_underdamped(omega, alpha, Gamma, w0)*(coth(beta*omega/2.)*(1-np.cos(omega*t)))

def integrand_Im(omega, t, alpha, beta, Gamma, w0):
    if omega ==0:
        return t*J_underdamped(1., alpha, Gamma, w0)
    else:
        return (1/omega**2)*J_underdamped(omega, alpha, Gamma, w0)*np.sin(omega*t)


def integral_converge(f, a):
    step = 60.+np.random.random()
    x = step
    I = 0
    der = 1000.
    while (abs(f(x))>1E-9):
        I += quad(f, a, x)[0]
        der = abs(f(a)-f(x))/step
        a+=step
        x+=step
    return I # Converged integral

def plot_integrand(beta, alpha, Gamma, w0):
    omega = np.linspace(0.,1000,1000)
    data = [[], [], [], []]
    T = [0.001, 0.01, 0.001, 0.01]
    for w in omega:
        for j, t in enumerate(T):
            if j <=1:
                data[j].append(integrand_Re(w, t, alpha, beta, Gamma, w0))
            else:
                data[j].append(integrand_Im(w, t, alpha, beta, Gamma, w0))
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', figsize=(11,8))
    axes = [eval("ax{}".format(i)) for i in range(1,5)]
    for i, dat in enumerate(data): # each vibronic level of interest
        axes[i].set_title(T[i])
        axes[i].plot(omega,data[i])
        axes[i].set(xlabel=r"decay rate",
               ylabel="time")
    plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=1.0)
    plt.show()

def exact_decay(t, alpha, beta, Gamma, w0):
    fRe = lambda x: integrand_Re(x, t, alpha, beta, Gamma, w0)
    fIm = lambda x: integrand_Im(x, t, alpha, beta, Gamma, w0)
    I = integral_converge(fRe, 0.)
    I += 1j*integral_converge(fIm, 0.)
    return -I

def exact_solution_at_t(t, eps, alpha, beta, Gamma, w0, rho_init):
    Decay = exact_decay(t, alpha, beta, Gamma, w0)
    sigma_x = (rho_init.matrix_element(E.dag(), G)+rho_init.matrix_element(E.dag(), G))*np.exp(Decay)*np.exp(-1j*eps*t)
    return sigma_x


def exact_dynamics(eps, alpha, wc, w0, Gamma, beta, rho_init, time_points, overdamped=False , silent=False):
    rho_t = []
    alpha = alpha
    print("eps: {}, alpha: {}, w0: {}, Gamma: {}, beta: {}".format(eps, alpha, w0, Gamma, beta))
    if overdamped:
        Gamma = (w0**2)/wc
    ti = time.time()
    for t in time_points:
        rho_t.append(exact_solution_at_t(t, eps, alpha, beta, Gamma, w0, rho_init))
        if not silent:
            if len(rho_t)%40 == 0:
                print("Exact solution {:0.3f} percent finished".format((float(len(rho_t))/len(time_points))*100))
    if not silent:
        print("Exact solution took {} to calculate.".format(time.time()- ti))
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

    plt.style.use('ggplot')
    plt.rcParams["axes.grid"] = True
    plt.rcParams["axes.edgecolor"] = "0.15"
    plt.rcParams["axes.linewidth"]  = 1.25
    G = ket([0])
    E = ket([1])
    rho_init = 0.5*(G+E)*(E.dag()+G.dag())
    #x = np.arange(1.,5,0.5)

    beta = 1/(0.695*300)
    eps = 8000.
    alpha = 1000/np.pi
    wc = 53.
    omega_0 = 500.
    Gamma = omega_0**2/wc
    #plot_integrand(beta, wc, Gamma, omega_0)
    timesteps = np.linspace(0, 0.02, 100)
    # defaults/API: exact_dynamics(eps, alpha, wc, w0, Gamma, beta, rho_init, time_points, overdamped=False , silent=False)
    rho_t = exact_dynamics(eps, alpha, wc, omega_0, Gamma, beta, rho_init, timesteps, overdamped=True)
    #reData = np.loadtxt("DATA/RealExactdata.dat")
    #imData = np.loadtxt("DATA/ImagExactdata.dat")
    #plt.plot(list(zip(*reData)[0]), list(zip(*reData)[1]), label='math real')
    #plt.plot(list(zip(*imData)[0]), list(zip(*imData)[1]), 'math imag')
    plt.plot(timesteps, np.array(rho_t).real, label='real')
    plt.plot(timesteps, np.array(rho_t).imag, label='imag')
    plt.legend()
    plt.show()

    """
    m = []
    frequencies = np.linspace(-1000,1000,100)
    for w in frequencies:
        m.append(absorption(w, eps, 200, alpha, beta, Gamma, omega_0, mu=1.))
    plt.plot(frequencies, m)
    #print exact_decay(1., alpha, beta, Gamma, omega_0)
    """
