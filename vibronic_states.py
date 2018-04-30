from qutip import ket, mesolve, qeye, tensor, thermal_dm, destroy, steadystate
import qutip as qt
import matplotlib.pyplot as plt
import numpy as np

import UD_liouv as RC
import driving_liouv as EM
import ME_checking as check
import exact_IB as exact
import scipy as sp
import phonon_weak_coupling as WC
from utils import J_overdamped, beta_f, J_underdamped, J_minimal_hard
reload(RC)
reload(EM)
reload(check)
reload(exact)
plt.style.use('ggplot')
def steadystate_coupling_dependence(prop_couplings, T_ph, eps, Gamma, w0, T_EM, Gamma_EM, overdamped=True):
    plt.figure()
    N = 15
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
    J = EM.J_minimal
    wc = []
    rc_nrwa = []
    rc_ns = []
    rc_s = []
    rc_naive = []
    count = 1

    data  = [[],[],[],[],[],[]]
    labels = ['g-0', 'e-0', 'g-1', 'e-1', 'g-2', 'e-2'] # This is true for specific Hamiltonian parameters
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for prop in prop_couplings:
        alpha_ph = prop*eps
        print count
        count+=1
        L_RC, H_RC, A_EM, A_nrwa, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps, T_ph, Gamma, w0, alpha_ph, N)
        lowest_6 = H_RC.eigenstates()[1][0:6]
        #L_wc = WC.L_phonon(alpha_ph, beta, Gamma, w0, overdamped=overdamped)
        #L_wc_EM = EM.L_EM_lindblad(eps, sigma, Gamma_EM, T_EM, J=J)
        # electromagnetic bath liouvillians
        #L_nrwa = EM.L_non_rwa(H_RC, A_nrwa, eps, Gamma_EM, T_EM, J) # Ignore this for now as it just doesn't work
        L_ns = EM.L_nonsecular(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
        #L_s = EM.L_vib_lindblad(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
        #L_naive = EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J)

        #print "NON-ROTATING-WAVE"


        for j, state in enumerate(lowest_6):
            data[j].append(steadystate(H_RC, [L_RC+L_ns]).matrix_element(state.dag(), state))

    for lab, dataset in zip(labels, data):
        plt.plot(prop_couplings, dataset, label=lab)

    plt.legend()
    plt.grid()
    plt.ylabel("Population")
    plt.xlabel(r"Phonon Coupling Strength: units of $\epsilon$")
    plt.show()

G = ket([0])
E = ket([1])
sigma = G*E.dag() # Definition of a sigma_- operator.
## ~~ Parameters ~~

eps = 0.1*8065.5 # TLS splitting
#eps = 2.*1519.3 # ps
T_EM = 6000. # Optical bath temperature
H_S = eps*E*E.dag()
#alpha_EM = 0.3 # System-bath strength (optical)
Gamma_EM = 0.1*5.309 #bare decay of electronic transition from inv. ps to in inv. cm
#Gamma_EM = 6.582E-7*1519.3
T_ph = 300. # Phonon bath temperature
overdamped = True
phonon_only = True
#fig1wc = 53.*0.188

#w0 = 200.*0.188
Gamma = 60. # Width of distribution
alpha_ph = 0.3*eps # Ind.-Boson frame coupling
beta = beta_f(T_ph)#1/(0.695*T_ph)
wc = 53. # Ind.-Boson frame phonon cutoff freq
shift = 0.5*np.pi*alpha_ph
w0 = eps+10 #-100 # underdamped SD parameter omega_0
optical_cutoff = 20000.
if overdamped:
    Gamma = w0**2/wc

couplings = np.linspace(0,0.1,30)
steadystate_coupling_dependence(couplings, T_ph, eps, Gamma, w0, T_EM, Gamma_EM, overdamped=True)
