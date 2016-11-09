from qutip import ket
from driving_liouv import L_vib_lindblad
from qutip import ket, mesolve, qeye, tensor, thermal_dm, destroy, steadystate
import matplotlib.pyplot as plt
import numpy as np
import UD_liouv as RC
import driving_liouv as diss

reload(RC)
reload(diss)

"""
def convergence_check(eps, T_EM, T_Ph, wc, alpha_ph, alpha_em, N):
    timelist = np.linspace(0,5,20000)
    coh_data = []
    G = ket([0])
    E = ket([1])
    for n in range(15,25):
        L_n, L_v, H, wRC = RC_function(eps, T_EM, T_Ph, wc, alpha_Ph, alpha_em, n)
        expects = [tensor(E*E.dag(), qeye(n)), tensor(G*E.dag(), qeye(n))]
        n_RC = Occupation(wRC, T_Ph)
        rho_0 = tensor(0.5*((E+G)*(E+G).dag()), thermal_dm(n, n_RC))
        DATA_v = mesolve(H, rho_0, timelist, [L_v], expects, progress_bar=True)
        coh_data.append(DATA_v.expect[0])
    plt.figure()
    for i in coh_data:
        plt.plot(i.expect[0])
    plt.legend()
"""

def SS_convergence_check(eps, T_EM, T_ph, wc, w0, alpha_ph, alpha_EM, expect_operator='excited', time_units='cm', start_n=14, end_n=20):
    # Only for Hamiltonians of rotating wave form
    G = ket([0])
    E = ket([1])
    ss_list = [] # steady states
    expect_operator = E*E.dag()
    if expect_operator == 'coherence':
        expect_operator = G*E.dag()
    else:
        pass
    for n in range(start_n,end_n):
        L_RC, H, A_EM, A_nrwa, wRC, kappa = RC.RC_function_UD(G*E.dag()+E*G.dag(), eps, T_ph, wc, w0, alpha_ph, n)
        L_s = L_vib_lindblad(H, A_EM, alpha_EM, T_EM)
        expects = []
        ss = steadystate(H, [L_RC+L_s]).ptrace(0)
        ss_E = (ss*expect_operator).tr()
        ss_list.append(ss_E)
        print "N=", n, "\n -----------------------------"
    return ss_list, range(start_n,end_n)
    plt.figure()
    plt.plot(range(15,25), ss_list)

def nonsec_check(eps, H, A, N):
    """
    Plots a scatter graph with a crude representation of how dominant non-secular terms are.
    The idea is that 'slowly' oscillating terms with large coefficients
        should contribute most to the non-secularity in the dynamics.
    """
    dipoles = []
    TD = []
    evals, evecs = H.eigenstates()
    for i in range(2*N):
        for j in range(2*N):
            for k in range(2*N):
                for l in range(2*N):
                    eps_ij = evals[i]-evals[j]
                    eps_kl = evals[k]-evals[l]
                    A_ij = A_em.matrix_element(evecs[i].dag(),evecs[j])
                    A_kl_conj = (A_em.dag()).matrix_element(evecs[l].dag(),evecs[k])
                    TD.append(eps_ij-eps_kl)
                    dipoles.append(A_ij*A_kl_conj)
    plt.figure()
    plt.scatter(TD, dipoles)
    plt.grid()
    plt.title(r"Non-secularity check: $\epsilon=$""%i"%eps)
    plt.xlabel(r"$\varphi_{ij}-\varphi_{kl}$")
    plt.ylabel(r"Dipole overlap of transitions $A_{ij}A_{kl}$")
    return TD, dipoles

def plot_spectrum(dyn, t):
    ss = dyn[-1]
    gg = dyn-ss
    spec = np.fft.fft(gg)
    freq = np.fft.fftfreq(t.shape[-1])
    plt.plot(freq, abs(spec.real))
