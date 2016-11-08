from qutip import ket, mesolve, qeye, tensor, thermal_dm, destroy, steadystate
import matplotlib.pyplot as plt
import numpy as np
import UD_liouv as RC
import driving_liouv as diss
import ME_checking as check
reload(RC)
reload(diss)
reload(check)



def plot_dynamics():
    if T_EM>0.0:
        #ss_ns = steadystate(H, [L_RC+L_ns]).ptrace(0)
        ss_v = steadystate(H, [L_RC+L_s]).ptrace(0)
        ss_n = steadystate(H, [L_RC+L_naive]).ptrace(0)
        #ss_g_ns = ss_ns.matrix_element(G.dag(), G)
        ss_g_v = ss_v.matrix_element(G.dag(), G)
        ss_g_n = ss_n.matrix_element(G.dag(), G)
        plt.axhline(1-ss_g_v, color='b', ls='--')
        #plt.axhline(1-ss_g_ns, color='g', ls='--')
        plt.axhline(1-ss_g_n, color='r', ls='--')

    plt.title(r"$\omega_0=$""%i"r"$cm^{-1}$, $\alpha_{ph}=$""%f"r"$cm^{-1}$, $T_{EM}=$""%i K" %(w0, alpha_ph, T_EM))
    #plt.plot(timelist, 1-DATA_nrwa.expect[0], label='nrwa', color='y')
    #plt.plot(timelist, 1-DATA_ns.expect[0].real, label='Non-secular', color='g')
    plt.plot(timelist, 1-DATA_s.expect[0].real, label='Vib. Lindblad', color='b')
    plt.plot(timelist, 1-DATA_naive.expect[0].real, label='Simple Lindblad', color='r')
    plt.ylabel("Excited state population")
    plt.xlabel("Time (cm)")
    plt.legend()
    #plt.savefig('dynamics_a'.format({alpha}))
    """

    plt.title(r"$\alpha_{ph}=$""%i"r"$cm^{-1}$, $T_{EM}=$""%i K" %(alpha_ph, T_EM))
    plt.plot(timelist, DATA_s.expect[1].real, label='Vib. Lindblad', color='b')
    #plt.plot(timelist, DATA_ns.expect[1].real, label='Non-secular', color='g',alpha=0.7)
    plt.plot(timelist, DATA_naive.expect[1].real, label='Simple Lindblad', color='r',alpha=0.4)
    plt.legend()
    plt.ylabel("Coherence")
    plt.xlabel("Time (cm)")
    """

def plot_spectrum(dyn, t):
    ss = dyn[-1]
    gg = dyn-ss
    spec = np.fft.fft(gg)
    freq = np.fft.fftfreq(t.shape[-1])
    plt.plot(freq, abs(spec.real))



if __name__ == "__main__":

    N = 5
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.

    eps = 18000. # TLS splitting

    T_EM = 6000.
    alpha_EM = 3.

    T_ph = 300.
    wc = 53. # Ind.-Boson frame cutoff freq
    w0 = 300.
    alpha_ph = 400. # 2./np.pi # Ind.-Boson frame coupling
    L_RC, H, A_EM, A_nrwa, wRC, kappa= RC.RC_function_UD(sigma, eps, T_ph, wc, w0, alpha_ph, N)

    # electromagnetic liouvillians
    L_nrwa = diss.L_nonrwa(H, A_nrwa, alpha_EM, T_EM) # Ignore this for now as it just doesn't work
    L_ns = diss.L_nonsecular(H, A_EM, alpha_EM, T_EM)
    L_s = diss.L_vib_lindblad(H, A_EM, alpha_EM, T_EM)
    L_naive = diss.L_EM_lindblad(eps, A_EM, alpha_EM, T_EM)

    # Set up the initial density matrix
    n_RC = diss.Occupation(wRC, T_ph)
    initial_sys = G*G.dag()
    #initial_sys = 0.5*(G+E)*(E.dag()+G.dag())
    rho_0 = tensor(initial_sys, thermal_dm(N, n_RC))

    # Expectation values and time increments needed to calculate the dynamics
    expects = [tensor(G*G.dag(), qeye(N)), tensor(G*E.dag(), qeye(N)), tensor(qeye(2), destroy(N).dag()*destroy(N))]
    timelist = np.linspace(0,10,20000)
    #nonsec_check(eps, H, A_em, N) # Plots a scatter graph representation of non-secularity.
    # Calculate dynamics
    #DATA_nrwa = mesolve(H, rho_0, timelist, [L_RC+L_nrwa], expects, progress_bar=True)
    #DATA_ns = mesolve(H, rho_0, timelist, [L_RC+L_ns], expects, progress_bar=True)
    DATA_s = mesolve(H, rho_0, timelist, [L_RC+L_s], expects, progress_bar=True)
    DATA_naive = mesolve(H, rho_0, timelist, [L_RC+L_naive], expects, progress_bar=True)
    #plt.figure()
    #plot_dynamics()
    #plt.show()
    check.SS_convergence_check(eps, T_EM, T_ph, wc, w0, alpha_ph, alpha_EM)





    #plot_spectrum(1-DATA_ns.expect[0], timelist)
    #np.savetxt('DATA_ns.txt', np.array([1- DATA_ns.expect[0], timelist]), delimiter = '\n', newline= ',')
    #plt.show()
    #raw_input("press return to close graph: ")
    #plt.close()
