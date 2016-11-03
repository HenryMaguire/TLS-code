from UD_liov import RC_function_UD
from driving_liouv import *



def plot_dynamics():
    plt.figure()
    if T_EM>0.0:
        ss_ns = steadystate(H, [L_RC+L_ns]).ptrace(0)
        ss_v = steadystate(H, [L_RC+L_v]).ptrace(0)
        ss_n = steadystate(H, [L_RC+L_n]).ptrace(0)
        ss_g_ns = ss_ns.matrix_element(G.dag(), G)
        ss_g_v = ss_v.matrix_element(G.dag(), G)
        ss_g_n = ss_n.matrix_element(G.dag(), G)
        plt.axhline(1-ss_g_v, color='b', ls='--')
        #plt.axhline(1-ss_g_ns, color='g', ls='--')
        plt.axhline(1-ss_g_n, color='r', ls='--')

    plt.title(r"$\omega_0=$""%i"r"$cm^{-1}$, $\alpha_{Ph}=$""%f"r"$cm^{-1}$, $T_{EM}=$""%i K" %(w0, alpha_Ph, T_EM))
    plt.plot(timelist, 1-DATA_ns.expect[0], label='Non-secular', color='g')
    plt.plot(timelist, 1-DATA_v.expect[0], label='Vib. Lindblad', color='b')
    plt.plot(timelist, 1-DATA_n.expect[0], label='Simple Lindblad', color='r')
    plt.ylabel("Excited state population")
    plt.xlabel("Time (cm)")
    plt.legend()

    plt.figure()
    plt.title(r"$\alpha_{Ph}=$""%i"r"$cm^{-1}$, $T_{EM}=$""%i K" %(alpha_Ph, T_EM))
    plt.plot(timelist, DATA_v.expect[1], label='Vib. Lindblad', color='b')
    plt.plot(timelist, DATA_ns.expect[1], label='Non-secular', color='g',alpha=0.7)
    plt.plot(timelist, DATA_n.expect[1], label='Simple Lindblad', color='r',alpha=0.4)
    plt.legend()
    plt.ylabel("Coherence")
    plt.xlabel("Time (cm)")
    plt.show()

if __name__ == "__main__":

    G = ket([0])
    E = ket([1])
    N = 6
    dipoles = []
    differences = []
    sigma = G*E.dag()

    eps = 500. # TLS splitting

    T_EM = 6000.
    alpha_em = 0.3

    T_Ph = 300.
    wc = 53. # Ind.-Boson frame cutoff freq
    w0 = 500.
    alpha_Ph = 200 # 2./np.pi # Ind.-Boson frame coupling


    L_RC, H, wRC, kappa = RC_function_UD(eps, T_EM, T_Ph, wc, w0, alpha_Ph, alpha_em, N)
    H, A_em, A_ph = Ham_RC(eps, wRC, kappa, N)
    L_ns = nonsecular(H, A_em, alpha_em, T_EM)

    n_RC = Occupation(wRC, T_Ph)
    rho_0 = tensor(0.5*((E+G)*(E+G).dag()), thermal_dm(N, n_RC))

    expects = [tensor(G*G.dag(), qeye(N)), tensor(G*E.dag(), qeye(N)), tensor(qeye(2), destroy(N).dag()*destroy(N))]
    timelist = np.linspace(0,15,20000)
    #nonsec_check(eps, H, A_em, N)

    DATA_ns = mesolve(H, rho_0, timelist, [L_RC+L_ns], expects, progress_bar=True)
    DATA_v = mesolve(H, rho_0, timelist, [L_RC+L_v], expects, progress_bar=True)
    DATA_n = mesolve(H, rho_0, timelist, [L_RC+L_n], expects, progress_bar=True)
