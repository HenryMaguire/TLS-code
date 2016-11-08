def convergence_check(eps, T_EM, T_Ph, wc, alpha_ph, alpha_em, N):
    expects = [tensor(G*E.dag(), qeye(N))]
    timelist = np.linspace(0,5,20000)
    coh_data = []

    for n in range(17,19):
        L_n, L_v, H, wRC = RC_function(eps, T_EM, T_Ph, wc, alpha_Ph, alpha_em, n)
        n_RC = Occupation(wRC, T_Ph)
        rho_0 = tensor(0.5*((E+G)*(E+G).dag()), thermal_dm(n, n_RC))
        DATA_v = mesolve(H, rho_0, timelist, [L_v], expects, progress_bar=True)
        coh_data.append(DATA_v.expect[0])
    plt.figure()
    for i in coh_data:
        plt.plot(i.expect[0])
    plt.legend()


def nonsec_check(eps, H, A, N):
    """
    Plots a scatter graph with a crude representation of how dominant non-secular terms are likely.
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
