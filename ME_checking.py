from qutip import ket
from qutip import ket, mesolve, qeye, tensor, thermal_dm, destroy, steadystate
import matplotlib.pyplot as plt
import numpy as np
import UD_liouv as RC
import driving_liouv as EM

reload(RC)
reload(EM)

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

def SS_convergence_check(sigma, eps, T_EM, T_ph, wc, w0, alpha_ph, alpha_EM, expect_op='excited', time_units='cm', start_n=14, end_n=20, method='direct'):
    """
    """
    # Only for Hamiltonians of rotating wave form
    G = ket([0])
    E = ket([1])
    ss_list_s,ss_list_ns,ss_list_naive  = [],[],[] # steady states
    r_vector = E # r_vector is the ket vector on the right in the .matrix_element operation. Default is E.
    l_vector = E.dag()
    N_values = range(start_n,end_n)
    if expect_op == 'coherence':
        l_vector = G.dag()
    else:
        pass
    for n in N_values:
        L_RC, H, A_EM, A_nrwa, wRC, kappa = RC.RC_function_UD(sigma, eps, T_ph, wc, w0, alpha_ph, n)
        L_s = EM.L_vib_lindblad(H, A_EM, alpha_EM, T_EM)
        L_ns = EM.L_nonsecular(H, A_EM, alpha_EM, T_EM)
        L_naive = EM.L_EM_lindblad(eps, A_EM, alpha_EM, T_EM)
        ss_s = steadystate(H, [L_RC+L_s], method=method).ptrace(0)
        ss_ns = steadystate(H, [L_RC+L_ns], method=method).ptrace(0)
        ss_naive = steadystate(H, [L_RC+L_naive], method=method).ptrace(0)
        ss_list_s.append(ss_s.matrix_element(l_vector, r_vector))
        ss_list_ns.append(ss_ns.matrix_element(l_vector, r_vector))
        ss_list_naive.append(ss_naive.matrix_element(l_vector, r_vector))
        print "N=", n, "\n -----------------------------"
    plt.figure()
    plt.ylim(0,0.4)
    plt.plot(N_values, ss_list_s, label='secular')
    plt.plot(N_values, ss_list_ns, label='non-secular')
    plt.plot(N_values, ss_list_naive, label='naive')
    plt.legend()
    plt.ylabel("Excited state population")
    plt.xlabel("RC Hilbert space dimension")
    p_file_name = "Notes/Images/Checks/SuperPop_convergence_a{:d}_Tem{:d}_w0{:d}_eps{:d}_{}.pdf".format(int(alpha_ph), int(T_EM), int(w0), int(eps), method)
    plt.savefig(p_file_name)
    return ss_list_s,ss_list_ns,ss_list_naive, p_file_name

def sec_nonsec_agreement(ss_list_s,ss_list_ns):
    diffs = np.array(ss_list_ns) - np.array(ss_list_s)
    ss_list_diff, p_file_name = diffs, None

    return ss_list_diff, p_file_name

def plot_SS_divergences(sigma, eps, T_EM, T_ph, wc, w0, alpha_ph, alpha_EM, N_values, eps_values, expect_op='excited', time_units='cm', method='direct'):
    # Set up a loop over different system splittings
    # Calculate all of the liouvillians and steady-states for each system
    G = ket([0])
    E = ket([1])
    ss_list_s,ss_list_ns,ss_list_naive  = [],[],[] # steady states
    r_vector = E # r_vector is the ket vector on the right in the .matrix_element operation. Default is E.
    l_vector = E.dag()
    if expect_op == 'coherence':
        l_vector = G.dag()
    else:
        pass
    i = 0
    for eps in eps_values:
        N = N_values[i]
        L_RC, H, A_EM, A_nrwa, wRC, kappa = RC.RC_function_UD(sigma, eps, T_ph, wc, w0, alpha_ph, N)
        L_s = EM.L_vib_lindblad(H, A_EM, alpha_EM, T_EM)
        L_ns = EM.L_nonsecular(H, A_EM, alpha_EM, T_EM)
        L_naive = EM.L_EM_lindblad(eps, A_EM, alpha_EM, T_EM)
        ss_s = steadystate(H, [L_RC+L_s], method=method).ptrace(0)
        ss_ns = steadystate(H, [L_RC+L_ns], method=method).ptrace(0)
        ss_naive = steadystate(H, [L_RC+L_naive], method=method).ptrace(0)
        ss_list_s.append(ss_s.matrix_element(l_vector, r_vector))
        ss_list_ns.append(ss_ns.matrix_element(l_vector, r_vector))
        ss_list_naive.append(ss_naive.matrix_element(l_vector, r_vector))
        print "Splitting =", eps, "\n -----------------------------"
        i+=1
    plt.ylim(0,0.5)
    plt.plot(eps_values, ss_list_s, label='secular')
    plt.plot(eps_values, ss_list_ns, label='non-secular')
    plt.plot(eps_values, ss_list_naive, label='naive')
    plt.legend()
    plt.ylabel("Excited state population")
    plt.xlabel(r"TLS splitting $(cm^{-1})$")
    p_file_name = "Notes/Images/Checks/Pop_SS_divergence_a{:d}_Tem{:d}_w0{:d}_{}.pdf".format(int(alpha_ph), int(T_EM), int(w0), method)
    return ss_list_s,ss_list_ns,ss_list_naive, p_file_name

def nonsec_check_H(H, A, N):
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

def nonsec_check_comb(H, A, alpha, T, N):
    """
    Plots a scatter graph with a crude representation of how dominant non-secular terms are.
    Ahsan's results suggest that it is in regimes where the occupation number is so large that Gamma(1 + 2N(w0)) becomes
    comparable to w0 then we are no longer able to approximate it looks like the secular
    approximation would break down in such cases that non-secularity becomes important.
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
                    A_ij = A.matrix_element(evecs[i].dag(),evecs[j])
                    A_kl_conj = (A.dag()).matrix_element(evecs[l].dag(),evecs[k])
                    N_occ = EM.Occupation(abs(eps_kl), T,)
                    TD.append(eps_ij-eps_kl)
                    dipoles.append(EM.rate_down(abs(eps_kl), N_occ, alpha)*A_ij*A_kl_conj.real)
    return TD, dipoles

def nonsec_check_A(H, A, alpha, T, N):
    """
    Plots a scatter graph with a crude representation of how dominant non-secular terms are.
    Ahsan's results suggest that it is in regimes where the occupation number is so large that Gamma(1 + 2N(w0)) becomes
    comparable to w0 then we are no longer able to approximate it looks like the secular
    approximation would break down in such cases that non-secularity becomes important.
    """
    rates = []
    TD = []
    evals, evecs = H.eigenstates()
    for i in range(2*N):
        for j in range(2*N):
                eps_ij = evals[i]-evals[j]
                A_ij = A.matrix_element(evecs[i].dag(),evecs[j])
                N_occ = EM.Occupation(abs(eps_ij), T,)
                TD.append(eps_ij)
                rates.append(EM.rate_down(abs(eps_ij), N_occ, alpha))
    return TD, rates

if __name__ == "__main__":
    N = 25
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.

    eps = 2000. # TLS splitting

    T_EM = 6000. # Optical bath temperature
    alpha_EM = 0.3 # System-bath strength (optical)

    T_ph = 300. # Phonon bath temperature
    wc = 53. # Ind.-Boson frame phonon cutoff freq
    w0 = 300. # underdamped SD parameter omega_0
    alpha_ph = 400. # Ind.-Boson frame coupling

    #Now we build all the operators
    """
    L_RC, H, A_EM, A_nrwa, wRC, kappa= RC.RC_function_UD(sigma, eps, T_ph, wc, w0, alpha_ph, N)
    L_s = EM.L_vib_lindblad(H, A_EM, alpha_EM, T_EM)
    L_ns = EM.L_nonsecular(H, A_EM, alpha_EM, T_EM)
    L_naive = EM.L_EM_lindblad(eps, A_EM, alpha_EM, T_EM)
    ss_naive = steadystate(H, [L_RC+L_naive]).ptrace(0)
    #TD, rates  = nonsec_check_A(H, A_EM, alpha_EM, T_EM, N)
    #plt.figure()
    #plt.scatter(TD, rates)
    #plt.show()
    """
    plt.figure()
    ss_list_s,ss_list_ns,ss_list_naive, p_file_name = SS_convergence_check(sigma, eps, T_EM, T_ph, wc, w0, alpha_ph, alpha_EM, start_n = 20, end_n=35)
    #eps_values = range(1000, 2000, 50)+range(2000, 4000, 500)+range(4000, 14000, 1000)
    #N_values = [30]*len(range(1000, 2000, 50)) + [20]*len(range(2000, 4000, 500)) + [12]*len(range(4000, 14000, 1000))
    #solver_method = 'power'
    #ss_list_s,ss_list_ns,ss_list_naive, p_file_name = plot_SS_divergences(sigma, eps, T_EM, T_ph, wc, w0, alpha_ph, alpha_EM, N_values, eps_values, method=solver_method)
    print "Plot saved: ",p_file_name
    plt.savefig(p_file_name)
    plt.close()
