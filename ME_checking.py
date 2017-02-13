
from qutip import ket, basis, mesolve, qeye, tensor, thermal_dm, destroy, steadystate, Qobj
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import UD_liouv as RC
import driving_liouv as EM
from math import factorial
from scipy.special import genlaguerre
reload(RC)
reload(EM)
def J_multipolar(omega, Gamma, omega_0):
    return Gamma*(omega**3)/(2*np.pi*(omega_0**3))

def J_minimal(omega, Gamma, omega_0):
    return Gamma*omega/(2*np.pi*omega_0)

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
        ss_list_s.append((ss_s*tensor(qeye(2), destroy(N).dag()*destroy(N))).tr())
        ss_list_ns.append((ss_ns*tensor(qeye(2), destroy(N).dag()*destroy(N))).tr())
        ss_list_naive.append((ss_naive*tensor(qeye(2), destroy(N).dag()*destroy(N))).tr())
        """
        ss_list_s.append(ss_s.matrix_element(l_vector, r_vector))
        ss_list_ns.append(ss_ns.matrix_element(l_vector, r_vector))
        ss_list_naive.append(ss_naive.matrix_element(l_vector, r_vector))
        """
        print "N=", n, "\n -----------------------------"
    plt.figure()
    plt.ylim(0,0.4)
    plt.plot(N_values, ss_list_s, label='secular')
    plt.plot(N_values, ss_list_ns, label='non-secular')
    plt.plot(N_values, ss_list_naive, label='naive')
    plt.legend()
    plt.ylabel("Excited state population")
    plt.xlabel("RC Hilbert space dimension")
    p_file_name = "Notes/Images/Checks/_convergence_a{:d}_Tem{:d}_w0{:d}_eps{:d}.pdf".format(int(alpha_ph), int(T_EM), int(w0), int(eps))
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
    p_file_name = "Notes/Images/Checks/1Pop_SS_divergence_a{:d}_Tem{:d}_w0{:d}.pdf".format(int(alpha_ph), int(T_EM), int(w0))
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
def secular_approx_check(H, A, Gamma, omega_0, T, N_max):
    """
    Checks whether frequency based or number state based arguments are equivalent in the TLS (no degeneracies.)
    """
    lazy_rates = []
    rig_rates = []
    lazy_freq = []
    rig_freq = []
    evals, evecs = H.eigenstates()
    for i in range(2*N):
        for j in range(2*N):
            for k in range(2*N):
                for l in range(2*N):
                    eps_ij = evals[i]-evals[j]
                    eps_kl = evals[k]-evals[l]
                    if abs(eps_kl)>0:
                        A_ij = A.matrix_element(evecs[i].dag(),evecs[j])
                        A_kl_conj = (A.dag()).matrix_element(evecs[l].dag(),evecs[k])
                        #N_occ = EM.Occupation(abs(eps_kl), T, time_units='ps')

                        if abs(eps_ij-eps_kl)==0:
                            rig_freq.append(eps_ij-eps_kl)
                            rig_rates.append(2*np.pi*J_minimal(abs(eps_kl), Gamma, omega_0)*A_ij*A_kl_conj)
                        else:
                            pass
                        if (i==k and j==l):
                            lazy_freq.append(eps_ij-eps_kl)
                            lazy_rates.append(2*np.pi*J_minimal(abs(eps_kl), Gamma, omega_0)*A_ij*A_kl_conj)

    return lazy_rates, rig_rates, lazy_freq, rig_freq

def rates(H, A, Gamma, eps, g, omega_0, T, N_max):
    """
    """
    multipolar_rates = []
    minimal_rates = []
    frequencies = []
    frequencies_zero=[]
    forbidden = []
    sec_rates = []
    ground, excited, evals, evecs = numerical_spectrum(H, eps, g, omega_0, N_max=N_max)
    print len(evals)
    for i in range(2*N_max):
        for j in range(2*N_max):
            for k in range(2*N_max):
                for l in range(2*N_max):
                    eps_ij = evals[i]-evals[j]
                    eps_kl = evals[k]-evals[l]
                    if abs(eps_kl)>0:
                        A_ij = A.matrix_element(evecs[i].dag(),evecs[j])
                        A_kl_conj = (A.dag()).matrix_element(evecs[l].dag(),evecs[k])
                        if A_ij==0.0 or A_kl_conj==0.0:
                            forbidden.append(eps_ij-eps_kl)
                        #N_occ = EM.Occupation(abs(eps_kl), T, time_units='ps')
                        frequencies.append(eps_ij-eps_kl)
                        multipolar_rates.append(2*np.pi*J_multipolar(abs(eps_kl), Gamma, omega_0)*A_ij*A_kl_conj)
                        minimal_rates.append(2*np.pi*J_minimal(abs(eps_kl), Gamma, omega_0)*A_ij*A_kl_conj)
                        if abs(eps_ij-eps_kl)==0:
                            sec_rates.append(2*np.pi*J_minimal(abs(eps_kl), Gamma, omega_0)*A_ij*A_kl_conj)
                    else:
                        frequencies_zero.append(eps_ij-eps_kl)

    return multipolar_rates, minimal_rates, sec_rates, frequencies, frequencies_zero, forbidden


def func(x, a, b, c, d):
    return a * np.exp(-b * x) + c*x + d

def ladder_spacing(N_check=4, N_max=40):
    fig = plt.figure()
    av_spacings= []
    n_list = []
    for n in range(4,N_max):
        L_RC, H, A_EM, A_nrwa, wRC, kappa= RC.RC_function_UD(sigma, eps, T_ph, wc, w0, alpha_ph, n, time_units='cm')
        n_list.append(n)
        spectrum= H.eigenenergies()

        sorted_spectrum = sorted(spectrum)[n:]
        print len(sorted_spectrum)
        spacings = []
        #diff_check = sorted_spectrum[N_check-1]-sorted_spectrum[N_check]
        for i in range(len(sorted_spectrum)-1):
            spacings.append(abs(sorted_spectrum[i]-sorted_spectrum[i+1]))
        #plt.scatter(n, diff_check, color='y')
        av_spacings.append(np.mean(spacings))
    n_list=np.arange(len(av_spacings))
    av_spacings=np.array(av_spacings)
    #plt.plot(n_list, func(n_list, *popt), 'r-', label="Fitted Curve")
    return n_list, av_spacings

def numerical_spectrum(H, epsilon, g, omega_0, N_max=6):
    energies, vecs = H.eigenstates()
    eps_prime = epsilon - (g**2)/omega_0
    ground = []
    excited = []
    enrgs = []
    vctrs = []
    # Creates a truncated spectrum so frequencies are all linearly spaced
    for state in zip(energies, vecs):
        E_i, v_i = state[0], state[1]
        if E_i > (eps_prime -10) and E_i<eps_prime+N_max*omega_0: # -1 is to so vib grnd is included
            excited.append((E_i,v_i))
            enrgs.append(E_i)
            vctrs.append(v_i)
        elif E_i< (N_max*omega_0):
            ground.append((E_i,v_i))
            enrgs.append(E_i)
            vctrs.append(v_i)
        else:
            pass
    return sorted(ground), sorted(excited), enrgs, vctrs

def analytic_spectrum(H, epsilon, g, omega_0, N_max=6):
    # Firstly just the center
    an_freqs = []
    num_freqs = []
    eps_prime = epsilon - (g**2)/omega_0
    an_freqs_calc = []
    for n in range(N_max-2):
        an_freqs_calc.append(n*omega_0)
        an_freqs_calc.append(eps_prime+n*omega_0)
        for m in range(N_max-2):
            for p in range(N_max-2):
                for q in range(N_max-2):
                    an_freqs.append((n-m)*omega_0 - (p-q)*omega_0)
                    #an_freqs.append((m-n)*omega_0 - (q-p)*omega_0)
    energies = H.eigenenergies()
    num_energies = []
    for E_i in energies:
        if E_i > (eps_prime -20) and E_i<eps_prime+N_max*omega_0:
            num_energies.append(np.round(E_i))
        elif E_i< (N_max*omega_0):
            num_energies.append(np.round(E_i))
        else:
            pass
    for E_i in num_energies:
        for E_j in num_energies:
            if abs(E_i-E_j)<epsilon/2:
                num_freqs.append(E_i-E_j)
    #print an_freqs_calc
    #print num_energies
    return np.array(an_freqs), np.array(num_freqs)

def associated_laguerre(n, k, x):
    L = 0
    for j in range(n+1):
        L+= ((-1)**j)*(factorial(n+k)/(factorial(n-j)*factorial(k+j)*factorial(j)))*(x**j)
    return L

def wavefunction_overlap(n, m, alpha):
    overlap = 0
    if m>=n:
        L_n_k = associated_laguerre(n, (m-n), abs(alpha)**2)
        return np.sqrt(factorial(n)/factorial(m))*(alpha**(m-n))*np.exp(-0.5*abs(alpha)**2)*L_n_k
    else:
        return ValueError("Only valid for m>=n.")

def num_wavefunction_overlap(excited, n, m, N):
    overlap = 0
    if m>=n:
        return tensor(qeye(2), qeye(N)).matrix_element(tensor(G, basis(N,n)).dag(), excited[m][1])
    else:
        return ValueError("Only valid for m>=n.")

def excited_in_uncoupled(excited, N):
    """
    Transform excited manifold vectors into the uncoupled basis
    """
    uncoupled = []
    for m in range(len(excited)):
        uncoup_vec = Qobj()
        for n in range(len(excited)):
            uncoup_vec+=tensor(qeye(2), qeye(N)).matrix_element(excited[n][1].dag(), tensor(E, basis(N,m)))*excited[n][1]
        uncoupled.append(uncoup_vec)
    return uncoupled
if __name__ == "__main__":
    N = 22
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.

    eps = 500.*8.066 # TLS splitting

    T_EM = 6000. # Optical bath temperature
    #alpha_EM = 0.3 # System-bath strength (optical)
    Gamma = 6.582E-4*8.066 #inv. cm
    T_ph = 300. # Phonon bath temperature
    wc = 53. # Ind.-Boson frame phonon cutoff freq
    w0 = 300. # overdamped SD parameter omega_0
    alpha_ph = 400 # Ind.-Boson frame coupling

    #Now we build all the operators

    L_RC, H, A_EM, A_nrwa, wRC, kappa= RC.RC_function_UD(sigma, eps, T_ph, wc, w0, alpha_ph, N, time_units='cm')

    """
    TD, rates  = nonsec_check_A(H, A_EM, alpha_EM, T_EM, N)
    plt.figure()
    plt.scatter(TD, rates)
    plt.show()

    plt.figure()
    #ss_list_s,ss_list_ns,ss_list_naive, p_file_name = SS_convergence_check(sigma, eps, T_EM, T_ph, wc, w0, alpha_ph, alpha_EM, start_n = 30, end_n=45)
    eps_values = range(1000, 2000, 50)+range(2000, 4000, 500)#+range(4000, 14000, 1000)
    N_values = [30]*len(range(1000, 2000, 50)) + [20]*len(range(2000, 4000, 500)) + [12]*len(range(4000, 14000, 1000))
    solver_method = 'eigen'
    ss_list_s,ss_list_ns,ss_list_naive, p_file_name = plot_SS_divergences(sigma, eps, T_EM, T_ph, wc, w0, alpha_ph, alpha_EM, N_values, eps_values, method=solver_method)
    print "Plot saved: ",p_file_name
    plt.savefig(p_file_name)
    plt.close()
    """

    multipolar_rates, minimal_rates, sec_rates, frequencies, frequencies_zero, forbidden = rates(H, A_EM, Gamma, eps, kappa, wRC, T_EM, 6)
    plt.figure()
    plt.scatter(frequencies, minimal_rates, color='r')
    plt.scatter(np.zeros(len(sec_rates)), sec_rates, color='b')
    plt.scatter(forbidden, np.zeros(len(forbidden)),color='y')
    plt.axvline(eps-(kappa**2/wRC))
    plt.axvline(-eps+(kappa**2/wRC))
    plt.axvline(2*(eps-(kappa**2/wRC)))
    plt.axvline(-2*(eps-(kappa**2/wRC)))
    plt.title("Vibronic Transition Rates at "r"$\alpha_{ph}=400cm^{-1}$")

    """
    plt.figure()
    lazy_rates, rig_rates, lazy_freq, rig_freq = secular_approx_check(H, A_EM, Gamma, eps, T_EM, N)
    plt.scatter(lazy_rates, lazy_freq, color='b')
    plt.scatter(rig_rates, rig_freq, color='r')

    n_list, av_spacings = ladder_spacing()
    plt.scatter(n_list, av_spacings)
    popt, pcov = curve_fit(func, n_list, av_spacings, p0=[300,0.3,0.001,500])
    plt.plot(np.arange(40), func(np.arange(40), *popt))
    plt.show()

    an_freqs, num_freqs = analytic_spectrum(H, eps, kappa, wRC)
    #print sorted(an_freqs)
    #print sorted(num_freqs)
    print wavefunction_overlap(0, 0, kappa/wRC)
    print wavefunction_overlap(0, 1, kappa/wRC)
    print wavefunction_overlap(0, 2, kappa/wRC)
    A = np.zeros((10, 10))
    for n in range(10):
        for m in range(10):
            if m>n:
                A[n][m] = abs(wavefunction_overlap(n, m, 5.))**2
            else:
                A[n][m] = abs(wavefunction_overlap(m, n, 5.))**2
    plt.imshow(A)
    plt.colorbar()
    plt.show()


    ground, excited, evals, evecs = numerical_spectrum(H, eps, kappa, wRC)
    #print num_wavefunction_overlap(excited, 0, 3, N)
    uncoupled_excited= excited_in_uncoupled(excited, N)
    for i in range(len(uncoupled_excited)):
        for j in range(len(uncoupled_excited)):
            print uncoupled_excited[i].dag()*ground[j][1]
    """
