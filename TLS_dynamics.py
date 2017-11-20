from qutip import ket, mesolve, qeye, tensor, thermal_dm, destroy, steadystate
import matplotlib.pyplot as plt
import numpy as np

import UD_liouv as RC
import driving_liouv as EM
import ME_checking as check
import exact_IB as exact
import scipy as sp
import OD_weak_coupling as WC
from utils import J_overdamped
reload(RC)
reload(EM)
reload(check)
reload(exact)
plt.style.use('ggplot')

if __name__ == "__main__":
    """
    Define all system and environment  parameters
    """
    N = 15
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
    eps = 0.1*8065.5 # TLS splitting
    #eps = 2.*1519.3 # ps
    T_EM = 6000. # Optical bath temperature
    H_S = eps*E*E.dag()
    #alpha_EM = 0.3 # System-bath strength (optical)
    Gamma_EM = 1.*5.309 #bare decay of electronic transition from inv. ps to in inv. cm
    #Gamma_EM = 6.582E-7*1519.3
    T_ph = 300. # Phonon bath temperature
    wc = 53. # Ind.-Boson frame phonon cutoff freq
    #wc = 53.*0.188
    w0 = 500. # underdamped SD parameter omega_0
    #w0 = 200.*0.188
    Gamma = w0**2/wc
    beta = 1/(0.695*T_ph)
    alpha_ph = 1/np.pi #0.05*eps/np.pi# Ind.-Boson frame coupling
    J = EM.J_minimal
    mu = 1.
    print "eps={:d}, T_EM={:d}, w_0={:d}, alpha_ph={:d}".format(int(eps), int(T_EM), int(w0), int(alpha_ph))
    """
    Now we build the sys-RC Hamiltonian and residual bath Liouvillian as well as generate mapped parameters
    """
    L_RC, H_RC, A_EM, A_nrwa, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps, T_ph, wc, w0, alpha_ph, N)
    L_ph = WC.L_phonon(J_overdamped, alpha_ph, beta, wc)
    # electromagnetic bath liouvillians

    #L_nrwa = EM.L_nonrwa(H, A_nrwa, eps, Gamma, T_EM) # Ignore this for now as it just doesn't work
    #L_ns = EM.L_nonsecular(H, A_EM, eps, Gamma_EM, T_EM, J=J)
    #L_s = EM.L_vib_lindblad(H, A_EM, eps, Gamma_EM, T_EM, J=J)
    #L_naive = EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J)

    # Set up the initial density matrix
    n_RC = EM.Occupation(wRC, T_ph)
    #initial_sys = G*G.dag()
    #initial_sys = E*E.dag()
    #initial_sys = 0.5*(G+E)*(E.dag()+G.dag())
    initial_sys = 0.5*(E+G)*(E+G).dag()

    rho_0 = tensor(initial_sys, thermal_dm(N, n_RC))
    #print (rho_0*tensor(qeye(2), destroy(N).dag()*destroy(N))).tr()

    # Expectation values and time increments needed to calculate the dynamics
    expects_wc = [G*G.dag(), E*G.dag()]
    expects = [tensor(G*G.dag(), qeye(N)), tensor(E*G.dag(), qeye(N)), tensor(qeye(2), destroy(N).dag()*destroy(N)), tensor(qeye(2), destroy(N).dag()+destroy(N))]
    a, b, points = 0, 0.1, 1000
    timelist = np.linspace(a, b, 8000)
    #nonsec_check(eps, H, A_em, N) # Plots a scatter graph representation of non-secularity. Could use nrwa instead.

    # Calculate dynamics

    DATA_wc = mesolve(H_S, initial_sys, timelist, [L_ph], expects_wc, progress_bar=True)
    #DATA_nrwa = mesolve(H, rho_0, timelist, [L_RC+L_nrwa], expects, progress_bar=True)

    DATA_rc = mesolve(H_RC, rho_0, timelist, [L_RC], expects, progress_bar=True)
    #T, rho_01 = exact.exact_dynamics(eps, alpha_ph, wc, w0, Gamma, beta, initial_sys, a, b, points)
    print "Exact dynamics calculated"
    with open('DATA/Exactdatawc.dat', 'rw+') as f:
        dat = f.read()
        times, points = zip(*[(float(i.split('\t')[0]), float(i.split('\t')[1])) for i in dat.split("\n")])
        plt.plot(times,points, label='mat', linestyle='dashed')
    #plt.plot(T, np.array(rho_01), label='exact')
    plt.plot(timelist, DATA_wc.expect[1], label='wc')
    plt.plot(timelist, DATA_rc.expect[1], label='sc', linestyle='dashed')
    plt.legend()
    plt.show()
    '''
    C = mu**2*(DATA.expect[1])
    plt.plot(timelist, C, label='RC')
    plt.plot(T, np.array(rho_01), label='exact')
    plt.legend()
    plt.figure()
    m = np.fft.fft(C)
    n = C.size
    freq = np.fft.fftfreq(n, 0.1)
    plt.scatter(freq, m)
    plt.show()
    '''
    #DATA_s = mesolve(H, rho_0, timelist, [L_RC+L_s], expects, progress_bar=True)
    #DATA_naive = mesolve(H, rho_0, timelist, [L_RC+L_naive], expects, progress_bar=True)

    #timelist*=0.188 # Convert time into ps

    #fig = plt.figure(figsize=(12, 6))
    #ax1 = fig.add_subplot(121)
    #ax2 = fig.add_subplot(122)
    #plot_dynamics(ax1)
    #plot_coherences(ax2)

    #p_file_name = "Notes/Images/Dynamics/Pop_min_a{:d}_N{:d}_Tem{:d}_w0{:d}_eps{:d}.pdf".format(int(alpha_ph), int(N), int(T_EM), int(w0), int(eps))
    #plt.savefig(p_file_name)
    #print "Figure saved: ", p_file_name
    """
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    plot_RC_pop(ax1)
    plot_RC_disp(ax2)
    p_file_name = "Notes/Images/Phonons/Pop_min_a{:d}_N{:d}_Tem{:d}_w0{:d}_eps{:d}.pdf".format(int(alpha_ph), int(N), int(T_EM), int(w0), int(eps))
    plt.savefig(p_file_name)
    print "Figure saved: ", p_file_name

    plt.figure()
    x = np.arange(0,400)
    plt.plot(x, RC.J_UD_SB(x, alpha_ph, w0, Gamma))



    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plot_RC_pop(ax1)
    J = EM.J_multipolar
    L_ns = EM.L_nonsecular(H, A_EM, eps, Gamma_EM, T_EM, J=J, time_units=time_units)
    #L_s = EM.L_vib_lindblad(H, A_EM, eps, Gamma_EM, T_EM, J=J, time_units=time_units)
    #L_naive = EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J, time_units=time_units)
    DATA_ns = mesolve(H, rho_0, timelist, [L_RC+L_ns], expects, progress_bar=True)
    plot_RC_pop(ax1, prefix='multi ')
    #DATA_s = mesolve(H, rho_0, timelist, [L_RC+L_s], expects, progress_bar=True)
    #DATA_naive = mesolve(H, rho_0, timelist, [L_RC+L_naive], expects, progress_bar=True)
    #plot_dynamics_spec(DATA_ns, DATA_s, DATA_naive, timelist)

    #np.savetxt('DATA/Dynamics/DATA_ns.txt', np.array([1- DATA_ns.expect[0], timelist]), delimiter = ',', newline= '\n')
    plt.show()
    """
