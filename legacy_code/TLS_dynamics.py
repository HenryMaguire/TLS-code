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
import imp
imp.reload(RC)
imp.reload(EM)
imp.reload(check)
imp.reload(exact)
plt.style.use('ggplot')
def steadystate_coupling_dependence(couplings, T_ph, eps, Gamma, w0, T_EM, Gamma_EM, overdamped=True):
    plt.figure()
    N = 30
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
    J = EM.J_minimal
    wc = []
    rc_nrwa = []
    rc_ns = []
    rc_s = []
    rc_naive = []
    i = 1
    for alpha_ph in couplings:
        alpha_ph = alpha_ph*eps
        print(i)
        i+=1
        L_RC, H_RC, A_EM, A_nrwa, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps, T_ph, Gamma, w0, alpha_ph, N)
        L_wc = WC.L_phonon(alpha_ph, beta, Gamma, w0, overdamped=overdamped)
        L_wc_EM = EM.L_EM_lindblad(eps, sigma, Gamma_EM, T_EM, J=J)
        # electromagnetic bath liouvillians
        L_nrwa = EM.L_non_rwa(H_RC, A_nrwa, eps, Gamma_EM, T_EM, J) # Ignore this for now as it just doesn't work
        L_ns = EM.L_nonsecular(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
        L_s = EM.L_vib_lindblad(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
        #L_naive = EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J)


        wc.append(steadystate(H_S, [L_wc+L_wc_EM]).matrix_element(E.dag(), E))
        #print "NON-ROTATING-WAVE"
        rc_nrwa.append(steadystate(H_RC, [L_RC+L_nrwa]).ptrace(0).matrix_element(E.dag(), E))
        #print "NON-SECULAR"
        rc_ns.append(steadystate(H_RC, [L_RC+L_ns]).ptrace(0).matrix_element(E.dag(), E))
        #print "SECULAR"
        rc_s.append(steadystate(H_RC, [L_RC+L_s]).ptrace(0).matrix_element(E.dag(), E))
        #print "NAIVE"
        #rc_naive.append(steadystate(H_RC, [L_RC+L_naive]).ptrace(0).matrix_element(E.dag(), E))
    plt.plot(couplings, wc, label='wc & naive')
    plt.plot(couplings, rc_nrwa, label='nrwa')
    plt.plot(couplings, rc_ns, label='ns')
    plt.plot(couplings, rc_s,label='s')
    plt.legend()
    plt.grid()
    plt.ylabel("Excited state population")
    plt.xlabel(r"Phonon Coupling Strength: units of $\epsilon$")
    plt.show()

def steadystate_bias_dependence(biases, T_ph, alpha_ph, Gamma, w0, cutoff, T_EM, Gamma_EM, overdamped=True):
    N = 15
    plt.figure()
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
    J = EM.J_minimal
    wc = []
    rc_nrwa = []
    rc_ns = []
    rc_s = []
    rc_naive = []
    i = 1
    print(w0, wc)
    for eps in biases:
        if overdamped:
            optical_cutoff = 20000.
            Gamma = w0**2/cutoff
        else:
            w0 = eps+0.5*np.pi*alpha_ph # underdamped SD parameter omega_0
        print(i)
        i+=1
        L_RC, H_RC, A_EM, A_nrwa, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps, T_ph, Gamma, w0, alpha_ph, N)
        L_wc = WC.L_phonon(alpha_ph, beta, Gamma, w0, overdamped=overdamped)
        L_wc_EM = EM.L_EM_lindblad(eps, sigma, Gamma_EM, T_EM, J=J)
        wc.append(steadystate(H_S, [L_wc+L_wc_EM]).matrix_element(E.dag(), E))
        del L_wc, L_wc_EM
        # electromagnetic bath liouvillians
        L_nrwa = EM.L_non_rwa(H_RC, A_nrwa, eps, Gamma_EM, T_EM, J, principal=False) # Ignore this for now as it just doesn't work
        L_ns = EM.L_nonsecular(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
        L_s = EM.L_vib_lindblad(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
        #L_naive = EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J)


        #print "NON-ROTATING-WAVE"
        rc_nrwa.append(steadystate(H_RC, [L_RC+L_nrwa]).ptrace(0).matrix_element(E.dag(), E))
        #print "NON-SECULAR"
        rc_ns.append(steadystate(H_RC, [L_RC+L_ns]).ptrace(0).matrix_element(E.dag(), E))
        #print "SECULAR"
        rc_s.append(steadystate(H_RC, [L_RC+L_s]).ptrace(0).matrix_element(E.dag(), E))
        #print "NAIVE"
        #rc_naive.append(steadystate(H_RC, [L_RC+L_naive]).ptrace(0).matrix_element(E.dag(), E))
    plt.plot(biases, wc, label='wc & naive')
    plt.plot(biases, rc_nrwa, label='nrwa')
    plt.plot(biases, rc_ns, label='ns')
    plt.plot(biases, rc_s,label='s')
    plt.title("Steady State Population (UD)")
    plt.legend()
    plt.grid()
    plt.ylabel("Excited state population")
    plt.xlabel(r"TLS Splitting $cm^{-1}$")
    plt.show()



def steadystate_mode_dependence(mode_freqs, eps, T_ph, alpha_ph, Gamma, T_EM, Gamma_EM):
    N = 15
    plt.figure()
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
    J = EM.J_minimal
    wc = []
    rc_nrwa = []
    rc_ns = []
    rc_s = []
    rc_naive = []
    i = 1
    for w0 in mode_freqs:
        print(i)
        i+=1
        L_RC, H_RC, A_EM, A_nrwa, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps, T_ph, Gamma, w0, alpha_ph, N)
        L_wc = WC.L_phonon(alpha_ph, beta, Gamma, w0, overdamped=False)
        L_wc_EM = EM.L_EM_lindblad(eps, sigma, Gamma_EM, T_EM, J=J)
        wc.append(steadystate(H_S, [L_wc+L_wc_EM]).matrix_element(E.dag(), E))
        del L_wc, L_wc_EM
        # electromagnetic bath liouvillians
        L_nrwa = EM.L_non_rwa(H_RC, A_nrwa, eps, Gamma_EM, T_EM, J, principal=False) # Ignore this for now as it just doesn't work
        L_ns = EM.L_nonsecular(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
        L_s = EM.L_vib_lindblad(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
        #L_naive = EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J)


        #print "NON-ROTATING-WAVE"
        rc_nrwa.append(steadystate(H_RC, [L_RC+L_nrwa]).ptrace(0).matrix_element(E.dag(), E))
        #print "NON-SECULAR"
        rc_ns.append(steadystate(H_RC, [L_RC+L_ns]).ptrace(0).matrix_element(E.dag(), E))
        #print "SECULAR"
        rc_s.append(steadystate(H_RC, [L_RC+L_s]).ptrace(0).matrix_element(E.dag(), E))
        #print "NAIVE"
        #rc_naive.append(steadystate(H_RC, [L_RC+L_naive]).ptrace(0).matrix_element(E.dag(), E))
    plt.plot(mode_freqs, wc, label='wc & naive')
    plt.plot(mode_freqs, rc_nrwa, label='nrwa')
    plt.plot(mode_freqs, rc_ns, label='ns')
    plt.plot(mode_freqs, rc_s,label='s')
    plt.title("Steady State Population (UD)")
    plt.legend()
    plt.grid()
    plt.ylabel("Excited state population")
    plt.xlabel(r"Peak position $cm^{-1}$")
    plt.show()
if __name__ == "__main__":
    """
    Define all system and environment  parameters
    """

    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
    ## ~~ Parameters ~~
    N = 30
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
    w0 = eps/0.9 #-100 # underdamped SD parameter omega_0
    optical_cutoff = 20000.
    if overdamped:
        Gamma = w0**2/wc
    #biases = np.linspace(100, 0.3*8065.5, 30)
    #steadystate_bias_dependence(biases, T_ph, alpha_ph, Gamma, w0, wc, T_EM, Gamma_EM, overdamped=overdamped)
    #couplings = np.linspace(0, 300, 40)/np.pi
    #steadystate_coupling_dependence(couplings, T_ph, eps, Gamma, w0, T_EM, Gamma_EM)
    #mode_freqs = np.linspace(100, 2*eps, 60)
    #steadystate_mode_dependence(mode_freqs, eps, T_ph, alpha_ph, Gamma, T_EM, Gamma_EM)
    #couplings = np.linspace(0, 10, 60)
    #steadystate_coupling_dependence(couplings, T_ph, eps, Gamma, w0, T_EM, Gamma_EM)

    if True:
        """
        w = np.linspace(0,w0+2*wc, 1000)
        plt.figure()
        plt.plot(w, J_underdamped(w, alpha_ph, Gamma, w0))
        plt.axvline(eps-shift, label='Splitting',color='k')
        plt.ylabel("Coupling Strength")
        plt.xlabel(r"Frequency $cm^{-1}$")
        plt.title("Phonon Spectral Density")
        plt.legend()
        J = EM.J_minimal
        J = lambda x,y,z : EM.J_minimal_hard(x,y,z, optical_cutoff)
        print "eps={:d}, T_EM={:d}, w_0={:d}, alpha_ph={:d}".format(int(eps), int(T_EM), int(w0), int(alpha_ph))
        """
        """
        Now we build the sys-RC Hamiltonian and residual bath Liouvillian as well as generate mapped parameters
        """
        #L_RC, H_RC, A_EM, A_nrwa, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps, T_ph, Gamma, w0, alpha_ph, N)
        #L_wc = WC.L_phonon(alpha_ph, beta, Gamma, w0, overdamped=overdamped)
        if not phonon_only:
            L_wc_EM = EM.L_EM_lindblad(eps, sigma, Gamma_EM, T_EM, J=J)
            # electromagnetic bath liouvillians
            #L_nrwa2 = EM.L_non_rwa(H_RC, A_nrwa, eps, Gamma_EM, T_EM, J, principal=True) # Ignore this for now as it just doesn't work
            L_nrwa = EM.L_non_rwa(H_RC, A_nrwa, eps, Gamma_EM, T_EM, J, principal=False)
            L_ns = EM.L_nonsecular(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
            L_s = EM.L_vib_lindblad(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
            L_naive = EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J)
        # Set up the initial density matrix
        n_RC = EM.Occupation(w0, T_ph)
        #initial_sys = G*G.dag()
        #initial_sys = E*E.dag()
        #initial_sys = 0.5*(G+E)*(E.dag()+G.dag())
        initial_sys = 0.5*(E+G)*(E+G).dag()

        rho_0 = tensor(initial_sys, thermal_dm(N, n_RC))
        #print (rho_0*tensor(qeye(2), destroy(N).dag()*destroy(N))).tr()

        # Expectation values and time increments needed to calculate the dynamics
        expects_wc = [G*G.dag(), E*G.dag()]
        expects = [tensor(G*G.dag(), qeye(N)), tensor(E*G.dag(), qeye(N)), tensor(E*E.dag(), qeye(N)), tensor(qeye(2), destroy(N).dag()*destroy(N)), tensor(qeye(2), destroy(N).dag()+destroy(N))]
        timelist = np.linspace(0, 0.07, 250)
        #nonsec_check(eps, H, A_em, N) # Plots a scatter graph representation of non-secularity. Could use nrwa instead.

        # Calculate dynamics
        opts = qt.Options(nsteps=15000)
        if not phonon_only:
            print("WEAK_COUPLING")
            #DATA_wc = mesolve(H_S, initial_sys, timelist, [L_wc+L_wc_EM], expects_wc, progress_bar=True)
            print("NON-ROTATING-WAVE")
            DATA_nrwa = mesolve(H_RC, rho_0, timelist, [L_RC+L_nrwa1], expects, progress_bar=True, options=opts)
            #DATA_nrwa2 = mesolve(H_RC, rho_0, timelist, [L_RC+L_nrwa2], expects, progress_bar=True, options=opts)
            print("NON-SECULAR")
            #DATA_ns = mesolve(H_RC, rho_0, timelist, [L_RC+L_ns], expects, progress_bar=True)
            print("SECULAR")
            #DATA_s = mesolve(H_RC, rho_0, timelist, [L_RC+L_s], expects, progress_bar=True)
            print("NAIVE")
            #DATA_naive = mesolve(H_RC, rho_0, timelist, [L_RC+L_naive], expects, progress_bar=True)
            plt.figure()
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            #plt.plot(timelist, DATA_wc.expect[1], label='wc', color=colors[0])
            #plt.plot(timelist, DATA_wc.expect[1].imag, color=colors[0], ls='dashed')
            plt.plot(timelist, DATA_nrwa.expect[1], label='nrwa', color=colors[4])
            #plt.plot(timelist, DATA_nrwa2.expect[1], label='nrwa', color=colors[5])
            #plt.plot(timelist, DATA_nrwa.expect[1].imag, linestyle='dashed', color=colors[4])
            #plt.plot(timelist, DATA_ns.expect[1], label='ns', color=colors[1])
            #plt.plot(timelist, DATA_ns.expect[1].imag, color=colors[1], ls='dashed')

            #plt.plot(timelist, DATA_s.expect[1], label='sec', color=colors[2])
            #plt.plot(timelist, DATA_s.expect[1].imag, label='sec', color=colors[2], linestyle='dashed')

            #plt.plot(timelist, DATA_naive.expect[1], label='naive', color=colors[3])
            #plt.plot(timelist, DATA_naive.expect[1].imag, color=colors[3], linestyle='dashed')
            plt.legend()
            plt.ylabel("Coherence")
            plt.xlabel("Time")


            plt.figure()
            #plt.plot(timelist, DATA_wc.expect[0], label='wc')
            plt.plot(timelist, DATA_nrwa.expect[0], label='nrwa')
            #plt.plot(timelist, DATA_nrwa2.expect[0], label='nrwa')
            #plt.plot(timelist, DATA_ns.expect[0], label='ns')
            #plt.plot(timelist, DATA_s.expect[0], label='sec', linestyle='dashed')
            #plt.plot(timelist, DATA_naive.expect[0], label='naive')
            plt.legend()
            plt.show()
        elif False:
            rho_01 = exact.exact_dynamics(eps-shift, alpha_ph, wc, w0, Gamma, beta, initial_sys, timelist, overdamped=False)
            #print "Exact dynamics calculated"

            #DATA_wc = mesolve(H_S, initial_sys, timelist, [L_wc], expects_wc, progress_bar=True)
            DATA_sc = mesolve(H_RC, rho_0, timelist, [L_RC], expects, progress_bar=True)
            plt.figure()
            plt.plot(timelist, DATA_wc.expect[1], label='wc')
            plt.plot(timelist, DATA_sc.expect[1], label='RC')
            plt.scatter(timelist, np.array(rho_01), label='exact')
            plt.ylabel("Excited State Population")
            plt.legend()

            import frequency_analysis as fa
            imp.reload(fa)
            plt.figure()
            fa.plot_frequency(DATA_wc, timelist, label='WC', QOBJ=True)
            fa.plot_frequency(DATA_sc, timelist, label='SC', QOBJ=True)
            fa.plot_frequency(np.array(rho_01), timelist, label='Exact', QOBJ=False)#
            plt.legend()
            plt.show()
        else:
            exact_sol = []
            rc_sol = []
            prop_of_eps = np.linspace(0.001,0.12,5)
            for prop in prop_of_eps:
                alpha = prop*eps
                shift = 0.5*np.pi*alpha
                rho_01 = np.array(exact.exact_dynamics(eps-shift, alpha, wc, w0, Gamma, beta, initial_sys, timelist, overdamped=overdamped))
                #print "Exact dynamics calculated"

                #DATA_wc = mesolve(H_S, initial_sys, timelist, [L_wc], expects_wc, progress_bar=True)
                L_RC, H_RC, A_EM, A_nrwa, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps, T_ph, Gamma, w0, alpha, N)
                DATA_sc = mesolve(H_RC, rho_0, timelist, [L_RC], expects, progress_bar=True)
                exact_sol.append(rho_01.real)
                rc_sol.append(DATA_sc.expect[1].real)
            plt.figure()
            for i, sols in enumerate(zip(exact_sol, rc_sol)):
                plt.figure()
                plt.title("prop_of_eps: {}".format(prop_of_eps[i]))
                plt.plot(timelist, sols[0], label="exact")
                plt.plot(timelist, sols[1], label="rc")
                plt.legend(loc='right')
            plt.show()
        """
        with open('DATA/Exactdatasc.dat', 'rw+') as f:
            dat = f.read()
            times, points = zip(*[(float(i.split('\t')[0]), float(i.split('\t')[1])) for i in dat.split("\n")])
            plt.plot(times,points, label='mat', linestyle='dashed')
        """
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
