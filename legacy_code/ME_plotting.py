from qutip import ket, mesolve, qeye, tensor, thermal_dm, destroy, steadystate
import qutip as qt
import matplotlib.pyplot as plt
import numpy as np
import UD_liouv as RC
import driving_liouv as EM
import ME_checking as check
import scipy as sp
reload(RC)
reload(EM)
reload(check)
plt.style.use('ggplot')


def plot_RC_pop(ax, prefix='min. '):

    #ax.title(r"$\omega_0=$""%i"r"$cm^{-1}$, $\alpha_{ph}=$""%f"r"$cm^{-1}$, $T_{EM}=$""%i K" %(w0, alpha_ph, T_EM))
    #ax.plot(timelist, 1-DATA_nrwa.expect[0], label='nrwa', color='y')
    ax.plot(timelist, DATA_ns.expect[2].real, label=prefix+'Non-secular', color='g')
    #ax.plot(timelist, DATA_s.expect[2].real, label=prefix+'Vib. Lindblad', color='b')
    #ax.plot(timelist, DATA_naive.expect[2].real, label=prefix+'Simple Lindblad', color='r')
    ax.set_ylabel("Reaction-Coordinate "r"$\langle n\rangle$")
    ax.set_xlabel("Time (ps)")
    ax.legend()
<<<<<<< HEAD

def plot_RC_disp(ax, prefix='min. '):

    #ax.title(r"$\omega_0=$""%i"r"$cm^{-1}$, $\alpha_{ph}=$""%f"r"$cm^{-1}$, $T_{EM}=$""%i K" %(w0, alpha_ph, T_EM))
    #ax.plot(timelist, 1-DATA_nrwa.expect[0], label='nrwa', color='y')
    ax.plot(timelist, DATA_ns.expect[3].real, label='Non-secular', color='g')
    ax.plot(timelist, DATA_s.expect[3].real, label='Vib. Lindblad', color='b')
    ax.plot(timelist, DATA_naive.expect[3].real, label='Simple Lindblad', color='r')
    ax.set_ylabel("Reaction-Coordinate displacement")
    ax.set_xlabel("Time (ps)")
    ax.legend()

def plot_dynamics(ax, prefix='min. '):

    #if T_EM>0.0: # No need to plot SS for T=0
    ss_ns = steadystate(H, [L_RC+L_ns]).ptrace(0)
    #ss_v = steadystate(H, [L_RC+L_s]).ptrace(0)
    #ss_n = steadystate(H, [L_RC+L_naive]).ptrace(0)
    ss_g_ns = ss_ns.matrix_element(G.dag(), G)
    #ss_g_v = ss_v.matrix_element(G.dag(), G)
    #ss_g_n = ss_n.matrix_element(G.dag(), G)
    #ax.axhline(1-ss_g_v.real, color='b', ls='--')
    ax.axhline(1-ss_g_ns.real, color='g', ls='--')
    #ax.axhline(1-ss_g_n.real, color='r', ls='--')

    #ax.title(r"$\omega_0=$""%i"r"$cm^{-1}$, $\alpha_{ph}=$""%f"r"$cm^{-1}$, $T_{EM}=$""%i K" %(w0, alpha_ph, T_EM))
    #ax.plot(timelist, 1-DATA_nrwa.expect[0], label='nrwa', color='y')
=======

def plot_RC_disp(ax, prefix='min. '):

    #ax.title(r"$\omega_0=$""%i"r"$cm^{-1}$, $\alpha_{ph}=$""%f"r"$cm^{-1}$, $T_{EM}=$""%i K" %(w0, alpha_ph, T_EM))
    #ax.plot(timelist, 1-DATA_nrwa.expect[0], label='nrwa', color='y')
    ax.plot(timelist, DATA_ns.expect[3].real, label='Non-secular', color='g')
    ax.plot(timelist, DATA_s.expect[3].real, label='Vib. Lindblad', color='b')
    ax.plot(timelist, DATA_naive.expect[3].real, label='Simple Lindblad', color='r')
    ax.set_ylabel("Reaction-Coordinate displacement")
    ax.set_xlabel("Time (ps)")
    ax.legend()

def plot_dynamics(ax, prefix='min. '):

    #if T_EM>0.0: # No need to plot SS for T=0
    ss_ns = steadystate(H, [L_RC+L_ns]).ptrace(0)
    #ss_v = steadystate(H, [L_RC+L_s]).ptrace(0)
    #ss_n = steadystate(H, [L_RC+L_naive]).ptrace(0)
    ss_g_ns = ss_ns.matrix_element(G.dag(), G)
    #ss_g_v = ss_v.matrix_element(G.dag(), G)
    #ss_g_n = ss_n.matrix_element(G.dag(), G)
    #ax.axhline(1-ss_g_v.real, color='b', ls='--')
    ax.axhline(1-ss_g_ns.real, color='g', ls='--')
    #ax.axhline(1-ss_g_n.real, color='r', ls='--')

    #ax.title(r"$\omega_0=$""%i"r"$cm^{-1}$, $\alpha_{ph}=$""%f"r"$cm^{-1}$, $T_{EM}=$""%i K" %(w0, alpha_ph, T_EM))
    #ax.plot(timelist, 1-DATA_nrwa.expect[0], label='nrwa', color='y')
>>>>>>> 2017-overhaul-branch
    ax.plot(timelist, 1-DATA_ns.expect[0].real, label=prefix+'Non-secular', color='g')
    #ax.plot(timelist, 1-DATA_s.expect[0].real, label=prefix+'Vib. Lindblad', color='b')
    #ax.plot(timelist, 1-DATA_naive.expect[0].real, label=prefix+'Simple Lindblad', color='r')
    ax.set_ylabel("Excited state population")
    ax.set_xlabel("Time (ps)")
    ax.legend()
    #p_file_name = "Notes/Images/Dynamics/Pop_a{:d}_Tph{:d}_Tem{:d}_w0{:d}.pdf".format(int(alpha_ph), int(T_ph), int(T_EM), int(w0))
    #pl.savefig(p_file_name)
    #plt.close()

def plot_coherences(ax):
    #plt.title(r"$\alpha_{ph}=$""%i"r"$cm^{-1}$, $T_{EM}=$""%i K" %(alpha_ph, T_EM))
    #ax.plot(timelist, DATA_s.expect[1].real, label='Vib. Lindblad', color='b')
    ax.plot(timelist, DATA_ns.expect[1].real, label='Non-secular', color='g',alpha=0.7)
    #ax.plot(timelist, DATA_naive.expect[1].real, label='Simple Lindblad', color='r',alpha=0.4)
    ax.legend()
    ax.set_ylabel("Coherence")
    ax.set_xlabel("Time (ps)")
    file_name = "Notes/Images/Dynamics/Coh_a{:d}_Tph{:d}_Tem{:d}_w0{:d}.pdf".format(int(alpha_ph), int(T_ph), int(T_EM), int(w0))
    #plt.savefig(file_name)
    #print "File saved: ", file_name
    #plt.close()

def plot_dynamics_spec(DAT_ns, DAT_s, DAT_n,  t):
    DAT_list, lab_list = [DAT_ns, DAT_s, DAT_n,], [r'$|S_ns|$', r'$|S_s|$', r'$|S_n|$']
    plt.figure()
    for DAT, lab in zip(DAT_list, lab_list):
        dyn = 1-DAT.expect[1]
        ss = dyn[-1]
        gg = dyn-ss
        spec = np.fft.fft(gg)
        freq = np.fft.fftfreq(t.shape[-1])
        plt.plot(freq, abs(spec), label=lab)
    #plt.plot(freq, abs(spec), linestyle='dotted')
    plt.title("Frequency spectrum of vibronic TLS coherence")
    plt.legend()
    plt.ylabel("Magnitude" )
    plt.xlabel(r"Frequency- $\epsilon$ ($cm^{-1}$)")
    p_file_name = "Notes/Images/Spectra/Coh_spec_a{:d}_Tph{:d}_Tem{:d}_w0{:d}.pdf".format(int(alpha_ph), int(T_ph), int(T_EM), int(w0))
    #plt.xlim(-0.2, 0.5)
    plt.savefig(p_file_name)
    d_file_name = "DATA/Spectra/Coh_spec_a{:d}_Tph{:d}_TEM{:d}_w0{:d}.txt".format(int(alpha_ph), int(T_ph), int(T_EM), int(w0))
    np.savetxt(d_file_name, np.array([spec, freq]), delimiter = ',', newline= '\n')
    plt.close()

def plot_manifolds(ax, H):
    eigs = H.eigenenergies()
    p_title = "Vibronic manifolds: " + r"$\alpha_{ph}$="+"{:d}".format(int(alpha_ph))+ "$\epsilon$={:d}, $\omega_0=${:d}".format(int(eps), int(w0))
    #ax.title(p_title)
    plt.ylabel(r"Energy ($cm^{-1}$)")
    j = 0
    for i in eigs:
        col = 'b'
        """
        if j<H.shape[0]/2:
            col = 'b'
        else:
            col = 'r'
        """
        ax.axhline(i, color = col)
        j+=1
    #p_file_name = "Notes/Images/Spectra/Manifolds_a{:d}_Tph{:d}_Tem{:d}_w0{:d}_eps{:d}.pdf".format(int(alpha_ph), int(T_ph), int(T_EM), int(w0), int(eps))
    #plt.savefig(p_file_name)
    #print "Plot saved: ", p_file_name
    #plt.show()


def plot_nonsec(ax, TD, dipoles):
    ax.scatter(TD, dipoles)
    ax.grid()
    #ax.title(r"Non-secularity check: $\epsilon=$""%i"%eps)
    ax.set_xlabel(r"$\varphi_{ij}-\varphi_{kl}$")
    ax.set_ylabel(r"Down rate weighted by dipole of transitions $A_{ij}A_{kl}$")
<<<<<<< HEAD




if __name__ == "__main__":
<<<<<<< HEAD

    N = 8
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.

    eps = 6000. # TLS splitting

    """
    Define all system and environment  parameters
    """
    N = 10
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
=======
    """
    Define all system and environment  parameters
    """
    N = 10
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
>>>>>>> 2017-overhaul-branch
    eps = 1.1*8065.5 # TLS splitting
    #eps = 2.*1519.3 # ps
    T_EM = 6000. # Optical bath temperature
    #alpha_EM = 0.3 # System-bath strength (optical)
    Gamma_EM = 1.*5.309 #bare decay of electronic transition from inv. ps to in inv. cm
    #Gamma_EM = 6.582E-7*1519.3
    T_ph = 300. # Phonon bath temperature
    wc = 53. # Ind.-Boson frame phonon cutoff freq
    #wc = 53.*0.188
    w0 = 300. # underdamped SD parameter omega_0
<<<<<<< HEAD

    alpha_ph = 200. # Ind.-Boson frame coupling

    #w0 = 200.*0.188
    alpha_ph = 700 #0.05*eps/np.pi# Ind.-Boson frame coupling
    J = EM.J_minimal

=======
    #w0 = 200.*0.188
    alpha_ph = 700 #0.05*eps/np.pi# Ind.-Boson frame coupling
    J = EM.J_minimal
>>>>>>> 2017-overhaul-branch
    print "eps={:d}, T_EM={:d}, w_0={:d}, alpha_ph={:d}".format(int(eps), int(T_EM), int(w0), int(alpha_ph))
    """
    Now we build the sys-RC Hamiltonian and residual bath Liouvillian as well as generate mapped parameters
    """
    L_RC, H, A_EM, A_nrwa, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps, T_ph, wc, w0, alpha_ph, N)

    # electromagnetic bath liouvillians

<<<<<<< HEAD
    #L_nrwa = EM.L_nonrwa(H, A_nrwa, alpha_EM, T_EM) # Ignore this for now as it just doesn't work
    L_ns = EM.L_nonsecular(H, A_EM, alpha_EM, T_EM)
    L_s = EM.L_vib_lindblad(H, A_EM, alpha_EM, T_EM)
    L_naive = EM.L_EM_lindblad(eps, A_EM, alpha_EM, T_EM)

    # Set up the initial density matrix
    n_RC = EM.Occupation(wRC, T_ph)
    initial_sys = E*E.dag()

    #L_nrwa = EM.L_nonrwa(H, A_nrwa, eps, Gamma, T_EM) # Ignore this for now as it just doesn't work
    L_ns = EM.L_nonsecular(H, A_EM, eps, Gamma_EM, T_EM, J=J)
    #L_s = EM.L_vib_lindblad(H, A_EM, eps, Gamma_EM, T_EM, J=J)
    #L_naive = EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J)

    # Set up the initial density matrix
    n_RC = EM.Occupation(wRC, T_ph)
    #initial_sys = G*G.dag()

=======
    #L_nrwa = EM.L_nonrwa(H, A_nrwa, eps, Gamma, T_EM) # Ignore this for now as it just doesn't work
    L_ns = EM.L_nonsecular(H, A_EM, eps, Gamma_EM, T_EM, J=J)
    #L_s = EM.L_vib_lindblad(H, A_EM, eps, Gamma_EM, T_EM, J=J)
    #L_naive = EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J)

    # Set up the initial density matrix
    n_RC = EM.Occupation(wRC, T_ph)
    #initial_sys = G*G.dag()
>>>>>>> 2017-overhaul-branch
    #initial_sys = E*E.dag()
    initial_sys = 0.5*(G+E)*(E.dag()+G.dag())
    rho_0 = tensor(initial_sys, thermal_dm(N, n_RC))
    print (rho_0*tensor(qeye(2), destroy(N).dag()*destroy(N))).tr()

    # Expectation values and time increments needed to calculate the dynamics
<<<<<<< HEAD

    expects = [tensor(G*G.dag(), qeye(N)), tensor(E*G.dag(), qeye(N)), tensor(qeye(2), destroy(N).dag()*destroy(N))]
    timelist = np.linspace(0,10,60000) # you need lots of points so that coherences are well defined -> spectra

    expects = [tensor(G*G.dag(), qeye(N)), tensor(E*G.dag(), qeye(N)), tensor(qeye(2), destroy(N).dag()*destroy(N)), tensor(qeye(2), destroy(N).dag()+destroy(N))]
    timelist = np.linspace(0,0.04,8000)

=======
    expects = [tensor(G*G.dag(), qeye(N)), tensor(E*G.dag(), qeye(N)), tensor(qeye(2), destroy(N).dag()*destroy(N)), tensor(qeye(2), destroy(N).dag()+destroy(N))]
    timelist = np.linspace(0,0.04,8000)
>>>>>>> 2017-overhaul-branch
    #nonsec_check(eps, H, A_em, N) # Plots a scatter graph representation of non-secularity. Could use nrwa instead.

    # Calculate dynamics
    #DATA_nrwa = mesolve(H, rho_0, timelist, [L_RC+L_nrwa], expects, progress_bar=True)
    DATA_ns = mesolve(H, rho_0, timelist, [L_RC], expects, progress_bar=True)

    #DATA_s = mesolve(H, rho_0, timelist, [L_RC+L_s], expects, progress_bar=True)
    #DATA_naive = mesolve(H, rho_0, timelist, [L_RC+L_naive], expects, progress_bar=True)

    #timelist*=0.188 # Convert time into ps
<<<<<<< HEAD

    #fig = plt.figure(figsize=(12, 6))
    #ax1 = fig.add_subplot(121)
    #ax2 = fig.add_subplot(122)
    #plot_dynamics(ax1)
    #plot_coherences(ax2)

=======

    #fig = plt.figure(figsize=(12, 6))
    #ax1 = fig.add_subplot(121)
    #ax2 = fig.add_subplot(122)
    #plot_dynamics(ax1)
    #plot_coherences(ax2)

>>>>>>> 2017-overhaul-branch
    #p_file_name = "Notes/Images/Dynamics/Pop_min_a{:d}_N{:d}_Tem{:d}_w0{:d}_eps{:d}.pdf".format(int(alpha_ph), int(N), int(T_EM), int(w0), int(eps))
    #plt.savefig(p_file_name)
    #print "Figure saved: ", p_file_name
    """
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
<<<<<<< HEAD

    plot_dynamics(ax1)
    plot_manifolds(ax2, H)

    p_file_name = "Notes/Images/Dynamics/Pop_a{:d}_N{:d}_Tem{:d}_w0{:d}_eps{:d}_flat.pdf".format(int(alpha_ph), int(N), int(T_EM), int(w0), int(eps))

=======
>>>>>>> 2017-overhaul-branch
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
<<<<<<< HEAD

=======
=======
>>>>>>> d91a73078751b16e271f9eb594f26a765226f8fc
>>>>>>> 2017-overhaul-branch
