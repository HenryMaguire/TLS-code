from qutip import ket, mesolve, qeye, tensor, thermal_dm, destroy, steadystate
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
