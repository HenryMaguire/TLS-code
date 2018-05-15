import time

from qutip import ket, mesolve, qeye, tensor, thermal_dm, destroy, steadystate
import qutip as qt
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt

import UD_liouv as RC
import driving_liouv as EM
import exact_IB as exact
import scipy as sp

import phonon_weak_coupling as WC
from utils import J_overdamped, beta_f, J_underdamped, J_multipolar
from utils import ground_and_excited_states, initialise_TLS

reload(RC)
reload(EM)
reload(exact)
plt.style.use('ggplot')
plt.rcParams["axes.grid"] = True
plt.rcParams["axes.edgecolor"] = "0.15"
plt.rcParams["axes.linewidth"]  = 1.25
plt.rcParams['axes.facecolor'] = 'white'



def emission_spectra(init_sys, init_RC, prop_coupling, eps, Gamma, w0_prop=2.1,
                    T_ph=300., T_EM=0.,Gamma_EM=1., overdamped=False, N=9,
                    end_T_mult=10, tau_f_mult=1., per_tau=1.):
    wc = 53.
    # Start system+RC in adiabatic eigenstate and then time-evolve
    plt.close('all')
    alpha_ph = prop_coupling*eps/pi
    w0 = eps*w0_prop
    if overdamped:
        Gamma = (w0**2)/wc
    else:
        #N=20
        pass
    w = np.linspace(0., eps*1.5, 1000)
    plt.figure()
    plt.plot(w, J_underdamped(w, alpha_ph, Gamma, w0))
    plt.axvline(eps, label='Splitting',color='k')
    plt.ylabel("Coupling Strength")
    plt.xlabel(r"Frequency $cm^{-1}$")
    plt.title("Phonon Spectral Density")
    plt.legend()
    plt.show()
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
    J = J_multipolar
    count = 1
    I_RC = qt.qeye(N)
    # This is true for specific Hamiltonian parameters
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    dyn_DATA = []
    L_RC, H_RC, A_EM, A_nrwa, Z, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps,
                                                                   T_ph, Gamma,
                                                                   w0, alpha_ph, N, rotating=False)
    evals, states = H_RC.eigenstates()
    ground_list, excited_list = ground_and_excited_states(states)
    # work out how to initialise system rho
    init_rho = initialise_TLS(init_sys, init_RC, states, w0, T_ph)
    label = ['Naive', 'Full']
    #print len(excited_list), 'excited'
    # electromagnetic bath liouvillians
    final_t = end_T_mult/Gamma_EM
    timelist = np.linspace(0,final_t,int(380*final_t))
    options = qt.Options(nsteps=6000, store_states=True)
    f1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,6))
    f2, (ax3, ax4) = plt.subplots(1, 2, figsize=(14,6))
    E_op = tensor(E*E.dag(), I_RC)
    for i, L_EM in enumerate([EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J, silent=True),
                        EM.L_non_rwa(H_RC, A_nrwa, eps, Gamma_EM, T_EM, J, silent=True)]):
        #L_EM = EM.L_nonsecular(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
        #L_s = EM.L_vib_lindblad(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J, silent=True)
        ti = time.time()
        sigma_RC = tensor(sigma, I_RC)

        if label[i]=='Full':
            sigma_plus_RC, sigma_RC, sigma_0_RC = EM.RWA_system_ops(H_RC, tensor(sigma+sigma.dag(), I_RC))
            ls = 'solid'
        else:
            ls = 'dashed'
        P = mesolve(H_RC, init_rho, timelist, [L_RC+L_EM], progress_bar=None,options=options).states

        ax1.plot(timelist, [((E_op*p).tr()).real for p in P], label=label[i], linestyle=ls)
        if overdamped:
            steps_per_tau, tau_f = 1500*per_tau,  1.1
        else:
            if abs((w0-eps)/eps)<0.1:
                steps_per_tau, tau_f = 2000*per_tau, 2.1
            else:
                steps_per_tau, tau_f = 2000*per_tau, 1.1
        tau_f*=tau_f_mult
        #norm = np.array([(sigma_RC.dag()*sigma_RC*rho_t).tr() for rho_t in P])[0:tau_f*steps_per_tau]
        R= sum(P)
        print "Completed initial dynamics calculations for {} in {} seconds.".format(label[i], time.time()-ti)
        ti = time.time()
        Lambda_0 = sigma_RC*R
        del P

        taulist = np.linspace(0, tau_f, int(tau_f*steps_per_tau))
        Lambda_t = mesolve(H_RC, Lambda_0, taulist, [L_RC+L_EM],options=options)
        #sigma_t = mesolve(H_RC, init_rho, taulist, [L_RC+L_EM],
        #                  [sigma_RC.dag()*sigma_RC], options=options).expect[0]
        g_1 = np.array([(sigma_RC.dag()*l).tr() for l in Lambda_t.states])
        g_1/=abs(g_1[0])
        #else:
        #    g_1/=np.sqrt(norm*norm[0])

        #if label[i] == 'Full':
        #    g_1=g_1.conjugate()

        ax3.plot(taulist[0:int(tau_f*steps_per_tau/2.)],
                 g_1.real[0:int(tau_f*steps_per_tau/2.)], label=label[i], linestyle=ls)
        ax4.plot(taulist[0:int(tau_f*steps_per_tau/2.)],
                 g_1.imag[0:int(tau_f*steps_per_tau/2.)], label=label[i], linestyle=ls)
        freq, spec = qt.spectrum_correlation_fft(taulist, g_1)
        spec-= min(spec)
        spec = spec/sum(spec)
        ax2.plot(freq, spec.real, label=label[i], linestyle=ls)
        print "Completed correlation function calculations for {} in {} seconds.".format(label[i], time.time()-ti)
    ax1.set_xlim(0,final_t)
    #ax3.set_xlim(0,taulist[int((len(taulist)-1)/5.)])
    #ax4.set_xlim(0,taulist[int((len(taulist)-1)/5.)])
    ax1.set_xlabel(r"Time")
    ax2.set_xlabel(r"Frequency $cm^{-1}$")
    #ax2.axvline(eps, linestyle='dashed',color='k', alpha=0.4)
    ax2.set_xlim(0,2*eps)

    ax1.set_ylabel(r"Excited state population")
    ax2.set_ylabel(r"Fluorescence intensity (arb. units)")
    ax3.set_ylabel(r"$Re[g_1(\tau)]$") # coherence
    ax4.set_ylabel(r"$Im[g_1(\tau)]$") # coherence
    ax3.set_xlabel(r"$\tau$")
    ax4.set_xlabel(r"$\tau$")
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    plt.show()
    return freq, spec



reload(RC)
def absorption_spectra(prop_coupling, eps, Gamma, w0, T_ph=300., T_EM=0.,
                                    Gamma_EM=1., overdamped=True, N=3):
    # Start system+RC in adiabatic eigenstate and then time-evolve
    plt.close('all')
    alpha_ph = prop_coupling*eps/np.pi

    if overdamped:
        w0 = eps*1.01
        Gamma = (w0**2)/wc
    else:
        #N=20
        pass
    w = np.linspace(0., eps*1.5, 1000)
    plt.figure()
    plt.plot(w, J_underdamped(w, alpha_ph, Gamma, w0))
    plt.axvline(eps, label='Splitting',color='k')
    plt.ylabel("Coupling Strength")
    plt.xlabel(r"Frequency $cm^{-1}$")
    plt.title("Phonon Spectral Density")
    plt.legend()
    plt.show()
    G = ket([0])
    E = ket([1])
    sigma = G*E.dag()
    J = J_multipolar
    I_RC = qt.qeye(N)
    # This is true for specific Hamiltonian parameters
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    dyn_DATA = []

    L_RC, H_RC, A_EM, A_nrwa, Z, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps,
                                                                   T_ph, Gamma,
                                                                   w0, alpha_ph, N)

    evals, states = H_RC.eigenstates()
    ground_list, excited_list = ground_and_excited_states(states)
    # work out how to initialise system rho

    label = ['Naive', 'Full']
    #print len(excited_list), 'excited'
    # electromagnetic bath liouvillians
    final_t = 2
    timelist = np.linspace(0,final_t,450*final_t)
    options = qt.Options(nsteps=1500)
    f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8,12))
    E_op = tensor(E*E.dag(), I_RC)
    mu_ESA = tensor(G*E.dag()+E*G.dag(), I_RC) # Definition of a sigma_- operator.
    sigma_p_RWA, sigma_RWA, _ = EM.RWA_system_ops(H_RC, tensor(sigma+sigma.dag(), I_RC))
    mu = [mu_ESA, sigma_p_RWA+sigma_RWA]
    for i, L_EM in enumerate([EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM, J=J, silent=True),
                        EM.L_non_rwa(H_RC, A_nrwa, eps, Gamma_EM, T_EM, J, silent=True)]):
        #L_EM = EM.L_nonsecular(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J)
        #L_s = EM.L_vib_lindblad(H_RC, A_EM, eps, Gamma_EM, T_EM, J=J, silent=True)
        ti = time.time()
        init_rho = initialise_TLS(0, 0, states, w0, T_ph)
        init_rho = mu[i]*init_rho-init_rho*mu[i]
        S1 = mesolve(H_RC, init_rho, timelist,  [L_RC+L_EM], [mu[i]],
                             progress_bar=None,options=options).expect[0]
        ax1.plot(timelist, S1.real, label=label[i]+' real')
        ax2.plot(timelist, S1.imag, label=label[i]+' imag')
        #norm = np.array([(sigma_RC.dag()*sigma_RC*rho_t).tr() for rho_t in P])[0:tau_f*steps_per_tau]
        print "Completed response function calculations for {} in {} seconds.".format(label[i], time.time()-ti)

        freq, spec = qt.spectrum_correlation_fft(timelist, S1)
        #spec-= min(spec)
        spec = spec/sum(spec)
        ax3.plot(freq, spec, label=label[i])
    ax1.set_xlim(0,timelist[int((len(timelist)-1)/2.)])
    ax2.set_xlim(0,timelist[int((len(timelist)-1)/2.)])
    ax1.legend()
    ax3.axvline(eps, linestyle='dashed',color='k')
    #ax3.set_xlim(0,2000)
    ax2.legend()
    ax3.legend()
    plt.show()
    return freq, spec
