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
                    end_T_mult=10, tau_f_mult=1., per_tau=1., rotating=False,
                    T_increments=9, nsteps=2000, method='adams', order=12):
    # Start system+RC in adiabatic eigenstate and then time-evolve
    alpha_ph = prop_coupling*eps/pi
    w0 = eps*w0_prop
    if overdamped:
        wc = 53.
        Gamma = (w0**2)/wc
    w = np.linspace(0., eps*1.5, 1000)
    plt.figure()
    plt.plot(w, J_underdamped(w, alpha_ph, Gamma, w0))
    plt.axvline(eps, label='Splitting',color='k')
    plt.ylabel("Coupling Strength")
    plt.xlabel(r"Frequency $cm^{-1}$")
    plt.title("Phonon Spectral Density")
    plt.legend()
    plt.show()
    G, E = ket([0]), ket([1])
    sigma = G*E.dag() # Definition of a sigma_- operator.
    J = J_multipolar
    count = 1
    I_RC = qt.qeye(N)
    # This is true for specific Hamiltonian parameters

    L_RC, H_RC, A_EM, A_nrwa, Z, wRC, kappa, Gamma= RC.RC_function_UD(sigma, eps,
                                                                   T_ph, Gamma,
                                                                   w0, alpha_ph,
                                                                   N, rotating=rotating)
    H_RC_opt = RC.Ham_RC(sigma, eps, w0, kappa, N, rotating=False)[0]
    evals, states = H_RC.eigenstates()
    ground_list, excited_list = ground_and_excited_states(states)
    # work out how to initialise system rho
    init_rho = initialise_TLS(init_sys, init_RC, states, w0, T_ph)

    # electromagnetic bath liouvillians
    final_t = end_T_mult/Gamma_EM
    print "final t: ", final_t
    timelist = np.linspace(0, final_t, int(T_increments*final_t))
    options = qt.Options(nsteps=nsteps, store_states=True, method=method, order=order)
    E_op = tensor(E*E.dag(), I_RC)
    #
    pop_list = []
    g1_data = []
    spectrum_data = []
    freq_data = []
    #EM.L_non_rwa(H_RC_opt, A_nrwa, eps, Gamma_EM, T_EM, J, silent=True)]
    label = [ 'Full','Naive']
    for i, L_EM in enumerate([EM.L_non_rwa(H_RC_opt, A_nrwa, eps, Gamma_EM,
                            T_EM, J, silent=True),
                             EM.L_EM_lindblad(eps, A_EM, Gamma_EM, T_EM,
                             J=J, silent=True)]):
        ti = time.time()
        sigma_RC = tensor(sigma, I_RC)

        if i==2:
            sigma_plus_RC, sigma_RC, sigma_0_RC = EM.RWA_system_ops(H_RC_opt,
                                                tensor(sigma+sigma.dag(), I_RC))
        P = mesolve(H_RC, init_rho, timelist, [L_RC+L_EM], progress_bar=True,
                                                        options=options).states
        pop = [((E_op*p).tr()).real for p in P]

        if overdamped:
            steps_per_tau, tau_f = 1500*per_tau,  1.1
        else:
            steps_per_tau, tau_f = 2000*per_tau, 2.1

        tau_f *= tau_f_mult # extend tau

        print "Completed initial dynamics calculations for {} in {}  seconds.".format(label[i], time.time()-ti)
        ti = time.time()
        Lambda_0 = sigma_RC*sum(P)
        del P

        taulist = np.linspace(0, tau_f, int(tau_f*steps_per_tau))
        Lambda_t = mesolve(H_RC, Lambda_0, taulist, [L_RC+L_EM],options=options)

        g_1 = np.array([(sigma_RC.dag()*l).tr() for l in Lambda_t.states])
        g_1/=abs(g_1[0])

        spec = sp.fftpack.fft(g_1)
        dt = taulist[1]-taulist[0]
        freq = 2 * pi * np.array(sp.fftpack.fftfreq(spec.size, dt))
        spec = 2 * dt* np.real(spec)
        #freq, spec = qt.spectrum_correlation_fft(taulist, g_1) # old method
        spec-= min(spec)
        spec = spec/sum(spec)
        freq, spec = zip(*sorted(zip(freq, np.array(spec).real)))

        pop_list.append(pop)
        g1_data.append(g_1)
        spectrum_data.append(spec)
        freq_data.append(freq)
        print "Completed correlation function calculations for {} in {} seconds.".format(label[i], time.time()-ti)
    return timelist, pop_list, taulist, g1_data, spectrum_data, freq_data



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
