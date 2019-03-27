import time
import numpy as np

import matplotlib.pyplot as plt
import qutip as qt
from scipy import integrate

import UD_liouv as RC
import driving_liouv as EM

import phonon_weak_coupling as WC

from qutip import ket, basis, sigmam, sigmap, sigmaz, spre, sprepost, spost, destroy, mesolve, tensor, qeye, Qobj
from numpy import pi, linspace, sqrt
from utils import J_overdamped, beta_f, J_underdamped, J_minimal_hard, J_multipolar
from utils import ground_and_excited_states, initialise_TLS

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

kB = 0.695
EE = basis(2,1)*basis(2,1).dag()
GG = basis(2,0)*basis(2,0).dag()
Zero = tensor(basis(2,0), basis(2,0))
sz = EE - GG
sm = destroy(2)
d2 = tensor(sm, sz) # swapped around from the calculations
d1 = tensor(qeye(2), sm)
I_sys = Qobj(qeye(4), dims= [[2, 2], [2, 2]])
d1dag, d2dag = d1.dag(), d2.dag()

z_ket = basis(4,0)
l_ket = basis(4,1)
lp_ket = basis(4,2)
d_ket = basis(4,3)

D = d1dag*d1*d2dag*d2
Z = Qobj(z_ket*z_ket.dag(), dims=D.dims)
LUMO = Qobj(l_ket*l_ket.dag(), dims=D.dims)
LUMOp = Qobj(lp_ket*lp_ket.dag(), dims=D.dims)
n = d1dag*d1 + d2dag*d2 # populatiion operator

J_leads = J_underdamped

def build_H(eps1, eps2, U):
    return eps1*d1dag*d1 + eps2*d2dag*d2 + U*D

def build_general_ME(d1, d2):
    L = 2. * spre(O) * spost(Od) - spre(Od * O) - spost(Od * O)
    return L

def commutator_term1(O1, O2):
    # [O1, O2*rho]
    return spre(O1*O2)-sprepost(O2, O1)

def commutator_term2(O1, O2):
    # [rho*O1, O2]
    return spost(O1*O2)-sprepost(O2, O1)

def fermi_occ(eps, T, mu):
    kB = 0.695
    exp_part = np.exp((eps-mu)/(kB*T))
    return 1/(exp_part+1)

def current_from_L(H, L_full, L_track, obs_ops, method='direct'):
    obs_out = []
    use_precond = True
    if method in ['direct','svd', 'eigen', 'power']:
        use_precond = False
    rho_ss = qt.steadystate(H, [L_full], method=method, use_precond=use_precond)
    obs_out.append(-(qt.vector_to_operator(L_track*qt.operator_to_vector(rho_ss))*obs_ops[0]).tr())
    for obs in obs_ops[1::]:
        obs_out.append((rho_ss*obs).tr())

    del rho_ss
    return obs_out

def cauchyIntegrands(eps, J, height, width, pos, T, mu, ver=1):
    # Function which will be called within another function where other inputs
    # are defined locally
    F = 0
    if ver == -1:
        F = J(eps, height, width, pos)*(1-fermi_occ(eps, T, mu))
    elif ver == 1:
        F = J(eps, height, width, pos)*(fermi_occ(eps, T, mu))
    return F


def Lamdba_complex_rate(eps, J, mu, T, height, width, pos, type='m', plot_integrands=False, real_only=False):
    F_p = (lambda x: (cauchyIntegrands(x, J, height, width, pos, T, mu, ver=1)))
    F_m = (lambda x: (cauchyIntegrands(x, J, height, width, pos, T, mu, ver=-1)))
    if plot_integrands:
        w = np.linspace(0,4*eps, 300)
        plt.plot(w, F_p(w), label='+')
        plt.plot(w, F_m(w), label='-')
        plt.legend()
    if type=='m':
        if real_only:
            Pm=0.
        else:
            Pm = integrate.quad(F_m, 0, 5*eps, weight='cauchy', wvar=eps)[0]

        return pi*F_m(eps) + 1j*Pm
    elif type=='p':
        if real_only:
            Pp=0.
        else:
            Pp = integrate.quad(F_p, 0, 5*eps, weight='cauchy', wvar=eps)[0] #integral_converge(F_p, 0, eps)
        return pi*F_p(eps) + 1j*Pp
    else:
        raise ValueError

def additive_liouvillian(eps1 =0., eps2=900., T_L=77., mu_L=1000.,
                         width_L=1000., pos_L=900., height_L=1., T_R=77.,
                         mu_R=0., width_R=1000., pos_R=900., height_R=1.,
                        secular=False, real_only=False):
    J = J_leads
    eps = [eps1, eps2]
    width = [width_L, width_R]
    pos = [pos_L, pos_R]
    T = [T_L, T_R]
    mu = [mu_L, mu_R]
    height= [height_L, height_R]
    L = []
    d = [d1, d2]
    ddag = [d1dag, d2dag]
    #_ = Lamdba_complex_rate(eps2, J, mu[0], T[0], height[0], width[0], pos[0], type='m', plot_integrands=True)
    #_ = Lamdba_complex_rate(eps2, J, mu[1], T[1], height[1], width[1], pos[1], type='m', plot_integrands=True)
    for j in range(2): # over leads

        L_j = 0
        for q in range(2): # over sites
            for p in range(2):
                if secular and p!=q:
                    pass
                else:
                    L_j += Lamdba_complex_rate(eps[p], J, mu[j],
                                             T[j], height[j],
                                             width[j], pos[j],
                                             type='p',real_only=real_only)*commutator_term1(d[q], ddag[p])

                    L_j += Lamdba_complex_rate(eps[p], J, mu[j],
                                             T[j], height[j],
                                             width[j], pos[j],
                                             type='m',real_only=real_only)*commutator_term2(ddag[p], d[q])
                    L_j += Lamdba_complex_rate(eps[p], J, mu[j],
                                             T[j], height[j],
                                             width[j], pos[j],
                                             type='m',real_only=real_only).conjugate()*commutator_term1(ddag[q], d[p])
                    L_j += Lamdba_complex_rate(eps[p], J, mu[j],
                                             T[j], height[j],
                                             width[j], pos[j],
                                             type='p',real_only=real_only).conjugate()*commutator_term2(d[p], ddag[q])
        L.append(L_j)
    return -L[0], -L[1]


def simple_current_voltage(eps1=100., eps2=900., U=0., T_L=77.,
                         width_L=1000., pos_L=900., T_R=77.,
                         mu_R=0., width_R=1000., pos_R=1000.,
                               secular=True, mu_L_max =2000):
    H = build_H(eps1, eps2, U)
    mu_Ls = np. linspace(0, mu_L_max, 100)
    currents = []
    for mu_L in mu_Ls:
        L_L, L_R = additive_liouvillian(eps1=eps1, eps2=eps2,  mu_L=mu_L,T_L=T_L,
                         width_L=width_L, pos_L=pos_L, T_R=T_R,
                         mu_R=mu_R, width_R=width_R, pos_R=pos_R,
                               secular=secular)
        currents.append(current_from_L(H, L_L+L_R, L_R, n))

def L_R_lead_dissipators(H, A, T_L=77., mu_L=1000.,
                     width_L=1000., pos_L=900., height_L=1., T_R=77.,
                     mu_R=0., width_R=1000., pos_R=900., height_R=1.,
                        real_only=False, silent=True):
    ti = time.time()
    L_leads = []
    T = [T_L, T_R]
    mu = [mu_L, mu_R]
    Gamma_0 = [height_L, height_R]
    width =[width_L, width_R]
    pos = [pos_L, pos_R]
    evals, estates = H.eigenstates()
    Z = []
    J = J_underdamped
    dim = len(evals)
    Adag = A.dag()
    for j in range(2): # for left and right lead
        Zp_1, Zp_2, Zm_1, Zm_2 = 0,0,0,0
        for l in range(dim):
            for k in range(dim):
                e_lk = abs(evals[l]- evals[k])
                A_kl = A.matrix_element(estates[k].dag(), estates[l])
                Adag_lk = Adag.matrix_element(estates[l].dag(), estates[k])
                LK = estates[l]*estates[k].dag()
                KL = estates[k]*estates[l].dag()
                if e_lk != 0:
                    rate_up = Lamdba_complex_rate(e_lk, J, mu[j],
                                             T[j], Gamma_0[j],
                                             width[j], pos[j],
                                             type='p',real_only=real_only)
                    rate_down = Lamdba_complex_rate(e_lk, J, mu[j],
                                             T[j], Gamma_0[j],
                                             width[j], pos[j],
                                             type='m',real_only=real_only)
                    Zp_1 += LK*Adag_lk*rate_up
                    Zp_2 += LK*Adag_lk*rate_down
                    Zm_1 += KL*A_kl*rate_up.conjugate()
                    Zm_2 += KL*A_kl*rate_down.conjugate()
                else:
                    pass
                    """rup, rdown = limit_fermi_flat(Gamma_0[j], T[j], mu[j])
                    Zp_1 += LK*Adag_lk*rup
                    Zp_2 += LK*Adag_lk*rdown
                    Zm_1 += KL*A_kl*rup
                    Zm_2 += KL*A_kl*rdown"""
        Z.append([Zp_1, Zp_2, Zm_1, Zm_2])

    #print Z_plus_1+Z_plus_2, Z_minus_1+Z_minus_2

    for j in range(2):
        Zp_1, Zp_2, Zm_1, Zm_2 = Z[j][0],Z[j][1],Z[j][2],Z[j][3]
        L=0
        L += spre(A*Zp_1)-sprepost(Zp_1, A)
        L += -sprepost(A, Zp_2)+spost(Zp_2*A)
        L += spre(Adag*Zm_2)-sprepost(Zm_2, Adag)
        L += -sprepost(Adag, Zm_1)+spost(Zm_1*Adag)
        L_leads.append(L)
    if not silent:
        print("Calculating the lead dissipators took {} seconds.".format(time.time()-ti))
    return -L_leads[0],-L_leads[1]

def nonadditive_current_voltage(eps1=500., eps2=900., U=0., T_L=77.,
                         width_L=1000., pos_L=900., height_L=1.,T_R=77.,
                         mu_R=0., width_R=1000., pos_R=1000.,height_R=1.,
                               secular=True, mu_L_max =2000):
    H = build_H(eps1, eps2, U)
    n = d1dag*d1 + d2dag*d2
    mu_Ls = np. linspace(0, mu_L_max, 100)
    currents = []
    for mu_L in mu_Ls:
        L_L, L_R = L_R_lead_dissipators(H, d1+d2, T_L=T_L, mu_L=mu_L,
                     width_L=width_L, pos_L=pos_L, height_L=height_L, T_R=T_R,
                     mu_R=mu_R, width_R=width_R, pos_R=pos_R, height_R=height_R)
        currents.append(current_from_L(H, L_L+L_R, L_R, n))
    return mu_Ls-mu_R, currents

def compare_add_and_nonadd(eps1=100., eps2=1000., U=0., T_L=77., mu_L=1000.,
                         width_L=1000., pos_L=900., height_L=1., T_R=77.,
                         mu_R=0., width_R=1000., pos_R=900., height_R=1.):

    L_L, L_R = additive_liouvillian(eps1=eps1, eps2=eps2,
                                    T_L=T_L, T_R=T_R,
                                    mu_L=mu_L, mu_R=mu_R,
                                    width_L=width_L, width_R=width_R,
                                    pos_L=pos_L, pos_R=pos_R,
                                    height_L=height_L, height_R=height_R,
                                    real_only=False)

    H = build_H(eps1, eps2, U)
    timelist = np.linspace(0,9,2000)
    rho0 = Z
    elist = [LUMOp, LUMO, D, LUMOp+LUMO+D]
    #print H.dims, rho0.dims, L_R.dims, L_L.dims, elist[0].dims
    data = mesolve(H, rho0, timelist, [L_L+L_R], elist)
    #rho_ss_num = (-((1/(kB*T))*H-mu*n)).expm()
    #rho_ss = rho_ss_num/rho_ss_num.tr()
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    plt.figure()
    plt.plot(timelist, data.expect[0], label='LUMO+', color=colors[0])
    plt.plot(timelist, data.expect[1], label='LUMO', color=colors[1])
    plt.plot(timelist, data.expect[2], label='Double', color=colors[2])
    plt.plot(timelist, data.expect[0]+data.expect[1]+data.expect[2], color=colors[3])
    plt.ylabel("Excited level population")
    plt.legend()
    #plt.axhline((rho_ss*elist[0]).tr(), ls='dashed', color='r')
    L_L, L_R = L_R_lead_dissipators(H, d1+d2,
                                    T_L=T_L, T_R=T_R,
                                    mu_L=mu_L, mu_R=mu_R,
                                    width_L=width_L, width_R=width_R,
                                    pos_L=pos_L, pos_R=pos_R,
                                    height_L=height_L, height_R=height_R, real_only=False)
    data = mesolve(H, rho0, timelist, [L_L+L_R], elist)
    plt.plot(timelist, data.expect[0], ls='dashed', color=colors[0])
    plt.plot(timelist, data.expect[1], ls='dashed', color=colors[1])
    plt.plot(timelist, data.expect[2], ls='dashed', color=colors[2])
    plt.plot(timelist, data.expect[0]+data.expect[1]+data.expect[2], ls='dashed', color=colors[3])
    plt.show()

def Ham_RC(H, n_op, sigma_op, Omega, kappa, N, rotating=False):
    """
    Input: System splitting, RC freq., system-RC coupling and Hilbert space dimension
    Output: Hamiltonian, sigma_- and sigma_z in the vibronic Hilbert space
    """
    if rotating:
        eps=0.
    a = destroy(N)
    shift = (kappa**2)/Omega
    I_sys = Qobj(qeye(sigma_op.shape[0]), dims=sigma_op.dims)
    I_RC = qeye(N)
    H_S = tensor(H, I_RC) + shift*tensor(n_op, I_RC)
    H_S += kappa*tensor(n_op, (a + a.dag())) + tensor(I_sys, Omega*a.dag()*a)
    A_em = tensor(sigma_op, qeye(N))
    A_nrwa = tensor(sigma_op+sigma_op.dag(), qeye(N))
    A_ph = tensor(I_sys, (a + a.dag()))
    return H_S, A_em, A_nrwa, A_ph



def RC_function_UD(H, n_op, sigma_op, T_ph, Gamma, wRC, alpha_ph, N, silent=False,
                                            residual_off=False, rotating=False):
    # we define all of the RC parameters by the underdamped spectral density
    gamma = Gamma / (2. * np.pi * wRC)  # coupling between RC and residual bath
    if residual_off:
        gamma=0
    kappa= np.sqrt(np.pi * alpha_ph * wRC / 2.)  # coupling strength between the TLS and RC

    if not silent:
        print("w_RC={} | RC-res. coupling={:0.2f} | TLS-RC coupling={:0.2f} | Gamma_RC={:0.2f} | alpha_ph={:0.2f} | N={} |".format(wRC, gamma,  kappa, Gamma, alpha_ph, N))
    H, A_em, A_nrwa, A_ph = Ham_RC(H, n_op, sigma_op, wRC, kappa, N, rotating=rotating)
    L_RC, Z =  RC.liouvillian_build(H, A_ph, gamma, wRC, T_ph)
    return L_RC, H, A_em, A_nrwa, Z, wRC, kappa, Gamma

def current_vs_voltage_with_phonons(eps1=1000., eps2=1000., U= 0., T_ph=77,
                                    Gamma=10., w0=50., alpha_ph=10., N=10,
                                    Gamma_L=0.5, T_L=77., width_L=1000., pos_L=1000.,
                                    Gamma_R=0.5, T_R=77., mu_R=0., width_R=1000., pos_R=1000.,
                                   real_only=False):
    T_L = T_ph
    T_R = T_ph
    ti = time.time()
    mu_Ls = np. linspace(0, 30.*w0, 70)
    print("RC would need {} states to fill electronic gap.".format(eps/w0))
    currents = []

    H = build_H(eps1, eps2, U)
    L_RC, H_RC, A_EM, A_nrwa, Z, _, _, _, = RC_function_UD(H, n, d1+d2, T_ph, Gamma, w0, alpha_ph,
                                                         N, silent=True)
    d_RC = tensor(d1+d2, qeye(N))
    E_RC = tensor(n, qeye(N))
    #timelist = np.linspace(0, 3/Gamma_EM, 360)
    for i, mu_L in enumerate(mu_Ls):

        L_L, L_R = L_R_lead_dissipators(H_RC, d_RC, T_L=T_L, mu_L=mu_L,
                         width_L=width_L, pos_L=pos_L, height_L=height_L, T_R=T_R,
                         mu_R=mu_R, width_R=width_R, pos_R=pos_R, height_R=height_R,
                         real_only=real_only)
        currents.append(current_from_L(H_RC, L_RC+L_L+L_R, L_R, E_RC))
        if (i%10)==0:
            print(100*(float(i)/len(mu_Ls)), "% complete")
    print("Took {} seconds.".format(time.time()-ti))
    return mu_Ls-mu_R, currents

def current_vs_phonon_coupling(eps1, eps2, T_ph=77., Gamma=30., w0=70., U=0., N=10,
                               gamma_L=1., T_L=77., mu_L=1000., width_L=1000., pos_L=1000.,
                               gamma_R=1., T_R=77., mu_R=0., width_R=1000., pos_R=1000.,
                               height_L=1., height_R=1., method='iterative-gmres'):
    ti = time.time()
    eps = abs(eps2-eps1)
    alpha_prop = np. linspace(0., 1., 15)
    print("RC would need {} states to fill electronic gap.".format(eps/w0))
    currents_nonadd = []
    #timelist = np.linspace(0, 3/Gamma_EM, 360)
    H = build_H(eps1, eps2, U)
    d_sub = d1+d2
    d_RC = tensor(d_sub, qeye(N))
    E_RC = tensor(n, qeye(N))
    RC_occ = destroy(N).dag()*destroy(N)
    RC_occ_ = tensor(I_sys, RC_occ)
    obs_ops = [E_RC, RC_occ_]
    #timelist = np.linspace(0, 3/Gamma_EM, 360)

    for i, alphap in enumerate(alpha_prop):
        alpha_ph = alphap*eps/pi

        L_RC, H_RC, A_EM, A_nrwa, Z, _, _, _ = RC_function_UD(H, n, d_sub, T_ph,
                                                            Gamma, w0, alpha_ph,
                                                            N, silent=True)
        L_Lfull, L_Rfull = L_R_lead_dissipators(H_RC, d_RC,
                                            T_L=T_L, mu_L=mu_L, width_L=width_L,
                                            pos_L=pos_L, height_L=height_L,
                                            T_R=T_R, mu_R=mu_R, width_R=width_R,
                                            pos_R=pos_R, height_R=height_R)
        #currents_add.append(current_from_L(H_RC, L_RC+L_Ladd+L_Radd, L_Radd, E_RC))
        currents_nonadd.append(current_from_L(H_RC, L_RC+L_Lfull+L_Rfull,
                                            L_Rfull, obs_ops, method=method))
        del L_Lfull, L_Rfull, L_RC, H_RC, A_EM, A_nrwa, Z
        if (i%10)==0:
            print(100*(float(i)/len(alpha_prop)), "% complete")
    print("Took {} seconds.".format(time.time()-ti))
    return alpha_prop,  np.array(currents_nonadd).T


"""

a_prop3,  pc3 = current_vs_phonon_coupling(1000., 200., T_ph=77., Gamma=30., w0=50.,
                                           N=15, T_L=77., T_R=77., mu_R=0., mu_L=1500.
                                           , method='iterative-gmres')
#plt.plot(a_prop1, pc1, label='5', color=colors[0])
#plt.plot(a_prop2, pc2, label='8', color=colors[1])
plt.figure()
plt.plot(a_prop3, pc3[0], label='current', color=colors[2])
plt.figure()
plt.plot(a_prop3, pc3[1], label='phonon occ.', color=colors[2])
#plt.plot(a_prop4, pc4, label='15', color=colors[3])
plt.legend(loc='best')
plt.show()"""
