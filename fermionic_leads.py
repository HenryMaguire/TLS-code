from qutip import ket, basis, sigmam, sigmap, spre, sprepost, spost, destroy, mesolve
import numpy as np
from numpy import pi, linspace
import matplotlib.pyplot as plt
import qutip as qt

kB = 0.695


def J_flat(eps, Gamma_0):
    return Gamma_0

def fermi_occ(eps, T, mu):
    eps, T, mu = float(eps), float(T), float(mu)
    exp_part = np.exp((eps-mu)/(kB*T))
    return 1/(exp_part+1)

def additive_lead_dissipator(eps, d, T, mu, Gamma):
    ddag = d.dag()
    L = 0
    L+= pi*J_flat(eps, Gamma)*fermi_occ(eps, T, mu)*(spost(d*ddag)+spre(d*ddag)-2*sprepost(ddag, d))
    L+= pi*J_flat(eps, Gamma)*(1-fermi_occ(eps, T, mu))*(spost(ddag*d)+spre(ddag*d)-2*sprepost(d, ddag))
    return -L

def rate_up(e_lk, T, mu, Gamma_0):
    return pi*J_flat(e_lk, Gamma_0)*fermi_occ(e_lk, T, mu)

def rate_down(e_lk, T, mu, Gamma_0):
    return pi*J_flat(e_lk, Gamma_0)*(1-fermi_occ(e_lk, T, mu))

def limit_fermi(Gamma_0, T, mu):
    # up, down
    return pi*Gamma_0*(1-fermi_occ(2*mu, T, mu)), pi*Gamma_0*fermi_occ(2*mu, T, mu)

def non_additive_lead_dissipator(H, A, eps, T, mu, Gamma_0):
    L=0
    evals, estates = H.eigenstates()
    Zp_1 = 0
    Zp_2 = 0
    Zm_1 = 0
    Zm_2 = 0
    dim = len(evals)
    Adag = A.dag()
    for l in range(dim):
        for k in range(dim):
            e_lk = abs(evals[l]- evals[k])
            A_kl = A.matrix_element(estates[k].dag(), estates[l])
            Adag_lk = Adag.matrix_element(estates[l].dag(), estates[k])
            LK = estates[l]*estates[k].dag()
            KL = estates[k]*estates[l].dag()
            if e_lk != 0:
                Zp_1 += LK*Adag_lk*rate_up(e_lk, T, mu, Gamma_0)
                Zp_2 += LK*Adag_lk*rate_down(e_lk, T, mu, Gamma_0)
                Zm_1 += KL*A_kl*rate_up(e_lk, T, mu, Gamma_0)
                Zm_2 += KL*A_kl*rate_down(e_lk, T, mu, Gamma_0)
            else:
                rup, rdown = limit_fermi(Gamma_0, T, mu)
                Zp_1 += LK*Adag_lk*rup
                Zp_2 += LK*Adag_lk*rdown
                Zm_1 += KL*A_kl*rup
                Zm_2 += KL*A_kl*rdown
    #print Z_plus_1+Z_plus_2, Z_minus_1+Z_minus_2
    L += spre(A*Zp_1)-sprepost(Zp_1, A)
    L += -sprepost(A, Zp_2)+spost(Zp_2*A)
    L += spre(Adag*Zm_2)-sprepost(Zm_2, Adag)
    L += -sprepost(Adag, Zm_1)+spost(Zm_1*Adag)
    return -L

def L_R_lead_dissipators(H, A, eps, TL, muL, Gamma_0L, TR, muR, Gamma_0R):
    L_R = 0
    L_L = 0
    evals, estates = H.eigenstates()
    Zp_1L = 0
    Zp_2L = 0
    Zm_1L = 0
    Zm_2L = 0
    Zp_1R = 0
    Zp_2R = 0
    Zm_1R = 0
    Zm_2R = 0
    dim = len(evals)
    Adag = A.dag()
    for l in range(dim):
        for k in range(dim):
            e_lk = abs(evals[l]- evals[k])
            A_kl = A.matrix_element(estates[k].dag(), estates[l])
            Adag_lk = Adag.matrix_element(estates[l].dag(), estates[k])
            LK = estates[l]*estates[k].dag()
            KL = estates[k]*estates[l].dag()
            if e_lk != 0:
                Zp_1L += LK*Adag_lk*rate_up(e_lk, TL, muL, Gamma_0L)
                Zp_2L += LK*Adag_lk*rate_down(e_lk, TL, muL, Gamma_0L)
                Zm_1L += KL*A_kl*rate_up(e_lk, TL, muL, Gamma_0L)
                Zm_2L += KL*A_kl*rate_down(e_lk, TL, muL, Gamma_0L)
                Zp_1R += LK*Adag_lk*rate_up(e_lk, TR, muR, Gamma_0R)
                Zp_2R += LK*Adag_lk*rate_down(e_lk, TR, muR, Gamma_0R)
                Zm_1R += KL*A_kl*rate_up(e_lk, TR, muR, Gamma_0R)
                Zm_2R += KL*A_kl*rate_down(e_lk, TR, muR, Gamma_0R)
            else:
                rup, rdown = limit_fermi(Gamma_0L, TL, muL)
                Zp_1L += LK*Adag_lk*rup
                Zp_2L += LK*Adag_lk*rdown
                Zm_1L += KL*A_kl*rup
                Zm_2L += KL*A_kl*rdown
                rup, rdown = limit_fermi(Gamma_0R, TR, muR)
                Zp_1R += LK*Adag_lk*rup
                Zp_2R += LK*Adag_lk*rdown
                Zm_1R += KL*A_kl*rup
                Zm_2R += KL*A_kl*rdown
    #print Z_plus_1+Z_plus_2, Z_minus_1+Z_minus_2
    L_L += spre(A*Zp_1L)-sprepost(Zp_1L, A)
    L_L += -sprepost(A, Zp_2L)+spost(Zp_2L*A)
    L_L += spre(Adag*Zm_2L)-sprepost(Zm_2L, Adag)
    L_L += -sprepost(Adag, Zm_1L)+spost(Zm_1L*Adag)
    L_R += spre(A*Zp_1R)-sprepost(Zp_1R, A)
    L_R += -sprepost(A, Zp_2R)+spost(Zp_2R*A)
    L_R += spre(Adag*Zm_2R)-sprepost(Zm_2R, Adag)
    L_R += -sprepost(Adag, Zm_1R)+spost(Zm_1R*Adag)
    return -L_L,-L_R

def analytic_additive_current(eps, Gamma_L=0.1, Gamma_R=0.1, T_L=300., T_R=300., mu_L=0, mu_R=500):
    return 2*pi*(fermi_occ(eps, T_L, mu_L)-fermi_occ(eps, T_R, mu_R))*(Gamma_L*Gamma_R)/(Gamma_R+Gamma_L)

def current_from_L(H, L_full, L_track, n):
    rho_ss = qt.steadystate(H, [L_full])
    time_dependence = (qt.vector_to_operator(L_track*qt.operator_to_vector(rho_ss))*n).tr()
    return -time_dependence
