"""
The four electromagnetic liouvillians I am studying:
- no rotating wave approximation
- no secular approximation
- a secular approximation
- an approximation which says that the enlarged system eigenstates are the same as the uncoupled system eigenstates

"""

import numpy as np
import scipy as sp
from scipy import integrate
import qutip as qt
from qutip import destroy, tensor, qeye, spost, spre, sprepost
import time
from utils import J_minimal, beta_f, J_minimal_hard

def Occupation(omega, T, time_units='cm'):
    conversion = 0.695
    if time_units == 'ev':
        conversion == 8.617E-5
    if time_units == 'ps':
        conversion == 0.131
    else:
        pass
    n =0.
    beta = 0.
    if T ==0.: # First calculate beta
        n = 0.
        beta = np.infty
    else:
        # no occupation yet, make sure it converges
        beta = 1. / (conversion*T)
        if sp.exp(omega*beta)-1 ==0.:
            n = 0.
        else:
            n = float(1./(sp.exp(omega*beta)-1))
    return n

"""
def J_multipolar(omega, Gamma, omega_0):
    return Gamma*(omega**3)/(2*np.pi*(omega_0**3))

def J_minimal(omega, Gamma, omega_0):
    return Gamma*omega/(2*np.pi*omega_0)

def J_flat(omega, Gamma, omega_0):
    return Gamma
"""
'''
def cauchyIntegrands(omega, beta, J, ver):
    # Function which will be called within another function where J, beta and the eta are defined locally.
    F = 0
    if ver == 1:
        F = J(omega)*(coth(beta*omega/2.)+1)
    elif ver == -1:
        F = J(omega)*(coth(beta*omega/2.)-1)
    elif ver == 0:
        F = J(omega)
    return F

def Gamma(omega, beta, J):
    G = 0
    # Here I define the functions which "dress" the integrands so they have only 1 free parameter for Quad.
    F_p = (lambda x: (cauchyIntegrands(x, beta, J, 1)))
    F_m = (lambda x: (cauchyIntegrands(x, beta, J, -1)))
    F_0 = (lambda x: (cauchyIntegrands(x, beta, J, 0)))
    n = 30
    print "Cauchy int. convergence checks: ", F_p(4*n), F_p(4*n), F_p(4*n)
    if omega>0.:
        # These bits do the Cauchy integrals too
        G = (np.pi/2)*(coth(beta*omega/2.)-1)*J(omega)+ (1j/2.)*((integrate.quad(F_m, 0, n, weight='cauchy', wvar=omega)[0]+integrate.quad(F_m, n, 2*n, weight='cauchy', wvar=omega)[0]+integrate.quad(F_m, 2*n, 3*n, weight='cauchy', wvar=omega)[0]+integrate.quad(F_m, 3*n, 4*n, weight='cauchy', wvar=omega)[0]) - (integrate.quad(F_p, 0, n, weight='cauchy', wvar=-omega)[0]+integrate.quad(F_p, n, 2*n, weight='cauchy', wvar=-omega)[0]+integrate.quad(F_p, 2*n, 3*n, weight='cauchy', wvar=-omega)[0]+integrate.quad(F_p, 3*n, 4*n, weight='cauchy', wvar=-omega)[0]))
        #print integrate.quad(F_m, 0, n, weight='cauchy', wvar=omega), integrate.quad(F_p, 0, n, weight='cauchy', wvar=-omega)
    elif omega==0.:
        # The limit as omega tends to zero is zero for superohmic case?
        G = -(1j)*(integrate.quad(F_0, -1e-12, n, weight='cauchy', wvar=0)[0]+integrate.quad(F_0, n, 2*n, weight='cauchy', wvar=0)[0]+integrate.quad(F_0, 2*n, 3*n, weight='cauchy', wvar=0)[0]+integrate.quad(F_0, 3*n, 4*n, weight='cauchy', wvar=0)[0])
        #print (integrate.quad(F_0, -1e-12, 20, weight='cauchy', wvar=0)[0])
    elif omega<0.:
        G = (np.pi/2)*(coth(beta*abs(omega)/2.)+1)*J(abs(omega))+ (1j/2.)*((integrate.quad(F_m, 0, n, weight='cauchy', wvar=-abs(omega))[0]+integrate.quad(F_m, n, 2*n, weight='cauchy', wvar=-abs(omega))[0]+integrate.quad(F_m, 2*n, 3*n, weight='cauchy', wvar=-abs(omega))[0]+integrate.quad(F_m, 3*n, 4*n, weight='cauchy', wvar=-abs(omega))[0]) - (integrate.quad(F_p, 0, n, weight='cauchy', wvar=abs(omega))[0]+integrate.quad(F_p, n, 2*n, weight='cauchy', wvar=abs(omega))[0]+integrate.quad(F_p, 2*n, 3*n, weight='cauchy', wvar=abs(omega))[0]+integrate.quad(F_p, 3*n, 4*n, weight='cauchy', wvar=abs(omega))[0]))
        #print integrate.quad(F_m, 0, n, weight='cauchy', wvar=-abs(omega)), integrate.quad(F_p, 0, n, weight='cauchy', wvar=abs(omega))
    return G

def decay_rate(J, omega, Gamma, omega_0, T, time_units):
    """
    Decay rate for non-secular master equation. In my notes, this is called Gamma.
    """
    conversion = 0.695
    if time_units == 'ps': # allows conversion to picoseconds
        conversion == 8.617E-5
    else:
        pass

    beta = 0
    if T==0:
        beta = 10E12 # Some very large number, this seems dodgy.
    else:
        beta = 1. / (conversion* T)

    Decay = 0 # Initialise decay rate variable
    coth = lambda x : (sp.exp(2*x)+1)/(sp.exp(2*x)-1)
    if omega>0:
        Decay = np.pi*0.5*J(omega, Gamma, omega_0)*(coth(beta*omega/2)-1)
    elif omega ==0:
        Decay = np.pi*J(1, Gamma, omega_0) #/beta I thought this was supposed to be just J(1), but mathematica says otherwise
    else:
        Decay = 0.5*J(omega, Gamma, omega_0)*(coth(beta*omega/2)+1)
    return Decay
'''
import sympy
def coth(x):
    return float(sympy.coth(x))

def cauchyIntegrands(omega, beta, J, alpha, wc, ver):
    # J_overdamped(omega, alpha, wc)
    # Function which will be called within another function where J, beta and
    # the eta are defined locally
    F = 0
    if ver == 1:
        F = J(omega, alpha, wc)*(coth(beta*omega/2.)+1)
    elif ver == -1:
        F = J(omega, alpha, wc)*(coth(beta*omega/2.)-1)
    elif ver == 0:
        F = J(omega, alpha, wc)
    return F
import time
def integral_converge(f, a, omega):
    inc = 5.
    x = inc
    I = 0
    while abs(f(x))>5E-10:
        #print x, f(x)
        I += integrate.quad(f, a, x, weight='cauchy', wvar=omega)[0]
        a+=inc
        x+=inc
        #time.sleep(0.1)
    print "Integral converged"
    return I # Converged integral

def Gamma(omega, beta, J, alpha, wc, imag_part=True):
    G = 0
    # Here I define the functions which "dress" the integrands so they
    # have only 1 free parameter for Quad.
    F_p = (lambda x: (cauchyIntegrands(x, beta, J, alpha, wc, 1)))
    F_m = (lambda x: (cauchyIntegrands(x, beta, J, alpha, wc, -1)))
    F_0 = (lambda x: (cauchyIntegrands(x, beta, J, alpha, wc, 0)))
    w='cauchy'
    if omega>0.:
        # These bits do the Cauchy integrals too
        G = (np.pi/2)*(coth(beta*omega/2.)-1)*J(omega,alpha, wc)
        if imag_part:
            G += (1j/2.)*(integral_converge(F_m, 0,omega))
            G -= (1j/2.)*(integral_converge(F_p, 0,-omega))

        #print integrate.quad(F_m, 0, n, weight='cauchy', wvar=omega), integrate.quad(F_p, 0, n, weight='cauchy', wvar=-omega)
    elif omega==0.:
        G = (np.pi/2)*(2*alpha/beta)
        # The limit as omega tends to zero is zero for superohmic case?
        if imag_part:
            G += -(1j)*integral_converge(F_0, -1e-12,0)
        #print (integrate.quad(F_0, -1e-12, 20, weight='cauchy', wvar=0)[0])
    elif omega<0.:
        G = (np.pi/2)*(coth(beta*abs(omega)/2.)+1)*J(abs(omega),alpha, wc)
        if imag_part:
            G += (1j/2.)*integral_converge(F_m, 0,-abs(omega))
            G -= (1j/2.)*integral_converge(F_p, 0,abs(omega))
        #print integrate.quad(F_m, 0, n, weight='cauchy', wvar=-abs(omega)), integrate.quad(F_p, 0, n, weight='cauchy', wvar=abs(omega))
    return G

def L_non_rwa(H_vib, A, w_0, alpha, T_EM, J, principal=True):
    ti = time.time()
    beta = beta_f(T_EM)

    eVals, eVecs = H_vib.eigenstates()
    #J=J_minimal # J_minimal(omega, Gamma, omega_0)
    d_dim = len(eVals)
    G = 0
    for i in xrange(d_dim):
        for j in xrange(d_dim):
            eta = eVals[i]-eVals[j]
            s = eVecs[i]*(eVecs[j].dag())
            #print A.matrix_element(eVecs[i].dag(), eVecs[j])
            s*= A.matrix_element(eVecs[i].dag(), eVecs[j])
            s*= Gamma(eta, beta, J, alpha, w_0, imag_part=principal)
            G+=s
    G_dag = G.dag()
    # Initialise liouvilliian
    L =  qt.spre(A*G) - qt.sprepost(G, A)
    L += qt.spost(G_dag*A) - qt.sprepost(A, G_dag)
    print "Calculating non-RWA Liouvilliian took {} seconds."
    return -L
"""
def L_nonrwa(H_vib, sig_x, omega_0, Gamma, T, J, time_units='cm'):
    ti = time.time()
    d = H_vib.shape[0]
    evals, evecs = H_vib.eigenstates()
    Z = 0 # initalise rate operator
    Z_dag = 0
    for i in range(int(d)):
        for j in range(int(d)):
            eps_ij = (evals[i]-evals[j])
            sig_ij = sig_x.matrix_element(evecs[i].dag(), evecs[j])
            sig_ji = (sig_x.dag()).matrix_element(evecs[j].dag(), evecs[i])
            IJ = evecs[i]*evecs[j].dag()
            JI = evecs[j]*evecs[i].dag()
            if abs(sig_ij) >0:
                Z+= decay_rate(J, eps_ij, Gamma, omega_0, T, time_units)*sig_ij*IJ
                Z_dag += decay_rate(J, eps_ij, Gamma, omega_0, T, time_units).conjugate()*sig_ji*JI
    L = spre(sig_x*Z) - sprepost(Z, sig_x) + spost(Z_dag*sig_x) - sprepost(sig_x, Z_dag)
    print "It took ", time.time()-ti, " seconds to build the non-RWA Liouvillian"
    return -L
"""

def L_nonsecular(H_vib, A, eps, Gamma, T, J, time_units='cm'):
    #Construct non-secular liouvillian
    ti = time.time()
    d = H_vib.shape[0]
    evals, evecs = H_vib.eigenstates()
    X1, X2, X3, X4 = 0, 0, 0, 0
    for i in range(int(d)):
        for j in range(int(d)):
            eps_ij = abs(evals[i]-evals[j])
            A_ij = A.matrix_element(evecs[i].dag(), evecs[j])
            A_ji = (A.dag()).matrix_element(evecs[j].dag(), evecs[i])
            Occ = Occupation(eps_ij, T, time_units)
            IJ = evecs[i]*evecs[j].dag()
            JI = evecs[j]*evecs[i].dag()
            # 0.5*np.pi*alpha*(N+1)
            if abs(A_ij)>0 or abs(A_ji)>0:
                r_up = 2*np.pi*J(eps_ij, Gamma, eps)*Occ
                r_down = 2*np.pi*J(eps_ij, Gamma, eps)*(Occ+1)
                X3+= r_down*A_ij*IJ
                X4+= r_up*A_ij*IJ
                X1+= r_up*A_ji*JI
                X2+= r_down*A_ji*JI

    L = spre(A*X1) -sprepost(X1,A)+spost(X2*A)-sprepost(A,X2)
    L+= spre(A.dag()*X3)-sprepost(X3, A.dag())+spost(X4*A.dag())-sprepost(A.dag(), X4)
    print "It took ", time.time()-ti, " seconds to build the Non-secular RWA Liouvillian"
    return -0.5*L

def L_vib_lindblad(H_vib, A, eps, Gamma, T, J, time_units='cm'):
    '''
    Initially assuming that the vibronic eigenstructure has no
    degeneracy and the secular approximation has been made
    '''
    ti = time.time()
    d = H_vib.shape[0]
    ti = time.time()
    L = 0
    eig = H_vib.eigenstates()
    eVals = eig[0]
    eVecs = eig[1] # come out like kets
    l = 0
    occs=[]
    for i in range(int(d)):
        l = 0
        for j in range(int(d)):
            t_0 = time.time() # initial time reference for tracking slow calculations
            lam_ij = A.matrix_element(eVecs[i].dag(), eVecs[j])
            #lam_mn = (A.dag()).matrix_element(eVecs[n].dag(), eVecs[m])
            lam_ij_sq = lam_ij*lam_ij.conjugate()
            eps_ij = abs(eVals[i]-eVals[j])
            if lam_ij_sq>0:
                IJ = eVecs[i]*eVecs[j].dag()
                JI = eVecs[j]*eVecs[i].dag()
                JJ = eVecs[j]*eVecs[j].dag()
                II = eVecs[i]*eVecs[i].dag()

                Occ = Occupation(eps_ij, T, time_units)
                r_up = 2*np.pi*J(eps_ij, Gamma, eps)*Occ
                r_down = 2*np.pi*J(eps_ij, Gamma, eps)*(Occ+1)

                T1 = r_up*spre(II)+r_down*spre(JJ)
                T2 = r_up.conjugate()*spost(II)+r_down.conjugate()*spost(JJ)
                T3 = (r_up*sprepost(JI, IJ)+r_down*sprepost(IJ,JI))
                L += lam_ij_sq*(0.5*(T1 + T2) - T3)
                l+=1

    print "It took ", time.time()-ti, " seconds to build the vibronic Lindblad Liouvillian"
    return -L

def L_EM_lindblad(splitting, col_em, Gamma, T, J, time_units='cm'):
    # col_em is collapse operator
    ti = time.time()
    L = 0
    EMnb = Occupation(splitting, T, time_units)
    L+= 2*np.pi*J(splitting, Gamma, splitting)*(EMnb+1)*(sprepost(col_em, col_em.dag())-0.5*(spre(col_em.dag()*col_em) +spost(col_em.dag()*col_em)))
    L+= 2*np.pi*J(splitting, Gamma, splitting)*EMnb*(sprepost(col_em.dag(), col_em)-0.5*(spre(col_em*col_em.dag())+ spost(col_em*col_em.dag())))
    print "It took ", time.time()-ti, " seconds to build the electronic-Lindblad Liouvillian"
    return L
