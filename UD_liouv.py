"""
In this script we have four methods.
1) Ham_RC builds the RC-frame Hamiltonian and system operators for both bath interactions.
    It takes in the system splitting, the RC frequency, system-RC coupling and the Hilbert space dimension.

2) RCME_operators builds the collapse operators for the RCME. It takes as input the
    Hamiltonian, system operator for phonon interaction, RC-residual coupling strength and beta (inv. T).

3) liouvillian_build... builds the RCME Liouvillian. Taking RC-residual coupling, RC-freq. and Temperature (for beta).
    It also takes a default parameter time_units which is set 'cm', giving frequencies in inverse cm.
    Can be set to 'ps' for inv. picoseconds.

4) RC_function_UD dresses up the Liouvillian with all the mapped RC frame parameters.
    It calculates these in accordance with Jake's initial paper.
"""
import numpy as np
import scipy as sp
from qutip import destroy, tensor, qeye, spre, spost, sprepost

def J_UD_SB(omega, alpha, omega_0, Gamma):
    return (alpha*Gamma*(omega_0**2))*omega/((omega_0**2 - omega**2)**2+((Gamma**2)*(omega**2)))

def Ham_RC(sigma, eps, Omega, kap, N):
    """
    Input: System splitting, RC freq., system-RC coupling and Hilbert space dimension
    Output: Hamiltonian, sigma_- and sigma_z in the vibronic Hilbert space
    """
    a = destroy(N)
    H_S = eps*tensor(sigma.dag()*sigma, qeye(N)) + kap*tensor(sigma.dag()*sigma, (a + a.dag()))+tensor(qeye(2),Omega*a.dag()*a)
    A_em = tensor(sigma, qeye(N))
    A_nrwa = tensor(sigma+sigma.dag(), qeye(N))
    A_ph = tensor(qeye(2), (a + a.dag()))
    return H_S, A_em, A_nrwa, A_ph

def RCME_operators(H_0, A, gamma, beta):
    # This function will be passed a TLS-RC hamiltonian, RC operator, spectral density and beta
    # outputs all of the operators needed for the RCME (underdamped)
    dim_ham = H_0.shape[0]
    Chi = 0 # Initiate the operators
    Xi = 0
    eVals, eVecs = H_0.eigenstates()
    #print H_0
    #ti = time.time()
    for j in range(dim_ham):
        for k in range(dim_ham):
            e_jk = eVals[j] - eVals[k] # eigenvalue difference
            A_jk = A.matrix_element(eVecs[j].dag(), eVecs[k])
            outer_eigen = eVecs[j] * (eVecs[k].dag())
            if sp.absolute(A_jk) > 0:
                if sp.absolute(e_jk) > 0:
                    #print e_jk
                    # If e_jk is zero, coth diverges but J goes to zero so limit taken seperately

                    Chi += 0.5*np.pi*e_jk*gamma * ((sp.cosh(e_jk * beta / 2)) / sp.sinh(e_jk * beta / 2))*A_jk*outer_eigen # e_jk*gamma is the spectral density
                    Xi += 0.5*np.pi*e_jk*gamma * A_jk * outer_eigen
                else:
                    Chi += (np.pi*gamma*A_jk/beta)*outer_eigen # Just return coefficients which are left over
                    #Xi += 0 #since J_RC goes to zero

    return H_0, A, Chi, Xi

def liouvillian_build(H_0, A, gamma, wRC, T_C, time_units='cm'):
    conversion = 0.695
    if time_units == 'ev':
        conversion == 8.617E-5
    if time_units == 'ps':
        conversion == 0.131
    else:
        pass

    beta_C = 0.
    if T_C == 0.0:
        beta_C = np.infty
        RCnb = 0
        print "Temperature is too low, this won't work"
    else:
        beta_C = 1./(conversion * T_C)
        RCnb = float(1. / (sp.exp( beta_C * wRC)-1))
    # Now this function has to construct the liouvillian so that it can be passed to mesolve
    H_0, A, Chi, Xi = RCME_operators(H_0, A, gamma, beta_C)
    L = 0
    L-=spre(A*Chi)
    L+=sprepost(A, Chi)
    L+=sprepost(Chi, A)
    L-=spost(Chi*A)

    L+=spre(A*Xi)
    L+=sprepost(A, Xi)
    L-=sprepost(Xi, A)
    L-=spost(Xi*A)

    return L

def RC_function_OD(sigma, eps, T_Ph, wc, wRC, alpha_ph, N, time_units='cm'):

    # we define all of the RC parameters by the underdamped spectral density
    Gamma = (wRC**2)/wc
    gamma = Gamma / (2. * np.pi * wRC)  # no longer a free parameter that we normally use to fix wRC to the system splitting
    kappa= np.sqrt(np.pi * alpha_ph * wRC / 2.)  # coupling strength between the TLS and RC
    print "SB cutoff= ",wc, "RC oscillator frequency=",wRC, " splitting =",eps, "residual bath coupling=", gamma, " N=",N, "TLS-RC coupling=", kappa
    H, A_em, A_nrwa, A_ph = Ham_RC(sigma, eps, wRC, kappa, N)
    L_RC =  liouvillian_build(H, A_ph, gamma, wRC, T_Ph, time_units)

    return L_RC, H, A_em, A_nrwa, wRC, kappa, Gamma
