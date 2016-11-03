"""
The four electromagnetic liouvillians I am studying:
- no rotating wave approximation
- no secular approximation
- a secular approximation
- an approximation which says that the enlarged system eigenstates are the same as the uncoupled system eigenstates

"""

import numpy as np
from qutip import destroy, tensor, qeye, spost, spre, sprepost

def Occupation(omega, T, time_units):
    conversion = 0.695
    if time_units == 'ps': # allows conversion to picoseconds
        conversion == 7.13
    else:
        pass
    n =0.
    if T ==0.:
        n = 0.
    else:
        beta = 1. / (conversion* T)
        n = float(1./(sp.exp(omega*beta)-1))
    return n

def Gamma_1(epsilon, N, alpha):
    return 0.5*np.pi*alpha*N

def Gamma_2(epsilon, N, alpha):
    return 0.5*np.pi*alpha*(N+1)

def decay_rate(omega, J, T, time_units):
    """
    Decay rate for non-secular master equation. In my notes, this is called Gamma.
    """
    conversion = 0.695
    if time_units == 'ps': # allows conversion to picoseconds
        conversion == 7.13
    else:
        pass

    beta = 0
    if T==0:
        beta = 10E12 # Some very large number, this seems dodgy.
    else:
        beta = 1. / (conversion* T)

    Gamma = 0
    coth = lambda x : np.cosh(x)/np.sinh(x)
    if omega>0:
        Gamma = np.pi*0.5*J(omega)*(coth(beta*omega/2)-1)
    elif omega ==0:
        Gamma = np.pi*J(1)/beta # I thought this was supposed to be just J(1), but mathematica says otherwise
    else:
        Gamma = 0.5*J(omega)*(coth(beta*omega/2)+1)
    return Gamma

def L_nonrwa(H_vib, sig_x, alpha, T, time_units='cm'):
    J = lambda x : alpha
    d = H_vib.shape[0]
    evals, evecs = H_vib.eigenstates()
    Z = 0 # initalise rate operator
    for i in range(int(d)):
        for j in range(int(d)):
            eps_ij = abs((evals[i]-evals[j]).real)
            sig_ij = sig_x.matrix_element(evecs[i].dag(), evecs[j])
            IJ = evecs[i]*evecs[j].dag()
            if sig_ij >0:
                Z+= decay_rate(eps_ij, J, T, time_units)*sig_ij
                print Z
    L = spre(sig_x*Z) - sprepost(Z, sig_x) + spost(Z.dag()*sig_x) - sprepost(sig_x, Z.dag())
    return -L

def L_nonsecular(H_vib, sig, alpha, T, time_units='cm'):
    d = H_vib.shape[0]
    evals, evecs = H.eigenstates()
    X1, X2, X3, X4 = 0, 0, 0, 0
    for i in range(int(d)):
        for j in range(int(d)):
            eps_ij = abs((evals[i]-evals[j]).real)
            sig_ij = sig.matrix_element(evecs[i].dag(), evecs[j])
            sig_ji = (sig.dag()).matrix_element(evecs[j].dag(), evecs[i])
            #print sig_ji == sig_ij.conjugate()
            Occ = Occupation(eps_ij, T, time_units)
            IJ = evecs[i]*evecs[j].dag()
            JI = evecs[j]*evecs[i].dag()
            if (sig_ij or sig_ji) >0:
                X1+= Gamma_1(eps_ij, Occ, alpha)*sig_ji*JI
                X2+= Gamma_2(eps_ij, Occ, alpha)*sig_ji*JI
                X3+= Gamma_2(eps_ij, Occ, alpha)*sig_ij*IJ
                X4+= Gamma_1(eps_ij, Occ, alpha)*sig_ij*IJ
    L = spre(sig*X1) -sprepost(X1,sig)+spost(X2*sig)-sprepost(sig,X2)
    L+= spre(sig.dag()*X3)-sprepost(X3, sig.dag())+spost(X4*sig.dag())-sprepost(sig.dag(), X4)
    return -L

def L_vib_lindblad(H_vib, A, T, alpha_em, time_units='cm'):
    '''
    Initially assuming that the vibronic eigenstructure has no
    degeneracy and the secular approximation has been made
    '''
    d = H_vib.shape[0]
    ti = time.time()
    L = 0
    eig = H_vib.eigenstates()
    eVals = eig[0]
    #print len(evals), len(evals[(d/2.)::])
    eVecs = eig[1] # come out like kets
    l = 0
    occs=[]
    for m in range(int(d)):
        #print " Liouvillian is ", (float(m)/H_vib.shape[0])*100, " percent complete after ", int(time.time()-ti), " seconds. (", int((time.time()-ti)/60.), " minutes)."
        #print "There were ", l, " non-zero contributions."
        l = 0
        for n in range(int(d)):
            t_0 = time.time() # initial time reference for tracking slow calculations
            lam_nm = A.matrix_element(eVecs[m].dag(), eVecs[n])
            #lam_mn = A.matrix_element(eVecs[n].dag(), eVecs[m])
            lam_nm_sq = lam_nm*lam_nm.conjugate()
            eps_mn = abs(eVals[m]-eVals[n])
            if lam_nm_sq>0:
                MN = eVecs[m]*eVecs[n].dag()
                NM = eVecs[n]*eVecs[m].dag()
                NN = eVecs[n]*eVecs[n].dag()
                MM = eVecs[m]*eVecs[m].dag()

                Occ = Occupation(eps_mn, T, time_units)
                g2 = Gamma_1(eps_mn, Occ, alpha_em)
                g1 = Gamma_2(eps_mn, Occ, alpha_em)

                #T1 = 0.5*Gamma_1(eps_mn, Occ)*(spre(NN) - 2*sprepost(MN, NM)) + 0.5*Gamma_2(eps_mn, Occ)*(spre(MM) - 2*sprepost(NM, MN))
                T1 = g2*lam_nm_sq*(spre(MM))+g1*lam_nm_sq*(spre(NN))
                T2 = g2.conjugate()*lam_nm_sq*(spost(MM))+g1.conjugate()*lam_nm_sq*(spost(NN))
                T3 = 2*(g2*lam_nm_sq*sprepost(NM, MN)+g1*lam_nm_sq*(sprepost(MN,NM)))
                L += (T1 + T2 - T3)
                l+=1

    print "It took ", time.time()-ti, " seconds to build the Liouvillian"
    #eMatrix = np.array(eMatrix).reshape((H_vib.shape[0], H_vib.shape[0]))
    #plt.imshow(eMatrix)
    return -L

def L_EM_lindblad(splitting, col_em, alpha, T, time_units='cm'):
    L = 0
    EMnb = Occupation(splitting, T, time_units)
    L+= np.pi*alpha*(EMnb+1)*(sprepost(col_em, col_em.dag())-0.5*(spre(col_em.dag()*col_em) +spost(col_em.dag()*col_em)))
    L+= np.pi*alpha*(EMnb)*(sprepost(col_em.dag(), col_em)-0.5*(spre(col_em*col_em.dag())+ spost(col_em*col_em.dag())))
    return L
