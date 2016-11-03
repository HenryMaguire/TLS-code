"""
Description:
The four liouvillians I am studying.
- no rotating wave approximation
- no secular approximation
- a secular approximation
- an approximation which says that the enlarged system eigenstates are the same as the uncoupled system eigenstates

"""

import numpy as np
from qutip import destroy, tensor, qeye, spost, spre, sprepost

def L_non_rwa(H_vib):
    L = 0

    return L

def L_nonsecular(H_vib, sig, alpha, T):
    d = H_vib.shape[0]
    evals, evecs = H.eigenstates()
    X1, X2, X3, X4 = 0, 0, 0, 0
    for i in range(int(d)):
        for j in range(int(d)):
            eps_ij = abs((evals[i]-evals[j]).real)
            sig_ij = sig.matrix_element(evecs[i].dag(), evecs[j])
            sig_ji = (sig.dag()).matrix_element(evecs[j].dag(), evecs[i])
            #print sig_ji == sig_ij.conjugate()
            Occ = Occupation(eps_ij, T)
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

def L_EM_vib(H_vib, A, T, alpha_em):
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

                Occ = Occupation(eps_mn, T)
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

def L_EM_naive(splitting, col_em, T):
    L = 0
    EMnb = Occupation(splitting, T)
    L+= np.pi*(EMnb+1)*(sprepost(col_em, col_em.dag())-0.5*(spre(col_em.dag()*col_em) +spost(col_em.dag()*col_em)))
    L+= np.pi*(EMnb)*(sprepost(col_em.dag(), col_em)-0.5*(spre(col_em*col_em.dag())+ spost(col_em*col_em.dag())))
    return L
