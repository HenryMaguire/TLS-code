"""
Resonant (fairly broad ) underdamped spectral density
"""
N = 10
eps = 0.1*8065.5 # TLS splitting
#eps = 2.*1519.3 # ps
T_EM = 6000. # Optical bath temperature
H_S = eps*E*E.dag()
#alpha_EM = 0.3 # System-bath strength (optical)
Gamma_EM = 1.*5.309 #bare decay of electronic transition from inv. ps to in inv. cm
#Gamma_EM = 6.582E-7*1519.3
T_ph = 300. # Phonon bath temperature
overdamped = False
#wc = 53.*0.188
w0 = eps # underdamped SD parameter omega_0
#w0 = 200.*0.188
Gamma = 1000. # Width of distribution
alpha_ph = 1000/np.pi #0.05*eps/np.pi# Ind.-Boson frame coupling
beta = 1/(0.695*T_ph)
wc = 53. # Ind.-Boson frame phonon cutoff freq
mu = 1.


N = 10
eps = 0.1*8065.5 # TLS splitting
#eps = 2.*1519.3 # ps
T_EM = 6000. # Optical bath temperature
H_S = eps*E*E.dag()
#alpha_EM = 0.3 # System-bath strength (optical)
Gamma_EM = 1.*5.309 #bare decay of electronic transition from inv. ps to in inv. cm
#Gamma_EM = 6.582E-7*1519.3
T_ph = 300. # Phonon bath temperature
overdamped = true
#wc = 53.*0.188
w0 = 600 #eps # underdamped SD parameter omega_0
#w0 = 200.*0.188
Gamma = 1000. # Width of distribution
alpha_ph = 100/np.pi #0.05*eps/np.pi# Ind.-Boson frame coupling
beta = 1/(0.695*T_ph)
wc = 53. # Ind.-Boson frame phonon cutoff freq
mu = 1.
