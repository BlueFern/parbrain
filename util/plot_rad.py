import matplotlib.pyplot as plt
import numpy as np

#t = []
#r_sig = []
#r_bor = []
#r_nosig = []
#
#t, r_sig = np.loadtxt('rad_signal_j04.csv', delimiter=',', unpack=True)
##t, r_bor = np.loadtxt('rad_border_j04.csv', delimiter=',', unpack=True)
#t, r_nosig = np.loadtxt('rad_nosignal_j04.csv', delimiter=',', unpack=True)
#
#
#plt.plot(t[50:250],r_sig[50:250], ':',label="Stimulation")
##plt.plot(t[50:250],r_bor[50:250], label="Boundary")
#plt.plot(t[50:250], r_nosig[50:250], '--r', label="Non-stimulation")
#
##plt.legend(loc=1) # ['Neuronal input', 'Borderline', 'No input'])
#
#plt.xlabel('t')
#plt.ylabel('Radius')
#plt.show()




t = []
r_sig_new = []
r_sig_old = []
r_nosig_new = []
r_nosig_old = []

t, r_sig_new = np.loadtxt('rad_signal_j04.csv', delimiter=',', unpack=True)
t, r_sig_old = np.loadtxt('old_rad_signal_j04.csv', delimiter=',', unpack=True)
t, r_nosig_new = np.loadtxt('rad_nosignal_j04.csv', delimiter=',', unpack=True)
t, r_nosig_old = np.loadtxt('old_rad_nosignal_j04.csv', delimiter=',', unpack=True)

plt.plot(t[50:300],r_sig_new[50:300],label="Stim new")
plt.plot(t[50:300],r_sig_old[50:300], label="Stim old")
plt.plot(t[50:250], r_nosig_new[50:250], label="Non-stim new")
plt.plot(t[50:250], r_nosig_old[50:250], label="Non-stim old")

plt.legend(loc=1) # ['Neuronal input', 'Borderline', 'No input'])

plt.xlabel('t')
plt.ylabel('Radius')
plt.show()





