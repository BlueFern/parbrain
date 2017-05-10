import matplotlib.pyplot as plt
import numpy as np


t = []
R2072 = []
R2518 = []
R2258 = []
R2254 = []
Kp2072 = []
Kp2518 = []
Kp2258 = []
Kp2254 = []

#
R2072, Kp2072, t = np.loadtxt('test/cell20720.0.csv', delimiter=',', unpack=True)
R2518, Kp2518, t = np.loadtxt('test/cell25180.0.csv', delimiter=',', unpack=True)
R2258, Kp2258, t = np.loadtxt('test/cell22580.0.csv', delimiter=',', unpack=True)
R2254, Kp2254, t = np.loadtxt('test/cell22540.0.csv', delimiter=',', unpack=True)
#
##NKs = NKs*64434963 - 2000
#
plt.plot(t, R2072, label="Cell 2072")
plt.plot(t, R2518, label="Cell 2518")
plt.plot(t, R2258, label="Cell 2258")
plt.plot(t, R2254, label="Cell 2254")

#plt.plot(t, Kp2072, label="Cell 2072")
#plt.plot(t, Kp2518, label="Cell 2518")
#plt.plot(t, Kp2258, label="Cell 2258")
#plt.plot(t, Kp2254, label="Cell 2254")

plt.legend(loc=1)

plt.xlabel('Time (sec)')
plt.ylabel('Radius ($\mu$m)')

#plt.ylabel('Perivascular [K+] ($\mu$M)')
#
#
fig = plt.gcf()
fig.set_size_inches(4,6)
plt.savefig("CSD_Rplots.svg")


#Kp2072 = []
#Kp2518 = []
#Kp2258 = []
#Kp2254 = []
#
#t, Kp2072 = np.loadtxt('csvFiles/Kpcell2072.csv', delimiter=',', unpack=True)
#t, Kp2518 = np.loadtxt('csvFiles/Kpcell2518.csv', delimiter=',', unpack=True)
#t, Kp2258 = np.loadtxt('csvFiles/Kpcell2258.csv', delimiter=',', unpack=True)
#t, Kp2254 = np.loadtxt('csvFiles/Kpcell2254.csv', delimiter=',', unpack=True)
#
##NKs = NKs*64434963 - 2000
#
#plt.plot(t, Kp2072/1e3, label="Cell 2072")
#plt.plot(t, Kp2518/1e3, label="Cell 2518")
#plt.plot(t, Kp2258/1e3, label="Cell 2258")
#plt.plot(t, Kp2254/1e3, label="Cell 2254")
#
#plt.legend(loc=1)
#
#plt.xlabel('Time (sec)')
#plt.ylabel('Perivascular K+ (mM)')
##plt.ylabel('ECS [K+] ($\mu$M)')
#
#fig = plt.gcf()
#fig.set_size_inches(4,6)
#plt.savefig("CSD_Kpplots.svg")

#t = []
#Ke2072 = []
#
#t, Ke2072 = np.loadtxt('csvFiles/KeCell2072.csv', delimiter=',', unpack=True)
#
##NKs = NKs*64434963 - 2000
#
#plt.plot(t, Ke2072)
#
#
#plt.legend(loc=1)
#
#plt.xlabel('Time (sec)')
#plt.ylabel('Extracellular K+ (mM)')
##plt.ylabel('ECS [K+] ($\mu$M)')
#
#
#fig = plt.gcf()
##fig.set_size_inches(4,6)
#plt.savefig("CSD_Keplot2.svg")


