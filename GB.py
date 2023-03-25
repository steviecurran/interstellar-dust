#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import os
import sys
from scipy.optimize import curve_fit

# universal constants and units
pi = 3.14159265359
c = 2.99792458e8        # m/s
h = 6.62606957e-34      # J s
k = 1.3806488e-23       # J/K
pc = 3.08567758e16      # m
Msun = 1.9891e30        # kg
W = 2.8977729E-3 # Wien's displacement constant

lognu, logL  = np.loadtxt('FIR_range.dat', unpack=True)
nu = 10**lognu # conversion to linear
L = 10**logL

ymax = np.max(L)
xmax = np.max(nu)
testnu = xmax 
T_est = (W/c)*xmax
testflux = 2.0*pi*h/(c**2) * (xmax**3)/(np.exp(h*testnu/(k*xmax)) -1)
norm = ymax/testflux # *testflux # approx scaling required

print('Initial guesses are T = %1.1f K nu = %1.2e Hz norm = %1.3e' %(T_est,testnu,norm))

def fit_func(nu,B,T): #  x DATA THEN GUESSES
   return norm*2.0*pi*h/(c**2) * nu**(3+B)/(np.exp(h*nu/(k*T)) -1)
 
xdata = nu
ydata = L #/norm

p0 = [1.5,T_est] # guess for B and T
parameters, covariance = curve_fit(fit_func, xdata, ydata,p0) #,method = 'dogbox') 

print('Fit is', parameters) 
data= parameters.T # does bugger all, because string?

#### WRITE TO FILE FOR C CODE ######
np.savetxt('GB.dat', data, fmt='%s') # 
f = open("GB.dat", "a")
f.write(str(norm))      # ADD norm TO FIT DATA
f.close()
############ PLOT #############
#print(xdata,ydata)
x=np.log10(xdata)
y=np.log10(ydata)
nor = np.log10(norm)
print(x,y)

B = parameters[0]
T = parameters[1]

#To calculate the standard error of the parameters from the covariance, you take the square root of the diagonal elements of the matrix. You can do this in one line using functions from numpy.
SE = np.sqrt(np.diag(covariance))
dB = SE[0]
dT = SE[1]

print('T = %1.2f +/- %1.2f and B = %1.2f +/- %1.2f' %(T, dT, B,dB))

xplot = np.linspace(1e11, 1e13,100) 
fit_y = fit_func(xplot,B,T) # BUT WANT IN LOG SPACE
fit_y  = np.log10(fit_y)
xplot = np.log10(xplot)
#print(fit_y)
#print(xplot)

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(6,4))
ax = plt.gca()
#ax.set_xscale('log'); ax.set_yscale('log')

ax.set_xlim([11.65, 13.15])
ax.set_ylim([23, 26.5])

ax.plot(x,y,'o', markersize=5, c = 'k') #points
ax.plot(xplot, fit_y, '-', c = 'r') #, label='fit')
text = "T$_{dust}$ = %1.0f $\pm$ %1.0f K, $\\beta$ = %1.2f $\pm$ %1.2f " % (T,dT,B,dB) # not actually latex
ax.text(11.7,26.3, text, fontsize = 14, horizontalalignment='left',verticalalignment='top') 
plt.xlabel("Frequency, log$_{10}\\nu$ [Hz]")
plt.ylabel("Luminosity, log$_{10}L_{\\nu}$ [W Hz$^{-1}$]")
png = "GB_T=%1.0f_B=%1.2f.png" % (T,B) #eps = "convert %s %s.eps" % (png, plot)
plt.tight_layout()
plt.savefig(png);print("Plot written to", png)
plt.show()

