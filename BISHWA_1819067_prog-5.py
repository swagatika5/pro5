#specific heat of solids with temperature by i) Dulong-petit's law ii) Einstein Theory iii) Debye Theory

import numpy as np
import matplotlib.pyplot as plt

#Dulong-Petit's law
def Dulong_petit(T,N):
 return 3*N*k    

def Einstein_distribution(T,N,thetaE):   
 Cv = 3*N*k*((thetaE/T)**2)*np.exp(thetaE/T)/((np.exp(thetaE/T)-1)**2)
 return Cv   

def integralD(x):
    integral = (np.exp(x)*(x**4))/((np.exp(x)-1)**2)
    return integral

def integration(T,thetaD):    #Using Simpson's Rule
    h=0.01    #step-size
    summation=0.0
    upperlim=thetaD/T
    n = round((upperlim-0.0)/(2*h))
    #x = hplanck*oscfrequency/kb*T
    x = np.linspace((10**-5),upperlim,n)
    for i in range(n-1): 
      summation=summation+(h/3)*(integralD(x[i])+4*integralD(x[i]+h)+integralD(x[i]+2*h))
    return summation

def Debye_distribution(T,N,thetaD): 
 Cv = 9*N*k*((T/thetaD)**3)*(integration(T,thetaD))   
 return Cv     

def plotgraph(x,y,xlabelstr,ylabelstr,titlestr):
        plt.plot(x, y, lw=3)
        plt.xlabel(xlabelstr, size=12)
        plt.ylabel(ylabelstr, size=12)
        plt.title(titlestr)    
    
N = 6.022e23  #No. of atoms(one mole)
p = 10.0
q = 2000.0
print("Plotting Molar specific heat vs Temperature(10-2000 Kelvin)")
points = 500
x = np.linspace(p,q,points)

pi= 4*np.arctan(1.0)
k = 8.617*(10**-5)    #eV.K^-1
hplanck = 4.136*(10**-15)   #eV⋅Hz−1   

y1 = np.zeros(points)
for i in range(points):
  y1[i] = Dulong_petit(x[i],N)
plotgraph(x,y1,"Temperature","Molar Specific heat","Dulong-Petit law") 

thetaE = float(input("Enter the Einstein temperature of the solid: "))
# thetaE = hbar*oscangfreq/kb where oscangfreq is the natural angular frequency of vibration of a single atom in the solid.  
y2 = np.zeros(points)
for i in range(points):
  y2[i] = Einstein_distribution(x[i],N,thetaE)  
plotgraph(x,y2,"Temperature","Molar Specific heat","Einstein distribution function") 

#oscfrequency = thetaE*kb/hplanck
thetaD = float(input("Enter the Debye temperature of the solid: "))
#Debye_freq = ((9/(4*pi))*atoms_per_unit_vol/(2/pow(vel_trans,3)+1/pow(vel_long,3)))**(1/3)
#thetaD = hplanck*Debye_freq/kb
y3 = np.zeros(points)
for i in range(points):
  y3[i] = Debye_distribution(x[i],N,thetaD)
plotgraph(x,y3,"Temperature","Molar Specific heat","Debye distribution function")   

plt.plot(x,y1)
plt.plot(x,y2)
plt.plot(x,y3)
plt.legend(["Dulong-Petit's Law","Einstein Theory","Debye Theory"]) 
plt.show()