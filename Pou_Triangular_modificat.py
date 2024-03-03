# ----------------------------
# Quantum Potential Well
# ATOMIC UNITS
# V=4 for x<-2 or x>2
# ----------------------------
# Finite differences method as developed by Truhlar JCP 10 (1972) 123-132
#
# code by Jordi Faraudo
#
#
import numpy as np
import matplotlib.pyplot as plt
print("Introdueix el valor de la mida del pou en unitats atòmiques(Valor recomanat: 4):")
Lpou=float(input())
print("Introdueix el valor del camp elèctric en unitats atòmiques(Valor recomanat: 2):")
campelectr=float(input())
#Potencial triangular
def Velectr(x):
    V = campelectr*x
    return V
#Potential as a function of position
def getV(x):
    if (abs(x)<(Lpou/2)):
       potvalue = 0.0 +Velectr(x)
    else:
       potvalue = 4.0 
    return potvalue

#Discretized Schrodinger equation in n points (FROM 0 to n-1)
def Eq(n,h,x):
    F = np.zeros([n,n])
    for i in range(0,n):
        F[i,i] = -2*((h**2)*getV(x[i]) + 1)
        if i > 0:
           F[i,i-1] = 1
           if i < n-1:
              F[i,i+1] = 1
    return F

#-------------------------
# Main program
#-------------------------
# Interval for calculating the wave function [-L/2,L/2]
L = 3*Lpou
xlower = -L/2.0
xupper = L/2.0

#Discretization options
print("Introdueix la distancia entre punts(Valor recomanat: 0.1):")
print("Cal ser coherent amb la distancia del pou!!!")
h = float(input())  #discretization in space

#Create coordinates at which the solution will be calculated
x = np.arange(xlower,xupper+h,h)
#grid size (how many discrete points to use in the range [-L/2,L/2])
npoints=len(x)

print("Using",npoints, "grid points.")

#Calculation of discrete form of Schrodinger Equation
print("Calculating matrix...")
F=Eq(npoints,h,x)

#diagonalize the matrix F
print("Diagonalizing...")
eigenValues, eigenVectors = np.linalg.eig(F)

#Order results by eigenvalue
# w ordered eigenvalues and vs ordered eigenvectors
idx = eigenValues.argsort()[::-1]
w = eigenValues[idx]
vs = eigenVectors[:,idx]

#Energy Level
E = - w/(2.0*h**2)

#Print Energy Values
print("RESULTS:")
for k in range(0,4):
  Erel=E[k]/2.0  #ratio between energy level and potential well
  print("State ",k,": Energy = %.4f" %E[k],', E/V pou quadrat = '+'{:.4f}'.format(Erel))

#Init Wavefunction (empty list with npoints elements)
psi = [None]*npoints

#Calculation of normalised Wave Functions
for k in range(0,len(w)):
	psi[k] = vs[:,k]
	integral = h*np.dot(psi[k],psi[k])
	psi[k] = psi[k]/integral**0.5

#Plot Wave functions
print("Plotting")

#v = int(input("\n Quantum Number (enter 0 for ground state):\n>"))
Pou=np.zeros(npoints)
for i in range(npoints):
    Pou[i]=getV(x[i])
for v in range(0,5):
    plt.plot(x,psi[v],label=r'$\psi_v(x)$, k = ' + str(v))
    plt.plot(x,0.5*Pou)
    plt.title(r'$n=$'+ str(v) + r', $E_n$=' + '{:.4f}'.format(E[v]))
    plt.legend()
    plt.xlabel(r'$x$ (a.u.)')
    plt.ylabel(r'$\psi(x)$')
    plt.show()
#Plot de totes les funcions d'ona amb shift corresponent a l'energia (Factor 0.5 per a que el plot es pugui veure millor), el plot és un anàlisi qualitatiu.
#No ens aporta cap informació extra, ja que el shift és arbitrari
for v in range(0,5):
    plt.plot(x,psi[v] +0.5*E[v],label=r'$\psi_v(x)$, k = ' + str(v))
    plt.plot(x,0.5*Pou)
    plt.title(r'$n=$'+ str(v) + r', $E_n$=' + '{:.4f}'.format(E[v]))
    plt.legend()
    plt.xlabel(r'$x$ (a.u.)')
    plt.ylabel(r'$\psi(x)$')
plt.show()

print("Bye")
