# Compute the total spectra as sum of R and NR.
#---------------------------------
import numpy as np

#-------------------------------------
# Read the file for GB and SE

R=np.genfromtxt("spectra_rephase_T0-zp.txt", delimiter=" ")
NR=np.genfromtxt("spectra_non_rephase_T0-zp.txt", delimiter=" ")

n_row,n_colm=R.shape
n_step=int(np.sqrt(n_row))
dt=0.5

# compute the pure absorption spectra
S_R=R[:,2]+R[:,3]*1j
S_NR=NR[:,2]+NR[:,3]*1j

pur_abs=S_R+S_NR
print(pur_abs.shape)

# Frequency in eV.

dw=3.335*10000/(n_step*dt)
spect=open("2DES_T0.txt","w+")
k=0
for i in range(n_step):
 for j in range(n_step):
  freq_1=i*dw*0.0001239
  freq_2=j*dw*0.0001239
  spect.write(str(freq_1)+" "+str(freq_2)+" "+str(np.real(pur_abs[k]))+" "+str(np.imag(pur_abs[k]))+" "+str(abs(pur_abs[k]))+"\n")
  k+=1
 spect.write("\n")

spect.close()
