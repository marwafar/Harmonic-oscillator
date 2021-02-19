import numpy as np

#-------------------------------------
# Read the file for GB and SE

bleach=np.genfromtxt("response_GB_NR.txt", delimiter=" ")
emission=np.genfromtxt("response_SE_NR.txt", delimiter=" ")
absorp=np.genfromtxt("response_ESA_NR.txt", delimiter=" ")

print(bleach.shape)

n_row,n_colm=bleach.shape
n_step=int(np.sqrt(n_row))

dephas=35
dt=0.5
# compute the response function
S_GB=(bleach[:,2]+bleach[:,3]*1j)*complex(0.0,1.0)*np.exp(-(bleach[:,0]+bleach[:,1])**2/(2*dephas**2))
S_SE=(emission[:,2]+emission[:,3]*1j)*complex(0.0,1.0)*np.exp(-(bleach[:,0]+bleach[:,1])**2/(2*dephas**2))
S_ESA=(absorp[:,2]+absorp[:,3]*1j)*complex(0.0,-1.0)*np.exp(-(bleach[:,0]+bleach[:,1])**2/(2*dephas**2))

response=S_GB+S_SE+S_ESA

#print(S_GB[1],bleach[1,2],bleach[1,3])
#print(S_GB.shape)

# Zero padding.

S=np.zeros((3*n_step*3*n_step), dtype=complex)
k=0
f=0
for i in range(3*n_step):
 for j in range(3*n_step):
  if i<n_step and j<n_step:
   S[k]=response[f]
   f+=1
  k+=1


n_step=3*n_step

# Compute the fourier transform.
temp=np.zeros((n_step),dtype=complex)
BW=np.zeros((n_step,n_step), dtype=complex)
k=0

for i in range(n_step):
 for j in range(n_step):
  temp[j]=S[k]
  k+=1
 BW[i,:]=np.fft.ifft(temp)

FW=np.zeros((n_step,n_step), dtype=complex)

for j in range(n_step):
 for i in range(n_step):
  temp[i]=BW[i,j]
 FW[:,j]=np.fft.ifft(temp)

spect=open("spectra_non_rephase_T0-zp.txt","w+")
for i in range(n_step):
 for j in range(n_step):
  spect.write(str(i)+" "+str(j)+" "+str(np.real(FW[i,j]))+" "+str(np.imag(FW[i,j]))+" "+str(abs(FW[i,j]))+"\n")
# spect.write("\n")

spect.close()

