import numpy as np

#-------------------------------------
resp=np.genfromtxt("response.txt", delimiter=" ")

n_row,n_colm=resp.shape
n_step=n_row

dephas=35
fstau=41.34137333656
print(n_row,n_colm)
S= np.zeros((3*n_row), dtype=complex)


S_n=(resp[:,1]+resp[:,2]*1j)*complex(0.0,1.0)*np.exp(-resp[:,0]**2/(2*dephas**2))
f=open("response_gauss.txt","w+")
for i in range(n_row):
# if i >= n_row:
#  S[i]=complex(0.0,0.0)
# else:
#  S[i]=S_n[i]
 f.write(str(i)+" "+str(np.real(S_n[i]))+" "+str(np.imag(S_n[i]))+ " "+str(abs(S_n[i]))+"\n")

print(S_n.shape)
# Compute the fourier transform.
#for i in range(n_row,2*n_row):
 

#total_n=S_n.shape
#S=np.zeros((2*n_row), dtype=complex)
S[:n_row]=S_n 

freq=np.fft.ifft(S)

dt=0.5
dw=(3.335*10000/(3*n_step*dt))
spect=open("linear_spectra-zp.txt","w+")
for i in range(3*n_step):
 freq_i=(i*dw)*0.0001239
 spect.write(str(freq_i)+str(np.real(freq[i]))+" "+str(np.imag(freq[i]))+" "+str(abs(freq[i]))+"\n")

spect.close()
f.close()
