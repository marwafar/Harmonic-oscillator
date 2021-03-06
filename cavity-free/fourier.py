import numpy as np

#-------------------------------------
resp=np.genfromtxt("response.txt", delimiter=" ")

n_row,n_colm=resp.shape
n_step=n_row

dephas=30
fstau=41.34137333656
print(n_row,n_colm)
S=(resp[:,1]+resp[:,2]*1j)*complex(0.0,1.0)*np.exp(-resp[:,0]**2/(2*dephas**2))
f=open("response_gauss.txt","w+")
for i in range(n_row):
 f.write(str(i)+" "+str(np.real(S[i]))+" "+str(np.imag(S[i]))+ " "+str(abs(S[i]))+"\n")

print(S.shape)
# Compute the fourier transform.
 
freq=np.fft.ifft(S)

dt=0.5
dw=(3.335*10000/(n_step*dt))
spect=open("linear_spectra.txt","w+")
for i in range(n_step):
 freq_i=(i*dw)*0.0001239
 spect.write(str(freq_i)+str(np.real(freq[i]))+" "+str(np.imag(freq[i]))+" "+str(abs(freq[i]))+"\n")

spect.close()
f.close()
