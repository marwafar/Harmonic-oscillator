import propg as pg
import numpy as np
#-----------------------------
def response(ket_coef,bra_coef):
 # Compute the response function.

 R_3_t=complex(0.0,1.0)*(np.dot(ket_coef,np.conj(bra_coef)))

 return R_3_t
#-----------------------------
if __name__== "__main__":
 # The parameters of the molecular Hamiltonian.
 tot_grid=21
 ev_to_h=27.21138386
 freq=0.185/ev_to_h
 w_eg=2.0/ev_to_h
 d=0.87
 reorg_erg=0.5*d*d*freq
 

 # The parameter of the cavity.
 n_fock=5
 freq_cavity=w_eg+reorg_erg
 g=0.08/ev_to_h

 # The parameters for the propagation.
 dt=0.01*41.34137333656
 dt_fs=0.01
 output=50

 # Parameter for the 2D
 t2_step=15000
#---------------------------------------
 # Compute the initial wavefunction.
 erg,coef=pg.init_wf(tot_grid,freq)

 bra_gs=np.zeros((tot_grid,n_fock),dtype=complex)
 bra_ex=np.zeros((tot_grid,n_fock),dtype=complex)
 ket_gs=np.zeros((tot_grid,n_fock),dtype=complex)
 ket_ex=np.zeros((tot_grid,n_fock),dtype=complex)

 temp=np.zeros((tot_grid,n_fock),dtype=complex)
#------------------------------------------------
# Run the propagation.

 GB_file=open("response_GB_R.txt", "w+")
 dipol=np.zeros((n_fock,n_fock))
 for i in range(n_fock):
  dipol[i,i]=1.0

 for t in range(0,t2_step,output):

  wait=0
  t1=0
  t2=dt_fs*t
  t3=wait+t2

  t_time=150
  total=t_time+t3

  n_step=int(total/dt_fs)

  for step in range(n_step):
   time_fs=dt_fs*step
   out=int(step/output)*output

   # 1- compute the WF after interaction with laser at time t=0.
   if time_fs == t1:
    bra_ex[:,0]=coef[:,0].astype(complex)
    ket_gs[:,0]=coef[:,0].astype(complex)
   
   # 2- compute the WF after interaction with laser at time=t2.
   if time_fs==t2:
    for i in range(tot_grid):
     bra_gs[i,:]=np.dot(dipol[:,:],bra_ex[i,:])
     bra_ex[i,:]=complex(0.0,0.0)

   # 3- compute the WF after interaction with laser at time t=t3.
   if time_fs==t3:
    for i in range(tot_grid):
     ket_ex[i,:]=np.dot(dipol[:,:],ket_gs[i,:])
     ket_gs[i,:]=complex(0.0,0.0)

   # 4- compute the response function.
   if (time_fs>=t3) and (step==out):
    for i in range(tot_grid):
     temp[i,:]=np.dot(dipol[:,:],ket_ex[i,:])
    correlation=response(temp[:,0],bra_gs[:,0])

    GB_file.write(str(t2)+" "+str(time_fs)+" "+str(np.real(correlation))+ " "+str(np.imag(correlation))+" "+ str(abs(correlation))+ "\n")

   # Compute the bra and ket at time t during the propagation.
   bra_gs,bra_ex=pg.EOF(dt,bra_gs,bra_ex,tot_grid,freq,w_eg,reorg_erg,d,g,n_fock,freq_cavity)
   ket_gs,ket_ex=pg.EOF(dt,ket_gs,ket_ex,tot_grid,freq,w_eg,reorg_erg,d,g,n_fock,freq_cavity)

 GB_file.close()
