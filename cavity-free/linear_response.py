import propg as pg
import numpy as np
#-----------------------------
def response(ket_coef,bra_coef):
 # Compute the response function.

 R_3_t=complex(0.0,1.0)*(np.dot(ket_coef,np.conj(bra_coef)))

 return R_3_t
#------------------------------

if __name__== "__main__":
 # The parameters of the molecular Hamiltonian.
 tot_grid=21
 ev_to_h=27.21138386
 freq=0.185/ev_to_h
 w_eg=2.0/ev_to_h
 d=0.87
 reorg_erg=0.5*d*d*freq


 # The parameters for the propagation.
 dt=0.01*41.34137333656
 dt_fs=0.01
 output=50
 t_step=25000

#------------------------------
 # Compute the initial wavefunction.
 erg,coef=pg.init_wf(tot_grid,freq)

 bra_gs=np.zeros((tot_grid),dtype=complex)
 bra_ex=np.zeros((tot_grid),dtype=complex)
 ket_gs=np.zeros((tot_grid),dtype=complex)
 ket_ex=np.zeros((tot_grid),dtype=complex)

 temp=np.zeros((tot_grid),dtype=complex)
#------------------------------------------------

#------------------------------------------------
# Run the propagation.

 l_file=open("response.txt", "w+")

 for step in range(t_step):
  time_fs=dt_fs*step
  out=int(step/output)*output

  if step==0:
   ket_ex[:]=coef[:,0].astype(complex)
   bra_gs[:]=coef[:,0].astype(complex)

  if step==out:
   for i in range(tot_grid):
    temp[i]=1.0*ket_ex[i]+0.0*ket_gs[i]
   correlation=response(temp[:],bra_gs[:])
   l_file.write(str(time_fs)+" "+str(np.real(correlation))+ " "+str(np.imag(correlation))+" "+ str(abs(correlation))+ "\n")

  # Compute the bra and ket at time t during the propagation.
  bra_gs,bra_ex=pg.EOF(dt,bra_gs,bra_ex,tot_grid,freq,w_eg,reorg_erg,d)
  ket_gs,ket_ex=pg.EOF(dt,ket_gs,ket_ex,tot_grid,freq,w_eg,reorg_erg,d)

 l_file.close()

