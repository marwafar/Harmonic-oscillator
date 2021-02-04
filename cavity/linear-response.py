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
 reorg_erg=0.5*0.14/ev_to_h
 d=0.87

 # The parameter of the cavity.
 n_fock=5
 freq_cavity=1.9/ev_to_h
 g=0.08/ev_to_h

 # The parameters for the propagation.
 dt=0.01*41.34137333656
 dt_fs=0.01
 output=50
 t_step=20000
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

 l_file=open("response.txt", "w+")
 dipol=np.zeros((n_fock,n_fock))
 for i in range(n_fock):
  dipol[i,i]=1.0
  
 for step in range(t_step):
  time_fs=dt_fs*step
  out=int(step/output)*output
  
  if step==0:
   ket_ex[:,0]=coef[:,0].astype(complex)
   bra_gs[:,0]=coef[:,0].astype(complex)

  if step==out:
   for i in range(tot_grid):
    temp[i,:]=np.dot(dipol[:,:],ket_ex[i,:])
   correlation=response(temp[:,0],bra_gs[:,0])
   l_file.write(str(time_fs)+" "+str(np.real(correlation))+ " "+str(np.imag(correlation))+" "+ str(abs(correlation))+ "\n")

  # Compute the bra and ket at time t during the propagation.
  bra_gs,bra_ex=pg.EOF(dt,bra_gs,bra_ex,tot_grid,freq,w_eg,reorg_erg,d,g,n_fock,freq_cavity)
  ket_gs,ket_ex=pg.EOF(dt,ket_gs,ket_ex,tot_grid,freq,w_eg,reorg_erg,d,g,n_fock,freq_cavity)
 
 l_file.close()
