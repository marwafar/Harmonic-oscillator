import numpy as np
from numpy import linalg as LA
#---------------------------------
def diag(M):
 
 # Diagonlaize the matrix
 E,V=LA.eigh(M)
 return E,V
#--------------------------------
def x_Hermit(tot_grid):
 # Compute the grid points using Hermit polynomial.

 x_eq=0.0
 x_ij=np.zeros((tot_grid,tot_grid))
 for row in range(tot_grid):
  for colm in range(tot_grid):
   x_ij[row,colm]=np.sqrt((row+1.0)/(2.0))*float(row==colm-1)\
                    +x_eq*float(row==colm)+np.sqrt(row/2.0)*float(row==colm+1)
 x_i,vect=diag(x_ij)
 return x_i,vect
#---------------------------------
def weight_Hermit(tot_grid):
 # Compute the weight for each grid point.

 x_eq=0
 x_i,vect=x_Hermit(tot_grid)
 w_i= ((1.0/np.pi)**(-0.25)*np.exp(0.5*x_i*x_i)*vect[0,:])**2

 return w_i
#----------------------------------
def second_derivavtive_Hermit(tot_grid):
 # Compute the second derivative

 x_i,vect=x_Hermit(tot_grid)
 k=np.array([i+0.5 for i in range(tot_grid)])
 
 dif2mat=-2.0*np.matmul(np.matmul(vect.T,np.diag(k)),vect)
 dif2mat+=(np.diag(x_i))**2


 return dif2mat
#-------------------------------------
def init_wf(tot_grid,freq):
 # Compute the ground state vibrational wavefunction.

 x_i,vect=x_Hermit(tot_grid)
 dif2mat=second_derivavtive_Hermit(tot_grid)

 pot=0.5*freq*x_i*x_i
 H_ij=-0.5*freq*dif2mat+np.diag(pot)
 energy,coef=diag(H_ij)

 return energy,coef
#--------------------------------------------
def KEO(tot_grid,coef_gs,coef_ex,freq):

 dif2mat=second_derivavtive_Hermit(tot_grid)
 
 KEO_gs=np.zeros(tot_grid,dtype=complex)
 KEO_ex=np.zeros(tot_grid,dtype=complex)

 KEO_gs=-0.5*freq*(dif2mat.dot(coef_gs))
 KEO_ex=-0.5*freq*(dif2mat.dot(coef_ex))

 return KEO_gs,KEO_ex
#------------------------------------------
def PEO(tot_grid,coef_gs,coef_ex,freq,w_eg,reorg_erg,d):
 
 PEO_gs=np.zeros(tot_grid,dtype=complex)
 PEO_ex=np.zeros(tot_grid,dtype=complex)

 x_i,vect=x_Hermit(tot_grid)

 H_g=0.5*freq*x_i*x_i
 H_e=w_eg+reorg_erg+0.5*freq*x_i*x_i+freq*d*x_i

 PEO_gs=coef_gs*H_g
 PEO_ex=coef_ex*H_e

 return PEO_gs,PEO_ex
#-------------------------------------
def EOF(dt,coef_gs,coef_ex,tot_grid,freq,w_eg,reorg_erg,d):

 old_coef_gs=np.copy(coef_gs)
 old_coef_ex=np.copy(coef_ex)

 #RK1:
 KEO_gs,KEO_ex=KEO(tot_grid,coef_gs,coef_ex,freq)
 PEO_gs,PEO_ex=PEO(tot_grid,coef_gs,coef_ex,freq,w_eg,reorg_erg,d)

 RK1_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs)
 RK1_ex=complex(0.0,-1.0)*(KEO_ex+PEO_ex)

 coef_gs=old_coef_gs+dt/2.0*RK1_gs
 coef_ex=old_coef_ex+dt/2.0*RK1_ex

 #RK2:
 KEO_gs,KEO_ex=KEO(tot_grid,coef_gs,coef_ex,freq)
 PEO_gs,PEO_ex=PEO(tot_grid,coef_gs,coef_ex,freq,w_eg,reorg_erg,d)

 RK2_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs)
 RK2_ex=complex(0.0,-1.0)*(KEO_ex+PEO_ex)
 
 coef_gs=old_coef_gs+dt/2.0*RK2_gs
 coef_ex=old_coef_ex+dt/2.0*RK2_ex
  
 #RK3:
 KEO_gs,KEO_ex=KEO(tot_grid,coef_gs,coef_ex,freq)
 PEO_gs,PEO_ex=PEO(tot_grid,coef_gs,coef_ex,freq,w_eg,reorg_erg,d)

 RK3_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs)
 RK3_ex=complex(0.0,-1.0)*(KEO_ex+PEO_ex)
 
 coef_gs=old_coef_gs+dt*RK3_gs
 coef_ex=old_coef_ex+dt*RK3_ex

 #RK4:
 KEO_gs,KEO_ex=KEO(tot_grid,coef_gs,coef_ex,freq)
 PEO_gs,PEO_ex=PEO(tot_grid,coef_gs,coef_ex,freq,w_eg,reorg_erg,d)

 RK4_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs)
 RK4_ex=complex(0.0,-1.0)*(KEO_ex+PEO_ex)

 coef_gs=old_coef_gs+dt/6.0*(RK1_gs+2.0*RK2_gs+2.0*RK3_gs+RK4_gs)
 coef_ex=old_coef_ex+dt/6.0*(RK1_ex+2.0*RK2_ex+2.0*RK3_ex+RK4_ex)

 return coef_gs,coef_ex
#-------------------------------------
if __name__== "__main___":

 # The parameters of the molecular Hamiltonian.
 tot_grid=21
 ev_to_h=27.21138386
 freq=0.185/ev_to_h
 w_eg=2.0/ev_to_h
 reorg_erg=0.5*0.14/ev_to_h
 d=0.87


 # The parameters for the propagation.
 dt=0.01*41.34137333656
 nsteps=5
#---------------------------------------
 # Compute the initial wavefunction.
 erg,coef=init_wf(tot_grid,freq)
 
 coef_gs=np.zeros(tot_grid,dtype=complex)
 coef_ex=np.zeros(tot_grid,dtype=complex)

 coef_ex=coef.astype(complex)
 
#------------------------------------------

    

  
