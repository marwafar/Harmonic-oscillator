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
def KEO(tot_grid,coef_gs,coef_ex,freq,n_fock):

 dif2mat=second_derivavtive_Hermit(tot_grid)
 
 KEO_gs=np.zeros((tot_grid,n_fock),dtype=complex)
 KEO_ex=np.zeros((tot_grid,n_fock),dtype=complex)

 for m in range(n_fock):
  KEO_gs[:,m]=-0.5*freq*(dif2mat.dot(coef_gs[:,m]))
  KEO_ex[:,m]=-0.5*freq*(dif2mat.dot(coef_ex[:,m]))

 return KEO_gs,KEO_ex
#------------------------------------------
def PEO(tot_grid,coef_gs,coef_ex,freq,w_eg,reorg_erg,d,n_fock,freq_cavity):
 
 PEO_gs=np.zeros((tot_grid,n_fock),dtype=complex)
 PEO_ex=np.zeros((tot_grid,n_fock),dtype=complex)

 x_i,vect=x_Hermit(tot_grid)

 H_g=0.5*freq*x_i*x_i
 H_e=w_eg+reorg_erg+0.5*freq*x_i*x_i+freq*d*x_i

 for m in range(n_fock):
  PEO_gs[:,m]=coef_gs[:,m]*(H_g+(m+0.5)*freq_cavity)
  PEO_ex[:,m]=coef_ex[:,m]*(H_e+(m+0.5)*freq_cavity)

 return PEO_gs,PEO_ex
#------------------------------------------------------------------
def interaction(coef_gs,coef_ex,tot_grid,n_fock,freq_cavity,g):

 Hcm_01=np.zeros((tot_grid,n_fock),dtype=complex)
 Hcm_10=np.zeros((tot_grid,n_fock),dtype=complex)

 for n in range(n_fock):
  for m in range(n_fock):
   Hcm_01[:,n]+=g*coef_ex[:,m]*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n))
   Hcm_10[:,n]+=g*coef_gs[:,m]*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n))

 return Hcm_01,Hcm_10
#-----------------------------------------------------------
def EOF(dt,coef_gs,coef_ex,tot_grid,freq,w_eg,reorg_erg,d,g,n_fock,freq_cavity):

 old_coef_gs=np.copy(coef_gs)
 old_coef_ex=np.copy(coef_ex)

 #RK1:
 KEO_gs,KEO_ex=KEO(tot_grid,coef_gs,coef_ex,freq,n_fock)
 PEO_gs,PEO_ex=PEO(tot_grid,coef_gs,coef_ex,freq,w_eg,reorg_erg,d,n_fock,freq_cavity)
 Hcm_01,Hcm_10=interaction(coef_gs,coef_ex,tot_grid,n_fock,freq_cavity,g)

 RK1_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs+Hcm_01)
 RK1_ex=complex(0.0,-1.0)*(KEO_ex+PEO_ex+Hcm_10)

 coef_gs=old_coef_gs+dt/2.0*RK1_gs
 coef_ex=old_coef_ex+dt/2.0*RK1_ex

 #RK2:
 KEO_gs,KEO_ex=KEO(tot_grid,coef_gs,coef_ex,freq,n_fock)
 PEO_gs,PEO_ex=PEO(tot_grid,coef_gs,coef_ex,freq,w_eg,reorg_erg,d,n_fock,freq_cavity)
 Hcm_01,Hcm_10=interaction(coef_gs,coef_ex,tot_grid,n_fock,freq_cavity,g)

 RK2_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs+Hcm_01)
 RK2_ex=complex(0.0,-1.0)*(KEO_ex+PEO_ex+Hcm_10)
 
 coef_gs=old_coef_gs+dt/2.0*RK2_gs
 coef_ex=old_coef_ex+dt/2.0*RK2_ex
  
 #RK3:
 KEO_gs,KEO_ex=KEO(tot_grid,coef_gs,coef_ex,freq,n_fock)
 PEO_gs,PEO_ex=PEO(tot_grid,coef_gs,coef_ex,freq,w_eg,reorg_erg,d,n_fock,freq_cavity)
 Hcm_01,Hcm_10=interaction(coef_gs,coef_ex,tot_grid,n_fock,freq_cavity,g)

 RK3_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs+Hcm_01)
 RK3_ex=complex(0.0,-1.0)*(KEO_ex+PEO_ex+Hcm_10)
 
 coef_gs=old_coef_gs+dt*RK3_gs
 coef_ex=old_coef_ex+dt*RK3_ex

 #RK4:
 KEO_gs,KEO_ex=KEO(tot_grid,coef_gs,coef_ex,freq,n_fock)
 PEO_gs,PEO_ex=PEO(tot_grid,coef_gs,coef_ex,freq,w_eg,reorg_erg,d,n_fock,freq_cavity)
 Hcm_01,Hcm_10=interaction(coef_gs,coef_ex,tot_grid,n_fock,freq_cavity,g)

 RK4_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs+Hcm_01)
 RK4_ex=complex(0.0,-1.0)*(KEO_ex+PEO_ex+Hcm_10)

 coef_gs=old_coef_gs+dt/6.0*(RK1_gs+2.0*RK2_gs+2.0*RK3_gs+RK4_gs)
 coef_ex=old_coef_ex+dt/6.0*(RK1_ex+2.0*RK2_ex+2.0*RK3_ex+RK4_ex)

 return coef_gs,coef_ex
#-------------------------------------
def PES(freq,w_eg,reorg_erg,d,n_fock,freq_cavity,g):

 x_i,vect=x_Hermit(tot_grid)

 E_ad=np.zeros((2,2))
 H_ij=np.zeros((2*n_fock,2*n_fock))
 PES=open("polariton-PES.txt", "w+")

 H_g=0.5*freq*x_i*x_i
 H_e=w_eg+reorg_erg+0.5*freq*x_i*x_i+freq*d*x_i
  
 for x in range(tot_grid):
  E_ad[0,0]=H_g[x]
  E_ad[1,1]=H_e[x]

  for row in range(2*n_fock):
   a=int(row/n_fock)
   m=row%n_fock
   for colm in range(2*n_fock):
    b=int(colm/n_fock)
    n=colm%n_fock

    H_ij[row,colm]=E_ad[a,b]*float(m==n) + (n+0.5)*freq_cavity*float(a==b)*float(m==n)
    H_ij[row,colm]+=g*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n))*\
      float((a==1 and b==0) or (a==0 and b==1))

  E_pol,U=diag(H_ij)
  PES.write(str(x_i[x]) + " "+ " ".join(E_pol.astype(str))+"\n")
 
 PES.close
 return
#-------------------------------------
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
 nsteps=20000
 output=50
#-------------------------------------
# Compute the PES
 PES(freq,w_eg,reorg_erg,d,n_fock,freq_cavity,g)

#---------------------------------------
 # Compute the initial wavefunction.
 erg,coef=init_wf(tot_grid,freq)
 
 coef_gs=np.zeros((tot_grid,n_fock),dtype=complex)
 coef_ex=np.zeros((tot_grid,n_fock),dtype=complex)

 coef_ex[:,0]=coef[:,0].astype(complex)
 
#------------------------------------------
 # Run the propagation.
 
 pop_gs=np.zeros((n_fock),dtype=complex)
 pop_ex=np.zeros((n_fock),dtype=complex)

 # file name
 pop_ad_gs=open("ad_pop_gs.txt","w+") 
 pop_ad_ex=open("ad_pop_ex.txt","w+")
  
  
 for step in range(nsteps):
  time=dt*step/41.34137333656
  out=int(step/output)*output
 
  if (step==out):
#   print(out)
   for n in range(n_fock):
    pop_gs[n]=np.dot(coef_gs[:,n],np.conj(coef_gs[:,n]))
    pop_ex[n]=np.dot(coef_ex[:,n],np.conj(coef_ex[:,n]))

   pop_ad_gs.write(str(time)+ " "+ " ".join(np.real(pop_gs).astype(str))+"\n")
   pop_ad_ex.write(str(time)+ " "+ " ".join(np.real(pop_ex).astype(str))+"\n")

  coef_gs,coef_ex=EOF(dt,coef_gs,coef_ex,tot_grid,freq,w_eg,reorg_erg,d,g,n_fock,freq_cavity)

 pop_ad_gs.close()
 pop_ad_ex.close()

