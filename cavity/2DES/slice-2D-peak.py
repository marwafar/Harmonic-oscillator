import numpy as np

#---------------------------
n_file=5
total_sum=np.zeros(5)
i=0
slice_data=open("UP_diagonal.txt","w+")

for k in range(0,25,5):
 data=np.loadtxt("2DES_T%02d.txt"%k)
 n_row,n_colm=data.shape

# for polaritons(lower and upper: w1=2.27244-2.29999; w3=1.7628-1.9281)
 for row in range(n_row):
   if (data[row,0]>=2.27) and (data[row,0]<=2.3):
    if (data[row,1]>=2.27) and(data[row,1]<=2.3):
     total_sum[i]+=data[row,3]

 slice_data.write(str(k)+ " "+ str(total_sum[i])+"\n")
 i+=1

slice_data.close()
