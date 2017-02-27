########################
#### Main script #######
#### Yu-Han Kao ########
#### kaoyh@umich.edu ###
#### 12/05/2016 ########
########################

#this script plots the surface of paired parameter using the .txt files created from pair_furface.py 
#The files stored RSS values of paried parameter are needed

import numpy as np
import matplotlib.pyplot as plt
import math
param_name=['beta', 'repH', 'x', 'pi_mua', 'beta_m', 'mu_m']

opttemp=[1.41530277e+01, 1.54673938e+03, 2.02896263e+00, 4.18336111e+00, 3.00024024e-02, 3.16770764e-01] #fitted param values
opt_beta, opt_repH, opt_x, opt_pi_mua, opt_beta_m, opt_mu_m=opttemp
opt_param=[abs(opt_beta), abs(opt_repH), abs(opt_x), abs(opt_pi_mua), abs(opt_beta_m), abs(opt_mu_m)]


temp_ee=[(a,b) for a in range(len(param_name)) for b in range(len(param_name))]
temp_ee2=[]

for j,k in temp_ee: #drawing pairs
	temp_ele1=(j,k)
	temp_ele2=(k,j)
	if temp_ele2 not in temp_ee2 and j!=k:
		temp_ee2.append(temp_ele1)
print temp_ee2

for q,r in temp_ee2: #iterate through all the .txt files from pair_surface.py

	temp=open('pf_surface_'+param_name[q]+'_'+param_name[r]+'_human_only_2017.txt','r')

	temp2=temp.read()


	varname=['var1','var2','var11','var22','error']

	temp3=temp2.split('\n')
	print 'pf_surface_'+param_name[q]+'_'+param_name[r]+'_2017.txt'

	#read the files
	readin={} 
	for i in range(len(temp3)):
		temp_len=len(temp3[i])
		temp4=temp3[i][1:temp_len-2]

		temp6=[]
		temp5=temp4.split(',')

		for j in temp5:
			temp6.append(float(j))
		readin[varname[i]]=temp6


	var1=np.array(readin['var1'])
	var2=np.array(readin['var2'])
	#error=np.array(readin['error'])

	error=readin['error']
	logerror=[]
	for i in error:
		j=math.log((i+1.0),1.2)
		logerror.append(j)
	logerror=np.array(logerror)

	logerror=np.reshape(logerror,[len(var1),len(var2)])
	
	#error=np.reshape(error,[len(var1),len(var2)])
	#print error[:,0]

	var2_me, var1_me=np.meshgrid(var2, var1)

	#print max(var1)
	#print max(var2)
	#print var1_me

	#making the heatmap
	fig, ax=plt.subplots()
	im=ax.pcolormesh(var2, var1, logerror)
	#im=ax.pcolormesh(var2, var1, error)
	ax.plot(opt_param[r],opt_param[q],'*',color='black') ####
	color=fig.colorbar(im)
	color.set_label('temp')
	ax.axis('tight')
	ax.set_xlim(xmin=0.0)
	ax.set_ylim(ymin=0.0)
	plt.title('figure_surface_'+param_name[q]+'_'+param_name[r])
	plt.xlabel(param_name[r])
	plt.ylabel(param_name[q])
	#plt.savefig('surface_log_'+param_name[q]+'_'+param_name[r]+'_vec_0530_2')
	plt.savefig('surface_'+param_name[q]+'_'+param_name[r]+'_human_2017')

#===end===#


