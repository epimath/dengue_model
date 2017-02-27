########################
#### Surface of ########
#### parameter pairs ###
#### Yu-Han Kao ########
#### kaoyh@umich.edu ###
#### 12/05/2016 ########
########################

#this script produces the likelihood surface for each parameter pair
#the files below are needed:
## dengue_model.py
## param_est.py
## profile.py



import scipy.integrate as spi
import scipy.optimize as spopt
import numpy as np
import pylab as pl
from scipy.interpolate import interp1d
import copy
import math
import pickle

from dengue_model import dengue_model
from param_est import param_est
from profile.py import profile




'''======================== data input ========================'''

ori_data_case=[1.0, 0.0, 1.0, 1.0, 3.0, 1.0, 0.0, 2.0, 1.0, 0.0, 
0.0, 0.0, 0.0, 1.0, 3.0, 1.0, 1.0, 1.0, 1.0, 0.0, 6.0, 4.0, 25.0, 
20.0, 29.0, 37.0, 71.0, 54.0, 77.0, 52.0, 70.0, 95.0, 70.0, 66.0, 
78.0, 53.0, 67.0, 52.0, 62.0, 33.0, 13.0, 14.0, 5.0, 0.0, 2.0, 1.0, 
1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0] #2010 kaohsung

sim_start=20
data_case=ori_data_case[sim_start::]

'''========================parameters========================'''

beta = 0.5                #transmission rate of humans* 
repH = 2*4.19559287e+02   #reporting rate of human cases
alpha = 0.14              #intrinsic incubation rate (humans)-fixed in model
eta = 0.2                 #human recovery rate-fixed in model
x = 1.054                 #oviposition-fertilization rate
pi_mua = 0.0528           #aquatic death rate
beta_m = 1.0              #transmission rate of mosquito
mu_m = 1.066              #mosquito death rate
gamma = 0.1               #extrinsic incubation rate (mosquito)-fixed in model

#human initial state
I0=data_case[0]/repH
E0=2*data_case[0]/repH#0.0#0.01
S0=1.0-I0-E0

#mosquito initial state
A0=1.0
S0_m=1.0/mu_m
E0_m=0.0
I0_m=0.0


ini=[S0, E0, I0, A0, S0_m, E0_m, I0_m]
param=[beta, repH, x, pi_mua, beta_m, mu_m] #alpha, eta, gamma are fixed and given in the model function; C is also fixed as 1 in model function
time_step=np.array(range(0,7*(53-sim_start),7))

##initiate functions
dengue=dengue_model()
dengue_est=param_est(dengue.model, ini, time_step)

'''========================opt simulation========================'''
# update parameters and ini status
opttemp=[1.41530277e+01, 1.54673938e+03, 2.02896263e+00, 4.18336111e+00, 3.00024024e-02, 3.16770764e-01]

opt_beta, opt_repH, opt_x, opt_pi_mua, opt_beta_m, opt_mu_m=opttemp
opt_param=[abs(opt_beta), abs(opt_repH), abs(opt_x), abs(opt_pi_mua), abs(opt_beta_m), abs(opt_mu_m)]


#human initial state
I0=data_case[0]/abs(opt_repH)
E0=2*data_case[0]/abs(opt_repH)
S0=1.0-I0-E0
#mosquito initial state
A0=1.0
S0_m=1.0/abs(opt_mu_m)
E0_m=0.0
I0_m=0.0

opt_ini=[S0, E0, I0, A0, S0_m, E0_m, I0_m]

#intigration restuls
opt_temp=dengue_model() #new model obj
opt_res=opt_temp.ode_run(opt_temp.model, opt_ini, time_step, opt_param)


'''========================plotting========================'''

#summing weekly data
res_inc_t1=np.append(np.array([0]),spi.cumtrapz(abs(opt_repH)*alpha*opt_res[:,1]))
res_inc_t0=np.append(np.array([0]),spi.cumtrapz(abs(opt_repH)*alpha*opt_res[0:-1,1]))
data_mw=7*(res_inc_t1+6.0/7.0) - 7*np.append(np.array([0]), (res_inc_t0+6.0/7.0))

'''
#plotting origianal data and fitted curve
pl.figure()
pl.plot(data_mw,'r', label='data_model_weekly')
pl.plot(data_case,'g', label='data_2010_tw')
pl.legend()


#plotting simulated epidemic curve
pl.figure()
#pl.plot(opt_res[:,0],'b', label='SH')
#pl.plot(opt_res[:,1],'r', label='EH')
#pl.plot(opt_res[:,2],'g', label='IH')
pl.plot(opt_res[:,3],'m', label='A')
pl.plot(opt_res[:,4],'y', label='Sm')
#pl.plot(opt_res[:,5],'k', label='Em')
#pl.plot(opt_res[:,6],'c', label='Im')
pl.legend()
pl.show()
#pl.savefig('temp_fig')
'''

'''==============================calculate paried param likelihood surface =============================='''

dengue_profile=profile(opt_temp.model, opt_param, opt_ini, time_step, 10, 60)
param_name=['beta', 'repH', 'x', 'pi_mua', 'beta_m', 'mu_m']
templs=dengue_profile.pert_level

temp_ee=[(a,b) for a in range(len(opt_param)) for b in range(len(opt_param))]
temp_ee2=[]

for j,k in temp_ee: #drawing pairs
	temp_ele1=(j,k)
	temp_ele2=(k,j)
	if temp_ele2 not in temp_ee2 and j!=k:
		temp_ee2.append(temp_ele1)
#print temp_ee2

for q,r in temp_ee2:

	var1=opt_param[q]
	var1_ls=[]
	var1_ls=dengue_profile.adj_param_list(var1) #making the list of chaning values for var1

	var2=opt_param[r]
	var2_ls=[]
	var2_ls=dengue_profile.adj_param_list(var2) #making the list of chaning values for var2

	var_rest_temp=copy.copy(opt_param)
	indices=(q,r)
	var_rest=[z for s,z in enumerate(var_rest_temp) if s not in indices]

	error_ls=[]
	var1_run=[]
	var2_run=[]
	for n in range(len(var1_ls)):
		var1_temp=var1_ls[n]
		for k in range(len(var2_ls)):
			var2_temp=var2_ls[k]
			var1_run.append(var1_temp)
			var2_run.append(var2_temp)
			newparam=copy.copy(var_rest)
			newparam.insert(q, var1_temp)
			newparam.insert(r, var2_temp) #create new param list
			est=param_est(dengue.model, ini, time_step)
			err=est.residues(newparam, data_mw) #calculate RSS
			error_ls.append(err)

	print len(var1_run),len(var2_run),len(error_ls)
	surf_data=open('pf_surface_'+param_name[q]+'_'+param_name[r]+'_human_only_2017.txt','w')
	surf_data.write(str(var1_ls)+'\n'+str(var2_ls)+'\n'+str(var1_run)+'\n'+str(var2_run)+'\n'+str(error_ls))
	surf_data.close()

#===end===#

