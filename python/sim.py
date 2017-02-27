########################
#### Intervention sim ##
#### Yu-Han Kao ########
#### kaoyh@umich.edu ###
#### 12/05/2016 ########
########################

#this script simulates the misleading results from identifiability while implementing interventions 
#the file below are needed:
## dengue_model_intervention.py
## dengue_model.py

import scipy.integrate as spi
import scipy.optimize as spopt
import numpy as np
import pylab as pl
import copy
import math

from dengue_model_intervention import dengue_model_itv
from dengue_model import dengue_model


ori_data_case=[1.0, 0.0, 1.0, 1.0, 3.0, 1.0, 0.0, 2.0, 1.0, 0.0, 
0.0, 0.0, 0.0, 1.0, 3.0, 1.0, 1.0, 1.0, 1.0, 0.0, 6.0, 4.0, 25.0, 
20.0, 29.0, 37.0, 71.0, 54.0, 77.0, 52.0, 70.0, 95.0, 70.0, 66.0, 
78.0, 53.0, 67.0, 52.0, 62.0, 33.0, 13.0, 14.0, 5.0, 0.0, 2.0, 1.0, 
1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0] #2010 kaohsung dengue case

sim_start=20 #started at week20
data_case=ori_data_case[sim_start::] #case data used in simulation

'''========================parameters and initial statust========================'''
beta=0.5                  #transmission rate of humans* 
repH=2*4.19559287e+02     #reporting rate of human cases
alpha=0.14                #intrinsic incubation rate (humans)-fixed in model
eta=0.2                   #human recovery rate-fixed in model
x=1.054                   #oviposition-fertilization rate
pi_mua=0.0528             #aquatic death rate
beta_m=1.0                #transmission rate of mosquito
mu_m=1.066                #mosquito death rate
gamma=0.1                 #extrinsic incubation rate (mosquito)-fixed in model


time_step=np.array(range(0,7*(53-sim_start),7))

opttemp1 = [1.41530277e+01, 1.54673938e+03, 2.02896263e+00, 4.18336111e+00, 3.00024024e-02, 3.16770764e-01]
opttemp2 = [38.102752439869036, 1625.421800668292, 0.12733483979297147, 0.14794345511043694, 0.02318660893926542,0.46390293357605983]

opt_beta1, opt_repH1, opt_x1, opt_pi_mua1, opt_beta_m1, opt_mu_m1=opttemp1
opt_param1=[abs(opt_beta1), abs(opt_repH1), abs(opt_x1), abs(opt_pi_mua1), abs(opt_beta_m1), abs(opt_mu_m1)]

opt_beta2, opt_repH2, opt_x2, opt_pi_mua2, opt_beta_m2, opt_mu_m2=opttemp2
opt_param2=[abs(opt_beta2), abs(opt_repH2), abs(opt_x2), abs(opt_pi_mua2), abs(opt_beta_m2), abs(opt_mu_m2)]

#human initial state for case 1
I01=data_case[0]/abs(opt_repH1)
E01=2*data_case[0]/abs(opt_repH1)
S01=1.0-I01-E01
#mosquito initial state for case 1
A01=1.0
S0_m1=1.0/abs(opt_mu_m1)
E0_m1=0.0
I0_m1=0.0
opt_ini1=[S01, E01, I01, A01, S0_m1, E0_m1, I0_m1]

#human initial state for case 2
I02=data_case[0]/abs(opt_repH2)
E02=2*data_case[0]/abs(opt_repH2)
S02=1.0-I02-E02
#mosquito initial state for case 2
A02=1.0
S0_m2=1.0/abs(opt_mu_m2)
E0_m2=0.0
I0_m2=0.0
opt_ini2=[S02, E02, I02, A02, S0_m2, E0_m2, I0_m2]

#intigration restuls
opt_temp1=dengue_model_itv() #new model obj
opt_res1=opt_temp1.ode_run(opt_temp1.model, opt_ini1, time_step, opt_param1)

opt_temp2=dengue_model_itv() #new model obj
opt_res2=opt_temp2.ode_run(opt_temp2.model, opt_ini2, time_step, opt_param2)

opt_temp_ori=dengue_model() #new model obj
opt_res_ori=opt_temp_ori.ode_run(opt_temp_ori.model, opt_ini1, time_step, opt_param1)


res_inc1_t1=np.append(np.array([0]),spi.cumtrapz(abs(opt_repH1)*alpha*opt_res1[:,1]))
res_inc1_t0=np.append(np.array([0]),spi.cumtrapz(abs(opt_repH1)*alpha*opt_res1[0:-1,1]))
data_mw1=7*(res_inc1_t1+6.0/7.0) - 7*np.append(np.array([0]), (res_inc1_t0+6.0/7.0))
print data_mw1


res_inc2_t1=np.append(np.array([0]),spi.cumtrapz(abs(opt_repH2)*alpha*opt_res2[:,1]))
res_inc2_t0=np.append(np.array([0]),spi.cumtrapz(abs(opt_repH2)*alpha*opt_res2[0:-1,1]))
data_mw2=7*(res_inc2_t1+6.0/7.0) - 7*np.append(np.array([0]), (res_inc2_t0+6.0/7.0))
print data_mw2

res_inc_ori_t1=np.append(np.array([0]),spi.cumtrapz(abs(opt_repH1)*alpha*opt_res_ori[:,1]))
res_inc_ori_t0=np.append(np.array([0]),spi.cumtrapz(abs(opt_repH1)*alpha*opt_res_ori[0:-1,1]))
data_mw_ori=7*(res_inc_ori_t1+6.0/7.0) - 7*np.append(np.array([0]), (res_inc_ori_t0+6.0/7.0))
print data_mw_ori

print opt_res1
print opt_res2
print opt_res_ori


##R0
Sm_DFE = (opt_param1[2] - opt_param1[3]*opt_param1[5])/(opt_param1[2]*opt_param1[5])
print Sm_DFE
Mu = 0;
R0 = (math.sqrt(Sm_DFE*alpha*opt_param1[4]*opt_param1[0])*gamma)/math.sqrt((alpha+Mu)*(eta+Mu)*opt_param1[5]*gamma*(gamma+opt_param1[5]))
print R0

Sm_DFE = (opt_param2[2] - opt_param2[3]*opt_param2[5])/(opt_param2[2]*opt_param2[5])
print Sm_DFE
Mu = 0;
R0 = (math.sqrt(Sm_DFE*alpha*opt_param2[4]*opt_param2[0])*gamma)/math.sqrt((alpha+Mu)*(eta+Mu)*opt_param2[5]*gamma*(gamma+opt_param2[5]))
print R0


#plotting origianal data and fitted curve
pl.figure()
pl.plot
pl.plot(data_mw1,'r',linewidth=3.5, label= 'exp1_intervention')
pl.plot(data_mw2,'k', linewidth=3.5, label= 'exp2_intervention')
#pl.plot(data_mw1,'r',linewidth=3.5)# label='data_model_weekly'
#pl.plot(data_mw2,'g',linewidth=3.5)
#pl.plot(data_case,'bo', data_case, 'k')#label='data_2010_tw'
pl.plot(data_mw_ori,'b', linewidth=3.5, label='fit_2010_tw')
pl.xticks(fontsize = 18)
pl.yticks(fontsize = 18)
pl.legend()
pl.show()
#pl.savefig('fitting_eta_2_0422')



#plotting simulated epidemic curve
pl.figure()
#pl.plot(opt_res[:,0],'b', label='SH')
#pl.plot(opt_res[:,1],'r', label='EH')
#pl.plot(opt_res[:,2],'g', label='IH')

pl.plot(opt_res1[:,3],'m', label='A_exp1')
#pl.plot(opt_res1[:,4]+opt_res1[:,5]+opt_res1[:,6], 'k', label='adult')
#pl.plot(opt_res1[:,4],'y', label='Sm')
#pl.plot(opt_res[:,5],'k', label='Em')
#pl.plot(opt_res[:,6],'c', label='Im')

pl.plot(opt_res2[:,3],'r', label='A_exp2')
#pl.plot(opt_res2[:,4]+opt_res2[:,5]+opt_res2[:,6], 'b', label='adult')
#pl.plot(opt_res2[:,4],'g', label='Sm')

pl.plot(opt_res_ori[:,3],'k', label='A_ori')

pl.legend()
pl.show()
#pl.savefig('temp_fig')

