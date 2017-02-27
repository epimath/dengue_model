########################
#### Run FIM ###########
#### Yu-Han Kao ########
#### kaoyh@umich.edu ###
#### 12/05/2016 ########
########################

#this script runs fim.py to calculate fisher information matrix

from fim import fim
from dengue_model import dengue_model
import scipy.integrate as spi
import scipy.optimize as spopt
import numpy as np
import math

#model setup:
#for initial state
ori_data_case=[1.0, 0.0, 1.0, 1.0, 3.0, 1.0, 0.0, 2.0, 1.0, 0.0, 
0.0, 0.0, 0.0, 1.0, 3.0, 1.0, 1.0, 1.0, 1.0, 0.0, 6.0, 4.0, 25.0, 
20.0, 29.0, 37.0, 71.0, 54.0, 77.0, 52.0, 70.0, 95.0, 70.0, 66.0, 
78.0, 53.0, 67.0, 52.0, 62.0, 33.0, 13.0, 14.0, 5.0, 0.0, 2.0, 1.0, 
1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0] #2010 kaohsung dengue case
sim_start=20 #started at week20
data_case=ori_data_case[sim_start::] #case data used in simulation

opttemp = [14.1530277, 1546.73938, 2.02896263, 4.18336111, 0.0300024024, 0.316770764] #from the fitting in dengue_run.py

opt_beta, opt_repH, opt_x, opt_pi_mua, opt_beta_m, opt_mu_m = opttemp
param = [abs(opt_beta), abs(opt_repH), abs(opt_x), abs(opt_pi_mua), abs(opt_beta_m), abs(opt_mu_m)]
alpha = 0.14

#human initial state
I0=data_case[0]/abs(opt_repH)
E0=2*data_case[0]/abs(opt_repH)
S0=1.0-I0-E0
#mosquito initial state
A0=1.0
S0_m=1.0/abs(opt_mu_m)
E0_m=0.0
I0_m=0.0
ini=[S0, E0, I0, A0, S0_m, E0_m, I0_m]
time_step=np.array(range(0,7*(53-sim_start),7))

#simulated case data
dengue=dengue_model()
opt_res=dengue.ode_run(dengue.model, ini, time_step, param)
res_inc_t1=np.append(np.array([0]),spi.cumtrapz(abs(opt_repH)*alpha*opt_res[:,1]))
res_inc_t0=np.append(np.array([0]),spi.cumtrapz(abs(opt_repH)*alpha*opt_res[0:-1,1]))
data_mw=7*(res_inc_t1+6.0/7.0) - 7*np.append(np.array([0]), (res_inc_t0+6.0/7.0))
#print 'data_mw', data_mw

#fisher information matrix 
fim_dengue=fim(dengue.model, time_step)
fim=fim_dengue.Fim_x(param, data_mw) #default pert = 0.01
fim_r=fim_dengue.Fim_rank(fim) #set the default tolerance to match matlab
print fim
print 'rank: ',fim_r





