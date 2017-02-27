########################
#### Main script #######
#### Yu-Han Kao ########
#### kaoyh@umich.edu ###
#### 12/05/2016 ########
########################

#this script runs parameter estimation and profile likelihood of the dengue model. 
#the files below are needed:
## dengue_model.py
## param_est.py
## profile.py
## residual_sigma.py


import scipy.integrate as spi
import scipy.optimize as spopt
import numpy as np
import pylab as pl
import copy
import math
import itertools
import time
import pickle

from dengue_model import dengue_model
from param_est import param_est
from profile import profile
from residual_sigma import sigma

'''
#plot formatting(optional)
from matplotlib import rc
matplotlib.style.use('fivethirtyeight')
matplotlib.rcParams.update({'text.usetex': True, 'font.family': 'sans-serif', 'font.sans-serif': 'Computer Modern Sans Serif', 'font.size': 24, 'figure.figsize': (8,6), 'font.weight': 'bold'})
rc('text.latex',preamble = '\usepackage{sfmath}')
'''

ori_data_case = [1.0, 0.0, 1.0, 1.0, 3.0, 1.0, 0.0, 2.0, 1.0, 0.0, 
0.0, 0.0, 0.0, 1.0, 3.0, 1.0, 1.0, 1.0, 1.0, 0.0, 6.0, 4.0, 25.0, 
20.0, 29.0, 37.0, 71.0, 54.0, 77.0, 52.0, 70.0, 95.0, 70.0, 66.0, 
78.0, 53.0, 67.0, 52.0, 62.0, 33.0, 13.0, 14.0, 5.0, 0.0, 2.0, 1.0, 
1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0] #2010 kaohsung dengue case

sim_start = 20 #fitting start at week20
data_case = ori_data_case[sim_start::] #case data used in simulation

'''========================parameters and initial conditions========================'''
#set up parameters (set to the fitted values)
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
I0 = data_case[0]/repH
E0 = 2*data_case[0]/repH
S0 = 1.0-I0-E0

#mosquito initial state
A0 = 1.0
S0_m = 1.0/mu_m
E0_m = 0.0
I0_m = 0.0


ini = [S0, E0, I0, A0, S0_m, E0_m, I0_m]
param = [beta, repH, x, pi_mua, beta_m, mu_m] 
#Note: alpha, eta, gamma are fixed and given in the model function; C is also fixed as 1 in model function
time_step = np.array(range(0,7*(53-sim_start),7))

##initiate functions
dengue = dengue_model()
dengue_est = param_est(dengue.model, ini, time_step)

'''========================optimizer========================'''

temp_param = param #placeholder
#temp_param=[1.41530277e+01, 1.54673938e+03, 2.02896263e+00, 4.18336111e+00, 3.00024024e-02, 3.16770764e-01]

optimizer = spopt.minimize(dengue_est.residuals, temp_param, args=(data_case, ), method='Nelder-Mead') #BFGS  Nelder-Mead
print "optimization restuls: "
print optimizer



# continue optimization from the previous fitted values if needed
count = 0
while count < 40:
	if optimizer.success==True and optimizer.fun < 93.75: ###105 #3276
		print 'done'
		print optimizer
		break
	else:
		optimizer=spopt.minimize(dengue_est.residuals, optimizer.x, args=(data_case, ), method='Nelder-Mead') #BFGS
		count+=1
		print count
		print 'not yet'
		print 'optimization results: ', optimizer


'''========================opt simulation========================'''

#update parameters and ini status
opttemp = optimizer.x
#opttemp = [1.41530277e+01, 1.54673938e+03, 2.02896263e+00, 4.18336111e+00, 3.00024024e-02, 3.16770764e-01]
opt_beta, opt_repH, opt_x, opt_pi_mua, opt_beta_m, opt_mu_m=opttemp
opt_param = [abs(opt_beta), abs(opt_repH), abs(opt_x), abs(opt_pi_mua), abs(opt_beta_m), abs(opt_mu_m)]


#human initial state
I0 = data_case[0]/abs(opt_repH)
E0 = 2*data_case[0]/abs(opt_repH)
S0 = 1.0-I0-E0
#mosquito initial state
A0 =1.0
S0_m = 1.0/abs(opt_mu_m)
E0_m = 0.0
I0_m = 0.0

opt_ini = [S0, E0, I0, A0, S0_m, E0_m, I0_m]

#intigration restuls
opt_dengue = dengue_model() #create new model obj
opt_res = opt_dengue.ode_run(opt_dengue.model, opt_ini, time_step, opt_param)
#print 'res',opt_res

#summing weekly data
res_inc_t1 = np.append(np.array([0]),spi.cumtrapz(abs(opt_repH)*alpha*opt_res[:,1]))
res_inc_t0 = np.append(np.array([0]),spi.cumtrapz(abs(opt_repH)*alpha*opt_res[0:-1,1]))
data_mw = 7*(res_inc_t1+6.0/7.0) - 7*np.append(np.array([0]), (res_inc_t0+6.0/7.0))
#print data_mw


#calculate sigma
#resisq=(data_mw-data_case)**2
#print 'sigma', ((sum(resisq))/float(len(data_mw)-len(param)))**(0.5)
dengue_sigma = sigma(data_case, data_mw, param)
print 'sigma: ', dengue_sigma


#calculate r0
Sm_DFE = (opt_param[2] - opt_param[3]*opt_param[5])/(opt_param[2]*opt_param[5]) #disease free equilibrium
Mu = 0; #placeholder in case we add it back in later
R0 = (math.sqrt(Sm_DFE*alpha*opt_param[4]*opt_param[0])*gamma)/math.sqrt((alpha+Mu)*(eta+Mu)*opt_param[5]*gamma*(gamma+opt_param[5]))
print R0

'''========================plotting========================'''

#plotting origianal data and fitted curve
pl.figure()
pl.plot(data_mw,'r', label='data_model_weekly')
pl.plot(data_case,'g', label='data_2010_tw')
pl.legend()
pl.show()
#pl.savefig('fitted epidemic')

#plotting simulated epidemic curve
pl.figure()
#pl.plot(opt_res[:,0],'b', label='SH')
#pl.plot(opt_res[:,1],'r', label='EH')
#pl.plot(opt_res[:,2],'g', label='IH')
pl.plot(opt_res[:,3],'m', label='A')
#pl.plot(opt_res[:,4]+opt_res[:,5]+opt_res[:,6], 'k', label='adult')
pl.plot(opt_res[:,4],'y', label='Sm')
#pl.plot(opt_res[:,5],'k', label='Em')
#pl.plot(opt_res[:,6],'c', label='Im')
pl.legend()
pl.show()


'''=====================================simulated data======================================='''

sim_case=data_mw
sim_noise=(np.random.randn(len(sim_case)))*0.2*sim_case
sim_case_ns=sim_case+sim_noise #mock noisy data for human incidence
#pl.figure()
#pl.plot(sim_noise,'g', label='noise')
#pl.plot(sim_case_ns,'r', label='data_model_noise_weekly')
#pl.legend()
#pl.savefig('sim_data')

sim_Am=opt_res[:,3]
sim_Sm=opt_res[:,4]
sim_Em=opt_res[:,5]
sim_Im=opt_res[:,6]

sim_noise_Am=(np.random.randn(len(sim_Am)))*0.2*sim_Am
sim_mo_a=sim_Am+sim_noise_Am #mock noisy data for collecting mosquito in Am
#pl.figure()
#pl.plot(sim_Am, '-k', label='Aquatic')
#pl.plot(sim_noise_Am,'g', label='noise')
#pl.plot(sim_mo_a,'r', label='data_model_noise_weekly')
#pl.legend()
#pl.show()

sim_noise_Sm=(np.random.randn(len(sim_Sm)))*0.2*sim_Sm
sim_mo_s=sim_Sm+sim_noise_Sm #mock noisy data for collecting mosquito in Sm
#pl.figure()
#pl.plot(sim_Am, '-y', label='Susceptible')
#pl.plot(sim_noise_Am,'g', label='noise')
#pl.plot(sim_mo_a,'r', label='data_model_noise_weekly')
#pl.legend()
#pl.show()

sim_noise_SEIm=(np.random.randn(len(sim_Am)))*0.2*(sim_Sm+sim_Em+sim_Im)
sim_mo=sim_Sm+sim_Em+sim_Im+sim_noise_SEIm #mock noisy data for collecting all adult mosquito 
#pl.figure()
#pl.plot(sim_Sm+sim_Em+sim_Im, '-g', label='adult_m')
#pl.plot(sim_noise_SEIm,'g', label='noise')
#pl.plot(sim_mo,'r', label='data_model_noise_weekly')
#pl.legend()
#pl.show()

sim_noise_EIm=(np.random.randn(len(sim_Am)))*0.2*(sim_Em+sim_Im)
sim_mo_i=sim_Em+sim_Im+sim_noise_EIm #mock noisy data for collecting infected mosquitoes
#pl.figure()
#pl.plot(sim_Em+sim_Im, '-b', label='infected_m')
#pl.plot(sim_noise_EIm,'g', label='noise')
#pl.plot(sim_mo_i,'r', label='data_model_noise_weekly')
#pl.legend()
#pl.show()

sim_adultm=sim_Sm+sim_Em+sim_Im #mock data for collecting adult mosquitoes
sim_iem=sim_Im+sim_Em #mock data for collecting infected mosquitoes


'''==============================calculate parmameter likelihood surface =============================='''

#making range list for parameter profiling
dengue_profile=profile(opt_dengue.model, opt_param, opt_ini, time_step, 10, 60) 
param_name=['beta', 'repH', 'x', 'pi_mua', 'beta_m', 'mu_m']


#profile likelihood:
########################

for i in range(len(opt_param)):
	proflog={}
	param_adj_ls, param_fix_temp, param_test_ind = dengue_profile.fit_fix(i,0)
	#iterate through fixed parameter list if there are any
	for fix_ind in param_fix_temp:
		proflog_temp={}
		param_fix=[dengue_profile.param[x] for x in fix_ind]
		fit_ind=list(set(param_test_ind)-set(fix_ind))
		param_fit=[dengue_profile.param[p_f] for p_f in fit_ind]

		for m in param_adj_ls:
			param_adj=[m]
			adj_ind=[i]

			print param_fit
			print param_adj
			
			optimizer_pl=spopt.minimize(dengue_profile.sse_proflike, param_fit, args=(sim_case, param_adj, fit_ind, adj_ind), method='Nelder-Mead')
			#optimizer_pl=spopt.minimize(dengue_profile.sse_vec_proflike, param_fit, args=(sim_case, sim_Am, param_adj, fit_ind, adj_ind), method='Nelder-Mead')
			#optimizer_pl=spopt.minimize(dengue_profile.sse_vec_proflike, param_fit, args=(sim_case, sim_Am, sim_adultm, param_adj, fit_ind, adj_ind), method='Nelder-Mead')
			#optimizer_pl=spopt.minimize(dengue_profile.sse_vec_proflike, param_fit, args=(sim_case, sim_Am, sim_Sm, sim_iem, param_adj, fit_ind, adj_ind), method='Nelder-Mead')

			count=0
			while count <5:
				if optimizer_pl.success==True:
					print optimizer_pl
					break
				else:
					optimizer_pl=spopt.minimize(dengue_profile.sse_proflike, optimizer_pl.x, args=(sim_case, param_adj, fit_ind, adj_ind), method='Nelder-Mead')#BFGS
					#optimizer_pl=spopt.minimize(dengue_profile.sse_vec_proflike, optimizer_pl.x, args=(sim_case, sim_Am, param_adj, fit_ind, adj_ind), method='Nelder-Mead')
					#optimizer_pl=spopt.minimize(dengue_profile.sse_vec_proflike, optimizer_pl.x, args=(sim_case, sim_Am, sim_adultm, param_adj, fit_ind, adj_ind), method='Nelder-Mead')
					#optimizer_pl=spopt.minimize(dengue_profile.sse_vec_proflike, optimizer_pl.x, args=(sim_case, sim_Am, sim_Sm, sim_iem, param_adj, fit_ind, adj_ind), method='Nelder-Mead')
					count+=1
					print count, 'not done yet...'
					print 'optimization results: ', optimizer_pl
					
			p_val_error=np.append(optimizer_pl.x,optimizer_pl.fun)
			temptemp=np.absolute(p_val_error)
			temptemp_abs=temptemp.tolist()
			temptemp_abs.append(optimizer_pl.success)
			proflog_temp[m]=temptemp_abs	
			
		fix_name=[param_name[z] for z in fix_ind]
		s="_"
		f_ns=str(fit_ind)[1:-1]
		fix_name_lb=str(s.join(fix_name)+'_['+f_ns+']')
		proflog[fix_name_lb]=proflog_temp
		#print proflog_temp

	pl_data=open('pf_data_human_01162017_'+param_name[i]+'.txt','w')
	pl_data.write (str(proflog))
	pl_data.close()
	pickle.dump(proflog,open('pf_data_human_01162017_'+param_name[i]+'.pickle','wb')) #save as pickle for plotting

#===end===#
