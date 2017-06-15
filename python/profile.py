########################
## Profile likelihood ##
#### Yu-Han Kao ########
#### kaoyh@umich.edu ###
#### 11/06/2016 ########
########################

#this script includes functions to calculate profile likelihood for parameters


import scipy.integrate as spi
import scipy.optimize as spopt
import numpy as np
import math
import pylab as pl
import copy
import itertools

class profile:

	def __init__(self, model, param, ini, time_step, pert_mag, pert_steps):
		self.param=param
		self.ini=ini
		self.time_step=time_step
		temptempls=np.linspace(1, pert_mag, pert_steps)
		temptempls.tolist
		templs=[]
		for k in temptempls:
			templs.append(1.0/float(k))
			templs.append(k)
		templs.pop(0)
		templs.sort()
		self.pert_level=templs
		self.model=model
		self.description = 'fit_fix: lists of values for the adjusted parameter and index for fixed parameters and fitted parameters;\n \
		sse_proflike: calculate profile likelihood with only human incidence data (scenario 1); \n \
		sse_vec_proflike: calculate profile likelihood with both human incidence data and mosquito population data (scenario 2-4)'


	#############################################################
	########## List for a range of values  ######################
	########## of the prameter we want to profile ###############
	#############################################################
	def adj_param_list (self, adj_var):
		param_adj_ls=[]
		for j in range(len(self.pert_level)): #create a list of values for the adjusted parameter
			param_adj_ls.append(adj_var*self.pert_level[j])
			if self.pert_level[j]<=1.0:
				param_adj_ls.append(adj_var+self.pert_level[j])
				param_adj_ls.append(adj_var-self.pert_level[j])
		param_adj_ls.sort()
		return param_adj_ls


	#############################################################
	####### lists for values of the parameter we want to ########
	####### profile and index of fixed and fitted parameters ####
	#############################################################
	
	def fit_fix(self, adj_ind=None, fixed_n=0):
		param_test=copy.copy(self.param)
		param_test_ind=list(set(range(len(self.param)))-set([adj_ind]))#the index of the parameters we want to fit while profiling
		param_adj_temp=param_test.pop(adj_ind) 
		param_test=np.array(param_test)
		param_adj_ls=[] #tthe list of values for the parameter we want to profile (adjusted parameter) across the range of perturbation
		param_adj_ls_up=[param_adj_temp]
		param_adj_ls_low=[param_adj_temp]
		for j in range(len(self.pert_level)): #create a list of perturbation values for the adjusted parameter
			param_adj_ls.append(param_adj_temp*self.pert_level[j])
			if self.pert_level[j]<=1.0:
				param_adj_ls.append(abs(param_adj_temp+self.pert_level[j]))
				param_adj_ls.append(abs(param_adj_temp-self.pert_level[j]))
		param_adj_ls.sort()
		for n in param_adj_ls:
			if n > param_adj_temp:
				param_adj_ls_up.append(n)
			if n < param_adj_temp:
				param_adj_ls_low.append(n)
		param_adj_ls_up.sort()
		param_adj_ls_low.sort(reverse=True)
		param_adj_ls_all=[param_adj_ls_up, param_adj_ls_low]
		param_fix_temp=itertools.combinations(param_test_ind, fixed_n) #choose the/how many parameters you want to fix
		fit_fix_list=[param_adj_ls_all, param_fix_temp, param_test_ind]
		return fit_fix_list

	#############################################################
	########## Profile likelihood function ######################
	########## with human data_generalized ######################
	#############################################################

	def sse_proflike(self, param_fit, case_data, param_adj=None, fit_ind=None, adj_ind=None):	
		#specify parameter sets: parameters for fitting and paramter for adjusting 
		param_fit_ls=list(param_fit)
		param_opt_ls=list(self.param)
		param_adj_ls=list(param_adj) #adjusted parameter list
		alpha = 0.14
		#print param_opt_ls

		#substitute fited parameter set
		try:
			for i in range(len(fit_ind)):
				param_opt_ls.pop(fit_ind[i])
				param_opt_ls.insert(fit_ind[i], param_fit_ls[i])
		except:
			pass
		print param_opt_ls

		#substitute adjusted parameter set
		try:
			for j in range(len(adj_ind)):
				param_opt_ls.pop(adj_ind[j])
				param_opt_ls.insert(adj_ind[j], param_adj_ls[j])
		except:
			pass
		print param_opt_ls

		param_opt_arr=np.array(param_opt_ls)
		print param_opt_arr
		beta, repH, x, pi_mua, beta_m, mu_m=param_opt_arr
		absparam=(abs(beta), abs(repH), abs(x), abs(pi_mua), abs(beta_m), abs(mu_m))
		S0, E0, I0, A0, S0_m, E0_m, I0_m=self.ini
		I0=case_data[0]/abs(repH)
		E0=2*case_data[0]/abs(repH)
		S0=1-I0-E0
		S0_m=1.0/abs(mu_m)
		temp_ini=(S0, E0, I0, A0, S0_m, E0_m, I0_m)
		

		res = spi.odeint(self.model, temp_ini, self.time_step, args=(absparam,), mxstep=2000)
		res_inc_t1=np.append(np.array([0]),spi.cumtrapz(abs(repH)*alpha*res[:,1]))
		res_inc_t0=np.append(np.array([0]),spi.cumtrapz(abs(repH)*alpha*res[0:-1,1]))
		data_mw=7*(res_inc_t1+6.0/7.0) - 7*np.append(np.array([0]), (res_inc_t0+6.0/7.0))
		#err=np.sum((case_data-data_mw)**2) #normal SSE
		err=np.sum((case_data-data_mw)**2/np.mean(case_data)) #weighted SSE changed on 02/06/2016
		#err=np.sum((case_data-data_mw)**2/(np.mean(data_mw))**0.5) #weighted SSE with squared root demoninator 12/17/2016 try
		print param_fit, param_adj
		print err
		return err


	#############################################################
	########## Profile likelihood function ######################
	########## with human and vec data_generalized ##############
	#############################################################

	def sse_vec_proflike(self, param_fit, case_data, vec_data1, param_adj=None, fit_ind=None, adj_ind=None): 
	#def sse_vec_proflike(self, param_fit, case_data, vec_data1, vec_data2, param_adj=None, fit_ind=None, adj_ind=None):
	#def sse_vec_proflike(self, param_fit, case_data, vec_data1, vec_data2, vec_data3, param_adj=None, fit_ind=None, adj_ind=None):
		#specify parameter sets: parameters for fitting and paramter for adjusting 
		alpha = 0.14
		param_fit_ls=list(param_fit)
		param_opt_ls=list(self.param)
		param_adj_ls=list(param_adj)
		#param_fix_ls=list(param_fix)
		print param_opt_ls

		#substitute fited parameter set
		try:
			for i in range(len(fit_ind)):
				param_opt_ls.pop(fit_ind[i])
				param_opt_ls.insert(fit_ind[i], param_fit_ls[i])
		except:
			pass
		#print param_opt_ls

		#substitute adjusted parameter set
		try:
			for j in range(len(adj_ind)):
				param_opt_ls.pop(adj_ind[j])
				param_opt_ls.insert(adj_ind[j], param_adj_ls[j])
		except:
			pass
		print param_opt_ls

		param_opt_arr=np.array(param_opt_ls)
		print param_opt_arr
		beta, repH, x, pi_mua, beta_m, mu_m=param_opt_arr
		absparam=(abs(beta), abs(repH), abs(x), abs(pi_mua), abs(beta_m), abs(mu_m))
		S0, E0, I0, A0, S0_m, E0_m, I0_m=self.ini
		I0=case_data[0]/abs(repH)
		E0=2*case_data[0]/abs(repH)
		S0=1-I0-E0
		S0_m=1.0/abs(mu_m)
		temp_ini=(S0, E0, I0, A0, S0_m, E0_m, I0_m)
		

		res = spi.odeint(self.model, temp_ini, self.time_step, args=(absparam,), mxstep=2000)
		res_inc_t1=np.append(np.array([0]),spi.cumtrapz(abs(repH)*alpha*res[:,1]))
		res_inc_t0=np.append(np.array([0]),spi.cumtrapz(abs(repH)*alpha*res[0:-1,1]))
		data_mw=7*(res_inc_t1+6.0/7.0) - 7*np.append(np.array([0]), (res_inc_t0+6.0/7.0))
		err1=np.sum((case_data-data_mw)**2/np.mean(case_data)) #weighted SSE changed on 02/06/2016

		data_m_vadult=res[:,4]+res[:,5]+res[:,6] #adult mosquitoes 
		data_m_vsus=res[:,4] #adult susceptiable
		data_m_vinf=res[:,5]+res[:,6] #infected adult mosquitoes
		data_m_vlarva=res[:,3] #larva

		err2=np.sum((vec_data1-data_m_vlarva)**2/np.mean(vec_data1)) ### larva data only
		###err3=np.sum((vec_data2-data_m_vadult)**2/np.mean(vec_data2)) ### larva+adult
		###err3=np.sum((vec_data2-data_m_vsus)**2/np.mean(vec_data2)) ### larva+sus+inf
		###err4=np.sum((vec_data3-data_m_vinf)**2/np.mean(vec_data3)) ###larva+sus+inf

		err=(err1+err2)/2.0
		###err=(err1+err2+err3)/3.0
		###err=(err1+err2+err3+err4)/4.0

		print param_fit, param_adj
		print err1, err2
		###print err1, err2, err3
		###print err1, err2, err3, err4
		print err
		return err
#===end===#
