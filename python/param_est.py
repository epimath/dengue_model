########################
# Parameter estimation #
#### Yu-Han Kao ########
#### kaoyh@umich.edu ###
#### 11/05/2016 ########
########################

#this script includes cost function for parameter estimation(residuals)
#(specific for model in dengue_model.py)


import scipy.integrate as spi
import scipy.optimize as spopt
import numpy as np
import math
import pylab as pl
import copy

class param_est:
	def __init__(self, model, ini, time_step):
		self.model = model
		self.ini = ini
		self.time_step = time_step
		self.description = 'residuals: calculate sum of squared errors, can be used for the built-in optimization function of scipy'
	
	#############################################################
	########## Residuals calculation function ###################
	########## with human data ################################## 
	#############################################################

	def residuals(self, param, case_data):
		
		alpha=0.14
		beta, repH, x, pi_mua, beta_m, mu_m=param
		absparam=(abs(beta), abs(repH), abs(x), abs(pi_mua), abs(beta_m), abs(mu_m))
		S0, E0, I0, A0, S0_m, E0_m, I0_m=self.ini
		
		I0=case_data[0]/abs(repH)
		E0=2*case_data[0]/abs(repH)
		S0=1-I0-E0
		S0_m=1.0/abs(mu_m)


		temp_ini=(S0, E0, I0, A0, S0_m, E0_m, I0_m)
		res = spi.odeint(self.model, temp_ini, self.time_step, args=(absparam,), mxstep=3000)
		res_inc_t1=np.append(np.array([0]),spi.cumtrapz(abs(repH)*alpha*res[:,1]))
		res_inc_t0=np.append(np.array([0]),spi.cumtrapz(abs(repH)*alpha*res[0:-1,1]))
		data_mw=7*(res_inc_t1+6.0/7.0) - 7*np.append(np.array([0]), (res_inc_t0+6.0/7.0))
		#err=math.log(np.sum((case_data-data_mw)**2)) #log of sse
		err=np.sum((case_data-data_mw)**2/np.mean(data_mw))
		#err=np.sum((case_data-data_mw)**2/(np.mean(data_mw))**0.5) #weighted SSE with squared root demoninator 12/17/2016 try
		print err, param
		return err
#===end===#
