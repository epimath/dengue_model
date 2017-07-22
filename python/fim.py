########################
######### FIM ##########
#### Yu-Han Kao ########
#### kaoyh@umich.edu ###
#### 12/05/2016 ########
########################
#this script constructs Fisher information matrix

import scipy.integrate as spi
import scipy.optimize as spopt
import numpy as np
import math


#pert: the level of perturbations, the default is 1%
class fim:
	def __init__(self, model, time_step):
		self.model = model
		self.time_step = time_step
		self.description = 'fisher information matirx using numerical approximation approach'

	def Fim_x(self, param, data_case, pert=0.01):
		#X
		alpha = 0.14
		listX = []

		for i in range(len(param)):
			param_1 = np.array(param)
			param_1[i] = abs(param[i]*(1+pert))
			param_2 = np.array(param)
			param_2[i] = abs(param[i]*(1-pert))
			repH_1 = param_1[1]
			repH_2 = param_2[1]
			mu_m_1 = param_1[5]
			mu_m_2 = param_2[5]

			#initial condition1
			#human initial state
			I0_1 = data_case[0]/abs(repH_1)
			E0_1 = 2*data_case[0]/abs(repH_1)
			S0_1 = 1.0-I0_1-E0_1
			#mosquito initial state
			A0_1 = 1.0
			S0_m_1 = 1.0/abs(mu_m_1)
			E0_m_1 = 0.0
			I0_m_1 = 0.0
			ini_1 = [S0_1, E0_1, I0_1, A0_1, S0_m_1, E0_m_1, I0_m_1]

			#initial condition2
			#human initial state
			I0_2 = data_case[0]/abs(repH_2)
			E0_2 = 2*data_case[0]/abs(repH_2)
			S0_2 = 1.0-I0_2-E0_2
			#mosquito initial state
			A0_2 = 1.0
			S0_m_2 = 1.0/abs(mu_m_2)
			E0_m_2 = 0.0
			I0_m_2 = 0.0
			ini_2 = [S0_2, E0_2, I0_2, A0_2, S0_m_2, E0_m_2, I0_m_2]

			#measurement
			res1 = spi.odeint(self.model, ini_1, self.time_step, args=(param_1,))
			res2 = spi.odeint(self.model, ini_2, self.time_step, args=(param_2,))
			res1_inc_t1 = np.append(np.array([0]),spi.cumtrapz(abs(repH_1)*alpha*res1[:,1]))
			res1_inc_t0 = np.append(np.array([0]),spi.cumtrapz(abs(repH_1)*alpha*res1[0:-1,1]))
			y1 = 7*(res1_inc_t1+6.0/7.0) - 7*np.append(np.array([0]), (res1_inc_t0+6.0/7.0))
			res2_inc_t1 = np.append(np.array([0]),spi.cumtrapz(abs(repH_2)*alpha*res2[:,1]))
			res2_inc_t0 = np.append(np.array([0]),spi.cumtrapz(abs(repH_2)*alpha*res2[0:-1,1]))
			y2 = 7*(res2_inc_t1+6.0/7.0) - 7*np.append(np.array([0]), (res2_inc_t0+6.0/7.0))
			#measurement
			subX = (y1-y2)/(2*pert*param[i])
			listX.append(subX.tolist()[1::]) #ignore t0 since it's fixed in our case
		X = np.matrix(listX)
		#print 'X',X
		FIM = np.dot(X,X.transpose())
		return FIM
	
	def Fim_rank(self, FIM):
		F_rank = np.linalg.matrix_rank(FIM,tol=3.5763e-07) #set the tol to match it from matlab
		#F_rank = np.linalg.matrix_rank(FIM)
		#print np.linalg.svd(FIM)[1].max()
		#print np.finfo(FIM.dtype).eps
		#print np.linalg.svd(FIM)[1].max()*max(FIM.shape)*np.finfo(FIM.dtype).eps #default tol value in python
		return F_rank
#===end===#
