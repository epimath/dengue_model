########################
#### Dengue model ######
#### Yu-Han Kao ########
#### kaoyh@umich.edu ###
#### 11/04/2016 ########
########################

#this script includes dengue model (model):
##human: SEIR
##vector: ASEI
#and model simulation function (ode_run)

import scipy.integrate as spi
import scipy.optimize as spopt
import numpy as np
import math
import pylab as pl
import copy

class dengue_model:

	def __init__(self):
		self.description = '.model: SEIR(human)-ASEI(vector) model; .ode_run: model simulation'


	#############################################################
	########## Model equations without climate forcing ##########   ###changed
	#############################################################

	def model(self, ini, time_step, param): #only for one year, population is scaled to 1, mosquito-->initial value 
		Y=np.zeros(7)
		V = ini
		beta, repH, x, pi_mua, beta_m, mu_m=param
		alpha=0.14
		eta=0.2
		gamma=0.1
		
		#human SEIR
		Y[0] =  - beta * V[0] * V[6]
		Y[1] = beta * V[0] * V[6] - alpha * V[1]
		Y[2] = alpha * V[1] - eta * V[2]
		#vector ASEI
		Y[3] = x * (V[4] + V[5] + V[6]) * (1 - V[3]) - pi_mua * V[3]
		Y[4] = V[3] - beta_m * V[4] * V[2] - mu_m * V[4]
		Y[5] = beta_m * V[4] * V[2] - gamma * V[5] - mu_m * V[5]
		Y[6] = gamma * V[5] - mu_m * V[6]

		return Y

	#############################################################
	########## ODE running function #############################
	#############################################################

	def ode_run(self, model, ini, time_step, param):
		RES = spi.odeint(model, ini, time_step, args=(param,))
		return RES

#===end===#

