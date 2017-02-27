########################
######## Sigma #########
#### Yu-Han Kao ########
#### kaoyh@umich.edu ###
#### 12/05/2016 ########
########################

# compute sigma square from residuals
def sigma (data, sim_data, param):
	resisq=(sim_data-data)**2
	sigmasq=(sum(resisq))/float(len(data)-len(param))
	sigma=sigmasq**(0.5)

	return sigma
#===end===#


