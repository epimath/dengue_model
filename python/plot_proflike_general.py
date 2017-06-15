import pylab as pl
import copy
import scipy.stats as stats
import numpy as np
import pickle
import math
from matplotlib import rc

rc('mathtext', default='regular')

#calculate threshold
sigmasq=119.476832507
threshold=sigmasq*stats.chi2.ppf(0.95,6)
########################

fname="pf_data_human_repH.pickle"#"pf_data_human_beta.pickle"
#ftemp=open(fname,'rb')

adj_param_name=fname[14:-7]
print adj_param_name

#temp=ast.literal_eval(fread) #not safe, if it's not for python structures: ast.literal_eval
temp = pickle.load(open(fname, "rb"))
#print temp 
#temp={'beta_repH_x_mu_m[3]': {0.5423457790386971: [0.24560618766389075, 30.0028143130626, True], -0.58265422096130293: [0.2523688902063602, 55.958229900817734, True], 0.052168222379837134: [4.440892098500626e-16, 15372.038982169945, True], 0.2923457790386971: [0.16158856230760343, 20.06431074640071, True], 1.4173457790386972: [0.30025688513928006, 5388.542523401168, True], 0.41734577903869707: [0.21573231786569233, 0.0, True], 3.3387662323095766: [0.3300873004318885, 169202.39526093786, True]}, 'beta_x_pi_mua_mu_m[1]': {0.5423457790386971: [1633.6112686409954, 349.9576420336353, True], -0.58265422096130293: [1798.253134969279, 582.2181980889633, True], 0.052168222379837134: [34154.70786506937, 490926.9375069752, True], 0.2923457790386971: [1559.8529708321164, 650.9957563362558, True], 1.4173457790386972: [9273.768429092856, 12928.0353367974, True], 0.41734577903869707: [1306.1073753366747, 0.0, True], 3.3387662323095766: [20701.801945942665, 35301.9013982237, True]}, 'beta_repH_x_pi_mua[5]': {0.5423457790386971: [1.181148226298753, 5.432804862424733, True], -0.58265422096130293: [1.2002567921438292, 9.54631093066374, True], 0.052168222379837134: [0.6196169618059533, 60.632132767937456, True], 0.2923457790386971: [1.0191054165045723, 5.510943031407855, True], 1.4173457790386972: [1.4265738688714467, 435.0338547121907, True], 0.41734577903869707: [1.111771103713506, 0.0, True], 3.3387662323095766: [1.6027310248657938, 5746.755328594114, True]}, 'beta_repH_pi_mua_mu_m[2]': {0.5423457790386971: [0.29864959749087305, 52.18475830100809, True], -0.58265422096130293: [0.2905275730912641, 99.55577105594526, True], 0.052168222379837134: [2505833.6086765565, 15372.012344275468, True], 0.2923457790386971: [0.46480396687690073, 27.337238537754626, True], 1.4173457790386972: [0.2436250558353356, 13231.419874216695, True], 0.41734577903869707: [0.34261392318020867, 0.0, True], 3.3387662323095766: [0.19887697004416538, 342276.7922081303, True]}, 'repH_x_pi_mua_mu_m[0]': {0.5423457790386971: [6.685721825692037, 0.02460543222632908, True], -0.58265422096130293: [6.231437494871134, 0.042819365143053095, True], 0.052168222379837134: [68.38632178105297, 0.22349087343012694, True], 0.2923457790386971: [12.301377491061697, 0.025397386049334152, True], 1.4173457790386972: [2.631120628438513, 1.4171457940985774, True], 0.41734577903869707: [8.652632499615974, 0.0, True], 3.3387662323095766: [1.1839277061547437, 9.795696760272092, True]}}
param_name=['beta', 'repH', 'x', 'pi_mua', 'beta_m', 'mu_m']
#param_name_test=param_name.remove(adj_param_name)
temp_fit_name=temp.keys()
print temp_fit_name
fit_name_list=[]
for n in temp_fit_name:
	mao=list(n.split('[')[1].rstrip(']').split(','))
	print mao
	meow=[int(m) for m in mao]
	fit_name_list.append([int(m) for m in mao])
print fit_name_list

temp_values=temp.values()
count=0
for i in temp_values: #iterate throught each fixed parameter set
	x_value=[abs(x) for x in i.keys()]
	y_value_ls=i.values()
	temp_yvalue=[] # fitted values for parameters
	temperror=[] # SSE values
	tempcov=[] #if the optimization is converged

	for y_ls in y_value_ls:
		temperror.append(y_ls[-2])
		tempcov.append(y_ls[-1])
		temp_yvalue.append(y_ls[:-2])

	small_error=[] #small sse
	small_x=[] #x value of small sse
	conv_error=[]
	conv_x=[]
	true_error=[]
	true_x=[]
	for k in range(len(temperror)):
		if temperror[k] <= 0.0000000005 and tempcov[k]==True:
		###if temperror[k] <= 0.1 and tempcov[k]==True:
			true_error.append(temperror[k])
			true_x.append(x_value[k])
		if tempcov[k]==True:
			conv_error.append(temperror[k])
			conv_x.append(x_value[k])
		if temperror[k] <= 0.1 and tempcov[k]==True:
			small_error.append(temperror[k])
			small_x.append(x_value[k])

	fit_name_all=[]
	for f in range(len(fit_name_list[count])):
		fit_name_all.append(param_name[fit_name_list[count][f]])
	fit_names='_'.join(fit_name_all)
	print fit_names
	
	print temp_yvalue
	#pl.rcParams['xtick.major.pad'] = 8
	#pl.rcParams['ytick.major.pad'] = 8

	pl.figure(figsize=(12,9))
	pl.axhline(y=threshold, ls='--', linewidth=5.5, label='Approx 95% CI', color='g')
	#pl.axhline(y=neg_threshold, ls='--', linewidth=5.5, label='Approx 95% CI', color='g')
	###pl.scatter(x_value,temperror, label='sum of squared', color='r' )
	###pl.scatter(conv_x,conv_error, label='SSE', color='r' ) #label='sum of squared'
	pl.scatter(small_x,small_error, label='SSE', color='r',s=70 ) #sse<0.1 #label='small sum of squared'
	pl.scatter(true_x,true_error, label='SSE of true value', color='b',s=70 )
	###pl.ylim(ymax=0.1, ymin=-0.02)
	pl.ylim(ymax=threshold*1.2)#, ymin=neg_threshold*1.2)
	###pl.title(adj_param_name)
	###pl.title('Profile Likelihood with Human data')
	pl.xlabel(r'$\beta_{mh}$',fontsize=30)
	###pl.xlabel(adj_param_name)
	###pl.ylabel(fit_names+'_SSE')
	pl.ylabel('Sum of Squared Error',fontsize=30)
	pl.xticks(fontsize = 26)
	pl.yticks(fontsize = 26)
	pl.show()


	#pl.xscale('log')
	#pl.legend()
	#pl.savefig('human_model_0131_error_'+param_name[i])

	###################
	'''
	#plots with break
	## adapt from broken_axis.py by Paul Ivanov 
	## https://github.com/matplotlib/matplotlib/blob/master/examples/pylab_examples/broken_axis.py
	pl.figure()
	# plot the same data on both axes
	ax.plot(small_x,small_error, marker='o',ls='None',color='red',label='SSE')
	ax.plot(true_x,true_error, marker='o',ls='None',color='blue',label='SSE of true value')
	ax.axhline(y=threshold,c='green',ls='--', label='Approx 95% CI')
	ax2.plot(small_x,small_error, marker='o',ls='None',color='red')
	ax2.plot(true_x,true_error, marker='o',ls='None',color='blue')
	#ax2.axhline(y=neg_threshold,c='green')
	#pl.ylabel('Sum of Squared Error', fontsize=20)
	pl.xlabel(r'$\beta_{hm}$',fontsize=20)

	# zoom in the data
	ax.set_ylim(120.0,160.0)  #outliers only
	ax2.set_ylim(-0.02, 0.2)  #most of the data

	# hide the spines between ax and ax2
	ax.spines['bottom'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax.xaxis.tick_top()
	ax.tick_params(labeltop='off')  #don't put tick labels at the top
	ax2.xaxis.tick_bottom()

	d = .015  # size of the diagonal lines
	kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
	ax.plot((-d, +d), (-d, +d), **kwargs)        #top-left diagonal
	ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  #top-right diagonal

	kwargs.update(transform=ax2.transAxes)  #switch to the bottom axes
	ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  #bottom-left diagonal
	ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  #bottom-right diagonal
	ax.legend()
	'''
	#pl.show()
	#pl.close()
	


	for f in range(len(fit_name_list[count])):
		fit_name=param_name[fit_name_list[count][f]]
		y_value=[]
		small_fit=[] #parameter values with small sse and converged
		for ys in range(len(temp_yvalue)):
			y_value.append(abs(temp_yvalue[ys][f]))
		for k in range(len(temperror)):
			if temperror[k] <= 0.1 and tempcov[k]==True:
				small_fit.append(y_value[k])
		
		pl.scatter(x_value, y_value,label=fit_name,color='r')
		print y_value
		pl.scatter(small_x,small_fit, label=fit_name+'_small sse', color='b')
		#print adj_param_name
		pl.title(adj_param_name)
		pl.xlabel(adj_param_name+'_param value')
		pl.ylabel(fit_name)
		#pl.xscale('log')
		#pl.legend()
		#pl.savefig('human_0118_fix_'+param_name[i]+'_test_'+temp_p_name[m])
		#pl.savefig('human_model_0131_error_'+param_name[i])
		pl.show()
		pl.close()
	
	count+=1

#===end===#

