import numpy as np
from sahara_work import Criticality as cr
from copy import deepcopy as cdc
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
import math
def pvaluenew(burst,boundary):
# Null hypothesis test. When P_value is very large, we could not reject
# the Null hypothesis that the distribution of burst follows power law.
# Usually, we use 0.05 as criteria.
	
	alpha, xmin, ks, Loglike = cr.tplfit(burst,boundary) #OG matlab hardcodes this to 40. I think that's correct. I'm going to do this. alpha, xmin, ks, Loglike = cr.tplfit(burst,boundary)
	
		
	print(f"xmin from pval test: {xmin}")

	# print(xmin)
	# print(alpha)
	xmax = np.max(burst)
	N   = np.size(burst)
	k   = 0
	z   = cdc(burst)
	z   = z[z>=xmin]
	n   = np.size(z)
	cdf = np.cumsum(np.histogram(z,np.arange(xmin,xmax+2))[0]/n) 

	idx = np.where(np.logical_and(xmin<=burst, burst<=xmax))[0]
	s = np.unique(burst[idx])
	smin = np.min(s)
	A = 1/np.sum(np.power(s, -alpha[0]))
	fit = np.cumsum(A*np.power(np.arange(xmin,xmax+1), -alpha[0]))
	KS = np.max(np.abs(cdf-fit))
	
	hfig, hax = plt.subplots(ncols = 1, nrows = 1)
	sns.despine()
	plt.yticks(fontsize = 13)
	plt.xticks(fontsize = 13)
	#plt.ion()
	hax.plot(np.arange(xmin,xmax+1), fit, zorder = 10000, label = 'Power law CDF', color = '#0504aa')
	hax.plot(np.arange(xmin,xmax+1), cdf, zorder = 10005, label = 'Experimental CDF', color = '#80013f')
	hax.legend()
	hax.set_title('Cumulative Distribution Function ',fontsize = 12)
	hax.set_xscale('log')
	shape, loc, scale = stats.lognorm.fit(z, floc=0)
	sig = shape
	mu = math.log(scale)
	sns.despine()
	# cfig, cax = plt.subplots(ncols = 1, nrows = 1)
	# hax.plot(np.arange(xmin,xmax), cdf)
	
	ks = []
	j = 1
	Niter = 1000

	while j<Niter:

		if not j % 400:
			print(str(j) + " loops completed")


		########################## Inverse Method  #############################################

		#N = 20*np.size(burst[burst>=xmin]) #N = 10*np.size(burst[burst>=xmin]) OG matlab has 20
		N = 10*np.size(burst[burst>=xmin])
		syn_data = np.floor((xmin-1/2)*np.power((1-np.random.uniform(0,1,N)), (1/(1-alpha[0]))) + 1/2)
		#syn_data = np.floor(np.heaviside(xmax-syn_data, 1/2) * syn_data) #syn_data = np.floor(np.heaviside(xmax-syn_data, 1/2) * syn_data) OG  matlab doesnt have a second parameter but i think its necessary for python
		syn_data = np.floor(np.heaviside(xmax-syn_data, 1/2) * syn_data)
		syn_data = np.delete(syn_data, np.where(syn_data == 0)[0])
		syn_data = syn_data[0:np.size(burst[burst>=xmin])]

		########################## Accept-Reject Method  #############################################

		# N = 2*np.size(burst[burst>=xmin+1])
		# beta = alpha[0] - 1
		# umax = np.power(xmin,-beta)  
		# u = umax * np.random.uniform(0,1,N)       
		# r = np.floor(np.power(u,-(np.power(beta,-1))))
		# syn_data = r * np.floor(np.heaviside((P(r,xmin,alpha) .* Q(xmin,xmin,beta))./(P(xmin,xmin,alpha) .* Q(r,xmin,beta)   ) - rand(1,N) ));
		# syn_data = floor(heaviside(xmax-syn_data)) .* syn_data;
		# % syn_data(syn_data==0)=[];
		# Ind = find(syn_data > 0);
		# syn_data = syn_data(Ind(1:length(burst(burst>=xmin))))

		############################################################################################

		## OKAY SO --- OG matlab doesnt index into the syn_data at all. So i'm commenting it out and trying it
		idx_syn = np.where(np.logical_and(xmin<=syn_data, syn_data<=xmax))[0]
		X = syn_data[idx_syn]

		# THIS IS BASED OFF THE MATLAB
		# X = syn_data
		alpha_syn, xmin_syn, ks_syn, Loglike_syn = cr.tplfit(syn_data,xmin) #calculate exponent for surrogated data 
		a = alpha_syn[0]   
		
		if np.abs(a-alpha[0])<=0.1 and a >1.0:   # Make sure we use the right surrogated data
			            
			n = np.size(X)
			cdf = np.cumsum(np.histogram(X,np.arange(xmin,xmax+2))[0]/n)
			
			s = np.unique(X)
			smin = min(s);      
			smax = max(s);
			A = 1/np.sum(np.power(s,-alpha[0]));
			fit = np.cumsum(A*np.power(np.arange(xmin,xmax+1), -alpha[0]));
			# A = 1/ sum(s.^-a );
			# fit = cumsum ( A*(xmin:xmax).^-a );
			# ks = np.max(np.abs(cdf - fit))

			######## ks is for surrogated and perfect power law with -alpha
			# as exponent #######################
			ks.append(np.max(np.abs(cdf - fit)))
			j = j + 1
			#hax.plot(np.arange(xmin,xmax), fit, )
			hax.plot(np.arange(xmin,xmax+1), cdf, color = '#647d8e', alpha  = 0.05, linewidth = 0.3)
	#plt.show()
	ks = np.asarray(ks)
	P_value = np.sum(np.sign(ks[ks>=KS]))/Niter

	ylim = hax.get_ylim()[-1]
	xlim = hax.get_xlim()[-1]
	txt_b = hax.text(xlim*0.5, ylim*0.5, 'p = ' + str(P_value), fontsize = 15)

	print('P_value = ' + str(P_value))
	print('KS = ' + str(KS))

	AA = [P_value, KS]

	return P_value, ks, hfig, xmin

























