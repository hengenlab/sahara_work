import numpy as np
from sahara_work import Criticality_new as cr
from sahara_work import Criticality as cr_old
import matplotlib.pyplot as plt
import seaborn as sns

def scaling_plots(Result, burst, burstMin, burstMax, alpha, T, tMin,tMax, beta, TT, Sm, sigma, fit_sigma, pltname, saveloc):
    # burst PDF
    fig1, ax1 = plt.subplots(nrows = 1, ncols = 3, figsize = [10, 6])
    pdf = np.histogram(burst, bins = np.arange(1, np.max(burst)+2))[0]
    ax1[0].plot(np.arange(1,np.max(burst)+1),pdf/np.sum(pdf), marker = 'o', markersize = 3, linestyle = 'None', color = '#2138ab', alpha = 0.75)
    ax1[0].set_yscale('log')
    ax1[0].set_xscale('log')

    x = np.arange(burstMin, burstMax+1)
    y = (np.size(np.where(burst == burstMin)[0])/np.power(burstMin, -alpha))*np.power(x, -alpha)
    y = y/np.sum(pdf)

    ax1[0].plot(x,y, color = '#c5c9c7');
    ax1[0].set_xlabel('AVsize')
    ax1[0].set_ylabel('PDF(D)')
    ax1[0].set_title('AVsize PDF, ' + str(np.round(alpha[0], 3)))

    # time pdf
    tdf = np.histogram(T,bins = np.arange(1, np.max(T)+2))[0]
    ax1[1].plot(np.arange(1,np.max(T)+1),tdf/np.sum(tdf), marker = 'o', markersize = 3, linestyle = 'None', color = '#7bfdc7', alpha = 0.75)
    ax1[1].set_yscale('log')
    ax1[1].set_xscale('log')
    sns.despine()    

    x = np.arange(tMin,tMax+1)
    y = np.size(np.where(T == tMin))/(np.power(tMin,-beta))*np.power(x,-beta)
    y = y/np.sum(tdf)
    ax1[1].plot(x,y, color = '#c5c9c7')
    ax1[1].set_xlabel('AVduration')
    ax1[1].set_ylabel('PDF(D)')
    ax1[1].set_title('AVdura PDF, ' + str(np.round(beta[0], 3))) 

    # figure out how to plot shuffled data

    #scaling relation 
    ax1[2].plot(TT, ((np.power(TT,sigma)/np.power(TT[7],sigma))*Sm[7]), label = 'pre', color = '#4b006e')
    ax1[2].plot(TT, (np.power(TT,fit_sigma[0])/np.power(TT[7],fit_sigma[0])*Sm[7]),'b', label = 'fit', linestyle = '--', color = '#826d8c')
    ax1[2].plot(TT, Sm,'o', color = '#fb7d07', alpha = 0.75)
    ax1[2].set_xscale('log')
    ax1[2].set_yscale('log')
    ax1[2].set_ylabel('<S>')
    ax1[2].set_title('Difference = ' + str(np.round(Result['df'][0], 3)))
    
    plt.tight_layout()
    plt.legend()
    plt.savefig(saveloc+"/"+pltname+'scaling_relations')

    return fig1

def AV_analysis_new(burst, T, bm, tm, pltname, saveloc, flag = 1, burst_shuffled=None, T_shuffled=None, plot_shuffled=False, plot=True):
    Result = {}
    burstMax, burstMin, alpha = cr.EXCLUDE(burst[burst < np.power(np.max(burst),.9)], bm)
    #burstMax, burstMin, alpha = cr.EXCLUDE(burst, bm)
    
    # for testing, exclude burst turned off
    # burstMax = np.max(burst)
    # burstMin = bm

    idx_burst = np.where(np.logical_and(burst<=burstMax, burst>=burstMin))[0]

    alpha, xmin, ks, L = cr.tplfit(burst[idx_burst], burstMin)

    print("burst min: ", burstMin)
    print("xmin: ", xmin)
    print("burst max:", burstMax)

    Result['burst'] = burst
    Result['alpha'] = alpha
    Result['xmin'] = burstMin
    Result['xmax'] = burstMax

    if flag == 2 :
        # pvalue test
        #Result['P_burst'], ks, hax_burst, ptest_bmin  = cr.pvaluenew(burst[idx_burst], alpha, xmin, ks, L)
        Result['P_burst'], ks, hax_burst, ptest_bmin  = cr_old.pvaluenew(burst[idx_burst], bm)




    #tMax, tMin, beta = cr.EXCLUDE(T[T < np.power(np.max(T),0.8)], tm)
    tMax, tMin, beta = cr.EXCLUDE(T, tm)

    # this is for testing, esentially canceling out the exclude function 
    # tMax = np.max(burst)
    # tMin=tm
    
    idx_time = np.where(np.logical_and(T >= tMin,T <= tMax + 1))[0]
    beta, new_tmin, tplfit_ks_time, L_t = cr.tplfit(T[idx_time], tMin)
    
    print(f'time min: {tMin}')
    print(f'new tmin: {new_tmin}')
    print(f'time max: {tMax}')

    Result['T'] = T
    Result['beta'] = beta
    Result['tmin'] = tMin
    Result['tmax'] = tMax

    if flag == 2:
        #pvalue for time
        #Result['P_t'], ks, hax_time, ptest_tmin  = cr.pvaluenew(T[idx_time], beta, tMin, tplfit_ks_time, L_t)
        Result['P_t'], ks, hax_time, ptest_tmin =  cr_old.pvaluenew(T[idx_time], tm)
    # scaling relation 
    TT = np.arange(1, np.max(T)+1)
    Sm = []
    for i in np.arange(0,np.size(TT)):
        Sm.append(np.mean(burst[np.where(T==TT[i])[0]]))
    Sm = np.asarray(Sm)
    Loc=np.where(Sm>0)[0]
    TT=TT[Loc]
    Sm=Sm[Loc]

    fit_sigma = np.polyfit(np.log(TT[np.where(np.logical_and(TT>tMin, TT<tMax))[0]]), np.log(Sm[np.where(np.logical_and(TT>tMin, TT<tMax))[0]]), 1)
    sigma = (beta - 1)/(alpha - 1)

    Result['pre'] = sigma
    Result['fit'] = fit_sigma
    Result['df'] = np.abs(sigma - fit_sigma[0])
    Result['TT'] = TT
    Result['Sm'] = Sm

    m = fit_sigma[0]
    c = fit_sigma[1]

    if plot:
        fig1 = scaling_plots(Result, burst, burstMin, burstMax, alpha, T, tMin, tMax, beta, TT, Sm, sigma, fit_sigma, pltname, saveloc)
        if flag == 2 :
            hax_burst.axes[0].set_xlabel('Size (S)', fontsize = 16)
            hax_burst.axes[0].set_ylabel('Prob(size < S)', fontsize = 16)
            hax_burst.savefig(saveloc+"/"+pltname+'pvalue_burst')

            hax_time.axes[0].set_xlabel('Duration (D)', fontsize = 16)
            hax_time.axes[0].set_ylabel('Prob(size < D)', fontsize = 16)
            hax_time.savefig(saveloc+"/"+pltname+'pvalue_time')
            Result['burst_cdf'] = hax_burst
            Result['time_cdf'] = hax_time
        Result['scaling_relation_plot'] = fig1

    return Result
