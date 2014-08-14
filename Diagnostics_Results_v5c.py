############################################################################
# This script correlates tower data with external meteorology for a variable.  Then calulates the correlation coefficient.  Then adjusts new time series and gap fills 
# Inputs Final DF
#
# Programmed by Jason (Dec 1, 2012)
############################################################################

import pandas as pd
import numpy as np
import os
import datetime as dt
import pylab as pl
import meteorologicalfunctions as metfuncs
from dateutil.relativedelta import relativedelta
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

from pylab import *
from numpy import array
import matplotlib.dates as mdates
import matplotlib.cbook as cbook

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

def Doplots_diurnal(mypathforResults,PlottingDF,variable_to_fill, Site_ID,units,item):
    ANN_label=str(item+"_NN")     #Do Monthly Plot
    print "Doing diurnal plot for month "
    #Do Diurnal Plots for all 12 months
    #create an X axis series for all 24 hours
    t = np.arange(1, 25, 1)
    NN_label='Fc'
    Plottemp = PlottingDF[[NN_label,item]][PlottingDF['day_night']!=3]
    #Plottemp = PlottingDF[[NN_label,item]].dropna(how='any')
	    
    figure(1)
    pl.subplot(321)
    pl.title('Diurnal '+item+' month = 1')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==1)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==1)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')    
 
    pl.subplot(322)
    pl.title('Diurnal '+item+' month = 3')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==3)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==3)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')      

    pl.subplot(323)
    pl.title('Diurnal '+item+' month = 5')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==5)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==5)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')      
    
    pl.subplot(324)
    pl.title('Diurnal '+item+' month = 7')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==7)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==7)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')  
    
    pl.subplot(325)
    pl.title('Diurnal '+item+' month = 9')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==9)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==9)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')  
    
    pl.subplot(326)
    pl.title('Diurnal '+item+' month = 11')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==11)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==11)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')  
    
    
    figure(1)
    pl.suptitle('ANN ensemble diurnal average for variable '+item+' at '+Site_ID)
    pl.subplots_adjust(top=0.85)
    pl.tight_layout()  
    pl.savefig(mypathforResults+'/ANN ensemble diurnal average for variable '+item+' at '+Site_ID)
    #pl.show() 
    pl.close()

def Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq):   
   
    print "Doing Mean plots for "+Site_ID
    #create a string of the items to be plotted to be used later to save filename
    list_string=''
    for z in list_in:
	#cretae string list
	list_string=list_string+' '+z

    #pdf = PdfPages(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq)+'.pdf')

    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    number_of_subplots=len(list_in)
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
      
    for i,v in enumerate(xrange(number_of_subplots)):
	item=list_in[i]
	item_Con=item+'_Con'
	item_Con_flag=item+'_Con_QCFlag'
	
	if item_Con=="Fc_ustar_Con" : item_Con="Fc_ustar"
	if item_Con=="Fc_ustar_Con" : item_Con_flag="Fc_Con_QCFlag"
	
	by = lambda x: lambda y: getattr(y, x)
	xdata1a=New_combined[item_Con].groupby([by('year'),by(plot_freq)]).mean()	
	xdata1b=New_combined[[item_Con,item_Con_flag]].groupby([by('year'),by(plot_freq)]).mean()
	#Here get a data series that is only periods where the percentage gap filled more than 30%
	xdata1b[xdata1b[item_Con_flag]>30.0]=np.nan
	xdata1b=xdata1b[item_Con]
	
	#ydata=np.arange(len(xdata1a))	
	
	startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
	enddate1=New_combined.index[len(New_combined)-1]
	endplot=dt.datetime(enddate1.year,enddate1.month,1)
	totalpoints=len(New_combined)
	if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
	if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
	if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
	if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])
	
	ax=pl.subplot(number_of_subplots,1,v)
	#ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
	ax.plot(ydata,xdata1a, color='#00BFFF',label=item) #DeepSkyBlue 
	ax.plot(ydata,xdata1b, color='#00008B',linewidth=2,label=item)	 #DarkBlue 
	
	ax.legend(loc='upper left')
	pl.ylabel(item+'\n'+ '('+get_units(item,Ws_label)+')')
	# format the ticks
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	ax.xaxis.set_minor_locator(months)
	ax.grid(True, which='both') 
	#horiz line
	ax.axhline(y=0,color='k')	
	
	#if v != number_of_subplots:
	#    pl.setp(ax.get_xticklabels(), visible=False)
	
	datemin = dt.date(ydata.min().year, 1, 1)
	datemax = dt.date(ydata.max().year+1, 1, 1)
	ax.set_xlim(datemin, datemax)
    
    pl.suptitle('Mean plots for for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq))  
    
    #pl.savefig(pdf, format='pdf')
    pl.savefig(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID + ' freq ' + str(plot_freq))
    #pl.show()
    pl.close()    
    
def Doplots_at_daily(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq):   
   
    print "Doing Mean plots for "+Site_ID
    #create a string of the items to be plotted to be used later to save filename
    running_freq=30    
    list_string=''
    for z in list_in:
	#cretae string list
	list_string=list_string+' '+z

    #pdf = PdfPages(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq)+'.pdf')

    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    number_of_subplots=len(list_in)
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
      
    for i,v in enumerate(xrange(number_of_subplots)):
	item=list_in[i]
	item_Con=item+'_Con'
	
	if item_Con=="Fc_ustar_Con" : item_Con="Fc_ustar"
	
	by = lambda x: lambda y: getattr(y, x)
	xdata1a=New_combined[item_Con].groupby([by('year'),by(plot_freq)]).mean()	
	
	#ydata=np.arange(len(xdata1a))	
	
	startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
	enddate1=New_combined.index[len(New_combined)-1]
	endplot=dt.datetime(enddate1.year,enddate1.month,1)
	totalpoints=len(New_combined)
	if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
	if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
	if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
	if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])
	
	ax=pl.subplot(number_of_subplots,1,v)
	#ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
	ax.plot(ydata,xdata1a, color='#90EE90',label=item)	#'LightGreen '
	#Plot running mean
	xdata1_smooth=smooth(xdata1a,running_freq,window='flat')[0:len(ydata)]
	ax.plot(ydata,xdata1_smooth, color='#008000',linewidth=2, label=item +' '+str(running_freq)+' day run mean')  #green
	
	ax.legend(loc='upper left')
	pl.ylabel(item+'\n'+ '('+get_units(item,Ws_label)+')')
	# format the ticks
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	ax.xaxis.set_minor_locator(months)
	ax.grid(True, which='both') 
	#horiz line
	ax.axhline(y=0,color='k')	
	
	#if v != number_of_subplots:
	#    pl.setp(ax.get_xticklabels(), visible=False)
	
	datemin = dt.date(ydata.min().year, 1, 1)
	datemax = dt.date(ydata.max().year+1, 1, 1)
	ax.set_xlim(datemin, datemax)
    
    pl.suptitle('Mean plots for for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq))  
    
    #pl.savefig(pdf, format='pdf')
    pl.savefig(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID + ' freq ' + str(plot_freq))
    #pl.show()
    pl.close()        


def Doplots_at_daily_other(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq):   
   
    print "Doing Other plots for "+Site_ID
    running_freq=30
    #Create new data and variables
    by = lambda x: lambda y: getattr(y, x)
    tempdata_other=New_combined.groupby([by('year'),by(plot_freq)]).mean()	    
    #Calculate Bowen ratio daily averages used
    tempdata_other['BR']=tempdata_other['Fh_Con']/tempdata_other['Fe_Con']
    #Calculate Evaporative fraction daily averages used    
    tempdata_other['EF']=tempdata_other['Fe_Con']/tempdata_other['Fn_Con']    
    #Calculate Evaporative fraction daily averages used    
    tempdata_other['WUE']=tempdata_other['Fc_Con']/tempdata_other['Fe_Con'] 
    #Calculate Evaporative fraction daily averages used    
    tempdata_other['RUE']=tempdata_other['Fc_Con']/tempdata_other['Fsd_Con'] 
    #Calculate Evaporative fraction daily averages used    
    tempdata_other['EBC']=(tempdata_other['Fh_Con']+tempdata_other['Fe_Con'])/(tempdata_other['Fn_Con'] -tempdata_other['Fg_Con'] )
    
    
    #create a string of the items to be plotted to be used later to save filename
    list_string=''
    for z in list_in:
	#cretae string list
	list_string=list_string+' '+z

    #pdf = PdfPages(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq)+'.pdf')

    fig=pl.figure(1, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='k')
    number_of_subplots=len(list_in)
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
      
    for i,v in enumerate(xrange(number_of_subplots)):
	item=list_in[i]
	xdata1a=tempdata_other[item]	
	
	startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
	enddate1=New_combined.index[len(New_combined)-1]
	endplot=dt.datetime(enddate1.year,enddate1.month,1)
	totalpoints=len(New_combined)
	if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
	if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
	if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
	if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])
	
	ax=pl.subplot(number_of_subplots,1,v)
	#ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
	ax.plot(ydata,xdata1a, color='#9370DB',label=(item + ' at '+ plot_freq))	#'LightGreen '
	#Plot running mean	
	xdata1_smooth=smooth(xdata1a,running_freq,window='flat')[0:len(ydata)]
	
	ax.plot(ydata,xdata1_smooth, color='#800080',linewidth=2, label=(item +' '+str(running_freq)+' day running mean'))
	#Calculate limits for plot.  Take the running mean and add some margins
	ylimit_max=max(xdata1_smooth)
	if ylimit_max>0:
	    ylimit_max=ylimit_max*1.2
	else:
	    ylimit_max=ylimit_max*0.8	
	    
	ylimit_min=min(xdata1_smooth)
	if ylimit_min<0:
	    ylimit_min=ylimit_min*1.2
	else:
	    ylimit_min=ylimit_min*0.8
	pl.ylim((ylimit_min,ylimit_max))
	ax.legend(loc='upper left')
	pl.ylabel(item,size=14)
	# format the ticks
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	ax.xaxis.set_minor_locator(months)
	ax.grid(True, which='both') 
	#horiz line
	ax.axhline(y=0,color='k')	
	
	#if v != number_of_subplots:
	#    pl.setp(ax.get_xticklabels(), visible=False)
	
	datemin = dt.date(ydata.min().year, 1, 1)
	datemax = dt.date(ydata.max().year+1, 1, 1)
	ax.set_xlim(datemin, datemax)
    
    pl.suptitle('Mean plots for for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq),size=20)  

    pl.savefig(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID + ' freq ' + str(plot_freq))
    #pl.show()
    pl.close()

def Doplots_at_daily_carbon_g(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq):   
   
    print "Doing Carbon plots in units mg for "+Site_ID
    #Create new data and variables
    running_freq=30
    by = lambda x: lambda y: getattr(y, x)
    tempdata_other=New_combined.groupby([by('year'),by(plot_freq)]).mean()	     

    fig=pl.figure(1, figsize=(17, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
    monthsFmt = mdates.DateFormatter('%M')
    
    xdata1a=tempdata_other['GPP_Con']*12/44/1000000*44.01*60*60*24
    #xdata2a=tempdata_other['Fc_Con']*12/44/1000000*44.01*60*60*24	
    xdata2a=tempdata_other['Fc_ustar']*12/44/1000000*44.01*60*60*24
    xdata3a=tempdata_other['Fre_Con']*12/44/1000000*44.01*60*60*24  
    
    startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
    enddate1=New_combined.index[len(New_combined)-1]
    endplot=dt.datetime(enddate1.year,enddate1.month,1)
    totalpoints=len(New_combined)
    if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
    if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
    if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
    if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])    
 
    xdata1_smooth=smooth(xdata1a,running_freq,window='flat')[0:len(ydata)]
    xdata2_smooth=smooth(xdata2a,running_freq,window='flat')[0:len(ydata)]
    xdata3_smooth=smooth(xdata3a,running_freq,window='flat')[0:len(ydata)]
	
    pl.plot(ydata,xdata1a, color='#90EE90',label=('GPP at '+ plot_freq))	
    pl.plot(ydata,xdata1_smooth, color='#006400',linewidth=2.5, label=('GPP '+str(running_freq)+' day running mean'))

    pl.plot(ydata,xdata2a, color='#FFD700',label=('Fc at '+ plot_freq))	
    pl.plot(ydata,xdata2_smooth, color='#DAA520',linewidth=2.5, label=('Fc ustar '+str(running_freq)+' day running mean'))
    
    pl.plot(ydata,xdata3a, color='#FF6347',label=('Fre at '+ plot_freq))	
    pl.plot(ydata,xdata3_smooth, color='#B22222',linewidth=2.5, label=('Fre '+str(running_freq)+' day running mean'))
 
    #Calculate limits for plot.  Take the running mean and add some margins
    ylimit_max=max(max(xdata1_smooth),max(xdata2_smooth),max(xdata3_smooth))
    if ylimit_max>0:
	ylimit_max=ylimit_max*1.2
    else:
	ylimit_max=ylimit_max*0.8	
	
    ylimit_min=min(min(xdata1_smooth),min(xdata2_smooth),min(xdata3_smooth))
    if ylimit_min<0:
	ylimit_min=ylimit_min*1.2
    else:
	ylimit_min=ylimit_min*0.8
    pl.ylim((ylimit_min,ylimit_max))

    pl.legend(loc='lower right')
    pl.ylabel('CO2 flux (g C m-2 d-1)',size=14)
    ## format the ticks
    #ax.xaxis.set_major_locator(years)
    #ax.xaxis.set_major_formatter(yearsFmt)
    #ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_minor_formatter(monthsFmt)
    #pl.grid(True, which='both') 
    ax.axhline(y=0,color='k')
    #datemin = dt.date(ydata.min().year, 1, 1)
    #datemax = dt.date(ydata.max().year+1, 1, 1)
    #ax.set_xlim(datemin, datemax)
    
    #pl.savefig(pdf, format='pdf')
    pl.suptitle('Timeseries Carbon plot for '+Site_ID + ' freq ' + str(plot_freq),size=24)
    pl.savefig(mypathforResults+'/Results carbon plots in g C m-2 d-1_'+Site_ID + ' freq ' + str(plot_freq))
    #pl.show()
    pl.close()

def Doplots_at_daily_EB(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq):   
   
    print "Doing Energy Balance plots for "+Site_ID
    #Create new data and variables
    running_freq=30
    by = lambda x: lambda y: getattr(y, x)
    tempdata_other=New_combined.groupby([by('year'),by(plot_freq)]).mean()	     

    fig=pl.figure(1, figsize=(17, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
    monthsFmt = mdates.DateFormatter('%M')
    
    xdata1a=tempdata_other['Fn_Con']
    xdata2a=tempdata_other['Fe_Con']	
    xdata3a=tempdata_other['Fh_Con'] 
    xdata4a=tempdata_other['Fg_Con']
    
    startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
    enddate1=New_combined.index[len(New_combined)-1]
    endplot=dt.datetime(enddate1.year,enddate1.month,1)
    totalpoints=len(New_combined)
    if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
    if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
    if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
    if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])    
 
    xdata1_smooth=smooth(xdata1a,running_freq,window='flat')[0:len(ydata)]
    xdata2_smooth=smooth(xdata2a,running_freq,window='flat')[0:len(ydata)]
    xdata3_smooth=smooth(xdata3a,running_freq,window='flat')[0:len(ydata)]
    xdata4_smooth=smooth(xdata4a,running_freq,window='flat')[0:len(ydata)]

    pl.plot(ydata,xdata1a, color='#848484')	
    pl.plot(ydata,xdata1_smooth, color='#000000',linewidth=2.5, label=('Fn '+str(running_freq)+' day running mean'))

    pl.plot(ydata,xdata2a, color='#81BEF7')	
    pl.plot(ydata,xdata2_smooth, color='#0101DF',linewidth=2.5, label=('Fe '+str(running_freq)+' day running mean'))
    
    pl.plot(ydata,xdata3a, color='#F5A9A9')	
    pl.plot(ydata,xdata3_smooth, color='#B40404',linewidth=2.5, label=('Fh '+str(running_freq)+' day running mean'))
    
    pl.plot(ydata,xdata4a, color='#58FA58')	
    pl.plot(ydata,xdata4_smooth, color='#31B404',linewidth=2.5, label=('Fg '+str(running_freq)+' day running mean'))    
 
    #Calculate limits for plot.  Take the running mean and add some margins
    ylimit_max=max(max(xdata1_smooth),max(xdata2_smooth),max(xdata3_smooth),max(xdata4_smooth))
    if ylimit_max>0:
	ylimit_max=ylimit_max*1.2
    else:
	ylimit_max=ylimit_max*0.8	
	
    ylimit_min=min(min(xdata1_smooth),min(xdata2_smooth),min(xdata3_smooth),min(xdata4_smooth))
    if ylimit_min<0:
	ylimit_min=ylimit_min*1.2
    else:
	ylimit_min=ylimit_min*0.8
    pl.ylim((ylimit_min,ylimit_max))

    pl.legend(loc='upper right')
    pl.ylabel('Energy Balance flux Daily average (W m-2)',size=14)
    ## format the ticks
    #ax.xaxis.set_major_locator(years)
    #ax.xaxis.set_major_formatter(yearsFmt)
    #ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_minor_formatter(monthsFmt)
    #pl.grid(True, which='both') 
    ax.axhline(y=0,color='k')
    #datemin = dt.date(ydata.min().year, 1, 1)
    #datemax = dt.date(ydata.max().year+1, 1, 1)
    #ax.set_xlim(datemin, datemax)
    
    #pl.savefig(pdf, format='pdf')
    pl.suptitle('Timeseries EB plot for '+Site_ID + ' freq ' + str(plot_freq),size=24)
    pl.savefig(mypathforResults+'/Results EB plots '+Site_ID + ' freq ' + str(plot_freq))
    #pl.show()
    pl.close()
    
def Doplots_at_daily_carbon_umol(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq):   
   
    print "Doing Carbon plots for "+Site_ID
    #Create new data and variables
    running_freq=30
    by = lambda x: lambda y: getattr(y, x)
    tempdata_other=New_combined.groupby([by('year'),by(plot_freq)]).mean()	     

    fig=pl.figure(1, figsize=(17, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
    monthsFmt = mdates.DateFormatter('%M')
    
    xdata1a=tempdata_other['GPP_Con']	
    #xdata2a=tempdata_other['Fc_Con']	
    xdata2a=tempdata_other['Fc_ustar']    
    xdata3a=tempdata_other['Fre_Con']  
    
    startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
    enddate1=New_combined.index[len(New_combined)-1]
    endplot=dt.datetime(enddate1.year,enddate1.month,1)
    totalpoints=len(New_combined)
    if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
    if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
    if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
    if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])    
 
    xdata1_smooth=smooth(xdata1a,running_freq,window='flat')[0:len(ydata)]
    xdata2_smooth=smooth(xdata2a,running_freq,window='flat')[0:len(ydata)]
    xdata3_smooth=smooth(xdata3a,running_freq,window='flat')[0:len(ydata)]
	
    pl.plot(ydata,xdata1a, color='#90EE90',label=('GPP at '+ plot_freq))	
    pl.plot(ydata,xdata1_smooth, color='#006400',linewidth=2.5, label=('GPP '+str(running_freq)+' day running mean'))

    pl.plot(ydata,xdata2a, color='#FFD700',label=('Fc ustar at '+ plot_freq))	
    pl.plot(ydata,xdata2_smooth, color='#DAA520',linewidth=2.5, label=('Fc '+str(running_freq)+' day running mean'))
    
    pl.plot(ydata,xdata3a, color='#FF6347',label=('Fre at '+ plot_freq))	
    pl.plot(ydata,xdata3_smooth, color='#B22222',linewidth=2.5, label=('Fre '+str(running_freq)+' day running mean'))
 
    #Calculate limits for plot.  Take the running mean and add some margins
    ylimit_max=max(max(xdata1_smooth),max(xdata2_smooth),max(xdata3_smooth))
    if ylimit_max>0:
	ylimit_max=ylimit_max*1.2
    else:
	ylimit_max=ylimit_max*0.8	
	
    ylimit_min=min(min(xdata1_smooth),min(xdata2_smooth),min(xdata3_smooth))
    if ylimit_min<0:
	ylimit_min=ylimit_min*1.2
    else:
	ylimit_min=ylimit_min*0.8
    pl.ylim((ylimit_min,ylimit_max))

    pl.legend(loc='lower right')
    pl.ylabel('CO2 flux (umol CO2 m-2 s-1)',size=14)
    ## format the ticks
    #ax.xaxis.set_major_locator(years)
    #ax.xaxis.set_major_formatter(yearsFmt)
    #ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_minor_formatter(monthsFmt)
    #pl.grid(True, which='both') 
    ax.axhline(y=0,color='k')
    #datemin = dt.date(ydata.min().year, 1, 1)
    #datemax = dt.date(ydata.max().year+1, 1, 1)
    #ax.set_xlim(datemin, datemax)
    
    #pl.savefig(pdf, format='pdf')
    pl.suptitle('Timeseries Carbon plot for '+Site_ID + ' freq ' + str(plot_freq),size=24)
    pl.savefig(mypathforResults+'/Results carbon plots in umol CO2 m-2 s-1_'+Site_ID + ' freq ' + str(plot_freq))
    #pl.show()
    pl.close()


def get_units(x,Ws_label):
    if x == 'Ta':
	units = 'oC'
    elif x == 'ps':
	units = 'hPa'
    elif x == 'Ah':
	units = 'g m-3'
    elif x == Ws_label:
	units = 'm s-1'	
    elif x == 'Fsd':
	units = 'W m-2'
    elif x == 'Fsu':
	units = 'W m-2'
    elif x == 'Precip':
	units = 'mm'
    elif x == 'Fld':
	units = 'W m-2'
    elif x == 'Flu':
	units = 'W m-2'	
    elif x == 'Fc':
	units = 'umol m-2 s-1'
    elif x == 'Fc_ustar':
	units = 'umol m-2 s-1'
    elif x == 'Fe':
	units = 'W m-2'
    elif x == 'Fh':
	units = 'W m-2'
    elif x == 'Fg':
	units = 'W m-2'	
    else:
	units = ''
    return units

def Doplots_monthly_diff(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label):   
    
    print "Doing Mean Monthly DIFF (Best Alt Met - Tower) plot for "+Site_ID
    #create a string of the items to be plotted to be used later to save filename
    list_string=''
    for z in list_in:
	#cretae string list
	list_string=list_string+' '+z

    #pdf = PdfPages(mypathforResults+'/Mean Monthly DIFF plot for '+Site_ID+'.pdf')

    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    number_of_subplots=len(list_in)
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
      
    for i,v in enumerate(xrange(number_of_subplots)):
	item=list_in[i]
	item_Con=item+'_Con'
	if item in ['Fc','Fe','Fh','Fg']:
	    item_Corr=item+'_NN'
	else:
	    item_Corr=item+'_Corr'	
	v = v+1
	
	xdata1a=New_combined[item_Corr].groupby([lambda x: x.year,lambda x: x.month]).mean()-New_combined[item].groupby([lambda x: x.year,lambda x: x.month]).mean()
	#ydata=np.arange(len(xdata1a))	
	
	startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
	enddate1=New_combined.index[len(New_combined)-1]
	endplot=dt.datetime(enddate1.year,enddate1.month,1)
	totalpoints=len(New_combined)
	ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])

	ax=pl.subplot(number_of_subplots,1,v)
	ax.bar(ydata,xdata1a, width=20, color='DarkMagenta',label=item)	
	ax.legend(loc='upper left')
	pl.ylabel(item+'\n'+ '('+get_units(item,Ws_label)+')')
	# format the ticks
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	ax.xaxis.set_minor_locator(months)
	ax.grid(True, which='both') 
	#horiz line
	ax.axhline(y=0,color='k')	
	
	if v != number_of_subplots:
	    pl.setp(ax.get_xticklabels(), visible=False)
	
	datemin = dt.date(ydata.min().year, 1, 1)
	datemax = dt.date(ydata.max().year+1, 1, 1)
	ax.set_xlim(datemin, datemax)
    
    pl.suptitle('Mean Monthly DIFF (Best Alt Met - Tower) plot for for '+list_string+ ' at '+Site_ID)  
    
    #pl.savefig(pdf, format='pdf')
    pl.savefig(mypathforResults+'/Mean Monthly DIFF plot for '+list_string+ ' at '+Site_ID)
    #pl.show()
    pl.close()
        
def Plot_Pct_nans(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in):   
        list_string=''
	for z in list_in:
	    #cretae string list
	    list_string=list_string+' '+z
	    
	print 'Doing Plots of missing data for '+list_string+ ' at '+Site_ID
	
	#pdf = PdfPages(mypathforResults+'/Plots of missing data for '+list_string+ ' at '+Site_ID+'.pdf')
    
	fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	number_of_subplots=len(list_in)
	years    = mdates.YearLocator()   # every year
	months   = mdates.MonthLocator()  # every month
	yearsFmt = mdates.DateFormatter('%Y') 
	  
	for i,v in enumerate(xrange(number_of_subplots)):
	    item=list_in[i]
	    v = v+1
	    
	    if item == 'Fc_ustar': item = 'Fc'
	    
	    Pct_nan_DF=pd.read_csv(mypathforResults+'/'+'Nan counts and Pct filled for '+item+' at ' +Site_ID+'.csv')
	    xdata1a=Pct_nan_DF['Pct_nan']
	    xdata1b=Pct_nan_DF['Pct_notfilled']
	    totalpoints=len(Pct_nan_DF)
	    startplot=dt.datetime(int(Pct_nan_DF.ix[0,0]),int(Pct_nan_DF.ix[0,1]),1)
	    endplot=dt.datetime(int(Pct_nan_DF.ix[(totalpoints-1),0]),int(Pct_nan_DF.ix[(totalpoints-1),1]),1)
	    ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
    
	    ax=pl.subplot(number_of_subplots,1,v)
	    ax.bar(ydata,xdata1a, width=20, color='orange',label='Percent Nan')	
	    ax.bar(ydata,xdata1b, width=20, color='g',label='Percent not filled')	
	    ax.legend(loc='upper left')
	    pl.ylabel(item+ ' % missing')
	    # format the ticks
	    ax.xaxis.set_major_locator(years)
	    ax.xaxis.set_major_formatter(yearsFmt)
	    ax.xaxis.set_minor_locator(months)
	    ax.grid(True, which='both') 
	    #horiz line
	    ax.axhline(y=0,color='k')	
	    
	    if v != number_of_subplots:
		pl.setp(ax.get_xticklabels(), visible=False)
	    
	    datemin = dt.date(ydata.min().year, 1, 1)
	    datemax = dt.date(ydata.max().year+1, 1, 1)
	    ax.set_xlim(datemin, datemax)
	    #except:
		#pass
	
	pl.suptitle('Plots of missing data for '+list_string+ ' at '+Site_ID)  

	#pl.savefig(pdf, format='pdf')
	pl.savefig(mypathforResults+'/Plots of missing data for '+list_string+ ' at '+Site_ID)
	#pl.show()
	pl.close()    
	print 'Closed'

def regressionEBclosure(mypathforResults,New_combined,Site_ID,startdate,enddate):    
    print "Doing Energy Balance closure plots for "+Site_ID

    #First plot use ALL data and years
    ###################################
    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    #Do first plot - All hours
    #============================
    #Get data for all hours - so dont select any
    tempdata=New_combined[['Fe','Fh','Fg','Fn']].dropna(axis=0,how='any')
    xdata=tempdata['Fn']-tempdata['Fg']
    ydata=tempdata['Fe']+tempdata['Fh']
    ax1=pl.subplot(2,2,1)
    pl.title=('Energy Balance closure ALL hours for '+Site_ID )
    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
    ax1.plot(ydata,xdata, 'o', color='#00BFFF') #DeepSkyBlue 
    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 

    #Get regression stats
    #In this function calculate the linear regression independantly rather than using the NN stats.
    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(ydata,xdata)

    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
    pl.figtext(0.36,0.57,graphtext1, bbox=dict())

    x = np.linspace(min(xdata),max(xdata))
    y = slope * x + intercept
    ax1.plot(x, y, linewidth = 2,label='ALL hours')    
    
    ax1.legend(loc='upper left')
    pl.ylabel('Fh + Fe (W m-2)')
    pl.xlabel('Fn - Fg (W m-2)')

    
    #Do second plot - Daytime hours
    #=================================
    #Get data for all hours - so dont select any
    tempdata=New_combined[['Fe','Fh','Fg','Fn']][New_combined['day_night']==1].dropna(axis=0,how='any')
    xdata=tempdata['Fn']-tempdata['Fg']
    ydata=tempdata['Fe']+tempdata['Fh']
    ax2=pl.subplot(2,2,2)
    pl.title=('Energy Balance closure DAYTIME hours for '+Site_ID )  
    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
    ax2.plot(ydata,xdata,  'o', color='#00BFFF') #DeepSkyBlue 
    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 

    #Get regression stats
    #In this function calculate the linear regression independantly rather than using the NN stats.
    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(ydata,xdata)

    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
    pl.figtext(0.8,0.57,graphtext1, bbox=dict())

    x = np.linspace(min(xdata),max(xdata))
    y = slope * x + intercept
    ax2.plot(x, y, linewidth = 2,label='DAYTIME hours')    
    
    ax2.legend(loc='upper left')
    pl.ylabel('Fh + Fe (W m-2)')
    pl.xlabel('Fn - Fg (W m-2)')
    
    #Do third plot - Nighttime hours
    #===============================
    #Get data for all hours - so dont select any
    tempdata=New_combined[['Fe','Fh','Fg','Fn']][New_combined['day_night']!=1].dropna(axis=0,how='any')
    xdata=tempdata['Fn']-tempdata['Fg']
    ydata=tempdata['Fe']+tempdata['Fh']
    ax3=pl.subplot(2,2,3)
    pl.title=('Energy Balance closure NIGHTTIME hours for '+Site_ID )    
    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
    ax3.plot(ydata,xdata,  'o', color='#00BFFF') #DeepSkyBlue 
    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 

    #Get regression stats
    #In this function calculate the linear regression independantly rather than using the NN stats.
    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(ydata,xdata)

    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
    pl.figtext(0.36,0.12,graphtext1, bbox=dict())

    x = np.linspace(min(xdata),max(xdata))
    y = slope * x + intercept
    ax3.plot(x, y, linewidth = 2,label='NIGHTTIME hours')    
    
    ax3.legend(loc='upper left')
    pl.ylabel('Fh + Fe (W m-2)')
    pl.xlabel('Fn - Fg (W m-2)')
    
    #Do forth plot - DAILY AVG hours
    #============================
    #Get data for all hours - so dont select any
    tempdata=New_combined[['Fe','Fh','Fg','Fn']].dropna(axis=0,how='any')
    
    by = lambda x: lambda y: getattr(y, x)
    tempdata=tempdata.groupby([by('year'),by('day')]).mean()    


    xdata=tempdata['Fn']-tempdata['Fg']
    ydata=tempdata['Fe']+tempdata['Fh']
    ax4=pl.subplot(2,2,4)
    pl.title=('Energy Balance closure DAILY average for '+Site_ID )  
    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
    ax4.plot(ydata,xdata,  'o', color='#00BFFF') #DeepSkyBlue 
    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 

    #Get regression stats
    #In this function calculate the linear regression independantly rather than using the NN stats.
    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(ydata,xdata)

    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
    pl.figtext(0.8,0.12,graphtext1, bbox=dict())

    x = np.linspace(min(xdata),max(xdata))
    y = slope * x + intercept
    ax4.plot(x, y, linewidth = 2,label='DAILY averages')    
    
    ax4.legend(loc='upper left')
    pl.ylabel('Fh + Fe (W m-2)')
    pl.xlabel('Fn - Fg (W m-2)')    
    pl.suptitle('Energy balance closure at '+Site_ID+ ' '+str(startdate.year)+ ' to '+str(enddate.year),size=20)  
    
    pl.savefig(mypathforResults+'/'+'Energy balance closure at ' +Site_ID + ' all years') 
    #pl.show()
    pl.close()    
    print 'Closed'

    #Now do PLOTS BY YEAR
    #===========================
    tempdata_grouped=New_combined[['Fe','Fh','Fg','Fn','day_night']].groupby([lambda x: x.year])
    
    for name, group in tempdata_grouped:
	try:
		#Get year for plot lables
	    plotyear=group.index[0].year
	    
	    #First plot use ALL data and years
	    ###################################
	    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	    #Do first plot - All hours
	    #============================
	    #Get data for all hours - so dont select any
	    tempdata=group.dropna(axis=0,how='any')
	    xdata=tempdata['Fn']-tempdata['Fg']
	    ydata=tempdata['Fe']+tempdata['Fh']
	    ax1=pl.subplot(2,2,1)
	    pl.title=('Energy Balance closure ALL hours for '+Site_ID )
	    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
	    ax1.plot(ydata,xdata, 'o', color='#00CED1') #DeepSkyBlue 
	    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(ydata,xdata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.36,0.57,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax1.plot(x, y, linewidth = 2,label='ALL hours')    
	    
	    ax1.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')
	
	    
	    #Do second plot - Daytime hours
	    #=================================
	    #Get data for all hours - so dont select any
	    tempdata=group[['Fe','Fh','Fg','Fn']][group['day_night']==1].dropna(axis=0,how='any')
	    xdata=tempdata['Fn']-tempdata['Fg']
	    ydata=tempdata['Fe']+tempdata['Fh']
	    ax2=pl.subplot(2,2,2)
	    pl.title=('Energy Balance closure DAYTIME hours for '+Site_ID )  
	    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
	    ax2.plot(ydata,xdata,  'o', color='#00CED1') #DeepSkyBlue 
	    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(ydata,xdata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.8,0.57,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax2.plot(x, y, linewidth = 2,label='DAYTIME hours')    
	    
	    ax2.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')
	    
	    #Do third plot - Nighttime hours
	    #===============================
	    #Get data for all hours - so dont select any
	    tempdata=group[['Fe','Fh','Fg','Fn']][group['day_night']!=1].dropna(axis=0,how='any')
	    xdata=tempdata['Fn']-tempdata['Fg']
	    ydata=tempdata['Fe']+tempdata['Fh']
	    ax3=pl.subplot(2,2,3)
	    pl.title=('Energy Balance closure NIGHTTIME hours for '+Site_ID )    
	    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
	    ax3.plot(ydata,xdata,  'o', color='#00CED1') #DeepSkyBlue 
	    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(ydata,xdata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.36,0.12,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax3.plot(x, y, linewidth = 2,label='NIGHTTIME hours')    
	    
	    ax3.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')
	    
	    #Do forth plot - DAILY AVG hours
	    #============================
	    #Get data for all hours - so dont select any
	    tempdata=group[['Fe','Fh','Fg','Fn']].dropna(axis=0,how='any')
	    
	    by = lambda x: lambda y: getattr(y, x)
	    tempdata=tempdata.groupby([by('day')]).mean()    
	
	
	    xdata=tempdata['Fn']-tempdata['Fg']
	    ydata=tempdata['Fe']+tempdata['Fh']
	    ax4=pl.subplot(2,2,4)
	    pl.title=('Energy Balance closure DAILY average for '+Site_ID )  
	    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
	    ax4.plot(ydata,xdata,  'o', color='#00CED1') #DeepSkyBlue 
	    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(ydata,xdata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.8,0.12,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax4.plot(x, y, linewidth = 2,label='DAILY averages')    
	    
	    ax4.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')    
	    pl.suptitle('Energy balance closure at '+Site_ID+ ' for year '+str(plotyear),size=20)  
	    
	    pl.savefig(mypathforResults+'/'+'Energy balance closure at '+Site_ID+ ' for year '+str(plotyear)) 
	    #pl.show()
	    pl.close()    
	except:
	    print "WARNING plot not completed for ebergy balance at " + Site_ID + " for "+str(plotyear)
	    pass
	    
def regressionLasslop(mypathforResults,New_combined,Site_ID,startdate,enddate,freq_list,list_in):    
    
    for plot_freq in freq_list: 
	for var_to_plot in list_in:
	    print "Doing Fre plots for "+Site_ID + " at freq " + plot_freq
		
	    std_variable=var_to_plot+"_Con"
	    Lasslop_variable=var_to_plot+"_Lasslop"
	    if var_to_plot=='Fc': std_variable='Fc_ustar'
	    
	    by = lambda x: lambda y: getattr(y, x)
	    tempdata1=New_combined[[std_variable,Lasslop_variable ]]
	    tempdata1=tempdata1.dropna(how='any')
	    tempdata=tempdata1.groupby([by('year'),by(plot_freq)]).mean()	     
	
	    #First plot use ALL data and years
	    ###################################
	    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	    #Do first plot - All hours
	    #============================
	    #Get data for all hours - so dont select any
	    xdata=tempdata[std_variable]
	    ydata=tempdata[Lasslop_variable]
	    
	    #Convert all units to g C m-2 period-1
	    if plot_freq=='dayofyear': ydata=tempdata[Lasslop_variable]*60*60*24/1000*12/44
	    if plot_freq=='dayofyear': xdata=tempdata[std_variable]*60*60*24/1000*12/44
	    if plot_freq=='week': ydata=tempdata[Lasslop_variable]*60*60*24/1000*12/44   *7
	    if plot_freq=='week': xdata= tempdata[std_variable]*60*60*24/1000*12/44   *7
	    if plot_freq=='month': ydata=tempdata[Lasslop_variable]*60*60*24/1000*12/44  *30
	    if plot_freq=='month': xdata=tempdata[std_variable]*60*60*24/1000*12/44     *30
	    
	    ax1=pl.subplot(1,1,1)
	    pl.title=('Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID )
	    
	    #Define colours
	    if var_to_plot=='Fre': plotcolour='#990033'
	    if var_to_plot=='GPP': plotcolour='#66CC00'
	    if var_to_plot=='Fc': plotcolour='#CC9900'
	    
	    ax1.plot(xdata,ydata, 'o', color=plotcolour) 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	    num_points=len(ydata)
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err)) +'\n' +
	                   'n points   ' + str("{0:.2f}".format(num_points))       )  
	    pl.figtext(0.7,0.20,graphtext1, bbox=dict(),size=16)
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax1.plot(x, y, linewidth = 2)    
	    
	    #plot linear regression 1:1
	    x2 = np.linspace(min(xdata),max(xdata))
	    y2 = 1 * x2 + 0
	    ax1.plot(x2, y2, linewidth = 1,label='1:1',color='k',linestyle='--') 
	    
	    print "Min X  " + str(min(xdata)) + "            Max X  " + str(max(xdata))
	    
	    
	    ax1.legend(loc='upper left')
	    pl.xlabel(std_variable + ' (g C m-2 '+plot_freq+'-1)',size=16)
	    pl.ylabel(Lasslop_variable +' (g C m-2 '+plot_freq+'-1)',size=16)
	    
	    pl.suptitle('Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID ,size=20)  
	    
	    pl.savefig(mypathforResults+'/Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID+'.png') 
	    #pl.show()
	    pl.close()    
	    print 'Closed Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID

def regressionFre(mypathforResults,New_combined,Site_ID,startdate,enddate,freq_list,list_in):    
    
    for plot_freq in freq_list: 
	for var_to_plot in list_in:
	    print "Doing Fre plots for "+Site_ID + " at freq " + plot_freq
		
	    std_variable=var_to_plot
	    Lasslop_variable="Fre_Lasslop"
	    
	    by = lambda x: lambda y: getattr(y, x)
	    tempdata1=New_combined[[std_variable,Lasslop_variable ]]
	    tempdata1=tempdata1.dropna(how='any')
	    tempdata=tempdata1.groupby([by('year'),by(plot_freq)]).mean()	     
	
	    #First plot use ALL data and years
	    ###################################
	    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	    #Do first plot - All hours
	    #============================
	    #Get data for all hours - so dont select any
	    xdata=tempdata[std_variable]
	    ydata=tempdata[Lasslop_variable]
	    
	    #Convert all units to g C m-2 period-1
	    if plot_freq=='dayofyear': ydata=tempdata[Lasslop_variable]*60*60*24/1000*12/44
	    if plot_freq=='dayofyear': xdata=tempdata[std_variable]*60*60*24/1000*12/44
	    if plot_freq=='week': ydata=tempdata[Lasslop_variable]*60*60*24/1000*12/44   *7
	    if plot_freq=='week': xdata= tempdata[std_variable]*60*60*24/1000*12/44   *7
	    if plot_freq=='month': ydata=tempdata[Lasslop_variable]*60*60*24/1000*12/44  *30
	    if plot_freq=='month': xdata=tempdata[std_variable]*60*60*24/1000*12/44     *30
	    
	    ax1=pl.subplot(1,1,1)
	    pl.title=('Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID )
	    
	    #Define colours
	    if var_to_plot=='Fre_noct': plotcolour='#FF9966'
	    if var_to_plot=='Fre_NN': plotcolour='#FF6666'
	    if var_to_plot=='Fre_Con': plotcolour='#FF3366'
	    
	    ax1.plot(xdata,ydata, 'o', color=plotcolour) 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	    num_points=len(ydata)
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err)) +'\n' +
	                   'n points   ' + str("{0:.2f}".format(num_points))       )  
	    pl.figtext(0.7,0.20,graphtext1, bbox=dict(),size=16)
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax1.plot(x, y, linewidth = 2)    
	    
	    #plot linear regression 1:1
	    x2 = np.linspace(min(xdata),max(xdata))
	    y2 = 1 * x2 + 0
	    ax1.plot(x2, y2, linewidth = 1,label='1:1',color='k',linestyle='--') 
	    
	    print "Min X  " + str(min(xdata)) + "            Max X  " + str(max(xdata))
	    
	    
	    ax1.legend(loc='upper left')
	    pl.xlabel(std_variable + ' (g C m-2 '+plot_freq+'-1)',size=16)
	    pl.ylabel(Lasslop_variable +' (g C m-2 '+plot_freq+'-1)',size=16)
	    
	    pl.suptitle('Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID ,size=20)  
	    
	    pl.savefig(mypathforResults+'/Plot2 '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID+'.png') 
	    #pl.show()
	    pl.close()    
	    print 'Closed Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID
 
	
def cummulative_CO2_H2O(mypathforResults,New_combined,Site_ID,startdate,enddate,Rain_Con_label_variable_to_fill):    
    print "Doing Cummulative CO2 and H2O plots for "+Site_ID
 
    #First plot use ALL data and years
    ###################################
    #Create the Cummulative variables
    New_combined['Fc_Cumm']=New_combined.groupby([lambda x: x.year])['Fc_ustar'].cumsum()
    New_combined['Fc_Cumm']=New_combined['Fc_Cumm']
    New_combined['Fe_Cumm']=New_combined.groupby([lambda x: x.year])['Fe_Con'].cumsum()  
    New_combined['Precip_Cumm']=New_combined.groupby([lambda x: x.year])[Rain_Con_label_variable_to_fill].cumsum()
    
    #get to units of g.CO2.-m-2. multiply by 44 (MW of CO2) to get g CO2
    New_combined['Fc_Cumm']=New_combined['Fc_Cumm']*60*30*10000/1000000*44.01*12/44/1000000
    #Get from W.m-2 to mm
    New_combined['Fe_Cumm']=New_combined['Fe_Cumm']*60*30/1000000/2.45	#Get year for plot lables    

    #Group all data by Year
    tempdata_grouped=New_combined.groupby([lambda x: x.year])

    Num_plots=len(tempdata_grouped)
    cm = get_cmap('Dark2')
    #First plot Fc
   
    fig = figure(1,figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    ax  = fig.add_subplot(111)
    
    index=0
    for name, group in tempdata_grouped:
	plotyear=group.index[0].year
	index=index+1
	xdata=group['Fc_Cumm']
	plotyear=group.index[0].year
	color_rgb = cm(1.*index/Num_plots)
	ax12=pl.plot(xdata, '-', linewidth=2,color=color_rgb,label='Fc ustar '+str(plotyear))

    index=0
    for name, group in tempdata_grouped:
	plotyear=group.index[0].year
	index=index+1
	xdata=group['Fc_Cumm']
	plotyear=group.index[0].year
	color_rgb = cm(1.*index/Num_plots)
	#ax12=pl.plot(xdata, '-', linewidth=2,color=color_rgb,label='Fc '+str(plotyear))
	x_coord_Fc=xdata[-1]
	y_coord_Fc=len(xdata)-1    
	ax.annotate('test', xy=(x_coord_Fc, y_coord_Fc),  xycoords='data',
	            xytext=(-100, -100), textcoords='offset points',
	            arrowprops=dict(arrowstyle="->")
	            )	
	
    # Set the ticks and labels...
    #ticks = np.linspace(0, len(xdata), 12)
    #labels = range(ticks.size)
    #pl.xticks(ticks, labels)  
    
    pl.legend(loc='upper left')
    pl.ylabel('Cummulative Carbon flux (t C ha-1 y-1)')
    pl.xlabel('Time (30 min intervals)')
    pl.suptitle('Cummulative CO2 ustar plot for '+Site_ID,size=20) 
    pl.savefig(mypathforResults+'/'+'Cummulative CO2 plot for '+Site_ID) 
    #pl.show()
    pl.close()    

    cm = get_cmap('winter')
    #Second plot Fe
    ax1=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    index=0
    for name, group in tempdata_grouped:
	plotyear=group.index[0].year
	index=index+1
	xdata1=group['Fe_Cumm']
	xdata2=group['Precip_Cumm']
	plotyear=group.index[0].year  
	color_rgb = cm(1.*index/Num_plots)
	ax1=pl.plot(xdata1, '-', linewidth=2,color=color_rgb,label='ET '+str(plotyear))
	ax1=pl.plot(xdata2, '--', linewidth=2,color=color_rgb,label='Precip '+str(plotyear))

	index=0
	for name, group in tempdata_grouped:
	    plotyear=group.index[0].year
	    index=index+1
	    xdata1=group['Fe_Cumm']
	    xdata2=group['Precip_Cumm']
	    plotyear=group.index[0].year  
	    color_rgb = cm(1.*index/Num_plots)
	    #ax1=pl.plot(xdata1, '-', linewidth=2,color=color_rgb,label='ET '+str(plotyear))
	    #ax1=pl.plot(xdata2, '--', linewidth=2,color=color_rgb,label='Precip '+str(plotyear))
	    x_coord_et=xdata1[-1]
	    y_coord_et=len(xdata1)-1
	    x_coord_Precip=xdata2[-1]
	    y_coord_Precip=len(xdata2)-1
	    pl.annotate(name, xy=(x_coord_et, y_coord_et),  xycoords='data',
		        xytext=(50, 0), textcoords='offset points',
		        arrowprops=dict(arrowstyle="->")
		        )	
	    pl.annotate(name, xy=(x_coord_Precip, y_coord_Precip),  xycoords='data',
		        xytext=(50, 0), textcoords='offset points',
		        arrowprops=dict(arrowstyle="->")
		        )	

    pl.legend(loc='upper left')
    pl.ylabel('Cummulative water flux (mm)')
    pl.xlabel('Time (30 minute periods)')
    pl.suptitle('Cummulative H2O plot for '+Site_ID,size=20) 
    pl.savefig(mypathforResults+'/'+'Cummulative H2O plot for '+Site_ID) 
    #pl.show()
    pl.close()   


def mintimeseries_plot(mypathforResults,predicted,observed,regress,variable_to_fill, Site_ID,units,targets,output,item):    
    ANN_label=str(item+"_NN")    
    pl.plot( targets, 'b--' )
    pl.plot( output, 'k-' )
    pl.legend(('targets', 'output'))
    pl.xlabel('Time'); 
    pl.title('Outputs vs. target of trained network for '+item)
    pl.grid(True)
    pl.legend()
    pl.title('Tower vs ANN 30 min timeseries for '+item+' at ' +Site_ID)
    pl.ylabel(item + '('+units+')')
    pl.savefig(mypathforResults+'/'+'Tower vs ANN 30 min timeseries for '+item+' at ' +Site_ID) 
    #pl.show()
    pl.close()
    


def count_numbers(frame):
    return (frame.count())

def count_total(frame):
    return (frame.shape[0])


def do_nan_stats(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label):

    for item in list_in:
	print "Doing Stats for Nans and pct filled"
	
	#Define the variable (item) for counts. 
	item_Con=item+'_Con'
	if item in ['Fc','Fe','Fh','Fg']:
	    item_Corr=item+'_NN'
	else:
	    item_Corr=item+'_Corr'
	
	#Do stats and counts for Sws
	Numbers_by_month= New_combined[[item,item_Con,item_Corr]].groupby([lambda x: x.year,lambda x: x.month]).apply(count_numbers)
	Numbers_by_month['total_n'] = New_combined[item_Con].groupby([lambda x: x.year,lambda x: x.month]).apply(count_total)  
	#Numbers_by_month['Pct_nan'] = (float(Numbers_by_month['Ta_EC'])/float(Numbers_by_month['total_n']))*100
	Numbers_by_month['Pct_nan']=Numbers_by_month.apply(lambda x: int(((np.float(x['total_n'])-np.float(x[item]))/np.float(x['total_n']))*100), axis=1)
	Numbers_by_month['Pct_notfilled']=Numbers_by_month.apply(lambda x: int(((np.float(x['total_n'])-np.float(x[item_Con]))/np.float(x['total_n']))*100), axis=1)
	Numbers_by_month.astype('int32')
	
	#Write out file
	Numbers_by_month.to_csv(mypathforResults+'/'+'Nan counts and Pct filled for '+item+' at ' +Site_ID+'.csv')	
	Numbers_by_month.to_pickle(mypathforResults+'/'+'Nan counts and Pct filled for '+item+' at ' +Site_ID+'.df') 	
    
###########################################################################################################
##                 START MAIN CODE
###########################################################################################################
def basic_diags(myBaseforResults,New_combined,Site_id,list_in,Ws_label,do_results,Rain_label_variable_to_fill):     
    global Site_ID
    Site_ID=Site_id
    
    print "Starting Diagnostics"
    Rain_Con_label_variable_to_fill= Rain_label_variable_to_fill+'_Con'
    
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
	os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/Diagnostics"):
	os.mkdir(myBaseforResults+"/Diagnostics")
    mypathforResults=myBaseforResults+"/Diagnostics"  
    
    number_of_inputs=len(list_in)
    
    startdate=New_combined.index[0]
    enddate=New_combined.index[len(New_combined)-1]
    
    list_in=['Fc','Fe','Fh','Fg','Ta','Ah',Rain_label_variable_to_fill,Ws_label,'Ts','Sws','Fsd','Fsu','Fld','Flu','Fn'] 
    do_nan_stats(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label)
    
    ###############################
    #Do plots
    ###############################
    #Plots at specified fequencies.  Pass freq to function.  Can be 'day of year', 'week', 'month', 'year'
    #Plot timeseries of daily over all periods
    plot_freq = 'dayofyear'
    
    list_in=['Fc','Fe','Fh','Fg']
    Doplots_at_daily(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)
    list_in=['Fsu','Fsd','Flu','Fld']
    Doplots_at_daily(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)
    list_in=['Ta','Ah',Ws_label,'VPD']
    Doplots_at_daily(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)    
    list_in=['Ts','Sws',Rain_label_variable_to_fill]
    Doplots_at_daily(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)       
    #Note here that there are NEW variable names that will be dealt with in the function
    list_in=['BR','WUE','RUE','EBC','EF']
    Doplots_at_daily_other(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq) 

    plot_freq = 'week'
    
    list_in=['Fc','Fe','Fh','Fg']
    Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)
    list_in=['Fsu','Fsd']
    Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)
    list_in=['Ta','Ah',Ws_label]
    Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)    
    list_in=['Ts','Sws',Rain_label_variable_to_fill]
    Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)   

    plot_freq = 'month'
    
    #list_in=['Fc','Fe','Fh','Fg']
    #Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)
    #list_in=['Fsu','Fsd']
    #Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)
    #list_in=['Ta','Ah',Ws_label]
    #Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)    
    #list_in=['Ts','Sws',Rain_label_variable_to_fill]
    #Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq)       
    
    #Plot timeseries of monthly over all periods
    list_in=['Ta','Ah',Ws_label]
    Doplots_monthly_diff(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label)
    list_in=['Fc','Fe','Fh','Fg']
    Doplots_monthly_diff(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label)
    list_in=['Fsd','Fsu','Fld','Flu','Fn']
    Doplots_monthly_diff(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label)    
    #Plot Nans and Pct data missing
    list_in=['Ta','Ah',Ws_label]	
    Plot_Pct_nans(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in)
    
    list_in=['Fc_ustar','Fe','Fh','Fg']	
    Plot_Pct_nans(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in)

    list_in=['Fsd','Fsu','Fld','Flu','Fn']	
    Plot_Pct_nans(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in)    
    
    list_in=['Ts','Sws',Rain_label_variable_to_fill]	
    Plot_Pct_nans(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in)

    
    #Plots of EB closure
    ####################################################################################
    regressionEBclosure(mypathforResults,New_combined,Site_ID,startdate,enddate)
    
    # Do plots of Fre (Reichstein ustar approach vs Lasslop daytime Fre approach
    ####################################################################################
    #Define list of plots to make to process
    freq_list=["dayofyear", "week", "month"]
    list_in=["Fre","GPP","Fc"]	
    regressionLasslop(mypathforResults,New_combined,Site_ID,startdate,enddate,freq_list,list_in)    

    #Define list of plots to make to process
    freq_list=["dayofyear", "week", "month"]
    list_in=["Fre_noct","Fre_NN","Fre_Con"]	
    regressionFre(mypathforResults,New_combined,Site_ID,startdate,enddate,freq_list,list_in)   

    ##Plot time series of all 30 minute data
    #mintimeseries_plot(mypathforResults,predicted,observed,regress,item, Site_ID,units,targets,output,ANN_label)
    ##Plot regression of Tower versus ANN
    #regressionANN(mypathforResults,predicted,observed,regress,item, Site_ID,units,ANN_label)
    ##Plot diurnals for every second month 6 graphs
    #Doplots_diurnal(mypathforResults,New_combined,item, Site_ID,units,ANN_label)
    ##Plot timeseries of monthly over all periods
    #Doplots_monthly(mypathforResults,New_combined,item, Site_ID,units,ANN_label)
    
    #Do results plots
    if do_results==True:
	print "Starting Results"
	
	#Check for place to put results - does it exist? If not create
	if not os.path.isdir(myBaseforResults):
	    os.mkdir(myBaseforResults)
	#Then subdirectories
	if not os.path.isdir(myBaseforResults+"/Results"):
	    os.mkdir(myBaseforResults+"/Results")
	mypathforResults=myBaseforResults+"/Results"  	#Create directory for results

	#Do cummulative water and carbon plots
	cummulative_CO2_H2O(mypathforResults,New_combined,Site_ID,startdate,enddate,Rain_Con_label_variable_to_fill)
	
	#Do time series CO2 (GPP,Re plots)
	list_in=['Fc_ustar','GPP_Con','Fre_Con']
	plot_freq = 'dayofyear'
	Doplots_at_daily_carbon_umol(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq) 	
	Doplots_at_daily_carbon_g(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq) 	

	#Do time series Energy Balnce (Fn,Fe,Fh,Fg plots)
	list_in=['Fn_Con','Fe_Con','Fh_Con','Fg_Con']
	plot_freq = 'dayofyear'	
	Doplots_at_daily_EB(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq) 
	

	
    ###################################################
    #  File stuff
    ###################################################
    
    print"Finished Diagnostics at "+ Site_ID    
    
    
    
    
    
