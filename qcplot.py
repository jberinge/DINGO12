import ast
import constants as c
import time
import matplotlib.dates as mdt
import matplotlib.pyplot as plt
import numpy
import qcutils
import logging

log = logging.getLogger('qc.plot')

def plottimeseries(cf,nFig,dsa,dsb,si,ei):
    SiteName = dsa.globalattributes['site_name']
    Level = dsb.globalattributes['nc_level']
    dt = int(dsa.globalattributes['time_step'])
    Month = dsa.series['Month']['Data'][0]
    p = plot_setup(cf,nFig)
    log.info(' Plotting series: '+str(p['SeriesList']))
    L1XArray = numpy.array(dsa.series['DateTime']['Data'][si:ei])
    L2XArray = numpy.array(dsb.series['DateTime']['Data'][si:ei])
    p['XAxMin'] = min(L2XArray)
    p['XAxMax'] = max(L2XArray)
    p['loc'],p['fmt'] = get_ticks(p['XAxMin'],p['XAxMax'])
    plt.ioff()
    fig = plt.figure(int(nFig),figsize=(p['PlotWidth'],p['PlotHeight']))
    fig.clf()
    plt.figtext(0.5,0.95,SiteName+': '+p['PlotDescription'],ha='center',size=16)
    for ThisOne, n in zip(p['SeriesList'],range(p['nGraphs'])):
        if ThisOne in dsa.series.keys():
            aflag = dsa.series[ThisOne]['Flag']
            p['Units'] = dsa.series[ThisOne]['Attr']['units']
            p['YAxOrg'] = p['ts_YAxOrg'] + n*p['yaxOrgOffset']
            L1YArray,p['nRecs'],p['nNotM'],p['nMskd'] = get_yarray(dsa,ThisOne,si=si,ei=ei)
            # check the control file to see if the Y axis minima have been specified
            nSer = p['SeriesList'].index(ThisOne)
            p['LYAxMax'],p['LYAxMin'] = get_yaxislimitsfromcf(cf,nFig,'YLMax','YLMin',nSer,L1YArray)
            plot_onetimeseries_left(fig,n,ThisOne,L1XArray,L1YArray,p)
        if ThisOne in dsb.series.keys():
            bflag = dsb.series[ThisOne]['Flag']
            p['Units'] = dsb.series[ThisOne]['Attr']['units']
            p['YAxOrg'] = p['ts_YAxOrg'] + n*p['yaxOrgOffset']
            #Plot the Level 2 data series on the same X axis but with the scale on the right Y axis.
            L2YArray,p['nRecs'],p['nNotM'],p['nMskd'] = get_yarray(dsb,ThisOne,si=si,ei=ei)
            # check the control file to see if the Y axis minima have been specified
            nSer = p['SeriesList'].index(ThisOne)
            p['RYAxMax'],p['RYAxMin'] = get_yaxislimitsfromcf(cf,nFig,'YRMax','YRMin',nSer,L2YArray)
            plot_onetimeseries_right(fig,n,ThisOne,L2XArray,L2YArray,p)

            #Plot the diurnal averages.
            Hr2,Av2,Sd2,Mx2,Mn2=get_diurnalstats(dsb.series['Hdh']['Data'][si:ei],
                                                dsb.series[ThisOne]['Data'][si:ei],dt)
            Av2 = numpy.ma.masked_where(Av2==-9999,Av2)
            Sd2 = numpy.ma.masked_where(Sd2==-9999,Sd2)
            Mx2 = numpy.ma.masked_where(Mx2==-9999,Mx2)
            Mn2 = numpy.ma.masked_where(Mn2==-9999,Mn2)
            hr2_ax = fig.add_axes([p['hr1_XAxOrg'],p['YAxOrg'],p['hr2_XAxLen'],p['ts_YAxLen']])
            hr2_ax.hold(True)
            hr2_ax.plot(Hr2,Av2,'y-',Hr2,Mx2,'r-',Hr2,Mn2,'b-')
            section = qcutils.get_cfsection(cf,series=ThisOne,mode='quiet')
            if len(section)!=0:
                if 'DiurnalCheck' in cf[section][ThisOne].keys():
                    NSdarr = numpy.array(eval(cf[section][ThisOne]['DiurnalCheck']['NumSd']),dtype=float)
                    nSd = NSdarr[Month-1]
                    hr2_ax.plot(Hr2,Av2+nSd*Sd2,'r.',Hr2,Av2-nSd*Sd2,'b.')
            plt.xlim(0,24)
            plt.xticks([0,6,12,18,24])
            if n==0:
                hr2_ax.set_xlabel('Hour',visible=True)
            else:
                hr2_ax.set_xlabel('',visible=False)
                plt.setp(hr2_ax.get_xticklabels(), visible=False)
            #if n > 0: plt.setp(hr2_ax.get_xticklabels(), visible=False)

            # vertical lines to show frequency distribution of flags
            bins = numpy.arange(0.5,23.5)
            ind = bins[:len(bins)-1]+0.5
            index = numpy.where(numpy.mod(bflag,10)==0)    # find the elements with flag = 0, 10, 20 etc
            bflag[index] = 0                               # set them all to 0
            hist, bin_edges = numpy.histogram(bflag, bins=bins)
            ymin = hist*0
            delta = 0.01*(numpy.max(hist)-numpy.min(hist))
            bar_ax = fig.add_axes([p['hr2_XAxOrg'],p['YAxOrg'],p['bar_XAxLen'],p['ts_YAxLen']])
            bar_ax.set_ylim(0,numpy.max(hist))
            bar_ax.vlines(ind,ymin,hist)
            for i,j in zip(ind,hist):
                if j>0.05*numpy.max(hist): bar_ax.text(i,j+delta,str(int(i)),ha='center',size='small')
            if n==0:
                bar_ax.set_xlabel('Flag',visible=True)
            else:
                bar_ax.set_xlabel('',visible=False)
                plt.setp(bar_ax.get_xticklabels(), visible=False)
            #if n > 0: plt.setp(bar_ax.get_xticklabels(), visible=False)
        else:
            log.error('  plttimeseries: series '+ThisOne+' not in data structure')
    fig.show()
    fname = 'plots/'+SiteName.replace(' ','')+'_'+Level+'_'+p['PlotDescription'].replace(' ','')+'.png'
    fig.savefig(fname,format='png')

def plot_setup(cf,nFig):
    p = {}
    p['PlotDescription'] = cf['Plots'][str(nFig)]['Title']
    p['SeriesList'] = ast.literal_eval(cf['Plots'][str(nFig)]['Variables'])
    p['nGraphs'] = len(p['SeriesList'])
    p['PlotWidth'] = 13
    p['PlotHeight'] = 9
    p['ts_YAxOrg'] = 0.08
    p['ts_XAxOrg'] = 0.06
    p['ts_XAxLen'] = 0.6
    p['hr_XAxLen'] = 0.1
    p['ts_YAxLen'] = (0.85 - (p['nGraphs'] - 1)*0.02)/p['nGraphs']
    p['yaxOrgOffset'] = (0.85 - p['ts_YAxLen'])/(p['nGraphs'] - 1)
    p['hr1_XAxOrg'] = p['ts_XAxOrg']+p['ts_XAxLen']+0.07
    p['hr1_XAxLen'] = p['hr_XAxLen']
    p['hr2_XAxOrg'] = p['hr1_XAxOrg']+p['hr1_XAxLen']+0.05
    p['hr2_XAxLen'] = p['hr_XAxLen']
    p['bar_XAxOrg'] = p['hr1_XAxOrg']+p['hr1_XAxLen']+0.05+p['hr1_XAxLen']+0.05
    p['bar_XAxLen'] = p['hr_XAxLen']
    return p

def plot_onetimeseries_left(fig,n,ThisOne,xarray,yarray,p):
    ts_ax = fig.add_axes([p['ts_XAxOrg'],p['YAxOrg'],p['ts_XAxLen'],p['ts_YAxLen']])
    ts_ax.hold(False)
    p['ts_ax'] = ts_ax
    ts_ax.plot(xarray,yarray,'b-')
    ts_ax.xaxis.set_major_locator(p['loc'])
    ts_ax.xaxis.set_major_formatter(p['fmt'])
    ts_ax.set_xlim(p['XAxMin'],p['XAxMax'])
    ts_ax.set_ylim(p['LYAxMin'],p['LYAxMax'])
    if n==0:
        ts_ax.set_xlabel('Date',visible=True)
    else:
        ts_ax.set_xlabel('',visible=False)
    TextStr = ThisOne+'('+p['Units']+')'+str(p['nRecs'])+' '+str(p['nNotM'])+' '+str(p['nMskd'])
    txtXLoc = p['ts_XAxOrg']+0.01
    txtYLoc = p['YAxOrg']+p['ts_YAxLen']-0.025
    plt.figtext(txtXLoc,txtYLoc,TextStr,color='b',horizontalalignment='left')
    if n > 0: plt.setp(ts_ax.get_xticklabels(),visible=False)

def plot_onetimeseries_right(fig,n,ThisOne,xarray,yarray,p):
    if not p.has_key('ts_ax'):
        ts_ax = fig.add_axes([p['ts_XAxOrg'],p['YAxOrg'],p['ts_XAxLen'],p['ts_YAxLen']])
        ts_ax.hold(False)
        ts_ax.yaxis.tick_right()
        TextStr = ThisOne+'('+p['Units']+')'
        txtXLoc = p['ts_XAxOrg']+0.01
        txtYLoc = p['YAxOrg']+p['ts_YAxLen']-0.025
        plt.figtext(txtXLoc,txtYLoc,TextStr,color='b',horizontalalignment='left')
    else:
        ts_ax = p['ts_ax'].twinx()
    colour = 'r'
    if p.has_key('ts_ax'): del p['ts_ax']
    ts_ax.plot(xarray,yarray,'r-')
    ts_ax.xaxis.set_major_locator(p['loc'])
    ts_ax.xaxis.set_major_formatter(p['fmt'])
    ts_ax.set_xlim(p['XAxMin'],p['XAxMax'])
    ts_ax.set_ylim(p['RYAxMin'],p['RYAxMax'])
    if n==0:
        ts_ax.set_xlabel('Date',visible=True)
    else:
        ts_ax.set_xlabel('',visible=False)
    TextStr = str(p['nNotM'])+' '+str(p['nMskd'])
    txtXLoc = p['ts_XAxOrg']+p['ts_XAxLen']-0.01
    txtYLoc = p['YAxOrg']+p['ts_YAxLen']-0.025
    plt.figtext(txtXLoc,txtYLoc,TextStr,color='r',horizontalalignment='right')
    if n > 0: plt.setp(ts_ax.get_xticklabels(),visible=False)

def get_yarray(ds,ThisOne,si=0,ei=-1):
    yarray = numpy.ma.masked_where(abs(ds.series[ThisOne]['Data'][si:ei]-float(-9999))<c.eps,
                                        ds.series[ThisOne]['Data'][si:ei])
    nRecs = numpy.ma.size(yarray)
    nNotM = numpy.ma.count(yarray)
    nMskd = numpy.ma.count_masked(yarray)
    if numpy.ma.count(yarray)==0:
        yarray = numpy.ma.zeros(numpy.size(yarray))
    return yarray,nRecs,nNotM,nMskd
    
def get_yaxislimitsfromcf(cf,nFig,maxkey,minkey,nSer,YArray):
    if maxkey in cf['Plots'][str(nFig)].keys():                               # Y axis minima specified
        maxlist = ast.literal_eval(cf['Plots'][str(nFig)][maxkey])     # Evaluate the minima list
        if str(maxlist[nSer])=='Auto':             # This entry is 'Auto' ...
            YAxMax = numpy.ma.maximum(YArray)                        # ... so take the array minimum value
        else:
            YAxMax = float(maxlist[nSer])         # Evaluate the entry for this series
    else:
        YAxMax = numpy.ma.maximum(YArray)                            # Y axis minima not given, use auto
    if minkey in cf['Plots'][str(nFig)].keys():                               # Y axis minima specified
        minlist = ast.literal_eval(cf['Plots'][str(nFig)][minkey])     # Evaluate the minima list
        if str(minlist[nSer])=='Auto':             # This entry is 'Auto' ...
            YAxMin = numpy.ma.minimum(YArray)                        # ... so take the array minimum value
        else:
            YAxMin = float(minlist[nSer])         # Evaluate the entry for this series
    else:
        YAxMin = numpy.ma.minimum(YArray)                            # Y axis minima not given, use auto
    return YAxMax,YAxMin

def plotxy(cf,nFig,plt_cf,dsa,dsb,si,ei):
    SiteName = dsa.globalattributes['site_name']
    PlotDescription = cf['Plots'][str(nFig)]['Title']
    fig = plt.figure(int(nFig))
    
    fig.clf()
    plt.figtext(0.5,0.95,SiteName+': '+PlotDescription,ha='center',size=16)
    XSeries = ast.literal_eval(plt_cf['XSeries'])
    YSeries = ast.literal_eval(plt_cf['YSeries'])
    log.info(' Plotting xy: '+str(XSeries)+' v '+str(YSeries))
    if dsa == dsb:
        for xname,yname in zip(XSeries,YSeries):
            xa,flag = qcutils.GetSeriesasMA(dsa,xname,si=si,ei=ei)
            ya,flag = qcutils.GetSeriesasMA(dsa,yname,si=si,ei=ei)
            xyplot(xa,ya,sub=[1,1,1],regr=1,xlabel=xname,ylabel=yname)
    else:
        for xname,yname in zip(XSeries,YSeries):
            xa,flag = qcutils.GetSeriesasMA(dsa,xname,si=si,ei=ei)
            ya,flag = qcutils.GetSeriesasMA(dsa,yname,si=si,ei=ei)
            xb,flag = qcutils.GetSeriesasMA(dsb,xname,si=si,ei=ei)
            yb,flag = qcutils.GetSeriesasMA(dsb,yname,si=si,ei=ei)
            xyplot(xa,ya,sub=[1,2,1],xlabel=xname,ylabel=yname)
            xyplot(xb,yb,sub=[1,2,2],regr=1,xlabel=xname,ylabel=yname)
    fig.show()

def xyplot(x,y,sub=[1,1,1],regr=0,title=None,xlabel=None,ylabel=None):
    '''Generic XY scatter plot routine'''
    plt.subplot(sub[0],sub[1],sub[2])
    plt.plot(x,y,'b.')
    ax = plt.gca()
    if xlabel!=None: plt.xlabel(xlabel)
    if ylabel!=None: plt.ylabel(ylabel)
    if title!=None: plt.title(title)
    if regr!=0:
        coefs = numpy.ma.polyfit(x,y,1)
        xfit = numpy.ma.array([numpy.ma.minimum(x),numpy.ma.maximum(x)])
        yfit = numpy.polyval(coefs,xfit)
        r = numpy.ma.corrcoef(x,y)
        eqnstr = 'y = %.3fx + %.3f, r = %.3f'%(coefs[0],coefs[1],r[0][1])
        plt.plot(xfit,yfit,'r--',linewidth=3)
        plt.text(0.5,0.9,eqnstr,fontsize=8,horizontalalignment='center',transform=ax.transAxes)

def tsplot(x,y,sub=[1,1,1],title=None,xlabel=None,ylabel=None,colours=None,lineat=None):
    plt.subplot(sub[0],sub[1],sub[2])
    MTFmt = mdt.DateFormatter('%m/%Y')
    if (y.all() is numpy.ma.masked):
        y = numpy.ma.zeros(len(y))
    if colours!=None:
        plt.scatter(x,y,c=colours)
    else:
        plt.scatter(x,y)
    if lineat!=None:
        plt.plot((x[0],x[-1]),(float(lineat),float(lineat)))
    plt.xlim((x[0],x[-1]))
    ax = plt.gca()
    ax.xaxis.set_major_formatter(MTFmt)
    if title!=None:
        plt.title(title)
    if ylabel!=None:
        ax.yaxis.set_label_text(ylabel)
    if xlabel!=None:
        ax.xaxis.set_label_text(xlabel)

def get_diurnalstats(DecHour,Data,dt):
    nInts = 24*int((60/dt)+0.5)
    Hr = numpy.array([-9999]*nInts,dtype=numpy.float64)
    Av = numpy.array([-9999]*nInts,dtype=numpy.float64)
    Sd = numpy.array([-9999]*nInts,dtype=numpy.float64)
    Mx = numpy.array([-9999]*nInts,dtype=numpy.float64)
    Mn = numpy.array([-9999]*nInts,dtype=numpy.float64)
    for i in range(nInts):
        Hr[i] = float(i)*dt/60.
        li = numpy.where((abs(DecHour-Hr[i])<c.eps)&(abs(Data-float(-9999))>c.eps))
        if numpy.size(li)!=0:
            Av[i] = numpy.mean(Data[li])
            Sd[i] = numpy.std(Data[li])
            Mx[i] = numpy.max(Data[li])
            Mn[i] = numpy.min(Data[li])
    return Hr, Av, Sd, Mx, Mn

def get_ticks(start, end):
    from datetime import timedelta as td
    delta = end - start

    if delta <= td(minutes=10):
        loc = mdt.MinuteLocator()
        fmt = mdt.DateFormatter('%H:%M')
    elif delta <= td(minutes=30):
        loc = mdt.MinuteLocator(byminute=range(0,60,5))
        fmt = mdt.DateFormatter('%H:%M')
    elif delta <= td(hours=1):
        loc = mdt.MinuteLocator(byminute=range(0,60,15))
        fmt = mdt.DateFormatter('%H:%M')
    elif delta <= td(hours=6):
        loc = mdt.HourLocator()
        fmt = mdt.DateFormatter('%H:%M')
    elif delta <= td(days=1):
        loc = mdt.HourLocator(byhour=range(0,24,3))
        fmt = mdt.DateFormatter('%H:%M')
    elif delta <= td(days=3):
        loc = mdt.HourLocator(byhour=range(0,24,12))
        fmt = mdt.DateFormatter('%d/%m %H')
    elif delta <= td(weeks=2):
        loc = mdt.DayLocator()
        fmt = mdt.DateFormatter('%d/%m')
    elif delta <= td(weeks=12):
        loc = mdt.WeekdayLocator()
        fmt = mdt.DateFormatter('%d/%m')
    elif delta <= td(weeks=104):
        loc = mdt.MonthLocator()
        fmt = mdt.DateFormatter('%d/%m')
    elif delta <= td(weeks=208):
        loc = mdt.MonthLocator(interval=3)
        fmt = mdt.DateFormatter('%d/%m/%y')
    else:
        loc = mdt.MonthLocator(interval=6)
        fmt = mdt.DateFormatter('%d/%m/%y')
    return loc,fmt