import ast
import copy
import constants as c
import datetime
import numpy
import time
import qcts
import qcutils
import logging

log = logging.getLogger('qc.ck')

def cliptorange(data, lower, upper):
    data = rangecheckserieslower(data,lower)
    data = rangecheckseriesupper(data,upper)
    return data

def rangecheckserieslower(data,lower):
    if numpy.ma.isMA(data):
        data = numpy.ma.masked_where(data<lower,data)
    else:
        index = numpy.where((abs(data-numpy.float64(-9999))>c.eps)&(data<lower))[0]
        data[index] = numpy.float64(-9999)
    return data

def rangecheckseriesupper(data,upper):
    if numpy.ma.isMA(data):
        data = numpy.ma.masked_where(data>upper,data)
    else:
        index = numpy.where((abs(data-numpy.float64(-9999))>c.eps)&(data>upper))[0]
        data[index] = numpy.float64(-9999)
    return data

def CoordinateFluxGaps(cf,ds,Fc_in='Fc',Fe_in='Fe',Fh_in='Fh'):
    if not qcutils.cfoptionskey(cf,Key='CoordinateFluxGaps'): return
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='gapsvars'):
        vars = ast.literal_eval(cf['FunctionArgs']['gapsvars'])
        Fc_in = vars[0]
        Fe_in = vars[1]
        Fh_in = vars[2]
    Fc,f = qcutils.GetSeriesasMA(ds,Fc_in)
    Fe,f = qcutils.GetSeriesasMA(ds,Fe_in)
    Fh,f = qcutils.GetSeriesasMA(ds,Fh_in)
    index = numpy.ma.where((Fc.mask==True) | (Fe.mask==True) | (Fh.mask==True))[0]
    # the following for ... in loop is not necessary
    for i in range(len(index)):
        j = index[i]
        if Fc.mask[j]==False:
            Fc.mask[j]=True
            Fc[j] = numpy.float64(-9999)
            ds.series[Fc_in]['Flag'][j] = numpy.int32(19)
        if Fe.mask[j]==False:
            Fe.mask[j]=True
            Fe[j] = numpy.float64(-9999)
            ds.series[Fe_in]['Flag'][j] = numpy.int32(19)           
        if Fh.mask[j]==False:
            Fh.mask[j]=True
            Fh[j] = numpy.float64(-9999)
            ds.series[Fh_in]['Flag'][j] = numpy.int32(19)
    ds.series[Fc_in]['Data']=numpy.ma.filled(Fc,float(-9999))
    ds.series[Fe_in]['Data']=numpy.ma.filled(Fe,float(-9999))
    ds.series[Fh_in]['Data']=numpy.ma.filled(Fh,float(-9999))
    log.info(' Finished gap co-ordination')

def CreateNewSeries(cf,ds):
    '''Create a new series using the MergeSeries or AverageSeries instructions.'''
    log.info(' Checking for new series to create')
    for ThisOne in cf['Variables'].keys():
        if 'MergeSeries' in cf['Variables'][ThisOne].keys():
            qcts.MergeSeries(cf,ds,ThisOne,[0,10])
        if 'AverageSeries' in cf['Variables'][ThisOne].keys():
            qcts.AverageSeriesByElements(cf,ds,ThisOne)

def do_7500check(cf,ds):
    '''Rejects data values for series specified in LI75List for times when the Diag_7500
       flag is non-zero.  If the Diag_7500 flag is not present in the data structure passed
       to this routine, it is constructed from the QC flags of the series specified in
       LI75Lisat.  Additional checks are done for AGC_7500 (the LI-7500 AGC value),
       Ah_7500_Sd (standard deviation of absolute humidity) and Cc_7500_Sd (standard
       deviation of CO2 concentration).'''
    log.info(' Doing the 7500 check')
    LI75List = ['Ah_7500_Av','Cc_7500_Av','AhAh','CcCc','UzA','UxA','UyA','UzC','UxC','UyC']
    if 'Diag_7500' not in cf['Variables'].keys():
        ds.series[unicode('Diag_7500')] = {}
        nRecs = numpy.size(ds.series['xlDateTime']['Data'])
        ds.series['Diag_7500']['Flag'] = numpy.zeros(nRecs,dtype=numpy.int32)
        for ThisOne in ['Ah_7500_Av','Cc_7500_Av']:
            if ThisOne in ds.series.keys():
                index = numpy.where(ds.series[ThisOne]['Flag']!=0)[0]
                log.info(' do_7500check: ', ThisOne, ' rejected ',len(index))
                ds.series['Diag_7500']['Flag'] = ds.series['Diag_7500']['Flag'] + ds.series[ThisOne]['Flag']
    index = numpy.where(ds.series['Diag_7500']['Flag']!=0)
    log.info('  7500Check: Diag_7500 ' + str(numpy.size(index)))
    for ThisOne in ['AGC_7500','Ah_7500_Sd','Cc_7500_Sd','AhAh','CcCc']:
        if ThisOne in ds.series.keys():
            index = numpy.where(ds.series[ThisOne]['Flag']!=0)
            log.info('  7500Check: ' + ThisOne + ' ' + str(numpy.size(index)))
            ds.series['Diag_7500']['Flag'] = ds.series['Diag_7500']['Flag'] + ds.series[ThisOne]['Flag']
    index = numpy.where((ds.series['Diag_7500']['Flag']!=0))
    log.info('  7500Check: Total ' + str(numpy.size(index)))
    for ThisOne in LI75List:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(-9999)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(4)
        else:
            log.error('  qcck.do_7500check: series '+str(ThisOne)+' in LI75List not found in ds.series')
    if '7500Check' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',7500Check'

def CoordinateAh7500AndFcGaps(cf,ds,Fcvar='Fc'):
    '''Cleans up Ah_7500_Av based upon Fc gaps to for QA check on Ah_7500_Av v Ah_HMP.'''
    if not qcutils.cfoptionskey(cf,Key='CoordinateAh7500&FcGaps'): return
    log.info(' Doing the Ah_7500 check')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='AhcheckFc'):
        Fclist = ast.literal_eval(cf['FunctionArgs']['AhcheckFc'])
        Fcvar = Fclist[0]
    
    # index1  Index of bad Ah_7500_Av observations
    index1 = numpy.where((ds.series['Ah_7500_Av']['Flag']!=0) & (ds.series['Ah_7500_Av']['Flag']!=10))
    
    # index2  Index of bad Fc observations
    index2 = numpy.where((ds.series[Fcvar]['Flag']!=0) & (ds.series[Fcvar]['Flag']!=10))
    
    ds.series['Ah_7500_Av']['Data'][index2] = numpy.float64(-9999)
    ds.series['Ah_7500_Av']['Flag'][index2] = ds.series[Fcvar]['Flag'][index2]
    ds.series['Ah_7500_Av']['Flag'][index1] = ds.series['Ah_7500_Av']['Flag'][index1]
    if 'CoordinateAh7500AndFcGaps' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',CoordinateAh7500AndFcGaps'

def do_CSATcheck(cf,ds):
    '''Rejects data values for series specified in CSATList for times when the Diag_CSAT
       flag is non-zero.  If the Diag_CSAT flag is not present in the data structure passed
       to this routine, it is constructed from the QC flags of the series specified in
       CSATList.'''
    log.info(' Doing the CSAT check')
    if 'Wd_CSAT_Compass' in ds.series.keys():
        Wd = 'Wd_CSAT_Compass'
    else:
        Wd = 'Wd_CSAT'
    CSATList = ['Ux','Uy','Uz','Ws_CSAT',Wd,'Tv_CSAT',
                'UzT','UxT','UyT','UzA','UxA','UyA','UzC','UxC','UyC',
                'UxUz','UyUz','UxUy','UxUx','UyUy']
    if 'Diag_CSAT' not in cf['Variables'].keys():
        ds.series['Diag_CSAT']= {}
        nRecs = numpy.size(ds.series['xlDateTime']['Data'])
        ds.series['Diag_CSAT']['Flag'] = numpy.zeros(nRecs,dtype=numpy.int32)
        for ThisOne in ['Ux','Uy','Uz','Tv_CSAT']:
            if ThisOne in ds.series.keys():
                index = numpy.where(ds.series[ThisOne]['Flag']!=0)[0]
                log.info(' do_CSATcheck: ', ThisOne, ' rejected ',len(index))
                ds.series['Diag_CSAT']['Flag'] = ds.series['Diag_CSAT']['Flag'] + ds.series[ThisOne]['Flag']
    index = numpy.where(ds.series['Diag_CSAT']['Flag']!=0)
    log.info('  CSATCheck: Diag_CSAT ' + str(numpy.size(index)))
    for ThisOne in CSATList:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(-9999)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(3)
        else:
            log.error('  qcck.do_CSATcheck: series '+str(ThisOne)+' in CSATList not found in ds.series')
    if 'CSATCheck' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',CSATCheck'

def do_diurnalcheck(cf,ds,section='',series='',code=5):
    if len(section)==0: return
    if len(series)==0: return
    if 'DiurnalCheck' not in cf[section][series].keys(): return
    if 'NumSd' not in cf[section][series]['DiurnalCheck'].keys(): return
    dt = float(ds.globalattributes['time_step'])
    n = int((60./dt) + 0.5)             #Number of timesteps per hour
    nInts = int((1440.0/dt)+0.5)        #Number of timesteps per day
    Av = numpy.array([-9999]*nInts,dtype=numpy.float64)
    Sd = numpy.array([-9999]*nInts,dtype=numpy.float64)
    NSd = numpy.array(eval(cf[section][series]['DiurnalCheck']['NumSd']),dtype=float)
    for m in range(1,13):
        mindex = numpy.where(ds.series['Month']['Data']==m)[0]
        if len(mindex)!=0:
            lHdh = ds.series['Hdh']['Data'][mindex]
            l2ds = ds.series[series]['Data'][mindex]
            for i in range(nInts):
                li = numpy.where(abs(lHdh-(float(i)/float(n))<c.eps)&(l2ds!=float(-9999)))
                if numpy.size(li)!=0:
                    Av[i] = numpy.mean(l2ds[li])
                    Sd[i] = numpy.std(l2ds[li])
                else:
                    Av[i] = float(-9999)
                    Sd[i] = float(-9999)
            Lwr = Av - NSd[m-1]*Sd
            Upr = Av + NSd[m-1]*Sd
            hindex = numpy.array(n*lHdh,int)
            index = numpy.where(((l2ds!=float(-9999))&(l2ds<Lwr[hindex]))|
                                ((l2ds!=float(-9999))&(l2ds>Upr[hindex])))[0] + mindex[0]
            ds.series[series]['Data'][index] = numpy.float64(-9999)
            ds.series[series]['Flag'][index] = numpy.int32(code)
            ds.series[series]['Attr']['diurnalcheck_numsd'] = cf[section][series]['DiurnalCheck']['NumSd']
    if 'DiurnalCheck' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',DiurnalCheck'

def do_excludedates(cf,ds,section='',series='',code=6):
    if len(section)==0: return
    if len(series)==0: return
    if 'ExcludeDates' not in cf[section][series].keys(): return
    ldt = ds.series['DateTime']['Data']
    ExcludeList = cf[section][series]['ExcludeDates'].keys()
    NumExclude = len(ExcludeList)
    for i in range(NumExclude):
        ExcludeDateList = ast.literal_eval(cf[section][series]['ExcludeDates'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ExcludeDateList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ExcludeDateList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        ds.series[series]['Data'][si:ei] = numpy.float64(-9999)
        ds.series[series]['Flag'][si:ei] = numpy.int32(code)
        ds.series[series]['Attr']['ExcludeDates_'+str(i)] = cf[section][series]['ExcludeDates'][str(i)]
    if 'ExcludeDates' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',ExcludeDates'

def do_excludehours(cf,ds,section='',series='',code=7):
    if len(section)==0: return
    if len(series)==0: return
    if 'ExcludeHours' not in cf[section][series].keys(): return
    ldt = ds.series['DateTime']['Data']
    ExcludeList = cf[section][series]['ExcludeHours'].keys()
    NumExclude = len(ExcludeList)
    for i in range(NumExclude):
        ExcludeHourList = ast.literal_eval(cf[section][series]['ExcludeHours'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ExcludeHourList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ExcludeHourList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        for j in range(len(ExcludeHourList[2])):
            ExHr = datetime.datetime.strptime(ExcludeHourList[2][j],'%H:%M').hour
            ExMn = datetime.datetime.strptime(ExcludeHourList[2][j],'%H:%M').minute
            index = numpy.where((ds.series['Hour']['Data'][si:ei]==ExHr)&
                                (ds.series['Minute']['Data'][si:ei]==ExMn))[0] + si
            ds.series[series]['Data'][index] = numpy.float64(-9999)
            ds.series[series]['Flag'][index] = numpy.int32(code)
            ds.series[series]['Attr']['ExcludeHours_'+str(i)] = cf[section][series]['ExcludeHours'][str(i)]
    if 'ExcludeHours' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',ExcludeHours'

def do_linear(cf,ds):
    log.info(' Applying linear corrections ...')
    level = ds.globalattributes['nc_level']
    for ThisOne in cf['Variables'].keys():
        if qcutils.haskey(cf,ThisOne,'Linear'):
            qcts.ApplyLinear(cf,ds,ThisOne)
        if qcutils.haskey(cf,ThisOne,'Drift'):
            qcts.ApplyLinearDrift(cf,ds,ThisOne)
        if qcutils.haskey(cf,ThisOne,'LocalDrift'):
            qcts.ApplyLinearDriftLocal(cf,ds,ThisOne)
    if 'do_linear' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',do_linear'

def do_rangecheck(cf,ds,section='',series='',code=2):
    '''Applies a range check to data series listed in the control file.  Data values that
       are less than the lower limit or greater than the upper limit are replaced with
       -9999 and the corresponding QC flag element is set to 2.'''
    if len(section)==0: return
    if len(series)==0: return
    if 'RangeCheck' not in cf[section][series].keys(): return
    if 'Lower' in cf[section][series]['RangeCheck'].keys():
        lwr = numpy.array(eval(cf[section][series]['RangeCheck']['Lower']))
        valid_lower = numpy.min(lwr)
        lwr = lwr[ds.series['Month']['Data']-1]
        index = numpy.where((abs(ds.series[series]['Data']-numpy.float64(-9999))>c.eps)&
                                (ds.series[series]['Data']<lwr))
        ds.series[series]['Data'][index] = numpy.float64(-9999)
        ds.series[series]['Flag'][index] = numpy.int32(code)
        ds.series[series]['Attr']['rangecheck_lower'] = cf[section][series]['RangeCheck']['Lower']
    if 'Upper' in cf[section][series]['RangeCheck'].keys():
        upr = numpy.array(eval(cf[section][series]['RangeCheck']['Upper']))
        valid_upper = numpy.min(upr)
        upr = upr[ds.series['Month']['Data']-1]
        index = numpy.where((abs(ds.series[series]['Data']-numpy.float64(-9999))>c.eps)&
                                (ds.series[series]['Data']>upr))
        ds.series[series]['Data'][index] = numpy.float64(-9999)
        ds.series[series]['Flag'][index] = numpy.int32(code)
        ds.series[series]['Attr']['rangecheck_upper'] = cf[section][series]['RangeCheck']['Upper']
        ds.series[series]['Attr']['valid_range'] = str(valid_lower)+','+str(valid_upper)
    if 'RangeCheck' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',RangeCheck'

def do_qcchecks(cf,ds):
    level = ds.globalattributes['nc_level']
    log.info(' Doing the QC checks at level '+str(level))
    for series in ds.series.keys():
            do_qcchecks_oneseries(cf,ds,series=series)

def do_qcchecks_oneseries(cf,ds,series=''):
    section = qcutils.get_cfsection(cf,series=series,mode='quiet')
    if len(section)==0: return
    level = ds.globalattributes['nc_level']
    if level == 'L2':
        range_code        = 2
        diurnal_code      = 3
        excludedates_code = 6
        excludehours_code = 7
    if level == 'L3':
        excludedates_code = 6
        excludehours_code = 7
        range_code        = 16
        diurnal_code      = 17
    if level == 'L4':
        excludedates_code = 6
        excludehours_code = 7
        range_code        = 38
        diurnal_code      = 39
    # do the range check
    do_rangecheck(cf,ds,section=section,series=series,code=range_code)
    # do the diurnal check
    do_diurnalcheck(cf,ds,section=section,series=series,code=diurnal_code)
    # do exclude dates
    do_excludedates(cf,ds,section=section,series=series,code=excludedates_code)
    # do exclude hours
    do_excludehours(cf,ds,section=section,series=series,code=excludehours_code)
    if 'do_qcchecks' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',do_qcchecks'
