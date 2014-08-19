"""
    QC Data Function Module
    Used to perform the tasks queued by qcls.py
    """

import sys
import ast
import constants as c
import datetime
from matplotlib.dates import date2num
import meteorologicalfunctions as mf
import numpy
import qcck
import qcio
import qcutils
from scipy import interpolate
import time
import xlrd
from matplotlib.mlab import griddata
import xlwt
import logging

log = logging.getLogger('qc.ts')

def albedo(cf,ds):
    """
        Filter albedo measurements to:
            high solar angle specified by periods between 10.00 and 14.00, inclusive
            and
            full sunlight in which Fsd > 290 W/m2
        
        Usage qcts.albedo(ds)
        ds: data structure
        """
    log.info(' Applying albedo constraints')
    if 'albedo' not in ds.series.keys():
        if 'Fsd' in ds.series.keys() and 'Fsu' in ds.series.keys():
            Fsd,f = qcutils.GetSeriesasMA(ds,'Fsd')
            Fsu,f = qcutils.GetSeriesasMA(ds,'Fsu')
            albedo = Fsu / Fsd
            attr = qcutils.MakeAttributeDictionary(long_name='solar albedo',units='none',standard_name='solar_albedo')
            qcutils.CreateSeries(ds,'albedo',albedo,FList=['Fsd','Fsu'],Attr=attr)
        else:
            log.warning('  Fsd or Fsu not in ds, albedo not calculated')
            return
    else:
        albedo,f = qcutils.GetSeriesasMA(ds,'albedo')
        if 'Fsd' in ds.series.keys():
            Fsd,f = qcutils.GetSeriesasMA(ds,'Fsd')
        else:
            Fsd,f = qcutils.GetSeriesasMA(ds,'Fn')
    
    if qcutils.cfkeycheck(cf,ThisOne='albedo',key='Threshold'):
        Fsdbase = float(cf['Variables']['albedo']['Threshold']['Fsd'])
        ds.series['albedo']['Attr']['FsdCutoff'] = Fsdbase
    else:
        Fsdbase = 290.
    index = numpy.ma.where((Fsd < Fsdbase) | (ds.series['Hdh']['Data'] < 10) | (ds.series['Hdh']['Data'] > 14))[0]
    index1 = numpy.ma.where(Fsd < Fsdbase)[0]
    index2 = numpy.ma.where((ds.series['Hdh']['Data'] < 10) | (ds.series['Hdh']['Data'] > 14))[0]
    albedo[index] = numpy.float64(-9999)
    ds.series['albedo']['Flag'][index1] = numpy.int32(51)     # bad Fsd flag only if bad time flag not set
    ds.series['albedo']['Flag'][index2] = numpy.int32(52)     # bad time flag
    ds.series['albedo']['Data']=numpy.ma.filled(albedo,float(-9999))

def ApplyLinear(cf,ds,ThisOne):
    """
        Applies a linear correction to variable passed from qcls. Time period
        to apply the correction, slope and offset are specified in the control
        file.
        
        Usage qcts.ApplyLinear(cf,ds,x)
        cf: control file
        ds: data structure
        x: input/output variable in ds.  Example: 'Cc_7500_Av'
        """
    log.info('  Applying linear correction to '+ThisOne)
    if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'Linear'):
        data = numpy.ma.masked_where(ds.series[ThisOne]['Data']==float(-9999),ds.series[ThisOne]['Data'])
        flag = ds.series[ThisOne]['Flag'].copy()
        ldt = ds.series['DateTime']['Data']
        LinearList = cf['Variables'][ThisOne]['Linear'].keys()
        for i in range(len(LinearList)):
            LinearItemList = ast.literal_eval(cf['Variables'][ThisOne]['Linear'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(LinearItemList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            try:
                ei = ldt.index(datetime.datetime.strptime(LinearItemList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            Slope = float(LinearItemList[2])
            Offset = float(LinearItemList[3])
            data[si:ei] = Slope * data[si:ei] + Offset
            index = numpy.where(flag[si:ei]==0)[0]
            flag[si:ei][index] = numpy.int32(10)
            ds.series[ThisOne]['Data'] = numpy.ma.filled(data,float(-9999))
            ds.series[ThisOne]['Flag'] = flag

def ApplyLinearDrift(cf,ds,ThisOne):
    """
        Applies a linear correction to variable passed from qcls. The slope is
        interpolated for each 30-min period between the starting value at time 0
        and the ending value at time 1.  Slope0, Slope1 and Offset are defined
        in the control file.  This function applies to a dataset in which the
        start and end times in the control file are matched by the time period
        in the dataset.
        
        Usage qcts.ApplyLinearDrift(cf,ds,x)
        cf: control file
        ds: data structure
        x: input/output variable in ds.  Example: 'Cc_7500_Av'
        """
    log.info('  Applying linear drift correction to '+ThisOne)
    if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'Drift'):
        data = numpy.ma.masked_where(ds.series[ThisOne]['Data']==float(-9999),ds.series[ThisOne]['Data'])
        flag = ds.series[ThisOne]['Flag']
        ldt = ds.series['DateTime']['Data']
        DriftList = cf['Variables'][ThisOne]['Drift'].keys()
        for i in range(len(DriftList)):
            DriftItemList = ast.literal_eval(cf['Variables'][ThisOne]['Drift'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(DriftItemList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            try:
                ei = ldt.index(datetime.datetime.strptime(DriftItemList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            Slope = numpy.zeros(len(data))
            Slope0 = float(DriftItemList[2])
            Slope1 = float(DriftItemList[3])
            Offset = float(DriftItemList[4])
            nRecs = len(Slope[si:ei])
            for i in range(nRecs):
                ssi = si + i
                Slope[ssi] = ((((Slope1 - Slope0) / nRecs) * i) + Slope0)
            data[si:ei] = Slope[si:ei] * data[si:ei] + Offset
            flag[si:ei] = 10
            ds.series[ThisOne]['Data'] = numpy.ma.filled(data,float(-9999))
            ds.series[ThisOne]['Flag'] = flag

def ApplyLinearDriftLocal(cf,ds,ThisOne):
    """
        Applies a linear correction to variable passed from qcls. The slope is
        interpolated since the starting value at time 0 using a known 30-min
        increment.  Slope0, SlopeIncrement and Offset are defined in the control
        file.  This function applies to a dataset in which the start time in the
        control file is matched by dataset start time, but in which the end time
        in the control file extends beyond the dataset end.
        
        Usage qcts.ApplyLinearDriftLocal(cf,ds,x)
        cf: control file
        ds: data structure
        x: input/output variable in ds.  Example: 'Cc_7500_Av'
        """
    log.info('  Applying linear drift correction to '+ThisOne)
    if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'LocalDrift'):
        data = numpy.ma.masked_where(ds.series[ThisOne]['Data']==float(-9999),ds.series[ThisOne]['Data'])
        flag = ds.series[ThisOne]['Flag']
        ldt = ds.series['DateTime']['Data']
        DriftList = cf['Variables'][ThisOne]['LocalDrift'].keys()
        for i in range(len(DriftList)):
            DriftItemList = ast.literal_eval(cf['Variables'][ThisOne]['LocalDrift'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(DriftItemList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            try:
                ei = ldt.index(datetime.datetime.strptime(DriftItemList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            Slope = numpy.zeros(len(data))
            Slope0 = float(DriftItemList[2])
            SlopeIncrement = float(DriftItemList[3])
            Offset = float(DriftItemList[4])
            nRecs = len(Slope[si:ei])
            for i in range(nRecs):
                ssi = si + i
                Slope[ssi] = (SlopeIncrement * i) + Slope0
            data[si:ei] = Slope[si:ei] * data[si:ei] + Offset
            flag[si:ei] = numpy.int32(10)
            ds.series[ThisOne]['Data'] = numpy.ma.filled(data,float(-9999))
            ds.series[ThisOne]['Flag'] = flag

def AverageSeriesByElements(cf,ds,Av_out):
    """
        Calculates the average of multiple time series.  Multiple time series
        are entered and a single time series representing the average at each
        observational period is returned.
        
        Usage qcts.AverageSeriesByElements(cf,ds,Av_out)
        cf: control file object (must contain an entry for Av_out)
        ds: data structure
        Av_out: output variable to ds.  Example: 'Fg'
        Series_in: input variable series in ds.  Example: ['Fg_8cma','Fg_8cmb']
        """
    if Av_out not in cf['Variables'].keys(): return
    if Av_out in ds.averageserieslist: return
    srclist, standardname = qcutils.GetAverageSeriesKeys(cf,Av_out)
    log.info(' Averaging series in '+str(srclist)+' into '+Av_out)
    nSeries = len(srclist)
    if nSeries==0:
        log.error('  AverageSeriesByElements: no input series specified for'+str(Av_out))
        return
    if nSeries==1:
        TmpArr_data = ds.series[srclist[0]]['Data'].copy()
        TmpArr_flag = ds.series[srclist[0]]['Flag'].copy()
        Av_data = numpy.ma.masked_where(TmpArr_data==float(-9999),TmpArr_data)
        Mx_flag = TmpArr_flag
        SeriesNameString = srclist[0]
        SeriesUnitString = ds.series[srclist[0]]['Attr']['units']
    else:
        TmpArr_data = ds.series[srclist[0]]['Data'].copy()
        TmpArr_flag = ds.series[srclist[0]]['Flag'].copy()
        SeriesNameString = srclist[0]
        srclist.remove(srclist[0])
        for ThisOne in srclist:
            SeriesNameString = SeriesNameString+', '+ThisOne
            TmpArr_data = numpy.vstack((TmpArr_data,ds.series[ThisOne]['Data'].copy()))
            TmpArr_flag = numpy.vstack((TmpArr_flag,ds.series[ThisOne]['Flag'].copy()))
        TmpArr_data = numpy.ma.masked_where(TmpArr_data==float(-9999),TmpArr_data)
        Av_data = numpy.ma.average(TmpArr_data,axis=0)
        Mx_flag = numpy.min(TmpArr_flag,axis=0)
    ds.averageserieslist.append(Av_out)
    attr = qcutils.MakeAttributeDictionary(long_name='Element-wise average of series '+SeriesNameString,
                                       standard_name=standardname,units=ds.series[srclist[0]]['Attr']['units'])
    qcutils.CreateSeries(ds,Av_out,Av_data,FList=srclist,Attr=attr)

def CalculateAvailableEnergy(ds,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg'):
    """
        Calculate the average energy as Fn - G.
        
        Usage qcts.CalculateAvailableEnergy(ds,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg')
        ds: data structure
        Fa_out: output available energy variable to ds.  Example: 'Fa'
        Fn_in: input net radiation in ds.  Example: 'Fn'
        Fg_in: input ground heat flux in ds.  Example: 'Fg'
        """
    log.info(' Calculating available energy from Fn and Fg')
    Fn,f = qcutils.GetSeriesasMA(ds,Fn_in)
    Fg,f = qcutils.GetSeriesasMA(ds,Fg_in)
    Fa = Fn - Fg
    attr = qcutils.MakeAttributeDictionary(long_name='Available energy using '+Fn_in+','+Fg_in,units='W/m2')
    qcutils.CreateSeries(ds,Fa_out,Fa,FList=[Fn_in,Fg_in],Attr=attr)

def CalculateFluxes(cf,ds,Ta_name='Ta',ps_name='ps',Ah_name='Ah',wT_in='wT',wA_in='wA',wC_in='wC',uw_in='uw',vw_in='vw',Fhv_out='Fhv',Fe_out='Fe',Fc_out='Fc',Fm_out='Fm',ustar_out='ustar'):
    """
        Calculate the fluxes from the rotated covariances.
        
        Usage qcts.CalculateFluxes(ds)
        ds: data structure
        
        Pre-requisite: CoordRotation2D
        
        Accepts meteorological constants or variables
        """
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='CF'):
        args = ast.literal_eval(cf['FunctionArgs']['CF'])
        Ta_name = args[0]
        Ah_name = args[1]
        ps_name = args[2]
        wT_in = args[3]
        wA_in = args[4]
        wC_in = args[5]
        uw_in = args[6]
        vw_in = args[7]
        Fh_out = args[8]
        Fe_out = args[9]
        Fc_out = args[10]
        Fm_out = args[11]
        ustar_out = args[12]
    long_name = ''
    if 'Massman' in ds.globalattributes['Functions']:
        long_name = ' and frequency response corrected'
    Ta,f = qcutils.GetSeriesasMA(ds,Ta_name)
    ps,f = qcutils.GetSeriesasMA(ds,ps_name)
    Ah,f = qcutils.GetSeriesasMA(ds,Ah_name)
    rhom,f = qcutils.GetSeriesasMA(ds,'rhom')
    RhoCp,f = qcutils.GetSeriesasMA(ds,'RhoCp')
    Lv,f = qcutils.GetSeriesasMA(ds,'Lv')
    
    log.info(' Calculating fluxes from covariances')
    if wT_in in ds.series.keys():
        wT,f = qcutils.GetSeriesasMA(ds,wT_in)
        Fhv = RhoCp * wT
        attr = qcutils.MakeAttributeDictionary(long_name='Virtual heat flux, rotated to natural wind coordinates'+long_name,
                                           units='W/m2')
        qcutils.CreateSeries(ds,Fhv_out,Fhv,FList=[wT_in],Attr=attr)
    else:
        log.error('  CalculateFluxes: '+wT_in+' not found in ds.series, Fh not calculated')
    if wA_in in ds.series.keys():
        wA,f = qcutils.GetSeriesasMA(ds,wA_in)
        Fe = Lv * wA / float(1000)
        attr = qcutils.MakeAttributeDictionary(long_name='Latent heat flux, rotated to natural wind coordinates'+long_name,
                                           standard_name='surface_upward_latent_heat_flux',units='W/m2')
        qcutils.CreateSeries(ds,Fe_out,Fe,FList=[wA_in],Attr=attr)
    else:
        log.error('  CalculateFluxes: '+wA_in+' not found in ds.series, Fe not calculated')
    if wC_in in ds.series.keys():
        wC,f = qcutils.GetSeriesasMA(ds,wC_in)
        Fc = wC
        attr = qcutils.MakeAttributeDictionary(long_name='CO2 flux, rotated to natural wind coordinates'+long_name,units='mg/m2/s')
        qcutils.CreateSeries(ds,Fc_out,Fc,FList=[wC_in],Attr=attr)
    else:
        log.error('  CalculateFluxes: '+wC_in+' not found in ds.series, Fc_raw not calculated')
    if uw_in in ds.series.keys():
        if vw_in in ds.series.keys():
            uw,f = qcutils.GetSeriesasMA(ds,uw_in)
            vw,f = qcutils.GetSeriesasMA(ds,vw_in)
            vs = uw*uw + vw*vw
            Fm = rhom * numpy.ma.sqrt(vs)
            us = numpy.ma.sqrt(numpy.ma.sqrt(vs))
            attr = qcutils.MakeAttributeDictionary(long_name='Momentum flux, rotated to natural wind coordinates'+long_name,units='kg/m/s2')
            qcutils.CreateSeries(ds,Fm_out,Fm,FList=[uw_in,vw_in],Attr=attr)
            attr = qcutils.MakeAttributeDictionary(long_name='Friction velocity, rotated to natural wind coordinates'+long_name,units='m/s')
            qcutils.CreateSeries(ds,ustar_out,us,FList=[uw_in,vw_in],Attr=attr)
        else:
            log.error('  CalculateFluxes: wy not found in ds.series, Fm and ustar not calculated')
    else:
        log.error('  CalculateFluxes: wx not found in ds.series, Fm and ustar not calculated')
    if 'CalculateFluxes' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', CalculateFluxes'

def CalculateLongwave(ds,Fl_out,Fl_in,Tbody_in):
    """
        Calculate the longwave radiation given the raw thermopile output and the
        sensor body temperature.
        
        Usage qcts.CalculateLongwave(ds,Fl_out,Fl_in,Tbody_in)
        ds: data structure
        Fl_out: output longwave variable to ds.  Example: 'Flu'
        Fl_in: input longwave in ds.  Example: 'Flu_raw'
        Tbody_in: input sensor body temperature in ds.  Example: 'Tbody'
        """
    log.info(' Calculating longwave radiation')
    Fl_raw,f = qcutils.GetSeriesasMA(ds,Fl_in)
    Tbody,f = qcutils.GetSeriesasMA(ds,Tbody_in)
    Fl = Fl_raw + c.sb*(Tbody + 273.15)**4
    attr = qcutils.MakeAttributeDictionary(long_name='Calculated longwave radiation using '+Fl_in+','+Tbody_in,units='W/m2')
    qcutils.CreateSeries(ds,Fl_out,Fl,FList=[Fl_in,Tbody_in],Attr=attr)

def CalculateMeteorologicalVariables(ds,Ta_name='Ta',Tv_name='Tv_CSAT',ps_name='ps',Ah_name='Ah',RH_name='RH',Cc_name='Cc'):
    """
        Add time series of meteorological variables based on fundamental
        relationships (Stull 1988)

        Usage qcts.CalculateMeteorologicalVariables(ds,Ta_name,ps_name,Ah_name)
        ds: data structure
        Ta_name: data series name for air temperature
        ps_name: data series name for pressure
        Ah_name: data series name for absolute humidity
        RH_name: data series for relative humidity

        Variables added:
            rhom: density of moist air, mf.densitymoistair(Ta,ps,Ah)
            Lv: latent heat of vapourisation, mf.Lv(Ta)
            q: specific humidity, mf.specifichumidity(mr)
                where mr (mixing ratio) = mf.mixingratio(ps,vp)
            Cpm: specific heat of moist air, mf.specificheatmoistair(q)
            VPD: vapour pressure deficit, VPD = esat - e
        """
    log.info(' Adding standard met variables to database')
    #if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='MetVars'):
        #vars = ast.literal_eval(cf['FunctionArgs']['MetVars'])
        #Ta_name = vars[0]
        #ps_name = vars[1]
        #Ah_name = vars[2]
    # get the required data series
    Ta,f = qcutils.GetSeriesasMA(ds,Ta_name)
    Tv,f = qcutils.GetSeriesasMA(ds,Tv_name)
    ps,f = qcutils.GetSeriesasMA(ds,ps_name)
    Ah,f = qcutils.GetSeriesasMA(ds,Ah_name)
    Cc,f = qcutils.GetSeriesasMA(ds,Cc_name)
    # calculate RH if it has not been read from the L1 spreadsheet
    if RH_name not in ds.series.keys():
        RH = mf.RHfromabsolutehumidity(Ah,Ta)     # relative humidity in units of percent
    else:
        RH,f = qcutils.GetSeriesasMA(ds,RH_name)
    # do the calculations
    e = mf.vapourpressure(Ah,Ta)                  # vapour pressure from absolute humidity and temperature
    esat = mf.es(Ta)                              # saturation vapour pressure
    rhod = mf.densitydryair(Ta,ps,e)              # partial density of dry air
    rhom = mf.densitymoistair(Ta,ps,e)            # density of moist air
    rhow = mf.densitywatervapour(Ta,e)            # partial density of water vapour
    Lv = mf.Lv(Ta)                                # latent heat of vapourisation
    mr = mf.mixingratio(ps,e)                     # mixing ratio
    mrsat = mf.mixingratio(ps,esat)               # saturation mixing ratio
    q = mf.specifichumidity(mr)                   # specific humidity from mixing ratio
    qsat = mf.specifichumidity(mrsat)             # saturation specific humidity from saturation mixing ratio
    Cpd = mf.specificheatcapacitydryair(Tv)
    Cpw = mf.specificheatcapacitywatervapour(Ta,RH)
    RhoCp = mf.densitytimesspecificheat(rhow,Cpw,rhod,Cpd)
    Cpm = mf.specificheatmoistair(q)              # specific heat of moist air
    VPD = esat - e                                # vapour pressure deficit
    SHD = qsat - q                                # specific humidity deficit
    c_ppm = mf.co2_ppmfrommgpm3(Cc,Ta,ps)         # CO2 concentration in units of umol/mol
    h_ppt = mf.h2o_mmolpmolfromgpm3(Ah,Ta,ps)     # H2O concentration in units of mmol/mol
    # write the meteorological series to the data structure
    if RH_name not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='not defined')
        qcutils.CreateSeries(ds,RH_name,RH,FList=[Ta_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Vapour pressure',units='kPa',standard_name='water_vapor_partial_pressure_in_air')
    qcutils.CreateSeries(ds,'e',e,FList=[Ta_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Saturation vapour pressure',units='kPa')
    qcutils.CreateSeries(ds,'esat',esat,FList=[Ta_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Density of dry air',units='kg/m3')
    qcutils.CreateSeries(ds,'rhod',rhod,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Density of moist air',units='kg/m3',standard_name='air_density')
    qcutils.CreateSeries(ds,'rhom',rhom,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Partial density of water vapour',units='kg/m3')
    qcutils.CreateSeries(ds,'rhow',rhow,FList=[Ta_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Latent heat of vapourisation',units='J/kg')
    qcutils.CreateSeries(ds,'Lv',Lv,FList=[Ta_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific heat capacity of dry air',units='J/kg-K')
    qcutils.CreateSeries(ds,'Cpd',Cpd,FList=[Tv_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific heat capacity of water vapour',units='J/kg-K')
    qcutils.CreateSeries(ds,'Cpw',Cpw,FList=[Ta_name,RH_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific heat capacity of moist air',units='J/kg-K')
    qcutils.CreateSeries(ds,'Cpm',Cpm,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Product of air density and specific heat capacity',units='J/m3-K')
    qcutils.CreateSeries(ds,'RhoCp',RhoCp,FList=[Ta_name,Tv_name,RH_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific humidity',units='kg/kg',standard_name='specific_humidity')
    qcutils.CreateSeries(ds,'q',q,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Vapour pressure deficit',units='kPa',standard_name='water_vapor_saturation_deficit_in_air')
    qcutils.CreateSeries(ds,'VPD',VPD,FList=[Ta_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific humidity deficit',units='kg/kg')
    qcutils.CreateSeries(ds,'SHD',SHD,FList=[Ta_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='CO2 concentration',units='umol/mol')
    qcutils.CreateSeries(ds,'C_ppm',c_ppm,FList=[Cc_name,Ta_name,ps_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='H2O concentration',units='mmol/mol')
    qcutils.CreateSeries(ds,'H_ppt',h_ppt,FList=[Ah_name,Ta_name,ps_name],Attr=attr)
    if 'CalculateMetVars' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', CalculateMetVars'

def CalculateNetRadiation(ds,Fn_out,Fsd_in,Fsu_in,Fld_in,Flu_in):
    """
        Calculate the net radiation from the 4 components of the surface
        radiation budget.
        
        Usage qcts.CalculateNetRadiation(ds,Fn_out,Fsd_in,Fsu_in,Fld_in,Flu_in)
        ds: data structure
        Fn_out: output net radiation variable to ds.  Example: 'Fn_KZ'
        Fsd_in: input downwelling solar radiation in ds.  Example: 'Fsd'
        Fsu_in: input upwelling solar radiation in ds.  Example: 'Fsu'
        Fld_in: input downwelling longwave radiation in ds.  Example: 'Fld'
        Flu_in: input upwelling longwave radiation in ds.  Example: 'Flu'
        """
    log.info(' Calculating net radiation from 4 components')
    if Fsd_in in ds.series.keys() and Fsu_in in ds.series.keys() and Fld_in in ds.series.keys() and Flu_in in ds.series.keys():
        Fsd,f = qcutils.GetSeriesasMA(ds,Fsd_in)
        Fsu,f = qcutils.GetSeriesasMA(ds,Fsu_in)
        Fld,f = qcutils.GetSeriesasMA(ds,Fld_in)
        Flu,f = qcutils.GetSeriesasMA(ds,Flu_in)
        Fn = (Fsd - Fsu) + (Fld - Flu)
        attr = qcutils.MakeAttributeDictionary(long_name='Calculated net radiation using '+Fsd_in+','+Fsu_in+','+Fld_in+','+Flu_in,
                             standard_name='surface_net_allwave_radiation',units='W/m2')
        qcutils.CreateSeries(ds,Fn_out,Fn,FList=[Fsd_in,Fsu_in,Fld_in,Flu_in],Attr=attr)
    else:
        nRecs = int(ds.globalattributes['nc_nrecs'])
        Fn = numpy.array([-9999]*nRecs,dtype=numpy.float64)
        flag = numpy.ones(nRecs,dtype=numpy.int32)
        attr = qcutils.MakeAttributeDictionary(long_name='Calculated net radiation (one or more components missing)',
                             standard_name='surface_net_allwave_radiation',units='W/m2')
        qcutils.CreateSeries(ds,Fn_out,Fn,Flag=flag,Attr=attr)

def ComputeDailySums(cf,ds,SumList,SubSumList,MinMaxList,MeanList,SoilList):
    """
        Computes daily sums, mininima and maxima on a collection variables in
        the L4 dataset containing gap filled fluxes.  Sums are computed only
        when the number of daily 30-min observations is equal to 48 (i.e., no
        missing data) to avoid biasing.  Output to an excel file that specified
        in the control file.
        
        Usage qcts.ComputeDailySums(cf,ds)
        cf: control file
        ds: data structure
        
        Parameters loaded from control file:
            M1st: dataset start month
            M2nd: dataset end month
            SumList: list of variables to be summed
            SubSumList: list of variables to sum positive and negative observations separately
            MinMaxList: list of variables to compute daily min & max
            SoilList: list of soil moisture measurements groups
            SW0, SW10, etc: list of soil moisture sensors at a common level (e.g., surface, 10cm, etc)
        
        Default List of sums:
            Rain, ET, Fe_MJ, Fh_MJ, Fg_MJ, Fld_MJ, Flu_MJ, Fnr_MJ, Fsd_MJ,
            Fsu_MJ, Fc_g, Fc_mmol
        Default List of sub-sums (sums split between positive and negative observations)
            Fe_MJ, Fh_MJ, Fg_MJ
        Default List of min/max:
            Ta_HMP, Vbat, Tpanel, Fc_mg, Fc_umol
        Default List of soil moisture measurements:
        """
    OutList = []
    SumOutList = []
    SubOutList = []
    MinMaxOutList = []
    MeanOutList = []
    SoilOutList = []
    
    for ThisOne in SubSumList:
        if ThisOne not in SumList:
            SumList.append(ThisOne)
    
    for ThisOne in SumList:
        if ThisOne == 'ET':
            if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='ETin'):
                Invar = ast.literal_eval(cf['Sums']['ETin'])
            else:
                Invar = ['Fe']
            Fe,f = qcutils.GetSeriesasMA(ds,Invar[0])
            if 'Lv' in ds.series.keys():
                Lv,f = qcutils.GetSeriesasMA(ds,'Lv')
            else:
                Lv = c.Lv
            ET = Fe * 60 * 30 * 1000 / (Lv * c.rho_water)  # mm/30min for summing
            attr = qcutils.MakeAttributeDictionary(long_name='Evapotranspiration Flux',units='mm')
            qcutils.CreateSeries(ds,'ET',ET,FList=Invar,Attr=attr)
            SumOutList.append('ET')
            OutList.append('ET')
            if ThisOne in SubSumList:
                SubOutList.append('ET')
        elif ThisOne == 'Energy':
            if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='Energyin'):
                EnergyIn = ast.literal_eval(cf['Sums']['Energyin'])
            else:
                EnergyIn = ['Fe', 'Fh', 'Fg']
            Fe,f = qcutils.GetSeriesasMA(ds,EnergyIn[0])
            Fh,f = qcutils.GetSeriesasMA(ds,EnergyIn[1])
            Fg,f = qcutils.GetSeriesasMA(ds,EnergyIn[2])
            EnergyOut = ['Fe_MJ','Fh_MJ','Fg_MJ']
            for index in range(0,3):
                convert_energy(ds,EnergyIn[index],EnergyOut[index])
                OutList.append(EnergyOut[index])
                SumOutList.append(EnergyOut[index])
                if ThisOne in SubSumList:
                    SubOutList.append(EnergyOut[index])
        elif ThisOne == 'Radiation':
            if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='Radin'):
                RadiationIn = ast.literal_eval(cf['Sums']['Radin'])
            else:
                RadiationIn = ['Fld','Flu','Fn','Fsd','Fsu']
            Fld,f = qcutils.GetSeriesasMA(ds,RadiationIn[0])
            Flu,f = qcutils.GetSeriesasMA(ds,RadiationIn[1])
            Fnr,f = qcutils.GetSeriesasMA(ds,RadiationIn[2])
            Fsd,f = qcutils.GetSeriesasMA(ds,RadiationIn[3])
            Fsu,f = qcutils.GetSeriesasMA(ds,RadiationIn[4])
            RadiationOut = ['Fld_MJ','Flu_MJ','Fnr_MJ','Fsd_MJ','Fsu_MJ']
            for index in range(0,5):
                convert_energy(ds,RadiationIn[index],RadiationOut[index])
                OutList.append(RadiationOut[index])
                SumOutList.append(RadiationOut[index])
                if ThisOne in SubSumList:
                    log.error('  Subsum: Negative radiation flux not defined')
        elif ThisOne == 'Carbon':
            if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='Cin'):
                CIn = ast.literal_eval(cf['Sums']['Cin'])
            else:
                CIn = ['Fc']
            
            if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='GPPin'):
                GPPIn = ast.literal_eval(cf['Sums']['GPPin'])
                GPP,f = qcutils.GetSeriesasMA(ds,GPPIn[0])
                Re,f = qcutils.GetSeriesasMA(ds,GPPIn[1])
                GPP_mmol = GPP * 1800 / 1000
                Re_mmol = Re * 1800 / 1000
                attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min GPP',units='mmol/m2',standard_name='gross_primary_productivity_of_carbon')
                qcutils.CreateSeries(ds,'GPP_mmol',GPP_mmol,FList=[GPPIn[0]],Attr=attr)
                attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Re',units='mmol/m2')
                qcutils.CreateSeries(ds,'Re_mmol',Re_mmol,FList=[GPPIn[1]],Attr=attr)
                GPPOut = ['GPP_mmol','Re_mmol']
                for listindex in range(0,2):
                    OutList.append(GPPOut[listindex])
                    SumOutList.append(GPPOut[listindex])
            
            Fc,f = qcutils.GetSeriesasMA(ds,CIn[0])
            Fc_umol = Fc * 1e6 / (1000 * 44)               # umol/m2-s for min/max
            Fc_mmol = Fc_umol * 1800 / 1000                # mmol/m2-30min for summing
            Fc_g = Fc * 1800 / 1000                        # g/m2-30min for summing
            attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Flux',units='mmol/m2',standard_name='surface_upward_mole_flux_of_carbon_dioxide')
            qcutils.CreateSeries(ds,'Fc_mmol',Fc_mmol,FList=CIn,Attr=attr)
            attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Flux',units='g/m2')
            qcutils.CreateSeries(ds,'Fc_g',Fc_g,FList=CIn,Attr=attr)
            COut = ['Fc_g','Fc_mmol']
            for listindex in range(0,2):
                OutList.append(COut[listindex])
                SumOutList.append(COut[listindex])
                if ThisOne in SubSumList:
                    SubOutList.append(COut[listindex])
        elif ThisOne == 'PM':
            if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='PMin'):
                if 'Gst' not in ds.series.keys():
                    SumList.remove('PM')
                    info.error('  Penman-Monteith Daily sum: input Source not located')
                else:
                    Gst_mmol,f = qcutils.GetSeriesasMA(ds,'Gst')   # mmol/m2-s
                    Gst_mol =  Gst_mmol * 1800 / 1000                 # mol/m2-30min for summing
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Bulk Stomatal Conductance',units='mol/m2')
                    qcutils.CreateSeries(ds,'Gst_mol',Gst_mol,FList=['Gst'],Attr=attr)
                    PMout = 'Gst_mol'
                    if PMout not in OutList:
                        OutList.append(PMout)
                    if ThisOne in SubSumList:
                        log.error('  Subsum: Negative bulk stomatal conductance not defined')
                    SumOutList.append(PMout)
            else:
                info.error('  Penman-Monteith Daily sums: input Source not defined')
        else:
            OutList.append(ThisOne)
            SumOutList.append(ThisOne)
    
    for ThisOne in MinMaxList:
        if ThisOne == 'Carbon':
            if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='Cin'):
                CIn = ast.literal_eval(cf['Sums']['Cin'])
            else:
                CIn = ['Fc']
            Fc,f = qcutils.GetSeriesasMA(ds,CIn[0])
            Fc_umol = Fc * 1e6 / (1000 * 44)               # umol/m2-s for min/max
            attr = qcutils.MakeAttributeDictionary(long_name='Average Flux',units='umol/(m2 s)',standard_name='surface_upward_mole_flux_of_carbon_dioxide')
            qcutils.CreateSeries(ds,'Fc_umol',Fc_umol,FList=CIn,Attr=attr)
            attr = qcutils.MakeAttributeDictionary(long_name='Average Flux',units='mg/(m2 s)')
            qcutils.CreateSeries(ds,'Fc_mg',Fc,FList=CIn,Attr=attr)
            COut = ['Fc_mg','Fc_umol']
            for listindex in range(0,2):
                OutList.append(COut[listindex])
                MinMaxOutList.append(COut[listindex])
        elif ThisOne == 'PM':
            if ThisOne not in SumList:
                if 'Gst' not in ds.series.keys() or 'rst' not in ds.series.keys():
                    MinMaxList.remove('PM')
                    PMout = []
                    info.error('  Penman-Monteith Daily min/max: input Source not located')
                else:
                    PMout = ['rst','Gst']
            else:
                PMout = ['rst','Gst']
            if len(PMout) > 0:
                for listindex in range(0,2):
                    if PMout[listindex] not in OutList:
                        OutList.append(PMout[listindex])
                    MinMaxOutList.append(PMout[listindex])
        else:
            if ThisOne not in OutList:
                OutList.append(ThisOne)
            MinMaxOutList.append(ThisOne)
    
    for ThisOne in MeanList:
        if ThisOne == 'Energy' or ThisOne == 'Carbon' or ThisOne == 'Radiation':
            log.error(' Mean error: '+ThisOne+' to be placed in SumList')
        elif ThisOne == 'PM':
            if ThisOne not in MinMaxList and ThisOne not in SumList:
                if 'Gst' not in ds.series.keys() or 'rst' not in ds.series.keys():
                    MeanList.remove('PM')
                    PMout = []
                    info.error('  Penman-Monteith Daily mean: input Source not located')
                else:
                    PMout = ['rst','Gst']
            else:
                PMout = ['rst','Gst']
            if len(PMout) > 0:
                for listindex in range(0,2):
                    if PMout[listindex] not in OutList:
                        OutList.append(PMout[listindex])
                    MeanOutList.append(PMout[listindex])
        else:
            MeanOutList.append(ThisOne)
            if ThisOne not in OutList:
                OutList.append(ThisOne)
    
    if len(SoilList) > 0:
        for ThisOne in SoilList:
            if qcutils.cfkeycheck(cf,Base='Sums',ThisOne=ThisOne):
                vars = ast.literal_eval(cf['Sums'][ThisOne])
                for index in range(0,len(vars)):
                    SoilOutList.append(vars[index])
                OutList.append(ThisOne)
    
    xlFileName = cf['Files']['L4']['xlSumFilePath']+cf['Files']['L4']['xlSumFileName']
    xlFile = xlwt.Workbook()
    
    for ThisOne in OutList:
        xlSheet = xlFile.add_sheet(ThisOne)
        xlCol = 0
        if ThisOne in SumOutList:
            if ThisOne in SubOutList:
                write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoSum='True',DoSubSum='True')
            else:
                write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoSum='True')
        
        if ThisOne in MinMaxOutList:
            if ThisOne in MeanOutList:
                write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoMinMax='True',DoMean='True')
            else:
                write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoMinMax='True')
        
        if ThisOne in MeanOutList:
            if ThisOne not in MinMaxOutList:
                write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoMean='True')
        
        if ThisOne in SoilList:
            soilvars = ast.literal_eval(cf['Sums'][ThisOne])
            for n in soilvars:
                if n == soilvars[0]:
                    xC,xS = write_sums(cf,ds,n,xlCol,xlSheet,DoSoil='True')
                else:
                    xC,xS = write_sums(cf,ds,n,xlCol,xS,DoSoil='True')
                xlCol = xC + 1
        
    log.info(' Saving Excel file '+xlFileName)
    xlFile.save(xlFileName)

    log.info(' Daily sums: All done')

def convert_energy(ds,InVar,OutVar):
    """
        Integrate energy flux over 30-min time period.
        Converts flux in W/m2 to MJ/(m2 30-min)
        
        Usage qcts.convert_energy(ds,InVar,OutVar)
        ds: data structure
        InVar: name of input variable.  Example: 'Fe_gapfilled'
        OutVar: name of output variable.  Example: 'Fe_MJ'
        """
    Wm2,f = qcutils.GetSeriesasMA(ds,InVar)
    MJ = Wm2 * 1800 / 1e6
    attr = qcutils.MakeAttributeDictionary(long_name=ds.series[InVar]['Attr']['long_name'],units='MJ/m2',standard_name=ds.series[InVar]['Attr']['standard_name'])
    qcutils.CreateSeries(ds,OutVar,MJ,FList=[InVar],Attr=attr)

def CoordRotation2D(cf,ds):
    """
        2D coordinate rotation to force v = w = 0.  Based on Lee et al, Chapter
        3 of Handbook of Micrometeorology.  This routine does not do the third
        rotation to force v'w' = 0.
        
        Usage qcts.CoordRotation2D(ds)
        ds: data structure
        """
    # get the raw wind velocity components
    Ux,f = qcutils.GetSeriesasMA(ds,'Ux')          # longitudinal component in CSAT coordinate system
    Uy,f = qcutils.GetSeriesasMA(ds,'Uy')          # lateral component in CSAT coordinate system
    Uz,f = qcutils.GetSeriesasMA(ds,'Uz')          # vertical component in CSAT coordinate system
    # get the raw covariances
    UxUz,f = qcutils.GetSeriesasMA(ds,'UxUz')      # covariance(Ux,Uz)
    UyUz,f = qcutils.GetSeriesasMA(ds,'UyUz')      # covariance(Uy,Uz)
    UxUy,f = qcutils.GetSeriesasMA(ds,'UxUy')      # covariance(Ux,Uy)
    UyUy,f = qcutils.GetSeriesasMA(ds,'UyUy')      # variance(Uy)
    UxUx,f = qcutils.GetSeriesasMA(ds,'UxUx')      # variance(Ux)
    UzUz,f = qcutils.GetSeriesasMA(ds,'UzUz')      # variance(Ux)
    UzC,f = qcutils.GetSeriesasMA(ds,'UzC')        # covariance(Uz,C)
    UzA,f = qcutils.GetSeriesasMA(ds,'UzA')        # covariance(Uz,A)
    UzT,f = qcutils.GetSeriesasMA(ds,'UzT')        # covariance(Uz,T)
    UxC,f = qcutils.GetSeriesasMA(ds,'UxC')        # covariance(Ux,C)
    UyC,f = qcutils.GetSeriesasMA(ds,'UyC')        # covariance(Uy,C)
    UxA,f = qcutils.GetSeriesasMA(ds,'UxA')        # covariance(Ux,A)
    UyA,f = qcutils.GetSeriesasMA(ds,'UyA')        # covariance(Ux,A)
    UxT,f = qcutils.GetSeriesasMA(ds,'UxT')        # covariance(Ux,T)
    UyT,f = qcutils.GetSeriesasMA(ds,'UyT')        # covariance(Uy,T)
    nRecs = int(ds.globalattributes['nc_nrecs'])   # number of records
    # apply 2D coordinate rotation unless otherwise specified in control file
    rotate = True
    if ('Options' in cf) and ('2DCoordRotation' in cf['Options'].keys()):
        if not cf['Options'].as_bool('2DCoordRotation'): rotate = False
    if rotate:
        log.info(' Applying 2D coordinate rotation to wind components and covariances')
        # get the 2D and 3D wind speeds
        ws2d = numpy.ma.sqrt(Ux**2 + Uy**2)
        ws3d = numpy.ma.sqrt(Ux**2 + Uy**2 + Uz**2)
        # get the sine and cosine of the angles through which to rotate
        #  - first we rotate about the Uz axis by eta to get v = 0
        #  - then we rotate about the v axis by theta to get w = 0
        ce = Ux/ws2d          # cos(eta)
        se = Uy/ws2d          # sin(eta)
        ct = ws2d/ws3d        # cos(theta)
        st = Uz/ws3d          # sin(theta)
        # get the rotation angles
        theta = numpy.rad2deg(numpy.arctan2(st,ct))
        eta = numpy.rad2deg(numpy.arctan2(se,ce))
        # do the wind velocity components first
        u = Ux*ct*ce + Uy*ct*se + Uz*st           # longitudinal component in natural wind coordinates
        v = Uy*ce - Ux*se                         # lateral component in natural wind coordinates
        w = Uz*ct - Ux*st*ce - Uy*st*se           # vertical component in natural wind coordinates
        # now do the scalar covariances
        wT = UzT*ct - UxT*st*ce - UyT*st*se       # covariance(w,T) in natural wind coordinate system
        wA = UzA*ct - UxA*st*ce - UyA*st*se       # covariance(w,A) in natural wind coordinate system
        wC = UzC*ct - UxC*st*ce - UyC*st*se       # covariance(w,C) in natural wind coordinate system
        # now do the momentum covariances
        # full equations, Wesely PhD thesis via James Cleverly and EddyPro
        uw = UxUz*ce*(ct*ct-st*st) - 2*UxUy*ct*st*ce*se + UyUz*se*(ct*ct-st*st) - \
             UxUx*ct*st*ce*ce - UyUy*ct*st*se*se + UzUz*ct*st # covariance(w,x) in natural wind coordinate system
        uv = UxUy*ct*(ce*ce-se*se) + UyUz*st*ce - UxUz*st*se - \
             UxUx*ct*ce*se + UyUy*ct*ce*se                    # covariance(x,y) in natural wind coordinate system
        vw = UyUz*ct*ce - UxUz*ct*se - UxUy*st*(ce*ce-se*se) + \
             UxUx*st*ce*se - UyUy*st*ce*se                    # covariance(w,y) in natural wind coordinate system
    else:
        log.info(' 2D coordinate rotation disabled, using unrotated components and covariances')
        # dummy series for rotation angles
        theta = numpy.zeros(nRecs)
        eta = numpy.zeros(nRecs)
        # unrotated wind components
        u = Ux           # unrotated x xomponent
        v = Uy           # unrotated y xomponent
        w = Uz           # unrotated z xomponent
        # unrotated covariances
        wT = UzT       # unrotated  wT covariance
        wA = UzA       # unrotated  wA covariance
        wC = UzC       # unrotated  wC covariance
        uw = UxUz      # unrotated  uw covariance
        vw = UyUz      # unrotated  vw covariance
    # store the rotated quantities in the nc object
    # default behaviour of CreateSeries is to use the maximum value of the QC flag for any series specified in FList
    attr = qcutils.MakeAttributeDictionary(long_name='Horizontal rotation angle',units='deg')
    qcutils.CreateSeries(ds,'eta',eta,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Vertical rotation angle',units='deg')
    qcutils.CreateSeries(ds,'theta',theta,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Longitudinal component of wind-speed in natural wind coordinates',units='m/s')
    qcutils.CreateSeries(ds,'u',u,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Lateral component of wind-speed in natural wind coordinates',units='m/s')
    qcutils.CreateSeries(ds,'v',v,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Vertical component of wind-speed in natural wind coordinates',units='m/s')
    qcutils.CreateSeries(ds,'w',w,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Kinematic heat flux, rotated to natural wind coordinates',units='mC/s')
    qcutils.CreateSeries(ds,'wT',wT,FList=['Ux','Uy','Uz','UxT','UyT','UzT'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Kinematic vapour flux, rotated to natural wind coordinates',units='g/m2/s')
    qcutils.CreateSeries(ds,'wA',wA,FList=['Ux','Uy','Uz','UxA','UyA','UzA'],Attr=attr)
    #ReplaceRotatedCovariance(cf,ds,'wA','UzA')
    attr = qcutils.MakeAttributeDictionary(long_name='Kinematic CO2 flux, rotated to natural wind coordinates',units='mg/m2/s')
    qcutils.CreateSeries(ds,'wC',wC,FList=['Ux','Uy','Uz','UxC','UyC','UzC'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Momentum flux X component, corrected to natural wind coordinates',units='m2/s2')
    qcutils.CreateSeries(ds,'uw',uw,FList=['Ux','Uy','Uz','UxUz','UxUx','UxUy'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Momentum flux Y component, corrected to natural wind coordinates',units='m2/s2')
    qcutils.CreateSeries(ds,'vw',vw,FList=['Ux','Uy','Uz','UyUz','UxUy','UyUy'],Attr=attr)
    # if RotateFlag is set, force the QC flag value from the maximum of the FList series to 11
    #if qcutils.cfkeycheck(cf,Base='General',ThisOne='RotateFlag') and cf['General']['RotateFlag'] == 'True':
        #keys = ['eta','theta','u','v','w','wT','wA','wC','uw','vw']
        #for ThisOne in keys:
            #testseries,f = qcutils.GetSeriesasMA(ds,ThisOne)
            #mask = numpy.ma.getmask(testseries)
            #index = numpy.where(mask.astype(int)==1)
            #ds.series[ThisOne]['Flag'][index] = numpy.int32(11)
    if 'CoordRotation2D' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', CoordRotation2D'
    if qcutils.cfoptionskey(cf,Key='RelaxRotation'):
        RotatedSeriesList = ['wT','wA','wC','uw','vw']
        NonRotatedSeriesList = ['UzT','UzA','UzC','UxUz','UyUz']
        for ThisOne, ThatOne in zip(RotatedSeriesList,NonRotatedSeriesList):
            ReplaceWhereMissing(ds.series[ThisOne],ds.series[ThisOne],ds.series[ThatOne],FlagValue=20)
        if 'RelaxRotation' not in ds.globalattributes['Functions']:
            ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', RelaxRotation'

def CalculateFcStorage(cf,ds,Fc_out='Fc_storage',CO2_in='Cc_7500_Av'):
    """
    Calculate CO2 flux storage term in the air column beneath the CO2 instrument.  This
    routine assumes the air column between the sensor and the surface is well mixed.
    
    Usage qcts.CalculateFcStorage(cf,ds,Fc_out,CO2_in)
    cf: control file object    
    ds: data structure
    Fc_out: series label of the CO2 flux storage term
    CO2_in: series label of the CO2 concentration
    
    Parameters loaded from control file:
        zms: measurement height from surface, m
    """
    if 'Fc_storage' not in ds.series.keys():
        if qcutils.cfkeycheck(cf,Base='General',ThisOne='zms'):
            log.info(' Calculating Fc storage from single height concentration time series')
            zms = float(cf['General']['zms'])
            Cc,Cc_flag = qcutils.GetSeriesasMA(ds,CO2_in,si=0,ei=-1)
            attr_in = qcutils.GetAttributeDictionary(ds,CO2_in)
            # calculate the storage
            dc = numpy.ma.ediff1d(Cc,to_begin=0)                      # CO2 concentration difference from timestep to timestep
            dt=86400*numpy.ediff1d(ds.series['xlDateTime']['Data'],to_begin=float(30)/1440)    # time step in seconds from the Excel datetime values
            Fc_storage = zms*dc/dt                                    # calculate the CO2 flux based on storage below the measurement height
            descr = 'Fc infered from CO2 storage using single point CO2 measurement'
            units_out = attr_in['units'].replace('m3','m2/s')
            attr_out = qcutils.MakeAttributeDictionary(long_name=descr,units=units_out)
            qcutils.CreateSeries(ds,Fc_out,Fc_storage,FList=[CO2_in],Attr=attr_out)
        else:
            log.error('CalculateFcStorage: zms expected in General section of control file but not found')
    else:
        log.info('CalculateFcStorage: Fc_storage found in data structure, not calculated')

def CorrectFcForStorage(cf,ds,Fc_out='Fc',Fc_in='Fc',Fc_storage_in='Fc_storage'):
    """
    Correct CO2 flux for storage in the air column beneath the CO2 instrument.
    
    Usage qcts.CorrectFcForStorage(cf,ds,Fc_out,Fc_in,Fc_storage)
    cf: control file object    
    ds: data structure
    Fc_out: series label of the corrected CO2 flux
    Fc_in: series label of the input CO2 flux
    Fc_storage: series label of the CO2 flux storage term

    """
    if not qcutils.cfoptionskey(cf,Key='ApplyFcStorage'): return
    if (Fc_in not in ds.series.keys()) or (Fc_storage not in ds.series.keys()): return
    Fc,Fc_flag = qcutils.GetSeriesasMA(ds,Fc_in)
    attr_in = qcutils.GetAttributeDictionary(ds,Fc_in)
    Fc_storage,Fc_storage_flag = qcutils.GetSeriesasMA(ds,Fc_storage_in)
    attr_storage = qcutils.GetAttributeDictionary(ds,Fc_storage)
    if attr_in['units']!=attr_storage['units']:
        log.error('CorrectFcForStorage: units of Fc do not match those of storage term, storage not applied')
        return
    log.info(' Applying storage correction to Fc')
    Fc = Fc + Fc_storage
    long_name = attr_in['long_name'] + 'corrected for storage using supplied storage term'
    attr_out = qcutils.MakeAttributeDictionary(long_name=long_name, units=attr_in['units'])
    qcutils.CreateSeries(ds,Fc_out,Fc,FList=[Fc_in,CO2_in],Attr=attr_out)
    if 'CorrectFcForStorage' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', CorrectFcForStorage'

def CorrectIndividualFgForStorage(cf,ds):
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='CFgArgs'):
        List = cf['FunctionArgs']['CFgArgs'].keys()
        for i in range(len(List)):
            CFgArgs = ast.literal_eval(cf['FunctionArgs']['CFgArgs'][str(i)])
            CorrectFgForStorage(cf,ds,Fg_out=CFgArgs[0],Fg_in=CFgArgs[1],Ts_in=CFgArgs[2],SWC_in=CFgArgs[3])
        return

def CorrectFgForStorage(cf,ds,Fg_out='Fg',Fg_in='Fg',Ts_in='Ts',SWC_in='Sws'):
    """
        Correct ground heat flux for storage in the layer above the heat flux plate
        
        Usage qcts.CorrectFgForStorage(cf,ds,Fg_out,Fg_in,Ts_in,Sws_in)
        ds: data structure
        Fg_out: output soil heat flux variable to ds.  Example: 'Fg'
        Fg_in: input soil heat flux in ds.  Example: 'Fg_Av'
        Ts_in: input soil temperature in ds.  Example: 'Ts'
        
        Parameters loaded from control file:
            FgDepth: Depth of soil heat flux plates, m
            BulkDensity: soil bulk density, kg/m3
            OrganicContent: soil organic content, fraction
            SwsDefault: default value of soil moisture content used when no sensors present
        """
    if 'Soil' not in cf.keys():
        log.error(' CorrectFgForStorage: [Soil] section not found in control file, Fg not corrected')
        return
    if Fg_in not in ds.series.keys() or Ts_in not in ds.series.keys():
        log.error(' CorrectFgForStorage: '+Fg_in+' or '+Ts_in+' not found in data structure, Fg not corrected')
        return
    log.info(' Correcting soil heat flux for storage')
    d = max(0.0,min(0.5,float(cf['Soil']['FgDepth'])))
    bd = max(1200.0,min(2500.0,float(cf['Soil']['BulkDensity'])))
    oc = max(0.0,min(1.0,float(cf['Soil']['OrganicContent'])))
    mc = 1.0 - oc
    Fg,Fg_flag = qcutils.GetSeriesasMA(ds,Fg_in)  # raw soil heat flux
    nRecs = len(Fg)                               # number of records in series
    Ts,f = qcutils.GetSeriesasMA(ds,Ts_in)        # soil temperature
    Sws_default = min(1.0,max(0.0,float(cf['Soil']['SwsDefault'])))
    if len(SWC_in) == 0:
        slist = []
        if qcutils.cfkeycheck(cf,Base='Soil',ThisOne='SwsSeries'):
            slist = ast.literal_eval(cf['Soil']['SwsSeries'])
        if len(slist)==0:
            log.info('  CorrectFgForStorage: Sws_default used for whole series')
            Sws = numpy.ones(nRecs)*Sws_default
        elif len(slist)==1:
            Sws,f = qcutils.GetSeriesasMA(ds,slist[0])
        else:
            MergeSeries(ds,'Sws',slist,[0,10])
            Sws,f = qcutils.GetSeriesasMA(ds,'Sws')
    else:
        slist = SWC_in
        Sws,f = qcutils.GetSeriesasMA(ds,SWC_in)
    iom = numpy.where(numpy.mod(f,10)!=0)[0]
    if len(iom)!=0:
        log.info('  CorrectFgForStorage: Sws_default used for '+str(len(iom))+' values')
        Sws[iom] = Sws_default
    dTs = numpy.ma.zeros(nRecs)
    dTs[1:] = numpy.diff(Ts)
    dt = numpy.ma.zeros(nRecs)
    dt[1:] = numpy.diff(date2num(ds.series['DateTime']['Data']))*float(86400)
    dt[0] = dt[1]
    Cs = mc*bd*c.Cd + oc*bd*c.Co + Sws*c.rho_water*c.Cw
    S = Cs*(dTs/dt)*d
    Fg_o = Fg + S
    attr= qcutils.MakeAttributeDictionary(long_name='Soil heat flux corrected for storage',units='W/m2',standard_name='downward_heat_flux_in_soil')
    qcutils.CreateSeries(ds,Fg_out,Fg_o,FList=[Fg_in],Attr=attr)
    # save the input (uncorrected) soil heat flux series, this will be used if the correction is relaxed
    attr = qcutils.MakeAttributeDictionary(long_name='Soil heat flux uncorrected for storage',units='W/m2')
    qcutils.CreateSeries(ds,'Fg_Av',Fg,Flag=Fg_flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Soil heat flux storage',units='W/m2')
    qcutils.CreateSeries(ds,'S',S,FList=[Fg_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific heat capacity',units='J/m3/K')
    qcutils.CreateSeries(ds,'Cs',Cs,FList=[Fg_in],Attr=attr)
    if 'CorrectFgForStorage' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', CorrectFgForStorage'
    if qcutils.cfoptionskey(cf,Key='RelaxFgStorage'):
        ReplaceWhereMissing(ds.series['Fg'],ds.series['Fg'],ds.series['Fg_Av'],FlagValue=20)
        if 'RelaxFgStorage' not in ds.globalattributes['Functions']:
            ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', RelaxFgStorage'

def CorrectSWC(cf,ds):
    """
        Correct soil moisture data using calibration curve developed from
        collected soil samples.  To avoid unrealistic or unphysical extreme
        values upon extrapolation, exponential and logarithmic using ln
        functions are applied to small and large values, respectively.
        Threshold values where one model replaces the other is determined where
        the functions cross.  The logarithmic curve is constrained at with a
        point at which the soil measurement = field porosity and the sensor
        measurement is maximised under saturation at field capacity.
        
        Usage qcts.CorrectSWC(cf,ds)
        cf: control file
        ds: data structure
        
        Parameters loaded from control file:
            SWCempList: list of raw CS616 variables
            SWCoutList: list of corrected CS616 variables
            SWCattr:  list of meta-data attributes for corrected CS616 variables
            SWC_a0: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            SWC_a1: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            SWC_b0: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            SWC_b1: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            SWC_t: threshold parameter for switching from exponential to logarithmic model
            TDRempList: list of raw CS610 variables
            TDRoutList: list of corrected CS610 variables
            TDRattr:  list of meta-data attributes for corrected CS610 variables
            TDRlinList: list of deep TDR probes requiring post-hoc linear correction to match empirical samples
            TDR_a0: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            TDR_a1: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            TDR_b0: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            TDR_b1: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            TDR_t: threshold parameter for switching from exponential to logarithmic model
        """
    if not qcutils.cfoptionskey(cf,Key='CorrectSWC'): return
    log.info(' Correcting soil moisture data ...')
    SWCempList = ast.literal_eval(cf['Soil']['empSWCin'])
    SWCoutList = ast.literal_eval(cf['Soil']['empSWCout'])
    SWCattr = ast.literal_eval(cf['Soil']['SWCattr'])
    if cf['Soil']['TDR']=='Yes':
        TDRempList = ast.literal_eval(cf['Soil']['empTDRin'])
        TDRoutList = ast.literal_eval(cf['Soil']['empTDRout'])
        TDRlinList = ast.literal_eval(cf['Soil']['linTDRin'])
        TDRattr = ast.literal_eval(cf['Soil']['TDRattr'])
        TDR_a0 = float(cf['Soil']['TDR_a0'])
        TDR_a1 = float(cf['Soil']['TDR_a1'])
        TDR_b0 = float(cf['Soil']['TDR_b0'])
        TDR_b1 = float(cf['Soil']['TDR_b1'])
        TDR_t = float(cf['Soil']['TDR_t'])
    SWC_a0 = float(cf['Soil']['SWC_a0'])
    SWC_a1 = float(cf['Soil']['SWC_a1'])
    SWC_b0 = float(cf['Soil']['SWC_b0'])
    SWC_b1 = float(cf['Soil']['SWC_b1'])
    SWC_t = float(cf['Soil']['SWC_t'])
    
    for i in range(len(SWCempList)):
        log.info('  Applying empirical correction to '+SWCempList[i])
        invar = SWCempList[i]
        outvar = SWCoutList[i]
        attr = SWCattr[i]
        Sws,f = qcutils.GetSeriesasMA(ds,invar)
        
        nRecs = len(Sws)
        
        Sws_out = numpy.ma.empty(nRecs,float)
        Sws_out.fill(-9999)
        Sws_out.mask = numpy.ma.empty(nRecs,bool)
        Sws_out.mask.fill(True)
        
        index_high = numpy.ma.where((Sws.mask == False) & (Sws > SWC_t))[0]
        index_low = numpy.ma.where((Sws.mask == False) & (Sws < SWC_t))[0]
        
        Sws_out[index_low] = SWC_b0 * numpy.exp(SWC_b1 * Sws[index_low])
        Sws_out[index_high] = (SWC_a1 * numpy.log(Sws[index_high])) + SWC_a0
        
        attr = qcutils.MakeAttributeDictionary(long_name=attr,units='cm3 water/cm3 soil',standard_name='soil_moisture_content')
        qcutils.CreateSeries(ds,outvar,Sws_out,FList=[invar],Attr=attr)
    if cf['Soil']['TDR']=='Yes':
        for i in range(len(TDRempList)):
            log.info('  Applying empirical correction to '+TDRempList[i])
            invar = TDRempList[i]
            outvar = TDRoutList[i]
            attr = TDRattr[i]
            Sws,f = qcutils.GetSeriesasMA(ds,invar)
            
            nRecs = len(Sws)
            
            Sws_out = numpy.ma.empty(nRecs,float)
            Sws_out.fill(-9999)
            Sws_out.mask = numpy.ma.empty(nRecs,bool)
            Sws_out.mask.fill(True)
            
            index_high = numpy.ma.where((Sws.mask == False) & (Sws > TDR_t))[0]
            index_low = numpy.ma.where((Sws.mask == False) & (Sws < TDR_t))[0]
            
            Sws_out[index_low] = TDR_b0 * numpy.exp(TDR_b1 * Sws[index_low])
            Sws_out[index_high] = (TDR_a1 * numpy.log(Sws[index_high])) + TDR_a0
            
            attr = qcutils.MakeAttributeDictionary(long_name=attr,units='cm3 water/cm3 soil',standard_name='soil_moisture_content')
            qcutils.CreateSeries(ds,outvar,Sws_out,FList=[invar],Attr=attr)

def CorrectWindDirection(cf,ds,Wd_in):
    """
        Correct wind direction for mis-aligned sensor direction.
        
        Usage qcts.CorrectWindDirection(cf,ds,Wd_in)
        cf: control file
        ds: data structure
        Wd_in: input/output wind direction variable in ds.  Example: 'Wd_CSAT'
        """
    log.info(' Correcting wind direction')
    Wd,f = qcutils.GetSeriesasMA(ds,Wd_in)
    ldt = ds.series['DateTime']['Data']
    KeyList = cf['Variables'][Wd_in]['Correction'].keys()
    for i in range(len(KeyList)):
        ItemList = ast.literal_eval(cf['Variables'][Wd_in]['Correction'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ItemList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ItemList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        Correction = float(ItemList[2])
        Wd[si:ei] = Wd[si:ei] + Correction
    Wd = numpy.mod(Wd,float(360))
    ds.series[Wd_in]['Data'] = numpy.ma.filled(Wd,float(-9999))

def DailyAverageSws_Interpolated(cf,ds,Sws_out='Sws_daily',Sws_in='Sws'):
    '''
    Create a series of daily averaged soil moisture data and then interpolate this
    back on to the time step of the data.  This result is a time series of soil
    moisture data that is less noisy than the data at the original time step but
    still resolves day-to-day changes and seasonal trends.
    '''
    Sws,f = qcutils.GetSeriesasMA(ds,Sws_in)
    Ddd,f = qcutils.GetSeriesasMA(ds,'Ddd')
    pass

def do_attributes(cf,ds):
    """
        Import attriubes in xl2nc control file to netCDF dataset.  Included
        global and variable attributes.  Also attach flag definitions to global
        meta-data for reference.
        
        Usage qcts.do_attributes(cf,ds)
        cf: control file
        ds: data structure
        """
    log.info(' Getting the attributes given in control file')
    if 'Global' in cf.keys():
        for gattr in cf['Global'].keys():
            ds.globalattributes[gattr] = cf['Global'][gattr]
        ds.globalattributes['Flag0'] = 'Good data'
        ds.globalattributes['Flag1'] = 'QA/QC: -9999 in level 1 dataset'
        ds.globalattributes['Flag2'] = 'QA/QC: L2 Range Check'
        ds.globalattributes['Flag3'] = 'QA/QC: CSAT Diagnostic'
        ds.globalattributes['Flag4'] = 'QA/QC: LI7500 Diagnostic'
        ds.globalattributes['Flag5'] = 'QA/QC: L2 Diurnal SD Check'
        ds.globalattributes['Flag6'] = 'QA/QC: Excluded Dates'
        ds.globalattributes['Flag7'] = 'QA/QC: Excluded Hours'
        ds.globalattributes['Flag10'] = 'Corrections: Apply Linear'
        ds.globalattributes['Flag11'] = 'Corrections/Combinations: Coordinate Rotation (Ux, Uy, Uz, UxT, UyT, UzT, UxA, UyA, UzA, UxC, UyC, UzC, UxUz, UxUx, UxUy, UyUz, UxUy, UyUy)'
        ds.globalattributes['Flag12'] = 'Corrections/Combinations: Massman Frequency Attenuation Correction (Coord Rotation, Tv_CSAT, Ah_HMP, ps)'
        ds.globalattributes['Flag13'] = 'Corrections/Combinations: Virtual to Actual Fh (Coord Rotation, Massman, Ta_HMP)'
        ds.globalattributes['Flag14'] = 'Corrections/Combinations: WPL correction for flux effects on density measurements (Coord Rotation, Massman, Fhv to Fh, Cc_7500_Av)'
        ds.globalattributes['Flag15'] = 'Corrections/Combinations: Ta from Tv'
        ds.globalattributes['Flag16'] = 'Corrections/Combinations: L3 Range Check'
        ds.globalattributes['Flag17'] = 'Corrections/Combinations: L3 Diurnal SD Check'
        ds.globalattributes['Flag18'] = 'Corrections/Combinations: u* filter'
        ds.globalattributes['Flag19'] = 'Corrections/Combinations: Gap coordination'
        ds.globalattributes['Flag30'] = 'GapFilling: Flux Gap Filled by ANN (SOLO)'
        ds.globalattributes['Flag31'] = 'GapFilling: Flux Gap not Filled by ANN'
        ds.globalattributes['Flag32'] = 'GapFilling: Met Gap Filled from Climatology'
        ds.globalattributes['Flag33'] = 'GapFilling: Gap Filled from Ratios'
        ds.globalattributes['Flag34'] = 'GapFilling: Gap Filled by Interpolation'
        ds.globalattributes['Flag35'] = 'GapFilling: Gap Filled by Replacement'
        ds.globalattributes['Flag36'] = 'GapFilling: u* from Fh'
        ds.globalattributes['Flag37'] = 'GapFilling: u* not from Fh'
        ds.globalattributes['Flag38'] = 'GapFilling: L4 Range Check'
        ds.globalattributes['Flag39'] = 'GapFilling: L4 Diurnal SD Check'
        ds.globalattributes['Flag51'] = 'albedo: bad Fsd < threshold (290 W/m2 default) only if bad time flag (31) not set'
        ds.globalattributes['Flag52'] = 'albedo: bad time flag (not midday 10.00 to 14.00)'
        ds.globalattributes['Flag61'] = 'Penman-Monteith: bad rst (rst < 0) only if bad Uavg (35), bad Fe (33) and bad Fsd (34) flags not set'
        ds.globalattributes['Flag62'] = 'Penman-Monteith: bad Fe < threshold (0 W/m2 default) only if bad Fsd (34) flag not set'
        ds.globalattributes['Flag63'] = 'Penman-Monteith: bad Fsd < threshold (10 W/m2 default)'
        ds.globalattributes['Flag64'] = 'Penman-Monteith: Uavg == 0 (undefined aerodynamic resistance under calm conditions) only if bad Fe (33) and bad Fsd (34) flags not set'
        ds.globalattributes['Flag70'] = 'Partitioning Night: Re computed from exponential temperature response curves'
        ds.globalattributes['Flag80'] = 'Partitioning Day: GPP/Re computed from light-response curves, GPP = Re - Fc'
        ds.globalattributes['Flag81'] = 'Partitioning Day: GPP night mask'
        ds.globalattributes['Flag82'] = 'Partitioning Day: Fc > Re, GPP = 0, Re = Fc'
    for ThisOne in ds.series.keys():
        if ThisOne in cf['Variables']:
            if 'Attr' in cf['Variables'][ThisOne].keys():
                ds.series[ThisOne]['Attr'] = {}
                for attr in cf['Variables'][ThisOne]['Attr'].keys():
                    ds.series[ThisOne]['Attr'][attr] = cf['Variables'][ThisOne]['Attr'][attr]

def do_functions(cf,ds):
    log.info(' Getting variances from standard deviations & vice versa')
    if 'AhAh' in ds.series.keys() and 'Ah_7500_Sd' not in ds.series.keys():
        AhAh, flag = qcutils.GetSeriesasMA(ds,'AhAh')
        Ah_7500_Sd = numpy.ma.sqrt(AhAh)
        attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity from Li-7500, standard deviation',units='g/m3')
        qcutils.CreateSeries(ds,'Ah_7500_Sd',Ah_7500_Sd,Flag=flag,Attr=attr)
    if 'Ah_7500_Sd' in ds.series.keys() and 'AhAh' not in ds.series.keys():
        Ah_7500_Sd, flag = qcutils.GetSeriesasMA(ds,'Ah_7500_Sd')
        AhAh = Ah_7500_Sd*Ah_7500_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity from Li-7500, variance',units='(g/m3)2')
        qcutils.CreateSeries(ds,'AhAh',AhAh,Flag=flag,Attr=attr)
    if 'CcCc' in ds.series.keys() and 'Cc_7500_Sd' not in ds.series.keys():
        CcCc, flag = qcutils.GetSeriesasMA(ds,'CcCc')
        Cc_7500_Sd = numpy.ma.sqrt(CcCc)
        attr = qcutils.MakeAttributeDictionary(long_name='CO2 concentration from Li-7500, standard deviation',units='mg/m3')
        qcutils.CreateSeries(ds,'Cc_7500_Sd',Cc_7500_Sd,Flag=flag,Attr=attr)
    if 'Cc_7500_Sd' in ds.series.keys() and 'CcCc' not in ds.series.keys():
        Cc_7500_Sd, flag = qcutils.GetSeriesasMA(ds,'Cc_7500_Sd')
        CcCc = Cc_7500_Sd*Cc_7500_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='CO2 concentration from Li-7500, variance',units='(mg/m3)2')
        qcutils.CreateSeries(ds,'CcCc',CcCc,Flag=flag,Attr=attr)
    if 'Ux_Sd' in ds.series.keys() and 'UxUx' not in ds.series.keys():
        Ux_Sd, flag = qcutils.GetSeriesasMA(ds,'Ux_Sd')
        UxUx = Ux_Sd*Ux_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Longitudinal velocity component from CSAT, variance',units='(m/s)2')
        qcutils.CreateSeries(ds,'UxUx',UxUx,Flag=flag,Attr=attr)
    if 'Uy_Sd' in ds.series.keys() and 'UyUy' not in ds.series.keys():
        Uy_Sd, flag = qcutils.GetSeriesasMA(ds,'Uy_Sd')
        UyUy = Uy_Sd*Uy_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Lateral velocity component from CSAT, variance',units='(m/s)2')
        qcutils.CreateSeries(ds,'UyUy',UyUy,Flag=flag,Attr=attr)
    if 'Uz_Sd' in ds.series.keys() and 'UzUz' not in ds.series.keys():
        Uz_Sd, flag = qcutils.GetSeriesasMA(ds,'Uz_Sd')
        UzUz = Uz_Sd*Uz_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Vertical velocity component from CSAT, variance',units='(m/s)2')
        qcutils.CreateSeries(ds,'UzUz',UzUz,Flag=flag,Attr=attr)

def do_solo(cf,ds4,Fc_in='Fc',Fe_in='Fe',Fh_in='Fh',Fc_out='Fc',Fe_out='Fe',Fh_out='Fh'):
    ''' duplicate gapfilled fluxes for graphing comparison'''
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='SOLOvars'):
        invars = ast.literal_eval(cf['FunctionArgs']['SOLOvars'])
        Fc_in = invars[0]
        Fe_in = invars[1]
        Fh_in = invars[2]
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='SOLOplot'):
        outvars = ast.literal_eval(cf['FunctionArgs']['SOLOplot'])
        Fc_out = outvars[0]
        Fe_out = outvars[1]
        Fh_out = outvars[2]
    # add relevant meteorological values to L3 data
    log.info(' Adding standard met variables to database')
    CalculateMeteorologicalVariables(ds4)
    ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions']+', CalculateMetVars'
    if Fe_in in ds4.series.keys():
        Fe,flag = qcutils.GetSeriesasMA(ds4,Fe_in)
        attr = qcutils.MakeAttributeDictionary(long_name='ANN gapfilled Latent Heat Flux',units='W/m2',standard_name='surface_upward_latent_heat_flux')
        qcutils.CreateSeries(ds4,Fe_out,Fe,Flag=flag,Attr=attr)
    if Fc_in in ds4.series.keys():
        Fc,flag = qcutils.GetSeriesasMA(ds4,Fc_in)
        attr = qcutils.MakeAttributeDictionary(long_name='ANN gapfilled Carbon Flux',units='mg/m2/s')
        qcutils.CreateSeries(ds4,Fc_out,Fc,Flag=flag,Attr=attr)
    if Fh_in in ds4.series.keys():
        Fh,flag = qcutils.GetSeriesasMA(ds4,Fh_in)
        attr = qcutils.MakeAttributeDictionary(long_name='ANN gapfilled Sensible Heat Flux',units='W/m2',standard_name='surface_upward_sensible_heat_flux')
        qcutils.CreateSeries(ds4,Fh_out,Fh,Flag=flag,Attr=attr)

def do_sums(cf,ds):
    if not qcutils.cfoptionskey(cf,Key='DoSums',default=False): return
    # compute daily statistics
    if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='SumList'):
        SumList = ast.literal_eval(cf['Sums']['SumList'])
    else:
        SumList = ['Rain','ET','Energy','Radiation','Carbon']
    
    if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='SubSumList'):
        SubSumList = ast.literal_eval(cf['Sums']['SubSumList'])
    else:
        SubSumList = []
    
    if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='MinMaxList'):
        MinMaxList = ast.literal_eval(cf['Sums']['MinMaxList'])
    else:
        MinMaxList = ['Ta','Vbat','Tpanel','Carbon']
    
    if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='MeanList'):
        MeanList = ast.literal_eval(cf['Sums']['MeanList'])
    else:
        MeanList = ['Ta','Tpanel']
    
    if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='SoilList'):
        SoilList = ast.literal_eval(cf['Sums']['SoilList'])
    else:
        SoilList = []
    
    StatsList = SumList + MinMaxList + MeanList + SoilList
    if len(StatsList) > 0:
        qcts.ComputeDailySums(cf,ds,SumList,SubSumList,MinMaxList,MeanList,SoilList)

def Fc_WPL(cf,ds,Fc_wpl_out='Fc',Fc_raw_in='Fc',Fh_in='Fh',Fe_in='Fe',Ta_in='Ta',Ah_in='Ah',Cc_in='Cc',ps_in='ps'):
    """
        Apply Webb, Pearman and Leuning correction to carbon flux.  This
        correction is necessary to account for flux effects on density
        measurements.  Original formulation: Campbell Scientific
        
        Usage qcts.Fc_WPL(ds,Fc_wpl_out,Fc_raw_in,Fh_in,Fe_raw_in,Ta_in,Ah_in,Cc_in,ps_in)
        ds: data structure
        Fc_wpl_out: output corrected carbon flux variable to ds.  Example: 'Fc'
        Fc_raw_in: input carbon flux in ds.  Example: 'Fc'
        Fh_in: input sensible heat flux in ds.  Example: 'Fh'
        Fe_raw_in: input uncorrected latent heat flux in ds.  Example: 'Fe_raw'
        Ta_in: input air temperature in ds.  Example: 'Ta'
        Ah_in: input absolute humidity in ds.  Example: 'Ah'
        Cc_in: input co2 density in ds.  Example: 'Cc'
        ps_in: input atmospheric pressure in ds.  Example: 'ps'
        
        Used for fluxes that are raw or rotated.
        
        Pre-requisite: CalculateFluxes, CalculateFluxes_Unrotated or CalculateFluxesRM
        Pre-requisite: FhvtoFh
        Pre-requisite: Fe_WPL
        
        Accepts meteorological constants or variables
        """
    if 'DisableFcWPL' in cf['Options'] and cf['Options'].as_bool('DisableFcWPL'):
        log.error(" qcts.Fc_WPL: WPL correction for Fc disabled in control file")
        return
    log.info(' Applying WPL correction to Fc')
    Fc_raw,Fc_raw_flag = qcutils.GetSeriesasMA(ds,Fc_raw_in)
    Fh,f = qcutils.GetSeriesasMA(ds,Fh_in)
    Fe,f = qcutils.GetSeriesasMA(ds,Fe_in)
    Ta,f = qcutils.GetSeriesasMA(ds,Ta_in)
    TaK = Ta+c.C2K                                # air temperature from C to K
    Ah,f = qcutils.GetSeriesasMA(ds,Ah_in)
    Ah = Ah*c.g2kg                                # absolute humidity from g/m3 to kg/m3
    Cc,f = qcutils.GetSeriesasMA(ds,Cc_in)
    ps,f = qcutils.GetSeriesasMA(ds,ps_in)
    rhod,f = qcutils.GetSeriesasMA(ds,'rhod')
    RhoCp, f = qcutils.GetSeriesasMA(ds,'RhoCp')
    Lv,f = qcutils.GetSeriesasMA(ds,'Lv')
    sigma = Ah/rhod
    co2_wpl_Fe = (c.mu/(1+c.mu*sigma))*(Cc/rhod)*(Fe/Lv)
    co2_wpl_Fh = (Cc/TaK)*(Fh/RhoCp)
    Fc_wpl_data = Fc_raw+co2_wpl_Fe+co2_wpl_Fh
    Fc_wpl_flag = numpy.zeros(len(Fc_wpl_data))
    index = numpy.where(Fc_wpl_data.mask==True)[0]
    Fc_wpl_flag[index] = numpy.int32(14)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL corrected Fc',units='mg/m2/s')
    qcutils.CreateSeries(ds,Fc_wpl_out,Fc_wpl_data,Flag=Fc_wpl_flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL correction to Fc due to Fe',units='mg/m2/s')
    qcutils.CreateSeries(ds,'co2_wpl_Fe',co2_wpl_Fe,Flag=Fc_wpl_flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL correction to Fc due to Fh',units='mg/m2/s')
    qcutils.CreateSeries(ds,'co2_wpl_Fh',co2_wpl_Fh,Flag=Fc_wpl_flag,Attr=attr)

def Fe_WPL(cf,ds,Fe_wpl_out='Fe',Fe_raw_in='Fe',Fh_in='Fh',Ta_in='Ta',Ah_in='Ah',ps_in='ps'):
    """
        Apply Webb, Pearman and Leuning correction to vapour flux.  This
        correction is necessary to account for flux effects on density
        measurements.  Original formulation: Campbell Scientific
        
        Usage qcts.Fe_WPL(ds,Fe_wpl_out,Fe_raw_in,Fh_in,Ta_in,Ah_in,ps_in)
        ds: data structure
        Fe_wpl_out: output corrected water vapour flux variable to ds.  Example: 'Fe'
        Fe_raw_in: input water vapour flux in ds.  Example: 'Fe'
        Fh_in: input sensible heat flux in ds.  Example: 'Fh'
        Ta_in: input air temperature in ds.  Example: 'Ta'
        Ah_in: input absolute humidity in ds.  Example: 'Ah'
        ps_in: input atmospheric pressure in ds.  Example: 'ps'
        
        Used for fluxes that are raw or rotated.
        
        Pre-requisite: CalculateFluxes, CalculateFluxes_Unrotated or CalculateFluxesRM
        Pre-requisite: FhvtoFh
        
        Accepts meteorological constants or variables
        """
    if 'DisableFeWPL' in cf['Options'] and cf['Options'].as_bool('DisableFeWPL'):
        log.error(" qcts.Fe_WPL: WPL correction for Fe disabled in control file")
        return
    log.info(' Applying WPL correction to Fe')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='EWPL'):
        Eargs = ast.literal_eval(cf['FunctionArgs']['EWPL'])
        Fe_wpl_out = Eargs[0]
        Fe_raw_in = Eargs[1]
        Fh_in = Eargs[2]
        Ta_in = Eargs[3]
        Ah_in = Eargs[4]
        ps_in = Eargs[5]
    Fe_raw,Fe_raw_flag = qcutils.GetSeriesasMA(ds,Fe_raw_in)
    Fh,f = qcutils.GetSeriesasMA(ds,Fh_in)
    Ta,f = qcutils.GetSeriesasMA(ds,Ta_in)
    TaK = Ta + c.C2K                              # air temperature from C to K
    Ah,f = qcutils.GetSeriesasMA(ds,Ah_in)
    ps,f = qcutils.GetSeriesasMA(ds,ps_in)
    rhod,f = qcutils.GetSeriesasMA(ds,'rhod')     # density dry air
    rhom,f = qcutils.GetSeriesasMA(ds,'rhom')     # density moist air
    RhoCp,f = qcutils.GetSeriesasMA(ds,'RhoCp')
    Lv,f = qcutils.GetSeriesasMA(ds,'Lv')
    Ah = Ah*c.g2kg                                # absolute humidity from g/m3 to kg/m3
    sigma = Ah/rhod
    h2o_wpl_Fe = c.mu*sigma*Fe_raw
    h2o_wpl_Fh = (1+c.mu*sigma)*Ah*Lv*(Fh/RhoCp)/TaK
    Fe_wpl_data = Fe_raw+h2o_wpl_Fe+h2o_wpl_Fh
    Fe_wpl_flag = numpy.zeros(len(Fe_wpl_data))
    mask = numpy.ma.getmask(Fe_wpl_data)
    index = numpy.where(Fe_wpl_data.mask==True)[0]
    Fe_wpl_flag[index] = numpy.int32(14)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL corrected Fe',units='W/m2',standard_name='surface_upward_latent_heat_flux')
    qcutils.CreateSeries(ds,Fe_wpl_out,Fe_wpl_data,Flag=Fe_wpl_flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Fe (uncorrected for WPL)',units='W/m2')
    qcutils.CreateSeries(ds,'Fe_raw',Fe_raw,Flag=Fe_raw_flag,Attr=attr)
    if qcutils.cfoptionskey(cf,Key='RelaxFeWPL'):
        ReplaceWhereMissing(ds.series['Fe'],ds.series['Fe'],ds.series['Fe_raw'],FlagValue=20)
        if 'RelaxFeWPL' not in ds.globalattributes['Functions']:
            ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', RelaxFeWPL'

def FhvtoFh(cf,ds,Fh_out='Fh',Fhv_in='Fhv',Tv_in='Tv_CSAT',q_in='q',wA_in='wA',wT_in='wT'):
    '''
    Convert the virtual heat flux to the sensible heat flux.
    USEAGE:
     qcts.FhvtoFh_EP(cf,ds,Fhv_in='Fhv',RhoCp_in='RhoCp',Tv_in='Tv_CSAT',wA_in='wA',rhom_in='rhom',q_in='q',wT_in='wT')
    INPUT:
     All inputs are read from the data structure.
      Fhv_in   - label of the virtual heat flux series, default is 'Fhv'
      RhoCp_in - label of the RhoCp series, default is 'RhoCp'
      Tv_in    - label of the virtual temperature series, default is 'Tv_CSAT'
      wA_in    - label of the wA covariance series, default is 'wA'
      rhom_in  - label of the moist air density series, default is 'rhom'
      q_in     - label of the specific humidity series, default is 'q'
      wT_in    - label of the wT covariance series, default is 'wT'
    OUTPUT:
     All outputs are written to the data structure.
      Fh_out   - label of sensible heat flux, default is 'Fh'
    '''
    log.info(' Converting virtual Fh to Fh')
    # get the input series
    Fhv,Fhv_f = qcutils.GetSeriesasMA(ds,Fhv_in)                      # get the virtual heat flux
    Tv,Tv_f = qcutils.GetSeriesasMA(ds,Tv_in)                         # get the virtual temperature, C
    TvK = Tv + c.C2K                                                  # convert from C to K
    wA,wA_f = qcutils.GetSeriesasMA(ds,wA_in)                         # get the wA covariance, g/m2/s
    wA = wA * c.g2kg                                                  # convert from g/m2/s to kg/m2/s
    q,q_f = qcutils.GetSeriesasMA(ds,q_in)                            # get the specific humidity, kg/kg
    wT,wT_f = qcutils.GetSeriesasMA(ds,wT_in)                         # get the wT covariance, mK/s
    # get the utility series
    RhoCp,RhoCp_f = qcutils.GetSeriesasMA(ds,'RhoCp')                 # get rho*Cp
    rhom,rhom_f = qcutils.GetSeriesasMA(ds,'rhom')                    # get the moist air density, kg/m3
    # define local constants
    alpha = 0.51
    # do the conversion
    Fh = Fhv - RhoCp*alpha*TvK*wA/rhom - RhoCp*alpha*q*wT
    # put the calculated sensible heat flux into the data structure
    attr = qcutils.MakeAttributeDictionary(long_name='Sensible heat flux from virtual heat flux',
                                           units='W/m2',standard_name='surface_upward_sensible_heat_flux')
    qcutils.CreateSeries(ds,Fh_out,Fh,FList=[Fhv_in,Tv_in,wA_in,q_in,wT_in],Attr=attr)
    if 'FhvtoFh' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', FhvtoFh'
    if qcutils.cfoptionskey(cf,Key='RelaxFhvtoFh'):
        ReplaceWhereMissing(ds.series['Fh'],ds.series['Fh'],ds.series['Fhv'],FlagValue=20)
        if 'RelaxFhvtoFh' not in ds.globalattributes['Functions']:
            ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', RelaxFhvtoFh'

def FilterUstar(cf,ds,ustar_in='ustar',ustar_out='ustar_filtered'):
    """
    Filter ustar for low turbulence periods.  The filtering is done by checking the
    friction velocity for each time period.  If ustar is less than or equal to the
    threshold specified in the control file then ustar is set to missing.  If
    the ustar is greater than the threshold, no action is taken.  Filtering is not
    done "in place", a new series is created with the label given in the control file.
    The QC flag is set to 18 to indicate the missing low ustar values.
    
    Usage: qcts.FilterUstar(cf,ds)
    cf: control file object
    ds: data structure object
    """
    if ustar_out not in cf['Variables'].keys(): return
    if 'ustar_threshold' in cf['Variables'][ustar_out].keys():
        log.info(' Filtering ustar to remove values below threshold')
        ustar_threshold = float(cf['Variables'][ustar_out]['ustar_threshold'])
        ustar, ustar_flag = qcutils.GetSeriesasMA(ds,ustar_in)
        index = numpy.ma.where(ustar<=ustar_threshold)[0]
        ustar = numpy.ma.masked_where(ustar<=ustar_threshold,ustar)
        ustar_flag[index] = 18
        descr = 'ustar filtered for low turbulence conditions (<'+str(ustar_threshold)+')'
        units = qcutils.GetUnitsFromds(ds, ustar_in)
        attr = qcutils.MakeAttributeDictionary(long_name=descr,units=units)
        qcutils.CreateSeries(ds,ustar_out,ustar,Flag=ustar_flag,Attr=attr)
    else:
        log.error(' ustar threshold (ustar_threshold) not found in '+ustar_out+' section of control file')

def get_averages(Data):
    """
        Get daily averages on days when no 30-min observations are missing.
        Days with missing observations return a value of -9999
        Values returned are sample size (Num) and average (Av)
        
        Usage qcts.get_averages(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-float(-9999))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Av = -9999
    elif Num == 48:
        Av = numpy.ma.mean(Data[li])
    else:
        x = 0
        index = numpy.ma.where(Data.mask == True)[0]
        if len(index) == 1:
            x = 1
        elif len(index) > 1:
            for i in range(len(Data)):
                if Data.mask[i] == True:
                    x = x + 1
        
        if x == 0:
            Av = numpy.ma.mean(Data[li])
        else:
            Av = -9999
    return Num, Av

def get_minmax(Data):
    """
        Get daily minima and maxima on days when no 30-min observations are missing.
        Days with missing observations return a value of -9999
        Values returned are sample size (Num), minimum (Min) and maximum (Max)
        
        Usage qcts.get_minmax(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-float(-9999))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Min = -9999
        Max = -9999
    elif Num == 48:
        Min = numpy.ma.min(Data[li])
        Max = numpy.ma.max(Data[li])
    else:
        x = 0
        index = numpy.ma.where(Data.mask == True)[0]
        if len(index) == 1:
            x = 1
        elif len(index) > 1:
            for i in range(len(Data)):
                if Data.mask[i] == True:
                    x = x + 1
        
        if x == 0:
            Min = numpy.ma.min(Data[li])
            Max = numpy.ma.max(Data[li])
        else:
            Min = -9999
            Max = -9999
    return Num, Min, Max

def get_nightsums(Data):
    """
        Get nightly sums and averages on nights when no 30-min observations are missing.
        Nights with missing observations return a value of -9999
        Values returned are sample size (Num), sums (Sum) and average (Av)
        
        Usage qcts.get_nightsums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(Data.mask == False)[0]
    Num = numpy.size(li)
    if Num == 0:
        Sum = -9999
        Av = -9999
    else:
        x = 0
        for i in range(len(Data)):
            if Data.mask[i] == True:
                x = x + 1
        
        if x == 0:
            Sum = numpy.ma.sum(Data[li])
            Av = numpy.ma.mean(Data[li])
        else:
            Sum = -9999
            Av = -9999
    
    return Num, Sum, Av

def get_soilaverages(Data):
    """
        Get daily averages of soil water content on days when 15 or fewer 30-min observations are missing.
        Days with 16 or more missing observations return a value of -9999
        Values returned are sample size (Num) and average (Av)
        
        Usage qcts.get_soilaverages(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-float(-9999))>c.eps)
    Num = numpy.size(li)
    if Num > 33:
        Av = numpy.ma.mean(Data[li])
    else:
        Av = -9999
    return Num, Av

def get_subsums(Data):
    """
        Get separate daily sums of positive and negative fluxes when no 30-min observations are missing.
        Days with missing observations return a value of -9999
        Values returned are positive and negative sample sizes (PosNum and NegNum) and sums (SumPos and SumNeg)
        
        Usage qcts.get_subsums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-float(-9999))>c.eps)
    Num = numpy.size(li)
    if Num == 48:
        pi = numpy.ma.where(Data[li]>0)
        ni = numpy.ma.where(Data[li]<0)
        PosNum = numpy.size(pi)
        NegNum = numpy.size(ni)
        if PosNum > 0:
            SumPos = numpy.ma.sum(Data[pi])
        else:
            SumPos = 0
        if NegNum > 0:
            SumNeg = numpy.ma.sum(Data[ni])
        else:
            SumNeg = 0
    else:
        pi = numpy.ma.where(Data[li]>0)
        ni = numpy.ma.where(Data[li]<0)
        PosNum = numpy.size(pi)
        NegNum = numpy.size(ni)
        SumPos = -9999
        SumNeg = -9999
    return PosNum, NegNum, SumPos, SumNeg

def get_sums(Data):
    """
        Get daily sums when no 30-min observations are missing.
        Days with missing observations return a value of -9999
        Values returned are sample size (Num) and sum (Sum)
        
        Usage qcts.get_sums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-float(-9999))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Sum = -9999
    elif Num == 48:
        Sum = numpy.ma.sum(Data[li])
    else:
        x = 0
        index = numpy.ma.where(Data.mask == True)[0]
        if len(index) == 1:
            x = 1
        elif len(index) > 1:
            for i in range(len(Data)):
                if Data.mask[i] == True:
                    x = x + 1
        
        if x == 0:
            Sum = numpy.ma.sum(Data[li])
        else:
            Sum = -9999
    return Num, Sum

def get_qcflag(ds):
    """
        Set up flags during ingest of L1 data.
        Identifies missing observations as -9999 and sets flag value 1
        
        Usage qcts.get_qcflag(ds)
        ds: data structure
        """
    log.info(' Setting up the QC flags')
    nRecs = len(ds.series['xlDateTime']['Data'])
    for ThisOne in ds.series.keys():
        ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,dtype=numpy.int32)
        index = numpy.where(ds.series[ThisOne]['Data']==-9999)[0]
        ds.series[ThisOne]['Flag'][index] = numpy.int32(1)

def InvertSign(ds,ThisOne):
    log.info(' Inverting sign of '+ThisOne)
    index = numpy.where(abs(ds.series[ThisOne]['Data']-float(-9999))>c.eps)[0]
    ds.series[ThisOne]['Data'][index] = float(-1)*ds.series[ThisOne]['Data'][index]

def InterpolateOverMissing(cf,ds,series='',maxlen=1000):
    if series not in ds.series.keys(): return
    section = qcutils.get_cfsection(cf,series=series,mode='quiet')
    if len(section)==0: return
    DateNum = date2num(ds.series['DateTime']['Data'])
    iog = numpy.where(ds.series[series]['Data']!=float(-9999))[0]            # index of good values
    if len(iog)<2:
        log.info(' InterpolateOverMissing: Less than 2 good points available for series '+str(series))
        return
    f = interpolate.interp1d(DateNum[iog],ds.series[series]['Data'][iog])    # linear interpolation function
    iom = numpy.where((ds.series[series]['Data']==float(-9999))&             # index of missing values
                      (DateNum>=DateNum[iog[0]])&                          # that occur between the first
                      (DateNum<=DateNum[iog[-1]]))[0]                      # and last dates used to define f
    # Now we step through the indices of the missing values and discard
    # contiguous blocks longer than maxlen.
    # !!! The following code is klunky and could be re-written to be
    # !!! neater and faster.
    # First, define 2 temporary arrays used and initialise 2 counters.
    tmp1 = numpy.zeros(len(iom),int)
    tmp2 = numpy.zeros(len(iom),int)
    k=0
    n=0
    # step through the array of idices for missing values
    for i in range(len(iom)-1):
        dn = iom[i+1]-iom[i]        # change in index number from one element of iom to the next
        if dn==1:                   # if the change is 1 then we are still in a contiguous block
            tmp1[n] = iom[i]        # save the index into a temporary array
            n = n + 1               # increment the contiguous block length counter
        elif dn>1:                  # if the change is greater than 1 then we have come to the end of a contiguous block
            if n<maxlen:            # if the contiguous block length is less then maxlen
                tmp1[n]=iom[i]      # save the last index of the contiguous block
                tmp2[k:k+n+1] = tmp1[0:n+1]   # concatenate the indices for this block to any previous block with length less than maxlen
                k=k+n+1             # update the pointer to the concatenating array
            n=0                     # reset the contiguous block length counter
    if k>0:                         # do the interpolation only if 1 gap is less than maxlen
        tmp2[k] = iom[-1]               # just accept the last missing value index regardless
        iom_new = tmp2[:k+1]            # the array of missing data indices with contiguous block lengths less than maxlen
        ds.series[series]['Data'][iom_new] = f(DateNum[iom_new]).astype(numpy.float32)        # fill missing values with linear interpolations
        ds.series[series]['Flag'][iom_new] = numpy.int32(50)
    if 'InterpolateOverMissing' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', InterpolateOverMissing'

def MassmanStandard(cf,ds,Ta_in='Ta',Ah_in='Ah',ps_in='ps',ustar_in='ustar',ustar_out='ustar',L_in='L',L_out ='L',uw_out='uw',vw_out='vw',wT_out='wT',wA_out='wA',wC_out='wC'):
    """
       Massman corrections.
       The steps involved are as follows:
        1) calculate ustar and L using rotated but otherwise uncorrected covariances
       """
    if not qcutils.cfoptionskey(cf,Key='MassmanCorrection'): return
    if 'Massman' not in cf:
        log.info(' Massman section not found in control file, no corrections applied')
        return
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='MassmanVars'):
        MArgs = ast.literal_eval(cf['FunctionArgs']['MassmanVars'])
        Ta_in = MArgs[0]
        Ah_in = MArgs[1]
        ps_in = MArgs[2]
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='MassmanOuts'):
        MOut = ast.literal_eval(cf['FunctionArgs']['MassmanOuts'])
        ustar_in = MOut[0]
        ustar_out = MOut[1]
        L_in = MOut[2]
        L_out = MOut[3]
        uw_out = MOut[4]
        vw_out = MOut[5]
        wT_out = MOut[6]
        wA_out = MOut[7]
        wC_out = MOut[8]
    log.info(' Correcting for flux loss from spectral attenuation')
    zmd = float(cf['Massman']['zmd'])             # z-d for site
    angle = float(cf['Massman']['angle'])         # CSAT3-IRGA separation angle
    CSATarm = float(cf['Massman']['CSATarm'])     # CSAT3 mounting distance
    IRGAarm = float(cf['Massman']['IRGAarm'])     # IRGA mounting distance
    lLat = numpy.ma.sin(numpy.deg2rad(angle)) * IRGAarm
    lLong = CSATarm - (numpy.ma.cos(numpy.deg2rad(angle)) * IRGAarm)
    # *** Massman_1stpass starts here ***
    #  The code for the first and second passes is very similar.  It would be useful to make them the
    #  same and put into a loop to reduce the nu,ber of lines in this function.
    # calculate ustar and Monin-Obukhov length from rotated but otherwise uncorrected covariances
    Ta,f = qcutils.GetSeriesasMA(ds,Ta_in)
    Ah,f = qcutils.GetSeriesasMA(ds,Ah_in)
    ps,f = qcutils.GetSeriesasMA(ds,ps_in)
    nRecs = numpy.size(Ta)
    u,f = qcutils.GetSeriesasMA(ds,'u')
    uw,f = qcutils.GetSeriesasMA(ds,'uw')
    vw,f = qcutils.GetSeriesasMA(ds,'vw')
    wT,f = qcutils.GetSeriesasMA(ds,'wT')
    wC,f = qcutils.GetSeriesasMA(ds,'wC')
    wA,f = qcutils.GetSeriesasMA(ds,'wA')
    if ustar_in not in ds.series.keys():
        ustarm = numpy.ma.sqrt(numpy.ma.sqrt(uw ** 2 + vw ** 2))
    else:
        ustarm,f = qcutils.GetSeriesasMA(ds,ustar_in)
    if L_in not in ds.series.keys():
        Lm = mf.molen(Ta, Ah, ps, ustarm, wT, fluxtype='kinematic')
    else:
        Lm,f = qcutils.GetSeriesasMA(ds,Lm_in)
    # now calculate z on L
    zoLm = zmd / Lm
    # start calculating the correction coefficients for approximate corrections
    #  create nxMom, nxScalar and alpha series with their unstable values by default
    nxMom, nxScalar, alpha = qcutils.nxMom_nxScalar_alpha(zoLm)
    # now calculate the fxMom and fxScalar coefficients
    fxMom = nxMom * u / zmd
    fxScalar = nxScalar * u / zmd
    # compute spectral filters
    tao_eMom = ((c.lwVert / (5.7 * u)) ** 2) + ((c.lwHor / (2.8 * u)) ** 2)
    tao_ewT = ((c.lwVert / (8.4 * u)) ** 2) + ((c.lTv / (4.0 * u)) ** 2)
    tao_ewIRGA = ((c.lwVert / (8.4 * u)) ** 2) + ((c.lIRGA / (4.0 * u)) ** 2) \
                 + ((lLat / (1.1 * u)) ** 2) + ((lLong / (1.05 * u)) ** 2)
    tao_b = c.Tb / 2.8
    # calculate coefficients
    bMom = qcutils.bp(fxMom,tao_b)
    bScalar = qcutils.bp(fxScalar,tao_b)
    pMom = qcutils.bp(fxMom,tao_eMom)
    pwT = qcutils.bp(fxScalar,tao_ewT)
    # calculate corrections for momentum and scalars
    rMom = qcutils.r(bMom, pMom, alpha)        # I suspect that rMom and rwT are the same functions
    rwT = qcutils.r(bScalar, pwT, alpha)
    # determine approximately-true Massman fluxes
    uwm = uw / rMom
    vwm = vw / rMom
    wTm = wT / rwT
    # *** Massman_1stpass ends here ***
    # *** Massman_2ndpass starts here ***
    # we have calculated the first pass corrected momentum and temperature covariances, now we use
    # these to calculate the final corrections
    #  first, get the 2nd pass corrected friction velocity and Monin-Obukhov length
    ustarm = numpy.ma.sqrt(numpy.ma.sqrt(uwm ** 2 + vwm ** 2))
    Lm = mf.molen(Ta, Ah, ps, ustarm, wTm, fluxtype='kinematic')
    zoLm = zmd / Lm
    nxMom, nxScalar, alpha = qcutils.nxMom_nxScalar_alpha(zoLm)
    fxMom = nxMom * (u / zmd)
    fxScalar = nxScalar * (u / zmd)
    # calculate coefficients
    bMom = qcutils.bp(fxMom,tao_b)
    bScalar = qcutils.bp(fxScalar,tao_b)
    pMom = qcutils.bp(fxMom,tao_eMom)
    pwT = qcutils.bp(fxScalar,tao_ewT)
    pwIRGA = qcutils.bp(fxScalar,tao_ewIRGA)
    # calculate corrections for momentum and scalars
    rMom = qcutils.r(bMom, pMom, alpha)
    rwT = qcutils.r(bScalar, pwT, alpha)
    rwIRGA = qcutils.r(bScalar, pwIRGA, alpha)
    # determine true fluxes
    uwM = uw / rMom
    vwM = vw / rMom
    wTM = wT / rwT
    wCM = wC / rwIRGA
    wAM = wA / rwIRGA
    ustarM = numpy.ma.sqrt(numpy.ma.sqrt(uwM ** 2 + vwM ** 2))
    LM = mf.molen(Ta, Ah, ps, ustarM, wTM, fluxtype='kinematic')
    # write the 2nd pass Massman corrected covariances to the data structure
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true ustar',units='m/s')
    qcutils.CreateSeries(ds,ustar_out,ustarM,FList=['uw','vw'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Obukhov Length',units='m')
    qcutils.CreateSeries(ds,L_out,LM,FList=[Ta_in,Ah_in,ps_in,'wT'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Cov(uw)',units='m2/s2')
    qcutils.CreateSeries(ds,uw_out,uwM,FList=['uw',L_out],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Cov(vw)',units='m2/s2')
    qcutils.CreateSeries(ds,vw_out,vwM,FList=['vw',L_out],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Cov(wT)',units='mC/s')
    qcutils.CreateSeries(ds,wT_out,wTM,FList=['wT',L_out],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Cov(wA)',units='g/m2/s')
    qcutils.CreateSeries(ds,wA_out,wAM,FList=['wA',L_out],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Cov(wC)',units='mg/m2/s')
    qcutils.CreateSeries(ds,wC_out,wCM,FList=['wC',L_out],Attr=attr)
    # *** Massman_2ndpass ends here ***
    
    if qcutils.cfkeycheck(cf,Base='General',ThisOne='MassmanFlag') and cf['General']['MassmanFlag'] == 'True':
        keys = [ustar_out,L_out,uw_out,vw_out,wT_out,wA_out,wC_out]
        for ThisOne in keys:
            testseries,f = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where(mask.astype(int)==1)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(12)

def MergeSeries(cf,ds,series,okflags):
    """
        Merge two series of data to produce one series containing the best data from both.
        Calling syntax is: MergeSeries(cf,ds,series,okflags)
         where ds is the data structure containing all series
               series (str) is the label of the destination series
               okflags (list) is a list of QC flag values for which the data is considered acceptable
        If the QC flag for Primary is in okflags, the value from Primary is placed in destination.
        If the QC flag for Primary is not in okflags but the QC flag for Secondary is, the value
        from Secondary is placed in Destination.
        """
    # check to see if the series is specified in the control file
    section = qcutils.get_cfsection(cf,series=series)
    if len(section)==0: return
    # check to see if the entry for series in the control file has the MergeSeries key
    if 'MergeSeries' not in cf[section][series].keys(): return
    # check to see if the series has already been merged
    if series in ds.mergeserieslist: return
    # now get the source list and the standard name
    srclist, standardname = qcutils.GetMergeSeriesKeys(cf,series,section=section)
    nSeries = len(srclist)
    if nSeries==0:
        log.info(' MergeSeries: no input series specified for '+str(series))
        return
    if nSeries==1:
        log.info(' Merging series '+str(srclist)+' into '+series)
        if srclist[0] not in ds.series.keys():
            log.error('  MergeSeries: primary input series'+srclist[0]+'not found for'+str(series))
            return
        data = ds.series[srclist[0]]['Data'].copy()
        flag = ds.series[srclist[0]]['Flag'].copy()
        SeriesNameString = srclist[0]
        SeriesUnitString = ds.series[srclist[0]]['Attr']['units']
    else:
        log.info(' Merging series '+str(srclist)+' into '+series)
        if srclist[0] not in ds.series.keys():
            log.error('  MergeSeries: primary input series'+srclist[0]+'not found')
            return
        data = ds.series[srclist[0]]['Data'].copy()
        flag = ds.series[srclist[0]]['Flag'].copy()
        SeriesNameString = srclist[0]
        SeriesUnitString = ds.series[srclist[0]]['Attr']['units']
        srclist.remove(srclist[0])
        for ThisOne in srclist:
            if ThisOne in ds.series.keys():
                SeriesNameString = SeriesNameString+', '+ThisOne
                indx1 = numpy.zeros(numpy.size(data),dtype=numpy.int)
                indx2 = numpy.zeros(numpy.size(data),dtype=numpy.int)
                for okflag in okflags:
                    index = numpy.where((flag==okflag))[0]                             # index of acceptable primary values
                    indx1[index] = 1                                                   # set primary index to 1 when primary good
                    index = numpy.where((ds.series[ThisOne]['Flag']==okflag))[0]       # same process for secondary
                    indx2[index] = 1
                index = numpy.where((indx1!=1)&(indx2==1))[0]           # index where primary bad but secondary good
                data[index] = ds.series[ThisOne]['Data'][index]         # replace bad primary with good secondary
                flag[index] = ds.series[ThisOne]['Flag'][index]
            else:
                log.error('  MergeSeries: secondary input series'+ThisOne+'not found')
    ds.mergeserieslist.append(series)
    attr = qcutils.MakeAttributeDictionary(long_name='Merged from '+SeriesNameString,
                             standard_name=standardname,units=SeriesUnitString)
    qcutils.CreateSeries(ds,series,data,Flag=flag,Attr=attr)

def PT100(ds,T_out,R_in,m):
    log.info(' Calculating temperature from PT100 resistance')
    R,f = qcutils.GetSeriesasMA(ds,R_in)
    R = m*R
    T = (-c.PT100_alpha+numpy.sqrt(c.PT100_alpha**2-4*c.PT100_beta*(-R/100+1)))/(2*c.PT100_beta)
    attr = qcutils.MakeAttributeDictionary(long_name='Calculated PT100 temperature using '+str(R_in),units='degC')
    qcutils.CreateSeries(ds,T_out,T,FList=[R_in],Attr=attr)

def ReplaceRotatedCovariance(cf,ds,rot_cov_label,non_cov_label):
    log.info(' Replacing missing '+rot_cov_label+' when '+non_cov_label+' is good')
    cr_data,cr_flag = qcutils.GetSeriesasMA(ds,rot_cov_label)
    cn_data,cn_flag = qcutils.GetSeriesasMA(ds,non_cov_label)
    index = numpy.ma.where((cr_data.mask==True)&(cn_data.mask==False))[0]
    if len(index)!=0:
        ds.series[rot_cov_label]['Data'][index] = cn_data[index]
        ds.series[rot_cov_label]['Flag'][index] = numpy.int32(20)
    return

def ReplaceOnDiff(cf,ds,series=''):
    # Gap fill using data from alternate sites specified in the control file
    ts = ds.globalattributes['time_step']
    if len(series)!=0:
        ds_alt = {}                     # create a dictionary for the data from alternate sites
        open_ncfiles = []               # create an empty list of open netCDF files
        for ThisOne in series:          # loop over variables in the series list
            # has ReplaceOnDiff been specified for this series?
            if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'ReplaceOnDiff'):
                # loop over all entries in the ReplaceOnDiff section
                for Alt in cf['Variables'][ThisOne]['ReplaceOnDiff'].keys():
                    if 'FileName' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                        alt_filename = cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['FileName']
                        if 'AltVarName' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            alt_varname = cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['AltVarName']
                        else:
                            alt_varname = ThisOne
                        if alt_filename not in open_ncfiles:
                            n = len(open_ncfiles)
                            open_ncfiles.append(alt_filename)
                            ds_alt[n] = qcio.nc_read_series_file(alt_filename)
                        else:
                            n = open_ncfiles.index(alt_filename)
                        if 'Transform' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            AltDateTime = ds_alt[n].series['DateTime']['Data']
                            AltSeriesData = ds_alt[n].series[alt_varname]['Data']
                            TList = ast.literal_eval(cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['Transform'])
                            for TListEntry in TList:
                                qcts.TransformAlternate(TListEntry,AltDateTime,AltSeriesData,ts=ts)
                        if 'Range' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            RList = ast.literal_eval(cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['Range'])
                            for RListEntry in RList:
                                qcts.ReplaceWhenDiffExceedsRange(ds.series['DateTime']['Data'],ds.series[ThisOne],
                                                                 ds.series[ThisOne],ds_alt[n].series[alt_varname],
                                                                 RListEntry)
                    elif 'AltVarName' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                        alt_varname = ThisOne
                        if 'Range' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            RList = ast.literal_eval(cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['Range'])
                            for RListEntry in RList:
                                qcts.ReplaceWhenDiffExceedsRange(ds.series['DateTime']['Data'],ds.series[ThisOne],
                                                                 ds.series[ThisOne],ds.series[alt_varname],
                                                                 RListEntry)
                    else:
                        log.error('ReplaceOnDiff: Neither AltFileName nor AltVarName given in control file')
    else:
        log.error('ReplaceOnDiff: No input series specified')

def ReplaceWhereMissing(Destination,Primary,Secondary,FlagOffset=None,FlagValue=None):
    #print time.strftime('%X')+' Merging series '+Primary+' and '+Secondary+' into '+Destination
    p_data = Primary['Data'].copy()
    p_flag = Primary['Flag'].copy()
    s_data = Secondary['Data'].copy()
    s_flag = Secondary['Flag'].copy()
    if numpy.size(p_data)>numpy.size(s_data):
        p_data = p_data[0:numpy.size(s_data)]
    if numpy.size(s_data)>numpy.size(p_data):
        s_data = s_data[0:numpy.size(p_data)]
    index = numpy.where((abs(p_data-float(-9999))<c.eps)&
                        (abs(s_data-float(-9999))>c.eps))[0]
    p_data[index] = s_data[index]
    if FlagValue==None and FlagOffset!=None:
        p_flag[index] = s_flag[index] + numpy.int32(FlagOffset)
    elif FlagValue!=None and FlagOffset==None:
        p_flag[index] = numpy.int32(FlagValue)
    else:
        p_flag[index] = s_flag[index]
    Destination['Data'] = Primary['Data'].copy()
    Destination['Flag'] = Primary['Flag'].copy()
    Destination['Data'][0:len(p_data)] = p_data
    Destination['Flag'][0:len(p_flag)] = p_flag
    Destination['Attr']['long_name'] = 'Merged from original and alternate'
    Destination['Attr']['units'] = Primary['Attr']['units']

def ReplaceWhenDiffExceedsRange(DateTime,Destination,Primary,Secondary,RList):
    #print time.strftime('%X')+' Replacing '+Primary+' with '+Secondary+' when difference exceeds threshold'
    # get the primary data series
    p_data = numpy.ma.array(Primary['Data'])
    p_flag = Primary['Flag'].copy()
    # get the secondary data series
    s_data = numpy.ma.array(Secondary['Data'])
    s_flag = Secondary['Flag'].copy()
    # truncate the longest series if the sizes do not match
    if numpy.size(p_data)!=numpy.size(s_data):
        log.warning(' ReplaceWhenDiffExceedsRange: Series lengths differ, longest will be truncated')
        if numpy.size(p_data)>numpy.size(s_data):
            p_data = p_data[0:numpy.size(s_data)]
        if numpy.size(s_data)>numpy.size(p_data):
            s_data = s_data[0:numpy.size(p_data)]
    # get the difference between the two data series
    d_data = p_data-s_data
    # normalise the difference if requested
    if RList[3]=='s':
        d_data = (p_data-s_data)/s_data
    elif RList[3]=='p':
        d_data = (p_data-s_data)/p_data
    #si = qcutils.GetDateIndex(DateTime,RList[0],0)
    #ei = qcutils.GetDateIndex(DateTime,RList[1],0)
    Range = RList[2]
    Upper = float(Range[0])
    Lower = float(Range[1])
    index = numpy.ma.where((abs(d_data)<Lower)|(abs(d_data)>Upper))
    p_data[index] = s_data[index]
    p_flag[index] = 35
    Destination['Data'] = numpy.ma.filled(p_data,float(-9999))
    Destination['Flag'] = p_flag.copy()
    Destination['Attr']['long_name'] = 'Replaced original with alternate when difference exceeded threshold'
    Destination['Attr']['units'] = Primary['Attr']['units']

def savitzky_golay(y, window_size, order, deriv=0):
    ''' Apply Savitsky-Golay low-pass filter to data.'''
    try:
        window_size = numpy.abs(numpy.int(window_size))
        order = numpy.abs(numpy.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = numpy.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = numpy.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - numpy.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + numpy.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = numpy.concatenate((firstvals, y, lastvals))
    return numpy.convolve( m, y, mode='valid')

def Square(Series):
    tmp = numpy.array([-9999]*numpy.size(Series),Series.dtype)
    index = numpy.where(Series!=float(-9999))[0]
    tmp[index] = Series[index] ** 2
    return tmp

def SquareRoot(Series):
    tmp = numpy.array([-9999]*numpy.size(Series),Series.dtype)
    index = numpy.where(Series!=float(-9999))[0]
    tmp[index] = Series[index] ** .5
    return tmp

def TaFromTv(cf,ds,Ta_out='Ta_CSAT',Tv_in='Tv_CSAT',Ah_in='Ah',ps_in='ps'):
    # Calculate the air temperature from the virtual temperature, the
    # absolute humidity and the pressure.
    # NOTE: the virtual temperature is used in place of the air temperature
    #       to calculate the vapour pressure from the absolute humidity, the
    #       approximation involved here is of the order of 1%.
    log.info(' Calculating Ta from Tv')
    Tv,f = qcutils.GetSeriesasMA(ds,Tv_in)
    Ah,f = qcutils.GetSeriesasMA(ds,Ah_in)
    ps,f = qcutils.GetSeriesasMA(ds,ps_in)
    nRecs = int(ds.globalattributes['nc_nrecs'])
    Ta_flag = numpy.zeros(nRecs,numpy.int32)
    vp = mf.vapourpressure(Ah,Tv)
    mr = mf.mixingratio(ps,vp)
    q = mf.specifichumidity(mr)
    Ta_data = mf.tafromtv(Tv,q)
    mask = numpy.ma.getmask(Ta_data)
    index = numpy.where(mask.astype(numpy.int32)==1)
    Ta_flag[index] = 15
    attr = qcutils.MakeAttributeDictionary(long_name='Ta calculated from Tv using '+Tv_in,units='C',standard_name='air_temperature')
    qcutils.CreateSeries(ds,Ta_out,Ta_data,Flag=Ta_flag,Attr=attr)
    if 'TaFromTv' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', TaFromTv'

def TransformAlternate(TList,DateTime,Series,ts=30):
    # Apply polynomial transform to data series being used as replacement data for gap filling
    #print time.strftime('%X')+' Applying polynomial transform to '+ThisOne
    si = qcutils.GetDateIndex(DateTime,TList[0],ts=ts,default=0,match='exact')
    ei = qcutils.GetDateIndex(DateTime,TList[1],ts=ts,default=-1,match='exact')
    Series = numpy.ma.masked_where(abs(Series-float(-9999))<c.eps,Series)
    Series[si:ei] = qcutils.polyval(TList[2],Series[si:ei])
    Series = numpy.ma.filled(Series,float(-9999))

def UstarFromFh(cf,ds,us_out='uscalc',T_in='Ta', Ah_in='Ah', p_in='ps', Fh_in='Fh', u_in='Ws_CSAT', us_in='ustar'):
    # Calculate ustar from sensible heat flux, wind speed and
    # roughness length using Wegstein's iterative method.
    #  T is the air temperature, C
    #  p is the atmospheric pressure, kPa
    #  H is the sensible heat flux, W/m^2
    #  u is the wind speed, m/s
    #  z is the measurement height minus the displacement height, m
    #  z0 is the momentum roughness length, m
    log.info(' Calculating ustar from (Fh,Ta,Ah,p,u)')
    # get z-d (measurement height minus displacement height) and z0 from the control file
    if qcutils.cfkeycheck(cf,Base='Params',ThisOne='zmd') and qcutils.cfkeycheck(cf,Base='Params',ThisOne='z0'):
        zmd = float(cf['Params']['zmd'])   # z-d for site
        z0 = float(cf['Params']['z0'])     # z0 for site
    else:
        log.error('Parameters zmd or z0 not found in control file.  u* not determined from Fh')
        return
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='ustarFh'):
        args = ast.literal_eval(cf['FunctionArgs']['ustarFh'])
        us_out = args[0]
        T_in = args[1]
        Ah_in = args[2]
        p_in = args[3]
        Fh_in = args[4]
        u_in = args[5]
        us_in = args[6]
    T,T_flag = qcutils.GetSeries(ds,T_in)
    Ah,Ah_flag = qcutils.GetSeries(ds,Ah_in)
    p,p_flag = qcutils.GetSeries(ds,p_in)
    Fh,Fh_flag = qcutils.GetSeries(ds,Fh_in)
    u,u_flag = qcutils.GetSeries(ds,u_in)
    nRecs = numpy.size(Fh)
    us = numpy.zeros(nRecs,dtype=numpy.float64) + numpy.float64(-9999)
    us_flag = numpy.zeros(nRecs,dtype=numpy.int)
    for i in range(nRecs):
        if((abs(T[i]-float(-9999))>c.eps)&(abs(Ah[i]-float(-9999))>c.eps)&
           (abs(p[i]-float(-9999))>c.eps)&(abs(Fh[i]-float(-9999))>c.eps)&
           (abs(u[i]-float(-9999))>c.eps)):
            #print ds.series['DateTime']['Data'][i],T[i]
            us[i] = qcutils.Wegstein(T[i], Ah[i], p[i], Fh[i], u[i], z, z0)
            us_flag[i] = 36
        else:
            us[i] = numpy.float64(-9999)
            us_flag[i] = 37
    attr = qcutils.MakeAttributeDictionary(long_name='ustar from (Fh,Ta,Ah,p,u)',units='m/s')
    qcutils.CreateSeries(ds,us_out,us,Flag=us_flag,Attr=attr)
    return us_in, us_out

def write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoSum='False',DoMinMax='False',DoMean='False',DoSubSum='False',DoSoil='False'):
    monthabr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    if qcutils.cfkeycheck(cf,Base='Params',ThisOne='firstMonth'):
        M1st = int(cf['Params']['firstMonth'])
    else:
        M1st = 1
    if qcutils.cfkeycheck(cf,Base='Params',ThisOne='secondMonth'):
        M2nd = int(cf['Params']['secondMonth'])
    else:
        M2nd = 12
    log.info(' Doing daily sums for '+ThisOne)
    Units = ds.series[ThisOne]['Attr']['units']
    
    xlRow = 1
    if xlCol == 0:
        xlSheet.write(xlRow,xlCol,'Month')
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,'Day')
        xlCol = xlCol + 1
    xlSheet.write(xlRow,xlCol,'n')
    xlCol = xlCol + 1
    if DoMinMax == 'True':
        xlSheet.write(xlRow,xlCol,ThisOne+'_min')
        xlSheet.write(xlRow-1,xlCol,Units)
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,ThisOne+'_max')
        if DoMean == 'True':
            xlSheet.write(xlRow-1,xlCol,Units)
            xlCol = xlCol + 1
            xlSheet.write(xlRow,xlCol,ThisOne+'_mean')
    elif DoMinMax == 'False' and DoMean == 'True':
        xlSheet.write(xlRow,xlCol,ThisOne+'_mean')
    elif DoMinMax == 'False' and DoMean == 'False':
        xlSheet.write(xlRow,xlCol,ThisOne)
        
    xlSheet.write(xlRow-1,xlCol,Units)

    if DoSubSum == 'True':
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,'Pos n')
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,ThisOne+'_pos')
        xlSheet.write(xlRow-1,xlCol,Units)
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,'Neg n')
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,ThisOne+'_neg')
        xlSheet.write(xlRow-1,xlCol,Units)
    
    data = numpy.ma.masked_where(abs(ds.series[ThisOne]['Data']-float(-9999))<c.eps,ds.series[ThisOne]['Data'])
    for month in range(M1st,M2nd+1):
        if month == 1 or month == 3 or month == 5 or month == 7 or month == 8 or month == 10 or month == 12:
            dRan = 31
        if month == 2:
            if ds.series['Year']['Data'][0] % 4 == 0:
                dRan = 29
            else:
                dRan = 28
        if month == 4 or month == 6 or month == 9 or month == 11:
            dRan = 30
            
        for day in range(1,dRan+1):
            xlRow = xlRow + 1
            if ThisOne == 'rst' or ThisOne == 'Gst' or ThisOne == 'Gst_mol':
                di = numpy.where((ds.series['Month']['Data']==month) & (ds.series['Day']['Data']==day) & (ds.series[ThisOne]['Flag'] < 61))[0]
                ti = numpy.where((ds.series['Month']['Data']==month) & (ds.series['Day']['Data']==day))[0]
                nRecs = len(ti)
                check = numpy.ma.empty(nRecs,str)
                for i in range(nRecs):
                    index = ti[i]
                    check[i] = ds.series['Day']['Data'][index]
                if len(check) < 48:
                    di = []
            elif ThisOne == 'GPP' or ThisOne == 'GPP_mmol':
                di = numpy.where((ds.series['Month']['Data']==month) & (ds.series['Day']['Data']==day) & ((ds.series[ThisOne]['Flag'] != 31) & (ds.series[ThisOne]['Flag'] != 81)))[0]
                ti = numpy.where((ds.series['Month']['Data']==month) & (ds.series['Day']['Data']==day))[0]
                nRecs = len(ti)
                check = numpy.ma.empty(nRecs,str)
                for i in range(nRecs):
                    index = ti[i]
                    check[i] = ds.series['Day']['Data'][index]
                if len(check) < 48:
                    di = []
            elif ThisOne == 'Re_mmol':
                di = numpy.where((ds.series['Month']['Data']==month) & (ds.series['Day']['Data']==day) & ((ds.series[ThisOne]['Flag'] == 0) | (ds.series[ThisOne]['Flag'] > 69)))[0]
                ti = numpy.where((ds.series['Month']['Data']==month) & (ds.series['Day']['Data']==day))[0]
                nRecs = len(ti)
                check = numpy.ma.empty(nRecs,str)
                for i in range(nRecs):
                    index = ti[i]
                    check[i] = ds.series['Day']['Data'][index]
                if len(check) < 48:
                    di = []
            else:
                di = numpy.where((ds.series['Month']['Data']==month) & (ds.series['Day']['Data']==day))[0]
                nRecs = len(di)
                check = numpy.ma.empty(nRecs,str)
                for i in range(nRecs):
                    index = di[i]
                    check[i] = ds.series['Day']['Data'][index]
                if len(check) < 48:
                    di = []
            
            if DoSoil == 'True':
                Num,Av = get_soilaverages(data[di])
                if xlCol == 3:
                    xlCol = 2
                    xlSheet.write(xlRow,xlCol-2,monthabr[month-1])
                    xlSheet.write(xlRow,xlCol-1,day)
                else:
                    xlCol = xlCol - 1
            else:
                if DoSum == 'True':
                    Num,Sum = get_sums(data[di])
                if DoMinMax == 'True':
                    Num,Min,Max = get_minmax(data[di])
                if DoMean == 'True':
                    if DoMinMax == 'True':
                        Num2,Av = get_averages(data[di])
                    else:
                        Num,Av = get_averages(data[di])
                if DoSubSum == 'True':
                    PosNum,NegNum,SumPos,SumNeg = get_subsums(data[di])
                xlCol = 2
                xlSheet.write(xlRow,xlCol-2,monthabr[month-1])
                xlSheet.write(xlRow,xlCol-1,day)
            
            xlSheet.write(xlRow,xlCol,Num)
            xlCol = xlCol + 1
            if DoSoil == 'True':
                xlSheet.write(xlRow,xlCol,Av)
            elif DoMinMax == 'True':
                xlSheet.write(xlRow,xlCol,Min)
                xlCol = xlCol + 1
                xlSheet.write(xlRow,xlCol,Max)
                if DoMean == 'True':
                    xlCol = xlCol + 1
                    xlSheet.write(xlRow,xlCol,Av)
            elif DoMinMax == 'False' and DoMean == 'True':
                xlSheet.write(xlRow,xlCol,Av)
            elif DoSum == 'True':
                xlSheet.write(xlRow,xlCol,Sum)
                if DoSubSum == 'True':
                    xlCol = xlCol + 1
                    xlSheet.write(xlRow,xlCol,PosNum)
                    xlCol = xlCol + 1
                    xlSheet.write(xlRow,xlCol,SumPos)
                    xlCol = xlCol + 1
                    xlSheet.write(xlRow,xlCol,NegNum)
                    xlCol = xlCol + 1
                    xlSheet.write(xlRow,xlCol,SumNeg)
    
    if DoSoil == 'True': 
        return xlCol,xlSheet
    else:
        return
