from configobj import ConfigObj
import ast
import cfg
import copy
import constants as c
import csv
import datetime
import numpy
import os
import sys
import time
import Tkinter, tkFileDialog
import xlrd
import xlwt
import netCDF4
import logging
import qcts
import qcutils

log = logging.getLogger('qc.io')

class DataStructure(object):
    def __init__(self):
        self.series = {}
        self.globalattributes = {}
        self.mergeserieslist = []
        self.averageserieslist = []
        self.soloserieslist = []

def copy_datastructure(cf,ds_in):
    '''
    Return a copy of a data structure based on the following rules:
     1) if the netCDF file at the "copy_to" level does not exist
        then copy the existing data structure at the "input" level
        to create a new data structure at the "output" level.
    '''
    # assumptions that need to be checked are:
    #  - the start datetime of the two sets of data are the same
    #  - the end datetime of the L3 data is the same or after the
    #    end datetime of the the L4 data
    #    - if the end datetimes are the same then we are just re-processing something
    #    - if the end datetime for the L3 data is after the end date of the L4 data
    #      then more data has been added to this year and the user wants to gap fill
    #      the new data
    # modificatons to be made:
    #  - check the modification datetime of the L3 and L4 files:
    #     - if the L3 file is newer than the L4 file the disregard the "UseExistingOutFile" setting
    # get the output (L4) file name
    ct_filename = cf['Files']['file_path']+cf['Files']['out_filename']
    # if the L4 file does not exist then create the L4 data structure as a copy
    # of the L3 data structure
    if not os.path.exists(ct_filename):
        ds_out = copy.deepcopy(ds_in)
    # if the L4 file does exist ...
    if os.path.exists(ct_filename):
        # check to see if the user wants to use it
        if cf['Options']['UseExistingOutFile']!='Yes':
            # if the user doesn't want to use the existing L4 data then create
            # the L4 data structure as a copy of the L3 data structure
            ds_out = copy.deepcopy(ds_in)
        else:
            # the user wants to use the data from an existing L4 file
            # get the netCDF file name at the "input" level
            outfilename = get_outfilename_from_cf(cf)
            # read the netCDF file at the "input" level
            ds_file = nc_read_series(outfilename)
            dt_file = ds_file.series['DateTime']['Data']
            sd_file = str(dt_file[0])
            ed_file = str(dt_file[-1])
            # create a copy of the data
            ds_out = copy.deepcopy(ds_in)
            dt_out = ds_out.series['DateTime']['Data']
            ts = ds_out.globalattributes['time_step']
            sd_out = str(dt_out[0])
            ed_out = str(dt_out[-1])
            # get the start and end indices based on the start and end dates
            si = qcutils.GetDateIndex(dt_out,sd_file,ts=ts,default=0,match='exact')
            ei = qcutils.GetDateIndex(dt_out,ed_file,ts=ts,default=-1,match='exact')
            # now replace parts of ds_out with the data read from file
            for ThisOne in ds_file.series.keys():
                # check to see if the L4 series exists in the L3 data
                if ThisOne in ds_out.series.keys():
                    # ds_out is the copy of the L3 data, now fill it with the L4 data read from file
                    ds_out.series[ThisOne]['Data'][si:ei+1] = ds_file.series[ThisOne]['Data']
                    ds_out.series[ThisOne]['Flag'][si:ei+1] = ds_file.series[ThisOne]['Flag']
                else:
                    # if it doesn't, create the series and put the data into it
                    ds_out.series[ThisOne] = {}
                    ds_out.series[ThisOne] = ds_file.series[ThisOne].copy()
                    # check to see if we have to append data to make the copy of the L4 data now
                    # in the L3 data structure the same length as the existing L3 data
                    nRecs_file = int(ds_file.globalattributes['nc_nrecs'])
                    nRecs_out = int(ds_out.globalattributes['nc_nrecs'])
                    if nRecs_file < nRecs_out:
                        # there is more data at L3 than at L4
                        # append missing data to make the series the same length
                        nRecs_append = nRecs_out - nRecs_file
                        data = numpy.array([-9999]*nRecs_append,dtype=numpy.float64)
                        flag = numpy.ones(nRecs_append,dtype=numpy.int32)
                        ds_out.series[ThisOne]['Data'] = numpy.concatenate((ds_out.series[ThisOne]['Data'],data))
                        ds_out.series[ThisOne]['Flag'] = numpy.concatenate((ds_out.series[ThisOne]['Flag'],flag))
                    elif nRecs_file > nRecs_out:
                        # tell the user something is wrong
                        log.error('copy_datastructure: L3 file contains less data than L4 file')
                        # return an empty dictionary
                        ds_out = {}
                    else:
                        # nRecs_file and nRecs_out are equal so we do not need to do anything
                        pass
    return ds_out

def nc2xl(cf):
    # get the filenames
    outfilename = get_outfilename_from_cf(cf)
    xlfilename= outfilename.replace('.nc','.xls')
    # get the output list
    outputlist = get_outputlist_from_cf(cf,'xl')
    # read the netCDF file
    ds = nc_read_series(outfilename)
    # write the variables to the Excel file
    xl_write_series(ds,xlfilename,outputlist=outputlist)

def read_eddypro_full(csvname):
    ds = DataStructure()
    csvfile = open(csvname,'rb')
    csvreader = csv.reader(csvfile)
    n = 0
    adatetime = []
    us_data_list = []
    us_flag_list = []
    Fh_data_list = []
    Fh_flag_list = []
    Fe_data_list = []
    Fe_flag_list = []
    Fc_data_list = []
    Fc_flag_list = []
    for row in csvreader:
        if n==0:
            header=row
        elif n==1:
            varlist=row
            us_data_col = varlist.index('u*')
            us_flag_col = varlist.index('qc_Tau')
            Fh_data_col = varlist.index('H')
            Fh_flag_col = varlist.index('qc_H')
            Fe_data_col = varlist.index('LE')
            Fe_flag_col = varlist.index('qc_LE')
            Fc_data_col = varlist.index('co2_flux')
            Fc_flag_col = varlist.index('qc_co2_flux')
        elif n==2:
            unitlist=row
        else:
            adatetime.append(datetime.datetime.strptime(row[1]+' '+row[2],'%Y-%m-%d %H:%M'))
            us_data_list.append(float(row[us_data_col]))
            us_flag_list.append(float(row[us_flag_col]))
            Fh_data_list.append(float(row[Fh_data_col]))
            Fh_flag_list.append(float(row[Fh_flag_col]))
            Fe_data_list.append(float(row[Fe_data_col]))
            Fe_flag_list.append(float(row[Fe_flag_col]))
            Fc_data_list.append(float(row[Fc_data_col]))
            Fc_flag_list.append(float(row[Fc_flag_col]))
        n = n + 1
    nRecs = len(adatetime)
    adatetime = qcutils.RoundDateTime(adatetime,dt=30)
    ds.series['DateTime'] = {}
    ds.series['DateTime']['Data'] = adatetime
    ds.series['ustar'] = {}
    ds.series['ustar']['Data'] = numpy.array(us_data_list,dtype=numpy.float64)
    ds.series['ustar']['Flag'] = numpy.array(us_flag_list,dtype=numpy.int32)
    ds.series['Fh'] = {}
    ds.series['Fh']['Data'] = numpy.array(Fh_data_list,dtype=numpy.float64)
    ds.series['Fh']['Flag'] = numpy.array(Fh_flag_list,dtype=numpy.int32)
    ds.series['Fe'] = {}
    ds.series['Fe']['Data'] = numpy.array(Fe_data_list,dtype=numpy.float64)
    ds.series['Fe']['Flag'] = numpy.array(Fe_flag_list,dtype=numpy.int32)
    ds.series['Fc'] = {}
    ds.series['Fc']['Data'] = numpy.array(Fc_data_list,dtype=numpy.float64)
    ds.series['Fc']['Flag'] = numpy.array(Fc_flag_list,dtype=numpy.int32)
    return ds
    
def xl2nc(cf,InLevel):
    # get the data series from the Excel file
    ds = xl_read_series(cf)
    if len(ds.series.keys())==0: return 1
    # get the netCDF attributes from the control file
    qcts.do_attributes(cf,ds)
    # get a series of Python datetime objects from the Excel datetime
    qcutils.get_datetimefromxldate(ds)
    #check for gaps in the Excel datetime series
    has_gaps = qcutils.CheckTimeStep(ds,fix='gaps')
    # write the processing level to a global attribute
    ds.globalattributes['nc_level'] = str(InLevel)
    # get the start and end date from the datetime series unless they were
    # given in the control file
    if 'start_date' not in ds.globalattributes.keys():
        ds.globalattributes['start_date'] = str(ds.series['DateTime']['Data'][0])
    if 'end_date' not in ds.globalattributes.keys():
        ds.globalattributes['end_date'] = str(ds.series['DateTime']['Data'][-1])
    # get the year, month, day, hour, minute and second from the xl date/time
    qcutils.get_ymdhmsfromxldate(ds)
    # do any functions to create new series
    qcts.do_functions(cf,ds)
    # write the data to the netCDF file
    outfilename = get_outfilename_from_cf(cf)
    ncFile = nc_open_write(outfilename)
    nc_write_series(ncFile,ds)
    return 0

def get_controlfilecontents(ControlFileName):
    log.info(' Processing the control file ')
    if len(ControlFileName)!=0:
        cf = ConfigObj(ControlFileName)
        cf['ControlFileName'] = ControlFileName
    else:
        cf = ConfigObj()
    return cf

def get_controlfilename(path=''):
    log.info(' Choosing the control file ')
    root = Tkinter.Tk(); root.withdraw()
    name = tkFileDialog.askopenfilename(initialdir=path)
    root.destroy()
    return name

def get_ncdtype(Series):
    sd = Series.dtype.name
    dt = 'f'
    if sd=='float64': dt = 'd'
    if sd=='int32': dt = 'i'
    if sd=='int64': dt = 'l'
    return dt

def get_filename_dialog(path='.',title='Choose a file'):
    '''Get an Excel file name'''
    root = Tkinter.Tk(); root.withdraw()
    FileName = tkFileDialog.askopenfilename(parent=root,initialdir=path,title=title)
    root.destroy()
    return str(FileName)

def get_infilename_from_cf(cf):
    filename = ""
    if "Files" in cf.keys():
        if "file_path" in cf['Files'].keys():
            if "in_filename" in cf['Files'].keys():
                filename = cf['Files']['file_path']+cf['Files']['in_filename']
            else:
                log.error("get_infilename_from_cf: 'in_filename' key not found in 'Files' section of control file")
        else:
            log.error("get_infilename_from_cf: 'file_path' key not found in 'Files' section of control file")
    else:
        log.error("get_infilename_from_cf: 'Files' section not found in control file")
    return str(filename)

def get_outfilename_from_cf(cf):
    try:
        filename = cf['Files']['file_path']+cf['Files']['out_filename']
    except:
        log.error('get_outfilename_from_cf: Error getting out_filename from control file')
        raise Exception('out_filename, file_path or Files section not in control path')
    return str(filename)

def get_keyvalue_from_cf(section,key):
    try:
        value = section[key]
    except:
        log.error('get_keyvalue_from_cf: '+str(key)+' not found in '+str(section.name)+' section of control file')
        raise Exception('key not found in section')
    return value

def get_outputlist_from_cf(cf,filetype):
    try:
        outputlist = ast.literal_eval(cf['Output'][filetype])
    except:
        #log.info('get_outputlist_from_cf: Unable to get output list from Output section in control file')
        outputlist = None
    return outputlist

def get_seriesstats(cf,ds):
    # open an Excel file for the flag statistics
    level = ds.globalattributes['nc_level']
    out_filename = get_outfilename_from_cf(cf)
    xl_filename = out_filename.replace('.nc','_FlagStats.xls')
    log.info(' Writing flag stats to Excel file '+xl_filename)
    xlFile = xlwt.Workbook()
    xlFlagSheet = xlFile.add_sheet('Flag')
    # get the flag statistics
    xlRow = 0
    xlCol = 0
    xlFlagSheet.write(xlRow,xlCol,'0:')
    xlFlagSheet.write(xlRow,xlCol+1,ds.globalattributes['Flag0'])
    xlFlagSheet.write(xlRow,xlCol+2,'1:')
    xlFlagSheet.write(xlRow,xlCol+3,ds.globalattributes['Flag1'])
    xlFlagSheet.write(xlRow,xlCol+4,'2:')
    xlFlagSheet.write(xlRow,xlCol+5,ds.globalattributes['Flag2'])
    xlFlagSheet.write(xlRow,xlCol+6,'3:')
    xlFlagSheet.write(xlRow,xlCol+7,ds.globalattributes['Flag3'])
    xlFlagSheet.write(xlRow,xlCol+8,'4:')
    xlFlagSheet.write(xlRow,xlCol+9,ds.globalattributes['Flag4'])
    xlFlagSheet.write(xlRow,xlCol+10,'5:')
    xlFlagSheet.write(xlRow,xlCol+11,ds.globalattributes['Flag5'])
    xlFlagSheet.write(xlRow,xlCol+12,'6:')
    xlFlagSheet.write(xlRow,xlCol+13,ds.globalattributes['Flag6'])
    xlFlagSheet.write(xlRow,xlCol+14,'7:')
    xlFlagSheet.write(xlRow,xlCol+15,ds.globalattributes['Flag7'])
    xlRow = xlRow + 1
    xlFlagSheet.write(xlRow,xlCol,'10:')
    xlFlagSheet.write(xlRow,xlCol+1,ds.globalattributes['Flag10'])
    xlFlagSheet.write(xlRow,xlCol+2,'11:')
    xlFlagSheet.write(xlRow,xlCol+3,ds.globalattributes['Flag11'])
    xlFlagSheet.write(xlRow,xlCol+4,'12:')
    xlFlagSheet.write(xlRow,xlCol+5,ds.globalattributes['Flag12'])
    xlFlagSheet.write(xlRow,xlCol+6,'13:')
    xlFlagSheet.write(xlRow,xlCol+7,ds.globalattributes['Flag13'])
    xlFlagSheet.write(xlRow,xlCol+8,'14:')
    xlFlagSheet.write(xlRow,xlCol+9,ds.globalattributes['Flag14'])
    xlFlagSheet.write(xlRow,xlCol+10,'15:')
    xlFlagSheet.write(xlRow,xlCol+11,ds.globalattributes['Flag15'])
    xlFlagSheet.write(xlRow,xlCol+12,'16:')
    xlFlagSheet.write(xlRow,xlCol+13,ds.globalattributes['Flag16'])
    xlFlagSheet.write(xlRow,xlCol+14,'17:')
    xlFlagSheet.write(xlRow,xlCol+15,ds.globalattributes['Flag17'])
    xlFlagSheet.write(xlRow,xlCol+16,'18:')
    xlFlagSheet.write(xlRow,xlCol+17,ds.globalattributes['Flag18'])
    xlFlagSheet.write(xlRow,xlCol+18,'19:')
    xlFlagSheet.write(xlRow,xlCol+19,ds.globalattributes['Flag19'])
    xlRow = xlRow + 1
    xlFlagSheet.write(xlRow,xlCol,'30:')
    xlFlagSheet.write(xlRow,xlCol+1,ds.globalattributes['Flag30'])
    xlFlagSheet.write(xlRow,xlCol+2,'31:')
    xlFlagSheet.write(xlRow,xlCol+3,ds.globalattributes['Flag31'])
    xlFlagSheet.write(xlRow,xlCol+4,'32:')
    xlFlagSheet.write(xlRow,xlCol+5,ds.globalattributes['Flag32'])
    xlFlagSheet.write(xlRow,xlCol+6,'33:')
    xlFlagSheet.write(xlRow,xlCol+7,ds.globalattributes['Flag33'])
    xlFlagSheet.write(xlRow,xlCol+8,'34:')
    xlFlagSheet.write(xlRow,xlCol+9,ds.globalattributes['Flag34'])
    xlFlagSheet.write(xlRow,xlCol+10,'35:')
    xlFlagSheet.write(xlRow,xlCol+11,ds.globalattributes['Flag35'])
    xlFlagSheet.write(xlRow,xlCol+12,'36:')
    xlFlagSheet.write(xlRow,xlCol+13,ds.globalattributes['Flag36'])
    xlFlagSheet.write(xlRow,xlCol+14,'37:')
    xlFlagSheet.write(xlRow,xlCol+15,ds.globalattributes['Flag37'])
    xlFlagSheet.write(xlRow,xlCol+16,'38:')
    xlFlagSheet.write(xlRow,xlCol+17,ds.globalattributes['Flag38'])
    xlFlagSheet.write(xlRow,xlCol+18,'39:')
    xlFlagSheet.write(xlRow,xlCol+19,ds.globalattributes['Flag39'])
    bins = numpy.arange(-0.5,23.5)
    xlRow = 5
    xlCol = 1
    for Value in bins[:len(bins)-1]:
        xlFlagSheet.write(xlRow,xlCol,int(Value+0.5))
        xlCol = xlCol + 1
    xlRow = xlRow + 1
    xlCol = 0
    dsVarNames = ds.series.keys()
    dsVarNames.sort(key=unicode.lower)
    for ThisOne in dsVarNames:
        data,flag = qcutils.GetSeries(ds, ThisOne)
        hist, bin_edges = numpy.histogram(flag, bins=bins)
        xlFlagSheet.write(xlRow,xlCol,ThisOne)
        xlCol = xlCol + 1
        for Value in hist:
            xlFlagSheet.write(xlRow,xlCol,float(Value))
            xlCol = xlCol + 1
        xlCol = 0
        xlRow = xlRow + 1
    xlFile.save(xl_filename)

def load_controlfile(path=''):
    ''' 
    Returns a control file object.
    USAGE: cf = load_controlfile([path=<some_path_to_a_controlfile>])
    The "path" keyword is optional.
    '''
    name = get_controlfilename(path=path)
    cf = get_controlfilecontents(name)
    return cf

def nc_read_series(ncFullName):
    ''' Read a netCDF file and put the data and meta-data into a DataStructure'''
    log.info(' Reading netCDF file '+ncFullName)
    netCDF4.default_encoding = 'latin-1'
    ncFile = netCDF4.Dataset(ncFullName,'r')
    ds = DataStructure()
    gattrlist = ncFile.ncattrs()
    if len(gattrlist)!=0:
        for gattr in gattrlist:
            ds.globalattributes[gattr] = getattr(ncFile,gattr)
            if 'time_step' in ds.globalattributes: c.ts = ds.globalattributes['time_step']
    for ThisOne in ncFile.variables.keys():
        if '_QCFlag' not in ThisOne:
            # create the series in the data structure
            ds.series[unicode(ThisOne)] = {}
            # get the data variable object
            ds.series[ThisOne]['Data'] = ncFile.variables[ThisOne][:]
            # check for a QC flag and if it exists, load it
            if ThisOne+'_QCFlag' in ncFile.variables.keys():
                ds.series[ThisOne]['Flag'] = ncFile.variables[ThisOne+'_QCFlag'][:]
            else:
                nRecs = numpy.size(ds.series[ThisOne]['Data'])
                ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,dtype=numpy.int32)
            # get the variable attributes
            vattrlist = ncFile.variables[ThisOne].ncattrs()
            if len(vattrlist)!=0:
                ds.series[ThisOne]['Attr'] = {}
                for vattr in vattrlist:
                    ds.series[ThisOne]['Attr'][vattr] = getattr(ncFile.variables[ThisOne],vattr)
    ncFile.close()
    # get a series of Python datetime objects
    qcutils.get_datetimefromymdhms(ds)
    return ds

def nc_open_write(ncFullName,nctype='NETCDF4'):
    log.info(' Opening netCDF file '+ncFullName+' for writing')
    ncFile = netCDF4.Dataset(ncFullName,'w',format=nctype)
    return ncFile

def nc_write_series(ncFile,ds,outputlist=None):
    ds.globalattributes['QC_version'] = str(cfg.version_name)+' '+str(cfg.version_number)
    for ThisOne in ds.globalattributes.keys():
        setattr(ncFile,ThisOne,ds.globalattributes[ThisOne])
    t = time.localtime()
    rundatetime = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    setattr(ncFile,'nc_rundatetime',rundatetime)
    # we specify the size of the Time dimension because netCDF4 is slow to write files
    # when the Time dimension is unlimited
    nRecs = int(ds.globalattributes['nc_nrecs'])
    ncFile.createDimension('Time',nRecs)
    SeriesList = ds.series.keys()
    if outputlist==None:
        outputlist = SeriesList
    else:
        for ThisOne in outputlist:
            if ThisOne not in SeriesList:
                log.info(' Requested series '+ThisOne+' not found in data structure')
                outputlist.remove(ThisOne)
        if len(outputlist)==0:
            outputlist = SeriesList
    # can't write an array of Python datetime objects to a netCDF file
    # actually, this could be written as characters
    if 'DateTime' in outputlist: SeriesList.remove('DateTime')
    # now make sure the date and time series are in outputlist
    for ThisOne in ['xlDateTime','Year','Month','Day','Hour','Minute','Second','Hdh']:
        if ThisOne not in outputlist: outputlist.append(ThisOne)
    # sort the output list into alphabetic order
    outputlist.sort()
    # write everything else to the netCDF file
    for ThisOne in outputlist:
        dt = get_ncdtype(ds.series[ThisOne]['Data'])
        ncVar = ncFile.createVariable(ThisOne,dt,('Time',))
        ncVar[:] = ds.series[ThisOne]['Data'].tolist()
        for attr in ds.series[ThisOne]['Attr']:
            setattr(ncVar,attr,ds.series[ThisOne]['Attr'][attr])
        dt = get_ncdtype(ds.series[ThisOne]['Flag'])
        ncVar = ncFile.createVariable(ThisOne+'_QCFlag',dt,('Time',))
        ncVar[:] = ds.series[ThisOne]['Flag'].tolist()
        setattr(ncVar,'long_name','QC flag')
        setattr(ncVar,'units','none')
    ncFile.close()

def xl_read_flags(cf,ds,level,VariablesInFile):
    # First data row in Excel worksheets.
    FirstDataRow = int(get_keyvalue_from_cf(cf['Files'][level],'first_data_row')) - 1
    HeaderRow = int(get_keyvalue_from_cf(cf['Files']['in'],'header_row')) - 1
    # Get the full name of the Excel file from the control file.
    xlFullName = get_filename_from_cf(cf,level)
    # Get the Excel workbook object.
    if os.path.isfile(xlFullName):
        xlBook = xlrd.open_workbook(xlFullName)
    else:
        log.error(' Excel file '+xlFullName+' not found, choose another')
        xlFullName = get_filename_dialog(path='.',title='Choose an Excel file')
        if len(xlFullName)==0:
            return
        xlBook = xlrd.open_workbook(xlFullName)
    ds.globalattributes['xlFullName'] = xlFullName

    for ThisOne in VariablesInFile:
        if 'xl' in cf['Variables'][ThisOne].keys():
            log.info(' Getting flags for '+ThisOne+' from spreadsheet')
            ActiveSheet = xlBook.sheet_by_name('Flag')
            LastDataRow = int(ActiveSheet.nrows)
            HeaderList = [x.lower() for x in ActiveSheet.row_values(HeaderRow)]
            if cf['Variables'][ThisOne]['xl']['name'] in HeaderList:
                xlCol = HeaderRow.index(cf['Variables'][ThisOne]['xl']['name'])
                Values = ActiveSheet.col_values(xlCol)[FirstDataRow:LastDataRow]
                Types = ActiveSheet.col_types(xlCol)[FirstDataRow:LastDataRow]
                ds.series[ThisOne]['Flag'] = numpy.array([-9999]*len(Values),numpy.int32)
                for i in range(len(Values)):
                    if Types[i]==2: #xlType=3 means a date/time value, xlType=2 means a number
                        ds.series[ThisOne]['Flag'][i] = numpy.int32(Values[i])
                    else:
                        log.error('  xl_read_flags: flags for '+ThisOne+' not found in xl file')
    return ds

def xl_read_series(cf):
    # Instance the data structure object.
    ds = DataStructure()
    # get the filename
    FileName = get_infilename_from_cf(cf)
    if len(FileName)==0:
        log.error(' xl_read_series: file name not found in control file')
        return ds
    if not os.path.exists(FileName):
        log.error(' xl_read_series: input file specified in control file not found')
        return ds
    # convert from Excel row number to xlrd row number
    FirstDataRow = int(get_keyvalue_from_cf(cf['Files'],'in_firstdatarow')) - 1
    HeaderRow = int(get_keyvalue_from_cf(cf['Files'],'in_headerrow')) - 1
    # get the Excel workbook object.
    log.info(' Opening and reading Excel file '+FileName)
    xlBook = xlrd.open_workbook(FileName)
    log.info(' Opened and read Excel file '+FileName)
    ds.globalattributes['featureType'] = 'timeseries'
    ds.globalattributes['xl_filename'] = FileName
    ds.globalattributes['xl_datemode'] = str(xlBook.datemode)
    xlsheet_names = [x.lower() for x in xlBook.sheet_names()]
    # Get the Excel file modification date and time, these will be
    # written to the netCDF file to uniquely identify the version
    # of the Excel file used to create this netCDF file.
    s = os.stat(FileName)
    t = time.localtime(s.st_mtime)
    ds.globalattributes['xl_moddatetime'] = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    # Loop over the variables defined in the 'Variables' section of the
    # configuration file.
    for ThisOne in cf['Variables'].keys():
        if 'xl' in cf['Variables'][ThisOne].keys():
            if 'sheet' in cf['Variables'][ThisOne]['xl'].keys():
                xlsheet_name = cf['Variables'][ThisOne]['xl']['sheet']
                if xlsheet_name.lower() in xlsheet_names:
                    log.info(' Getting data for '+ThisOne+' from spreadsheet')
                    xlsheet_index = xlsheet_names.index(xlsheet_name.lower())
                    ActiveSheet = xlBook.sheet_by_index(xlsheet_index)
                    HeaderList = [x.lower() for x in ActiveSheet.row_values(HeaderRow)]
                    if cf['Variables'][ThisOne]['xl']['name'].lower() in HeaderList:
                        LastDataRow = int(ActiveSheet.nrows)
                        ds.series[unicode(ThisOne)] = {}
                        xlCol = HeaderList.index(cf['Variables'][ThisOne]['xl']['name'].lower())
                        Values = ActiveSheet.col_values(xlCol)[FirstDataRow:LastDataRow]
                        Types = ActiveSheet.col_types(xlCol)[FirstDataRow:LastDataRow]
                        ds.series[ThisOne]['Data'] = numpy.array([-9999]*len(Values),numpy.float64)
                        ds.series[ThisOne]['Flag'] = numpy.ones(len(Values),dtype=numpy.int32)
                        # we could use "where" and get rid of this for loop
                        for i in range(len(Values)):
                            if (Types[i]==3) or (Types[i]==2): #xlType=3 means a date/time value, xlType=2 means a number
                                ds.series[ThisOne]['Data'][i] = numpy.float64(Values[i])
                                ds.series[ThisOne]['Flag'][i] = numpy.int32(0)
                    else:
                        log.error('  xl_read_series: series '+ThisOne+' not found in xl file')
                else:
                    log.error('  xl_read_series: sheet '+xlsheet_name+' not found in xl file')
            else:
                log.error('  xl_read_series: key "sheet" not found in control file entry for '+ThisOne)
        else:
            log.error('  xl_read_series: key "xl" not found in control file entry for '+ThisOne)
    ds.globalattributes['nc_nrecs'] = str(len(ds.series['xlDateTime']['Data']))
    return ds

def xl_write_series(ds, xlfullname, outputlist=None):
    # open the Excel file
    log.info(' Opening and writing Excel file '+xlfullname)
    xlfile = xlwt.Workbook()
    # add sheets to the Excel file
    xlAttrSheet = xlfile.add_sheet('Attr')
    xlDataSheet = xlfile.add_sheet('Data')
    xlFlagSheet = xlfile.add_sheet('Flag')
    # write the global attributes
    log.info(' Writing the global attributes to Excel file '+xlfullname)
    xlcol = 0
    xlrow = 0
    xlAttrSheet.write(xlrow,xlcol,'Global attributes')
    xlrow = xlrow + 1
    globalattrlist = ds.globalattributes.keys()
    globalattrlist.sort()
    for ThisOne in [x for x in globalattrlist if 'Flag' not in x]:
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
        xlrow = xlrow + 1
    for ThisOne in [x for x in globalattrlist if 'Flag' in x]:
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
        xlrow = xlrow + 1
    # write the variable attributes
    log.info(' Writing the variable attributes to Excel file '+xlfullname)
    xlrow = xlrow + 1
    xlAttrSheet.write(xlrow,xlcol,'Variable attributes')
    xlrow = xlrow + 1
    xlcol_varname = 0
    xlcol_attrname = 1
    xlcol_attrvalue = 2
    variablelist = ds.series.keys()
    if outputlist==None:
        outputlist = variablelist
    else:
        for ThisOne in outputlist:
            if ThisOne not in variablelist:
                log.info(' Requested series '+ThisOne+' not found in data structure')
                outputlist.remove(ThisOne)
        if len(outputlist)==0:
            outputlist = variablelist
    outputlist.sort()
    if 'DateTime' in outputlist:
        outputlist.remove('DateTime')
    for ThisOne in outputlist:
        xlAttrSheet.write(xlrow,xlcol_varname,ThisOne)
        attributelist = ds.series[ThisOne]['Attr'].keys()
        attributelist.sort()
        for Attr in attributelist:
            xlAttrSheet.write(xlrow,xlcol_attrname,Attr)
            xlAttrSheet.write(xlrow,xlcol_attrvalue,ds.series[ThisOne]['Attr'][Attr])
            xlrow = xlrow + 1
    # write the Excel date/time to the data and the QC flags as the first column
    log.info(' Writing the datetime to Excel file '+xlfullname)
    nRecs = len(ds.series['xlDateTime']['Data'])
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    xlDataSheet.write(2,xlcol,'xlDateTime')
    for j in range(nRecs):
        xlDataSheet.write(j+3,xlcol,ds.series['xlDateTime']['Data'][j],d_xf)
        xlFlagSheet.write(j+3,xlcol,ds.series['xlDateTime']['Data'][j],d_xf)
    # remove xlDateTime from the list of variables to be written to the Excel file
    if 'xlDateTime' in outputlist:
        outputlist.remove('xlDateTime')
    # now start looping over the other variables in the xl file
    xlcol = xlcol + 1
    # loop over variables to be output to xl file
    for ThisOne in outputlist:
        # put up a progress message
        log.info(' Writing '+ThisOne+' into column '+str(xlcol)+' of the Excel file')
        # write the units and the variable name to the header rows in the xl file
        attrlist = ds.series[ThisOne]['Attr'].keys()
        if 'long_name' in attrlist:
            longname = ds.series[ThisOne]['Attr']['long_name']
        elif 'Description' in attrlist:
            longname = ds.series[ThisOne]['Attr']['Description']
        else:
            longname = None
        if 'units' in attrlist:
            units = ds.series[ThisOne]['Attr']['units']
        elif 'Units' in attrlist:
            units = ds.series[ThisOne]['Attr']['Units']
        else:
            units = None
        xlDataSheet.write(0,xlcol,longname)
        xlDataSheet.write(1,xlcol,units)
        xlDataSheet.write(2,xlcol,ThisOne)
        # loop over the values in the variable series (array writes don't seem to work)
        for j in range(nRecs):
            xlDataSheet.write(j+3,xlcol,float(ds.series[ThisOne]['Data'][j]))
        # check to see if this variable has a quality control flag
        if 'Flag' in ds.series[ThisOne].keys():
            # write the QC flag name to the xk file
            xlFlagSheet.write(2,xlcol,ThisOne)
            # specify the format of the QC flag (integer)
            d_xf = xlwt.easyxf(num_format_str='0')
            # loop over QV flag values and write to xl file
            for j in range(nRecs):
                xlFlagSheet.write(j+3,xlcol,int(ds.series[ThisOne]['Flag'][j]),d_xf)
        # increment the column pointer
        xlcol = xlcol + 1
    
    xlfile.save(xlfullname)
    