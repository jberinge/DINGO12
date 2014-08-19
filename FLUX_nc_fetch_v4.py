############################################################################
# This script correlates tower data with external meteorology for a variable.  Then calulates the correlation coefficient.  Then adjusts new time series and gap fills 
# Inputs:
#    base_dir_netcdf : Path where to save the netcdf dataset (e.g. '/.')
#    latlim : Latitude limits (in degrees) of the domain to examine (e.g [-40.0, -35.0])
#    lonlim : Longitude limits (in degrees) of the domain to examine (e.g [140.0, 150.0])
#    excel_file: Path to the Excel file containing the indices to correlate with (e.g 'CIall_mon_new_cropyears.xls')
#    index_name: Name of the column containing the index to correlate with the rainfall data (e.g 'Nino1+2_ANOM')
#    months_offset : Offset (in months) between the climate indices and the rainfall data 
#    output_prefix: Prefix of the output NetCDF file that will be created (e.g. '' will create a file "$(season)_$(index_name).nc)"
#
# Notes:
#   The excel file is suppose to have the dates on the first column.
#   The index name is expected to be found on the first row
#
# Programmed by Jason (Dec 1, 2012)
############################################################################

import pandas as pd
import datetime as dt
import xlrd
import numpy as np
import netCDF4

def excel_to_pydate(exceldate,Site_ID):
    # datemode: 0 for 1900-based, 1 for 1904-based
    if Site_ID=="AliceSprings":
        datemode=1
    else:
        datemode=0
        
    pyear, pmonth, pday, phour, pminute, psecond = xlrd.xldate_as_tuple(exceldate, datemode)
    py_date = dt.datetime(pyear, pmonth, pday, phour, pminute, psecond)
    return(py_date)

def fetch_flux(FLUXfilename,FluxFreq,Site_ID):
   
    #Define the Excel date time variable here
    Exceldatetime_var='xlDateTime'

    #Next Fetch the Flux Tower Data
    #start by Opening the NetCDF file
    nc_file = netCDF4.Dataset(FLUXfilename) 
    
    #List all variable names
    nc_variableNames = nc_file.variables.keys() 
      
    #Create some lists of size equal to length of vnames list.
    temp=list(xrange(len(nc_variableNames)))
    vartemp=list(xrange(len(nc_variableNames)))
    
    # call the function to convert to datetime from excel. Assume datemode: 0
    # Need to know the name of variable Exceldatetime_var, defined above
    times = [excel_to_pydate(elem,Site_ID) for elem in nc_file.variables[Exceldatetime_var]]
    vartemp[0] = times    
    
    Flux_Variable_Units = {}
    
    #Create series for each variable then loop through other variables in the list    
    for index, variable in enumerate(nc_variableNames):              
        temp[index]= nc_file.variables[variable]
        vartemp[index] = pd.Series(temp[index][:],name=variable)
        
        #Create a dictionary of variable attributes to use later
        try:
            instrument1=temp[index].instrument
        except:
            instrument1=""
        try:
            height1=temp[index].height
        except:
            height1=""            
        try:
            long_name1=temp[index].long_name
        except:
            long_name1=""        
        try:
            standard_name1=temp[index].standard_name
        except:
            standard_name1=""            
        try:
            serial_number1=temp[index].serial_number
        except:
            serial_number1=""   
        try:
            units1=temp[index].units
        except:
            units1=""            
        try:
            valid_range1=temp[index].valid_range
        except:
            valid_range1=""  
        try:
            ancillary_variables1=temp[index].ancillary_variables
        except:
            ancillary_variables1=""   
        try:
            rangecheck_upper1=temp[index].rangecheck_upper
        except:
            rangecheck_upper1=""            
        try:
            rangecheck_lower1=temp[index].rangecheck_lower
        except:
            rangecheck_lower1=""  
    
        Flux_Variable_Units.update({variable: units1})
                                    
       
    #Concatenate all the series together
    theDataFrame=pd.concat((vartemp[0:]),axis=1)
    #assign index (excel dfate time) to the dataframe
    theDataFrame.index=times
    
    #Define missing data value and apply to DataFrame
    missing=-9999
    FLUXDataFrame=theDataFrame.replace(missing,np.nan)
    
    #If File has CO2 flux units in mg m-2 s-1 then convert to umol m-2 s-1
    if Flux_Variable_Units['Fc']=="mg/m2/s":
        FLUXDataFrame["Fc"] =FLUXDataFrame["Fc"] /1000 / 44.01 * 1000000
        Flux_Variable_Units['Fc']="umol/m2/s"
        if "Fc_storage" in nc_variableNames:
            FLUXDataFrame["Fc_storage"] =FLUXDataFrame["Fc_storage"] /1000 / 44.01 * 1000000
            Flux_Variable_Units['Fc_storage']="umol/m2/s" 
            
           
    #Check dataframe for duplicates, pad as necessary and sort
    FLUXDataFrame.sort(inplace=True)
    FLUXDataFrame["index"] = FLUXDataFrame.index
    FLUXDataFrame.drop_duplicates(cols='index', take_last=True, inplace=True)
    del FLUXDataFrame["index"]
    #FLUXDataFrame=FLUXDataFrame.asfreq(FluxFreq, method=None)
    
    print FLUXDataFrame
    print "FINISHED fetching the Flux tower nc file"
    return(FLUXDataFrame,Flux_Variable_Units)
