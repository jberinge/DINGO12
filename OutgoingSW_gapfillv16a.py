import numpy as np
import pandas as pd
import datetime as dt
import os
from scipy import stats

#Import custom code modules required for adv processing
import Gap_Fill_climatology_v2 as gap_fill_climatology

###### User-set options ######

write_to_DF=True

###### Peripheral functions ######
def linear_trans(DF,proc,year):

    if proc==True:
        DF=DF.dropna(how='any',axis=0)
        temp_df=pd.DataFrame(np.array(stats.linregress(DF['MOD_Refl'],DF['Albedo'])).reshape(1,5),index=[year])
        if (temp_df.iloc[0,3]<0.05) & (abs(temp_df.iloc[0,4]/temp_df.iloc[0,0])<0.5):
            print str(year)+': data processed, p='+str(round(temp_df[3],3))
        else:
            temp_df[0]=1
            temp_df[1]=0
            print str(year)+': data processed, failed statistical QC, transform not applied... using untransformed MODIS values!'
    else:
        print str(year)+': insufficient data to calculate regression parameters... using untransformed MODIS data!'
        temp_df=pd.DataFrame(np.ones(5).reshape(1,5),index=[year])
    return temp_df


###### Main function ######                
def Fsu_gapfill(FLUXDataframe,VarToProcess,fluxfreq):
    

    global rec_length
    rec_length=int(fluxfreq[:2])
    
    print "Starting Fsu gapfilling..."
    print "Calculate daily averages for obs and add to MODIS_DF"
    
    ###### Calculate daily averages for obs and add to MODIS_DF ######
    
    #Output daily summed incoming and outgoing solar radiation from obs
    by = lambda x: lambda y: getattr(y, x)
    ObsDailySum_DF=pd.DataFrame({'Fsd_mean':FLUXDataframe['Fsd'].groupby([by('year'),by('dayofyear')]).mean(),
                            'Fsd_count':FLUXDataframe['Fsd'].groupby([by('year'),by('dayofyear')]).count(),
                            VarToProcess+'_mean':FLUXDataframe[VarToProcess].groupby([by('year'),by('dayofyear')]).mean(),
                            VarToProcess+'_count':FLUXDataframe[VarToProcess].groupby([by('year'),by('dayofyear')]).count(),
                            'MOD_Refl':FLUXDataframe['sur_refl_b07_new_interp'].groupby([by('year'),by('dayofyear')]).mean()}).reset_index()
    
    #Generate dates from group names (year and day of year)
    ObsDailySum_DF.index=(ObsDailySum_DF.apply(lambda x: dt.datetime.strptime(str(int(x['level_0']))+
                    '-'+str(int(x['level_1'])),'%Y-%j'), axis=1))
    
    ObsDailySum_DF=ObsDailySum_DF.drop(['level_0','level_1'],axis=1)
    
    #Throw out data that don't have all records, reindex and add to MODIS DF
    ObsDailySum_DF['Fsd_mean']=np.where(ObsDailySum_DF['Fsd_count']==1440/rec_length,ObsDailySum_DF['Fsd_mean'],np.nan)
    ObsDailySum_DF[VarToProcess+'_mean']=np.where(ObsDailySum_DF[VarToProcess+'_count']==1440/rec_length,ObsDailySum_DF[VarToProcess+'_mean'],np.nan)
    ObsDailySum_DF['Albedo']=ObsDailySum_DF[VarToProcess+'_mean']/ObsDailySum_DF['Fsd_mean']

        
    ###### Linear transform of MODIS albedo data - year by year ######
    
    #Check the number of years with sufficient data for regression
    temp_DF=ObsDailySum_DF.groupby(lambda x: x.year).count()
    temp_DF['Test']=np.logical_and(temp_DF['MOD_Refl'].map(lambda x: x>10),temp_DF['Albedo'].map(lambda x: x>10))
    
    #Pass data to the regression function, return regression parameters and statistics and add to DF
    stats_DF=pd.concat([linear_trans(ObsDailySum_DF[['MOD_Refl','Albedo']].loc[str(i)],temp_DF['Test'].loc[i],i) for i in temp_DF.index])
    stats_DF.columns=['Slope','Intercept','r-squared','p-value','standard_error']
          
    #Linear transform MODIS data using regression parameters
    transform_DF=pd.concat([pd.DataFrame({'mod_Trans':ObsDailySum_DF['MOD_Refl'].loc[str(i)]*stats_DF['Slope'].loc[i]+stats_DF['Intercept'].loc[i]}) for i in stats_DF.index])
    
    #Merge
    ObsDailySum_DF=pd.merge(ObsDailySum_DF,transform_DF,left_index=True,right_index=True,how='left')
    
    #Combine adjusted and non-adjusted MODIS data
    ObsDailySum_DF['MODIS_combined']=np.where(pd.isnull(ObsDailySum_DF['mod_Trans']),ObsDailySum_DF['MOD_Refl'],ObsDailySum_DF['mod_Trans'])
    
    #eturn ObsDailySum_DF
    
    
    ###### Apply gap-filled MODIS-derived albedo to incoming shortwave radiation ######
    
    #Reset both dataframes to hierarchical index (year,day of year)
    ObsDailySum_DF['DT']=ObsDailySum_DF.index
    ObsDailySum_DF.index=pd.MultiIndex.from_tuples(zip(*[np.array(ObsDailySum_DF['DT'].apply(lambda x: int(dt.datetime.strftime(x,'%Y')))),
                                                np.array(ObsDailySum_DF['DT'].apply(lambda x: int(dt.datetime.strftime(x,'%j'))))]),names=['Year', 'DOY'])
    
    FLUXDataframe['DT']=FLUXDataframe.index
    FLUXDataframe.index=pd.MultiIndex.from_tuples(zip(*[np.array(FLUXDataframe['DT'].apply(lambda x: int(dt.datetime.strftime(x,'%Y')))),
                                                        np.array(FLUXDataframe['DT'].apply(lambda x: int(dt.datetime.strftime(x,'%j'))))]),names=['Year', 'DOY'])
    
    #Calculate outgoing shortwave radiation using MODIS derived albedo and incoming shortwave radiation
    output_Fsu=pd.concat([pd.DataFrame(FLUXDataframe['Fsd_Corr'].ix[i]*ObsDailySum_DF['MODIS_combined'].ix[i]).set_index(FLUXDataframe['DT'].ix[i]) for i in ObsDailySum_DF.index])
    
    #Reset the main dataframe index to original datetimeindex
    FLUXDataframe.index=FLUXDataframe['DT']
    FLUXDataframe.drop('DT',axis=1)
            
    ###########################################    
    # Generate output and return new variables
    ############################################
    if write_to_DF==True:
   	
        #Create a new label for pandas df column for the contructed variable (and the QC flag) and column to fill label
        construct_label=str(VarToProcess+"_Con")
        FLUXDataframe[construct_label]=np.nan
        #add a column for the constructed QC flag
        #This will be 1 if valid data from the tower else 99
        construct_flag_label=str(VarToProcess+"_Con_QCFlag")
        FLUXDataframe[construct_flag_label]=np.nan
        #Create a series that is just the correlated output of CABLE adjusted to Tower
        corr_label=str(VarToProcess+"_Corr")
        FLUXDataframe[corr_label]=np.nan	
        
        #Start by filling with valid  values from tower
        FLUXDataframe[construct_flag_label][FLUXDataframe[VarToProcess].notnull()]=1
        FLUXDataframe[construct_label][FLUXDataframe[VarToProcess].notnull()]=FLUXDataframe[VarToProcess]
        #Fill with series calculated here
        FLUXDataframe[construct_flag_label][FLUXDataframe[construct_flag_label].isnull()]=99
        
        FLUXDataframe[corr_label]=output_Fsu
        FLUXDataframe[construct_label][FLUXDataframe[construct_flag_label]==99]=FLUXDataframe[corr_label]
                        
        #FLUXDataframe[construct_label][FLUXDataframe[construct_flag_label]==99]=output_Fsu
        #FLUXDataframe[corr_label]=output_Fsu
	
    ###########################################
    # Call function to do climatology gap fill#
    ###########################################
    #Do climatological gap filling if stillmissing data
    FLUXDataframe=gap_fill_climatology.climatology_monthly_diurnal(FLUXDataframe,VarToProcess)    

    print '\n Finished Fsu gapfilling returning DF...'
    return FLUXDataframe
