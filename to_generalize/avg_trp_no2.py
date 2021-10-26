
# This code is written to process monthly averaged map for both wrfchem and TROPOMI
# --- Meng Li, 2019. 5. 9
# --- Contact: meng.li@noaa.gov; meng.li.atm@gmail.com

import os
import numpy as np
from netCDF4 import Dataset
import wrf
from os import listdir
from os.path import isfile

from OutputPlot_Config import output_config

# baseline + lightning
Basedir_wrfoutput = '/scratch1/BMC/rcm2/mli/nyc18_lightning/run_12km_five18_bmcdVCP_fog_wofire_BEIS_0.5ISO/Output/'
Baseoutdir = '/scratch1/BMC/rcm2/mli/outdir_12km_noPM_baseline_bmc_cams/'

year = 2018
seasoname = 'm07'

'''
====================================================================
Comparison Between TROPOMI and WRF-Chem NO2 Trop. Columns Seasonaly
====================================================================
'''

Month = [7]
MonthStartDay = [1] # start day for each month, 1-based
MonthEndDay = [15]  # end day for each month, 1-based


# Use wind speed as a criterion?
UseWPD = False  # False: not use wind speed

#================Preparation Codes======================
#---
class file_management:
    def __init__(self):
        pass
    def subdirlist(self, indir, keyword=''):
        subdirlist = []
        print(indir)
        subdirlist_org = [x[0] for x in os.walk(indir)]
        for sd in subdirlist_org:
            if keyword in sd:
                subdirlist.append(sd)
        return subdirlist
                
        return subdirlist
    def filelist(self, indir, keyword=''):
        filelist = []
        filelist_org = [os.path.join(indir, f) for f in listdir(indir) if isfile(os.path.join(indir,f))]
        for f in filelist_org:
            if keyword in f:
                filelist.append(f)
        return filelist   
    
#---
def extractwrfcoord(lats='', lons=''):
    # extract one wrfdata
    fm = file_management()
    subdirlist = fm.subdirlist(Basedir_wrfoutput)
    ff = fm.filelist(subdirlist[1])[0]
    wrfin = Dataset(ff,'r',format = 'NETCDF4_CLASSIC')
        
    # get some attributes of the wrf domain 
    latdata = wrf.getvar(wrfin, 'XLAT', timeidx=0)[:,:]            # latitude
    londata = wrf.getvar(wrfin, 'XLONG', timeidx=0)[:,:]           # longitude
    wrflonlat = {'lon':londata, 'lat':latdata}

    if (lats == '') & (lons == ''):
        return wrflonlat
    else:
        xyinds = wrf.ll_to_xy(wrfin, lats, lons)
        return xyinds
    wrfin.close()

#==============MAIN PROGRAM STARTS HERE==============
    
# GET THE WRF COORDIATE INFORMATION OF WRF-CHEM
wrfcoord = extractwrfcoord()
wrflonlat = wrfcoord
# extract locations
wrflon = wrflonlat['lon']
wrflat = wrflonlat['lat']
xy = np.shape(wrflonlat['lon'])
 
# MAIN PROGRAM
def main():
    m = model_validation(year)
    m.evaluatedata()
        
class model_validation():
    
    def __init__(self, year):        
        self.year = year
        
    def evaluatedata(self):
        year = self.year
        
        # define the outdir based on use wind speed or not
        if UseWPD == False:
            Outdir = Baseoutdir + seasoname + '/'
        else:
            Outdir = Baseoutdir + seasoname + '/' + 'wpduvle'+str(maxwd) + '/'
        if os.path.isdir(Outdir):
            pass
        else:
            os.mkdir(Outdir)
        print('***Evaluation starts here: ', year, seasoname)
      
        # initialize data array
        no2_tomi = np.zeros([xy[0], xy[1]], dtype = np.float32) #TROPOMI NO2 columns for further sum
        no2_tomi[:,:] = 0.0       
        no2_wrfchem = np.zeros([xy[0], xy[1]], dtype = np.float32) #WRF-Chem NO2 columns for further sum
        no2_wrfchem[:,:] = 0.0
        
        num_tomi = np.zeros([xy[0], xy[1]], dtype = np.float32) #Number of observations
        num_tomi[:,:] = 0.0        
        num_wrfchem = np.zeros([xy[0], xy[1]], dtype = np.float32) #Number of observations or WRF-Chem, should be the same
        num_wrfchem[:,:] = 0.0

        no2_tomi_avg = np.zeros([xy[0], xy[1]], dtype = np.float32) #seasonal averaged TROPOMI NO2 columns
        no2_tomi_avg[:,:] = np.nan       
        no2_wrfchem_avg = np.zeros([xy[0], xy[1]], dtype = np.float32) #seasonal averaged WRF-Chem NO2 columns
        no2_wrfchem_avg[:,:] = np.nan
        
        # summerize each day
        for mind in range(len(Month)):
            month = Month[mind]
            daymin = MonthStartDay[mind]
            daymax = MonthEndDay[mind]

            for day in range(daymin,daymax+1):

                Indir = Baseoutdir + '{:02d}'.format(month) + '{:02d}'.format(day)+'/'
                fn = Indir+ 'no2_wrfchem_'+str(year)+'_'+str(month)+'_'+str(day)+'.nc'
                if isfile(fn):
                    # wrfchem daily no2 column 
                    fn = Indir+ 'no2_wrfchem_'+str(year)+'_'+str(month)+'_'+str(day)+'.nc'
                    print('--> reading wrfchem datafile of : ', fn)
                    ds = Dataset(fn,"r")
                    variable_wc = ds.variables['NO2'][:,:]
                    ds.close()

                    # tropomi daily no2 column
                    fn = Indir+ 'no2_tropomi_wchamf_'+str(year)+'_'+str(month)+'_'+str(day)+'.nc'
                    print('--> reading tropomi datafile of : ', fn)
                    ds = Dataset(fn,"r")
                    variable_tp = ds.variables['NO2'][:,:]
                    ds.close() 

                    if UseWPD == False:
                        # add to summary array
                        ind = np.where((variable_wc >= 0.0) & (variable_tp >= 0.0) & (variable_tp != np.nan))
                        print('check', ind)
                        print('wc',np.nanmin(variable_wc), np.nanmax(variable_wc))
                        print('wc',np.nanmin(variable_tp), np.nanmax(variable_tp))
                        no2_wrfchem[ind] += variable_wc[ind]
                        num_wrfchem[ind] += 1.0
                
                        #ind = np.where(variable_tp >= 0.0)
                        no2_tomi[ind] += variable_tp[ind]
                        num_tomi[ind] += 1.0

                    else:
			            # read the surface wind speed file
                        fnwpd_u = Indir + 'wrfchem_u10_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'.nc'
                        ds = Dataset(fnwpd_u,"r")
                        variable_wpdu = ds.variables['u10'][:,:]
                        print('--> reading wrfchem u10 of :', fnwpd_u)
                        print('    ',np.nanmin(variable_wpdu), np.nanmax(variable_wpdu))
                        ds.close()
                
                        fnwpd_v = Indir + 'wrfchem_v10_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'.nc'
                        ds = Dataset(fnwpd_v,"r")
                        variable_wpdv = ds.variables['v10'][:,:]
                        print('--> reading wrfchem u10 of :', fnwpd_v)
                        print('    ', np.nanmin(variable_wpdv), np.nanmax(variable_wpdv))
                        ds.close()
                        # the surface wind speed
                        wd = (variable_wpdu**2 + variable_wpdv**2)**0.5

                        # add to summary array
                        ind = np.where((variable_wc >= 0.0) & (np.absolute(wd) <= maxwd) & (variable_tp >= 0.0))
                        no2_wrfchem[ind] += variable_wc[ind]
                        num_wrfchem[ind] += 1.0
                
                        #ind = np.where((variable_tp >= 0.0) & (np.absolute(wd) <= maxwd))
                        no2_tomi[ind] += variable_tp[ind]
                        num_tomi[ind] += 1.0
                else:
                    print('--> NO wrfchem / tropomi file found: ', day, fn)
                    pass

        # calculate seasonal average
        no2_tomi_avg = no2_tomi / num_tomi
        no2_wrfchem_avg = no2_wrfchem / num_wrfchem

        #----Output and Plotting ---
        ot = output_config()
        pmin = 0.0	# pcolormap, min
        pmax = 1e16	# pcolormap, max
     
        # TROPOMI NO2 for comparison 
        fnnc = Outdir+ 'no2_tropomi_wchamf_'+str(year)+'_'+ seasoname +'_'+'mavg'+'.nc'
        ot.outputnc_2d(fnnc, no2_tomi_avg, 'NO2', 'molec cm-2')

        fn = Outdir+ 'no2_tropomi_wchamf_'+str(year)+'_'+ seasoname +'_'+'mavg'
        ot.plot_2dmap(fn, no2_tomi_avg,'NO2','molec cm-2', mindata=pmin, maxdata = pmax)

        print('TROPOMI NO2 after AMF revision:')
        print('    x1e15: min ',np.nanmin(no2_tomi_avg)/1e15, ' max ',np.nanmax(no2_tomi_avg)/1e15, ' median ', np.nanmedian(no2_tomi_avg)/1e15, ' mean ',np.nanmean(no2_tomi_avg)/1e15)
              
        # WRF-Chem NO2 for comparison
        fnnc = Outdir+ 'no2_wrfchem_'+str(year)+'_'+ seasoname +'_'+'mavg'+'.nc'
        ot.outputnc_2d(fnnc, no2_wrfchem_avg, 'NO2', 'molec cm-2')

        fn = Outdir+ 'no2_wrfchem_'+str(year)+'_'+ seasoname +'_'+ 'mavg'
        ot.plot_2dmap(fn, no2_wrfchem_avg,'NO2','molec cm-2', mindata=pmin, maxdata = pmax)

        print('WRF-Chem NO2:')
        print('    x1e15: min ',np.nanmin(no2_wrfchem_avg)/1e15, ' max ',np.nanmax(no2_wrfchem_avg)/1e15, ' median ', np.nanmedian(no2_wrfchem_avg)/1e15, ' mean ',np.nanmean(no2_wrfchem_avg)/1e15               )

        # Observation numbers
        fnnc = Outdir+ 'num_tomi_no2_'+str(year)+'_'+ seasoname +'_'+'mavg'+'.nc'
        ot.outputnc_2d(fnnc, num_tomi, 'number', 'unitless')

        # Correlations
        cal = calstatis()
        ind = np.where((no2_tomi_avg > 0.0) & (no2_wrfchem_avg > 0.0))
        x = no2_tomi_avg[ind].flatten()
        y = no2_wrfchem_avg[ind].flatten()
        mb = cal.calmb(y,x)
        nmb = cal.calnmb(y,x)
        fn = Outdir + 'corrl_trop_wrfchem_no2_'+str(year)+'_'+ seasoname +'_'+'mavg'
        ot.plot_scatter(fn, 'NO2 columns', x, y, mindata=0.0, maxdata=3e16, mb=mb, nmb=nmb)

        # minus
        minarr = np.zeros([xy[0], xy[1]], dtype=np.float)
        minarr[:,:] = np.nan
        minarr[:,:] = (no2_wrfchem_avg[:,:] - no2_tomi_avg[:,:])/1e15
        fn = Outdir + 'minus_wrfchem-tropomi_'+str(year)+'_'+seasoname + '_'+'mavg'
        ot.plot_minus(fn,  minarr, 'NO2', '$\mathregular{10^{15}}$ molec c$\mathregular{m^{-2}}$',  mindata=-3.0, maxdata = 3.0)
        


class calstatis:
    def __init__(self):
        pass
    def calcorrelation(self,modellist, obslist):
        #r = np.corrcoef(x=modellist, y=obslist)
        slope, intercept, r, p_value, stderr=stats.linregress(obslist, modellist)
        return r
    def calmb(self, modellist, obslist):
        minuslist = [(modellist[n] - obslist[n]) for n in range(len(modellist))]
        mb = np.sum(minuslist) / len(modellist)
        return mb
    def calrmse(self, modellist, obslist):
        minuslist = [pow((modellist[n] - obslist[n]),2) for n in range(len(modellist))]
        rmse = pow((np.sum(minuslist) / len(modellist)),0.5)
        return rmse
    def calnmb(self, modellist, obslist):
        minuslist = [(modellist[n] - obslist[n]) for n in range(len(modellist))]
        nmb = np.sum(minuslist) / np.sum(obslist)
        return nmb
    def calnme(self, modellist, obslist):
        minuslist = [(abs(modellist[n] - obslist[n])) for n in range(len(modellist))]
        nme = np.sum(minuslist) / np.sum(obslist)
        return nme


#---
if __name__ == '__main__':
        main()
    
 



