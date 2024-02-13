
# This pro is written to read the wrf-chem output 
# and process the TROPOMI data, to compare the observation with the model
# --- Meng Li
# --- 2019.04.11
# --- Contact: meng.li@noaa.gov; meng.li.atm@gmail.com

'''
===================================================
Import necessary packages and set the environment
===================================================
'''

from os import listdir
from os.path import isfile
import csv, wrf,os,sys
import numpy as np
from netCDF4 import Dataset
import multiprocessing
from datetime import datetime
import pytz
from timezonefinder import TimezoneFinder
import ESMF
import xesmf as xe
import math
ESMF.Manager(debug=True)

# import outputplot_config file
from OutputPlot_Config import output_config


# Get Basedir_tropomi, Baseoutdir, Basedir_wrfoutput, and Geofile from environment variables
Basedir_tropomi = os.environ.get('Basedir_tropomi')
Baseoutdir = os.environ.get('Baseoutdir')
Basedir_wrfoutput = os.environ.get('Basedir_wrfoutput')
Geofile = os.environ.get('Geofile')


'''
=============================================================
Generate Daily Trop. NO2 Columns for TROPOMI and WRF-Chem
=============================================================
'''

#==============Preparation Codes==================
#---
class file_management:
    def __init__(self):
        pass
    def subdirlist(self, indir, keyword=''):
        subdirlist = []
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

def extractwrfcorners():
    ds = Dataset(Geofile, "r")
    variable_latc = ds['XLAT_C'][0,:,:]
    variable_longc = ds['XLONG_C'][0,:,:] #1,285,441
    cornerdic = {'lon_c': variable_longc, "lat_c": variable_latc}
    #print(variable_latc)
    ds.close()
    return cornerdic

def extractwrfcoord(lats=[''], lons=['']):
    # extract one wrfdata
    fm = file_management()
    subdirlist = fm.subdirlist(Basedir_wrfoutput)
    ff = fm.filelist(subdirlist[1])[0]
    wrfin = Dataset(ff,'r',format = 'NETCDF4_CLASSIC')
        
    # get some attributes of the wrf domain 
    latdata = wrf.getvar(wrfin, 'XLAT', timeidx=0)[:,:]            # latitude
    londata = wrf.getvar(wrfin, 'XLONG', timeidx=0)[:,:]           # longitude
    wrflonlat = {'lon':londata, 'lat':latdata}

    if (len(lats) == 1) & (lats[0] == '') & (len(lons) == 1) & (lons[0] == ''):
        cornerdic = extractwrfcorners()
        wrflonlat.update(cornerdic)
        #print(wrflonlat)
        return wrflonlat
    else:
        xyinds = wrf.ll_to_xy(wrfin, lats, lons)
        return xyinds
    wrfin.close()


#=======MAIN PROGRAM STARTS HERE=============
    
# GET THE WRF COORDIATE INFORMATION
wrfcoord = extractwrfcoord()
wrflon = wrfcoord['lon']
wrflat = wrfcoord['lat']
wrflon_c = wrfcoord['lon_c']
wrflat_c = wrfcoord['lat_c']
xy = np.shape(wrflon)
tf = TimezoneFinder()

print(wrflon_c)
print(wrflat_c)

# MAIN PROGRAM
def main(year, month, day):
    m = model_validation(year, month, day)
    m.evaluatedata()
    pass
#---
class model_validation():
    
    def __init__(self, year, month, day):        
        self.year = year
        self.month = month
        self.day = day

    
    def findxy(self,coord_lons, coord_lats, lons, lats, iswrf):        
        arrayout = extractwrfcoord(lats=lats, lons=lons)
        arrayout = wrf.to_np(arrayout) # EDGE: coord_lons[0][0], coord_lats[0],[0], POINT OF (0,0)          
        return arrayout

    def extractloc(self, wrflon, wrflat, lons,lats):
        londata = wrflon
        latdata = wrflat
        lats_1d = lats.flatten() # change 2-d to 1-d, dims of lons and lats should be the same
        lons_1d  = lons.flatten() # change 2-d to 1-d       
        idx = self.findxy(londata, latdata, lons_1d, lats_1d, 1)
        idx = np.reshape(idx, [2, np.shape(lons)[0], np.shape(lons)[1]])
        return idx
        
    def evaluatedata(self):
        ot = output_config()
        wrflonlat = wrfcoord
        year = self.year
        month = self.month
        day = self.day

        Outdir = Baseoutdir + '{:02d}'.format(month) + '{:02d}'.format(day)+'/'
        # check if Outdir exsits, if not, create a new one.
        if os.path.isdir(Outdir):
            pass
        else:
            os.mkdir(Outdir)
        print('*** Evaluation starts here: ', year, month, day)

        # extract wrf-chem data directionary for each day
        w = wrf_chem_process()
        w.extractwrfdata(year, month, day)
        wrf_omi_avg =  w.wrf_omi_avg

        # no2 column and pressure of wrf-chem
        wrfpres = wrf_omi_avg['pres']
        wrfno2_omi_avg = wrf_omi_avg['no2']
        wrfno2_omi_ppm = wrf_omi_avg['no2ppm']

        # get the wrf grid cell center and boundaries
        lat_wrf = wrflat
        lon_wrf = wrflon
        lat_wrf_b = wrflat_c
        lon_wrf_b = wrflon_c
        # extract the tropomi data dictionary including all swath tracks covering US
        t = tropomi_process(year, month, day)
        t.avgtropomi()
        omidata_alltrack = t.omidata_alltrack
        
        print('--> total track number of tropomi is: '+str(len(omidata_alltrack)))
        
        # initialize arrays
        no2colomi_avg = np.zeros([xy[0], xy[1], len(omidata_alltrack)], dtype = float) # NO2 column for sum in standard product 
        no2colomi_avg[:,:,:] = np.nan

        no2numomi_avg = np.zeros([xy[0], xy[1], len(omidata_alltrack)], dtype = float) # NO2 number for sum in standard product
        no2numomi_avg[:,:,:] = 0.0

        ratio = np.zeros([xy[0], xy[1], len(omidata_alltrack)], dtype = float) # ratio of revised column / standard column
        ratio[:,:,:] = np.nan         

        tamf_omi = np.zeros([xy[0], xy[1], len(omidata_alltrack)], dtype = float) # trop. AMF in standard product
        tamf_omi[:,:,:] = np.nan         
        
        amf_model = np.zeros([xy[0], xy[1], len(omidata_alltrack)], dtype = float) # revised trop. AMF 
        amf_model[:,:,:] = np.nan  

        amf_total = np.zeros([xy[0], xy[1], len(omidata_alltrack)], dtype = float) # total AMF in standard product
        amf_total[:,:,:] = np.nan

        slantcol = np.zeros([xy[0], xy[1], len(omidata_alltrack)], dtype = float) # slant column density for sum in standard product
        slantcol[:,:,:] = 0.0
        
        zmax_model = np.zeros([xy[0], xy[1]], dtype = np.int32) # vertical layer
        zmax_model[:,:] = 0.0

        no2_wrf_forcmp = np.zeros([xy[0], xy[1]], dtype = np.float) # final wrf-chem no2 column for comparison
        no2_wrf_forcmp[:,:] = np.nan
        
        # processing each track into wrf domain and average 
        for trk in range(len(omidata_alltrack)):
            print('   -----> prosessing the track of ', trk)
            omidata = omidata_alltrack[trk]
            
            # cell centers and boundaries
            lon = np.squeeze(omidata['longitude'][:,:,:], axis=0)
            lat = np.squeeze(omidata['latitude'][:,:,:], axis=0)     
            lon_b = np.squeeze(omidata['longitude_bounds'][:,:,:,:], axis=0) # 3245x450x4
            lat_b = np.squeeze(omidata['latitude_bounds'][:,:,:,:], axis=0)

            locind = self.extractloc(wrflon, wrflat, lon, lat)  # locind: the value of x_wrf and y_wrf, omi domain shape

            # extract data in standard NO2 product
            no2colorg =np.squeeze(omidata['nitrogendioxide_tropospheric_column'][:,:,:], axis=0)# Trop. NO2 column in standard product
            tnx = np.shape(lat)[0]
            tny = np.shape(lat)[1]
            lat_trp_b = np.zeros([tnx+1, tny+1], dtype = np.float)
            lon_trp_b = np.zeros([tnx+1, tny+1], dtype = np.float)

            # Revised NO2 Tropomi product at original dimension
            no2colrev = np.zeros([tnx, tny],dtype=np.float)
            no2colrev[:,:] = np.nan

            qaorg = np.squeeze(omidata['qa_value'][:,:,:], axis=0)                              # quality flag
            tamf_tm5 = np.squeeze(omidata['air_mass_factor_troposphere'][:,:,:], axis=0)        # TM5 trop. AMF 
            amf_tm5 = np.squeeze(omidata['air_mass_factor_total'][:,:,:], axis=0)               # TM5 total AMF
            cldfrc = np.squeeze(omidata['cloud_fraction_crb'][:,:,:],axis=0)                    # cloud fraction
            tpreslev_tm5 = np.squeeze(omidata['preslev'][:,:,:], axis=0)                        # TM5 pressure level
            trplayer_tm5 = np.squeeze(omidata['tm5_tropopause_layer_index'][:,:,:], axis=0)     # TM5 tropopause  
            slant_col_single = np.squeeze(omidata['nitrogendioxide_slant_column_density'][:,:,:],axis=0)# NO2 Slant column density in standard product
            scatwts = np.squeeze(omidata['averaging_kernel'][:,:,:,:],axis=0)                   # averaging kernel


            for x_tomi in range(tnx):
                for y_tomi in range(tny):
                    x_wrf = locind[1,x_tomi, y_tomi] # MAX: 299, SOUTH-NORTH
                    y_wrf = locind[0,x_tomi, y_tomi] # MIN: 239, WEST-EAST
                    # get the cell boundaries of each swath
                    lat_0 = lat_b[x_tomi,y_tomi,0]
                    lat_1 = lat_b[x_tomi,y_tomi,1]
                    lat_2 = lat_b[x_tomi,y_tomi,2]
                    lat_3 = lat_b[x_tomi,y_tomi,3]

                    lon_0 = lon_b[x_tomi,y_tomi,0]
                    lon_1 = lon_b[x_tomi,y_tomi,1]
                    lon_2 = lon_b[x_tomi,y_tomi,2]
                    lon_3 = lon_b[x_tomi,y_tomi,3]

                    lat_trp_b[x_tomi,y_tomi] = lat_0
                    lat_trp_b[x_tomi,y_tomi+1] = lat_1
                    lat_trp_b[x_tomi+1,y_tomi+1] = lat_2
                    lat_trp_b[x_tomi+1,y_tomi] = lat_3

                    lon_trp_b[x_tomi,y_tomi] = lon_0
                    lon_trp_b[x_tomi,y_tomi+1] = lon_1
                    lon_trp_b[x_tomi+1,y_tomi+1] = lon_2
                    lon_trp_b[x_tomi+1,y_tomi] = lon_3
   
                    if (x_wrf < 0.0) or (x_wrf > xy[0]-1) or (y_wrf < 0.0) or (y_wrf > xy[1]-1):
                        pass
                    else:

                        # extract the no2 trop column, quality, cloud radiation fraction                            
                        value_tno2 = no2colorg[x_tomi, y_tomi]
                        value_slnt = slant_col_single[x_tomi, y_tomi]
                                            
                        # screen data
                        if qaorg[x_tomi, y_tomi] >= 0.75 and value_tno2 > 0.0e-30 and cldfrc[x_tomi, y_tomi] <= 0.5: # QUALITY CONTROL
                            value_tno2 = self.converunit(value_tno2, 'mole m-2', 'molec cm-2') # CONVERT THE UNIT
                            value_slnt = self.converunit(value_slnt, 'mole m-2', 'molec cm-2') # CONVERT THE UNIT

                            # add slant NO2 column for further averaging
                            no2numomi_avg[x_wrf, y_wrf, trk] += 1.0
                            slantcol[x_wrf, y_wrf, trk] += value_slnt
                            if value_slnt < 0.0:
                                print('Error for slant column here!' )
        
                            # revise amf using wrfchem NO2 vertical profile, start here
                            scatwts_vertical = scatwts[x_tomi, y_tomi, :]
                            tpreslev = tpreslev_tm5[x_tomi, y_tomi,:]
                            trplayer = trplayer_tm5[x_tomi, y_tomi]

                            wrfpreslayer = wrfpres[:,x_wrf, y_wrf]
                            wrfno2layer_molec = wrfno2_omi_avg[:,x_wrf, y_wrf] # mole cm^-2 by WRF layers
                            wrfno2layer = wrfno2_omi_ppm[:,x_wrf, y_wrf] # use unit of ppm to derive NO2 profile

                            # trop. AMF and total AMF in standard product
                            tamf_org = tamf_tm5[x_tomi, y_tomi]    # trop. amf
                            tamf_omi[x_wrf, y_wrf, trk] = tamf_org # add trop. amf to array
                            amf_total[x_wrf, y_wrf, trk] = amf_tm5[x_tomi, y_tomi] # add total amf to array

                         
                            # find the vertical index of wrf-chem corresponding to the tropomi tropopause
                            if type(trplayer) == np.int32:
                                X = abs(wrfpreslayer - tpreslev[trplayer])
                                zm_wrf = np.where(X == np.min(X))
                                zmax_model[x_wrf, y_wrf] = zm_wrf[0][0]

                                # calculate the revised trop. AMF, amf_model
                                scatwts_vertical = scatwts_vertical * amf_total[x_wrf,y_wrf,trk] / tamf_org # converting from AKs to tropospheric AKs
                                amf_model[x_wrf, y_wrf,trk] = self.calamfwrfchem(scatwts_vertical, wrfpreslayer, wrfno2layer, tpreslev, trplayer, zm_wrf[0][0], wrfno2layer_molec)*tamf_org

                                # summarize all columns in WRF-Chem from surface to the tropopause
                                ratio[x_wrf,y_wrf, trk] = tamf_org/amf_model[x_wrf, y_wrf, trk]
                                
                            else:
                                #no2wrf_forcmp[x_wrf, y_wrf, trk] = np.nansum(wrfno2layer[:])
                                ratio[x_wrf, y_wrf, trk] = 1.0


                            no2colrev[x_tomi, y_tomi] = value_tno2*ratio[x_wrf,y_wrf,trk]
            
                        else:
                            no2colrev[x_tomi, y_tomi] = np.nan

            
            # Regrid from revised TROPOMI to WRF-Chem grid, conservative method
            # Refs: https://xesmf.readthedocs.io/en/latest/notebooks/Pure_numpy.html?highlight=conservative#Regridding
            lon_wrf_value = lon_wrf.values
            lat_wrf_value = lat_wrf.values
            grid_in={'lon': lon, 'lat': lat, 'lon_b': lon_trp_b, 'lat_b': lat_trp_b }
            grid_out ={'lon': lon_wrf_value, 'lat': lat_wrf_value, 'lon_b': lon_wrf_b, 'lat_b': lat_wrf_b}
            regridder = xe.Regridder(grid_in, grid_out, 'conservative', ignore_degenerate=True)
            no2_trp_regrid = regridder(no2colrev)

            ind = np.where(no2_trp_regrid <= 0.0e-30)
            if (ind != []):
                no2_trp_regrid[ind] = np.nan

            no2colomi_avg[:,:,trk] = no2_trp_regrid[:,:]
            print('regridded NO2 column', np.nanmin(no2_trp_regrid), np.nanmax(no2_trp_regrid))


        #-- final averaged data ---
        # averaged slant NO2 columns in standard product, WRF-Chem domain
        slantcol2 = slantcol / no2numomi_avg

        # WRF-Chem arrays for comparison
        zmax_default = np.nanmax(zmax_model)
        print('zmax_default', zmax_default, np.nanmin(zmax_model))
        no2_wrf_forcmp[:, :] = np.nansum(wrfno2_omi_avg[0:zmax_default+1,:,:],axis=0)

        no2_omi_forcmp = np.nanmean(no2colomi_avg, axis=2)
        slantcol_forcmp = np.nanmean(slantcol2, axis=2)        

        print('OMI NO2 after AMF replacement:', np.nanmin(no2_omi_forcmp)/1e15, np.nanmax(no2_omi_forcmp)/1e15 )
        print('WRF-Chem NO2:', np.nanmin(no2_wrf_forcmp)/1e15, np.nanmax(no2_wrf_forcmp)/1e15)
         
      
       #----------------Output and Plotting-------------------
       # Configured in OutputPlot_Config.py     

        fnnc = Outdir+ 'no2_tropomi_wchamf_'+str(year)+'_'+str(month)+'_'+str(day)+'.nc'
        ot.outputnc_2d(fnnc, no2_omi_forcmp, 'NO2', 'molec cm-2')

        fnnc = Outdir+ 'no2_wrfchem_'+str(year)+'_'+str(month)+'_'+str(day)+'.nc'
        ot.outputnc_2d(fnnc, no2_wrf_forcmp, 'NO2', 'molec cm-2')
       
        outdata = np.nanmean(amf_model, axis=2)
        fnnc = Outdir + 'amf_model_'+str(year)+'_'+str(month)+'_'+str(day) + '.nc'
        ot.outputnc_2d(fnnc, outdata, 'amf_model', '1.0') 

        outdata = np.nanmean(tamf_omi, axis=2)
        fnnc = Outdir + 'amf_omi_'+str(year)+'_'+str(month)+'_'+str(day)
        ot.outputnc_2d(fnnc, outdata, 'amf_omi', '1.0')                
      
        outdata = np.nanmean(amf_total, axis=2)
        fnnc = Outdir + 'totalamf_omi_'+str(year)+'_'+str(month)+'_'+str(day)+'.nc'
        ot.outputnc_2d(fnnc, outdata, 'totalamf_omi', '1.0')

        fnnc = Outdir+ 'no2_total_slantcol_'+str(year)+'_'+str(month)+'_'+str(day)+'.nc'
        ot.outputnc_2d(fnnc,slantcol_forcmp , 'NO2', 'molec cm-2')
      
        fnnc = Outdir+ 'zmax_'+str(year)+'_'+str(month)+'_'+str(day)+'.nc'
        ot.outputnc_2d(fnnc, zmax_model, 'zmax', '1')


    def converunit(self, invalue, inunit, outunit):
        if (inunit == 'mole m-2' and outunit == 'molec cm-2'):
            #outvalue = invalue * 1.0*6.02e23/100.0/100.0
            outvalue = invalue * 6.02214e+19 # provided by the TROPOMI
        else:
            print('no unit is changed', invalue)
        return outvalue
    
    def calamfwrfchem_tm5(self, scatw, wrfpreslayer, wrfno2layer, tpreslev, trplayer):
        from scipy import interpolate 
        nume = 0.0
        deno = 0.0
        amf_wrfchem = np.nan
        
        f = interpolate.interp1d(np.log10(wrfpreslayer),wrfno2layer, fill_value="extrapolate")      
        wrfno2_omi_avg_tm5 = f(np.log10(tpreslev[:]))
                   
        if type(trplayer) != np.int32:
            amf_wrfchem = np.nan
        else:
            for l in range(trplayer+1): # add all tropospheric layers
                if l == 0:
                    deltapres = tpreslev[l]-tpreslev[l+1]
                else:
                    deltapres = tpreslev[l-1]-tpreslev[l]
                nume += scatw[l] *wrfno2_omi_avg_tm5[l]*deltapres
                deno += wrfno2_omi_avg_tm5[l]*deltapres
                
            amf_wrfchem = nume / deno
        return amf_wrfchem
 
    def calamfwrfchem(self, scatw, wrfpreslayer, wrfno2layer, tpreslev, trplayer, zwrftrop, wrfno2layer_molec):
        from scipy import interpolate 
        nume = 0.0
        deno = 0.0
        amf_wrfchem = np.nan
        
        tpreslev[0] = wrfpreslayer[0] # set the surface pressure to wrf one
        f = interpolate.interp1d(np.log10(tpreslev),scatw, fill_value="extrapolate")# relationship between pressure to avk
        wrfavk = f(np.log10(wrfpreslayer[:])) #wrf-chem averaging kernel
        
        # add all tropospheric layers in WRF
        for l in range(zwrftrop+1):
            if (np.isinf(wrfavk[l]) == True) | (wrfavk[l] <= 0.0):
                nume += wrfavk[l+1]*wrfno2layer_molec[l]
                deno += wrfno2layer_molec[l]
                print('error of ak here', l, wrfavk[l], wrfavk[l+1], wrfno2layer_molec[l])
            else:
                nume += wrfavk[l]*wrfno2layer_molec[l]
                deno += wrfno2layer_molec[l]
                
        amf_wrfchem = nume / deno
        return amf_wrfchem
   
#---    
class wrf_chem_process:
    def __init__(self):
        self.wrf_omi_avg = {}
         
    def extractwrfdata(self, year, month ,day):
        # extract the NO2 column and Pressure from WRF-Chem
        # average between 13:00 and 14:00 localtime

        fm = file_management()
        keywords = '{:02d}'.format(month) + '{:02d}'.format(day)
        subdirlist = fm.subdirlist(Basedir_wrfoutput, keyword = keywords)
        subdir = subdirlist[0]        

        if len(subdirlist) > 1:
            print('Warning: more than one directories for day of ' + keywords)
            
        infile_wrf_12 = os.path.join(subdir, 'wrfout_d01_'+str(year) + '-'+ '{:02d}'.format(month) + '-'+'{:02d}'.format(day)+'_'+'12:00:00')
        infile_wrf_18 = os.path.join(subdir, 'wrfout_d01_'+str(year) + '-'+ '{:02d}'.format(month) + '-'+'{:02d}'.format(day)+'_'+'18:00:00')
        wrfdata_perfile_12 = self.readwrfoutput(infile_wrf_12)
        wrfdata_perfile_18 = self.readwrfoutput(infile_wrf_18)

        no2_perfile_12 = wrfdata_perfile_12['no2']
        no2_perfile_18 = wrfdata_perfile_18['no2']

        pres_perfile_12 = wrfdata_perfile_12['pb2']
        pres_perfile_18 = wrfdata_perfile_18['pb2']

        ph_perfile_12 = wrfdata_perfile_12['ph']
        ph_perfile_18 = wrfdata_perfile_18['ph']

        phb_perfile_12 = wrfdata_perfile_12['phb']
        phb_perfile_18 = wrfdata_perfile_18['phb']

        t_perfile_12 = wrfdata_perfile_12['t2']
        t_perfile_18 = wrfdata_perfile_18['t2']
        
        layers = np.shape(no2_perfile_12)[1]
        
        wrf_omi_all_no2 = np.zeros([layers,xy[0],xy[1],2], dtype = np.float)
        wrf_omi_ppm_no2 = np.zeros([layers,xy[0],xy[1],2], dtype = np.float)
        wrf_omi_all_pres = np.zeros([layers,xy[0],xy[1],2], dtype = np.float)
        wrf_omi_all_no2[:,:,:,:] = np.nan
        wrf_omi_all_pres[:,:,:,:] = np.nan

        no2_u2  = np.zeros([layers,xy[0],xy[1],2], dtype = np.float)
        #pres_u2 = np.zeros([6, layers,xy[0],xy[1],2], dtype = np.float)
        ph_u2   = np.zeros([layers+1,xy[0],xy[1],2], dtype = np.float)
        phb_u2  = np.zeros([layers+1,xy[0],xy[1],2], dtype = np.float)
        t_u2  = np.zeros([layers,xy[0],xy[1],2], dtype = np.float)
        
        no2_u2[:,:,:,:] = np.nan
        #pres_u2[:,:,:,:] = np.nan
        ph_u2[:,:,:,:] = np.nan
        phb_u2[:,:,:,:] = np.nan
        t_u2[:,:,:,:] = np.nan

        for x in range(xy[0]):
            for y in range(xy[1]):
                # Get UTC time at local time 13:00 and 14:00 for each wrf grid
                lat = wrfcoord['lat'][x,y]
                lon = wrfcoord['lon'][x,y]
                utc_13 = 13 - self.utchourfromlatlon(lat, lon, year, month, day) # Local - UTC
                utc_14 = utc_13 + 1

                # read the corresponding file, 12 or 18
                if utc_13 >= 12.0 and utc_13 < 18.0:
                    hour_u13 = int(utc_13-12)

                    no2_u2[:,x,y,0]  = no2_perfile_12[hour_u13,:,x,y]
                    wrf_omi_all_pres[:,x,y,0] = pres_perfile_12[hour_u13,:,x,y]

                    ph_u2[:,x,y,0]   = ph_perfile_12[hour_u13,:,x,y]
                    phb_u2[:,x,y,0]  = phb_perfile_12[hour_u13,:,x,y]
                    t_u2[:,x,y,0]    = t_perfile_12[hour_u13,:,x,y]
                    
                    if utc_14 >= 12.0 and utc_14 < 18.0:
                        hour_u14 = int(utc_14 - 12)
                        no2_u2[:,x,y,1] = no2_perfile_12[hour_u14,:,x,y]
                        wrf_omi_all_pres[:,x,y,1] = pres_perfile_12[hour_u14,:,x,y]

                        ph_u2[:,x,y,1]  = ph_perfile_12[hour_u14,:,x,y]
                        phb_u2[:,x,y,1] = phb_perfile_12[hour_u14,:,x,y]
                        t_u2[:,x,y,1] = t_perfile_12[hour_u14,:,x,y]

                    else:
                        hour_u14 = int(utc_14 - 18)
                        no2_u2[:,x,y,1]  = no2_perfile_18[hour_u14,:,x,y]
                        wrf_omi_all_pres[:,x,y,1] = pres_perfile_18[hour_u14,:,x,y]

                        ph_u2[:,x,y,1]  = ph_perfile_18[hour_u14,:,x,y]
                        phb_u2[:,x,y,1] = phb_perfile_18[hour_u14,:,x,y]

                        t_u2[:,x,y,1] = t_perfile_18[hour_u14,:,x,y]
                else:
                    hour_u13 = int(utc_13-18)
                    no2_u2[:,x,y,0] = no2_perfile_18[hour_u13,:,x,y]
                    wrf_omi_all_pres[:,x,y,0] = pres_perfile_18[hour_u13,:,x,y]

                    ph_u2[:,x,y,0] = ph_perfile_18[hour_u13,:,x,y]
                    phb_u2[:,x,y,0] = phb_perfile_18[hour_u13,:,x,y]

                    t_u2[:,x,y,0] = t_perfile_18[hour_u13,:,x,y]
                    
                    hour_u14 = int(utc_14-18)
                    no2_u2[:,x,y,1]  = no2_perfile_18[hour_u14,:,x,y]
                    wrf_omi_all_pres[:,x,y,1] = pres_perfile_18[hour_u14,:,x,y]

                    ph_u2[:,x,y,1] = ph_perfile_18[hour_u14,:,x,y]
                    phb_u2[:,x,y,1] = phb_perfile_18[hour_u14,:,x,y]

                    t_u2[:,x,y,1] = t_perfile_18[hour_u14,:,x,y]

        # calculatethe NO2 column, convert the unit
        for l in range(layers):
            ad = wrf_omi_all_pres[l,:,:,0]*(28.97e-3)/(8.314*t_u2[l,:,:,0])
            #print('ad', np.nanmin(ad), np.nanmax(ad))
            zh = ((ph_u2[l+1,:,:,0] + phb_u2[l+1,:,:,0]) - (ph_u2[l,:,:,0]+phb_u2[l,:,:,0]))/9.81
            wrf_omi_all_no2[l,:,:,0] = no2_u2[l,:,:,0]*zh*6.022e23/(28.97e-3)*1e-10*ad[:,:]
            wrf_omi_ppm_no2[l,:,:,0] = no2_u2[l,:,:,0]
            ad2 = wrf_omi_all_pres[l,:,:,1]*(28.97e-3)/(8.314*t_u2[l,:,:,1])
            zh2 = ((ph_u2[l+1,:,:,1] + phb_u2[l+1,:,:,1]) - (ph_u2[l,:,:,1]+phb_u2[l,:,:,1]))/9.81
            wrf_omi_all_no2[l,:,:,1] = no2_u2[l,:,:,1]*zh2*6.022e23/(28.97e-3)*1e-10*ad2[:,:]
            wrf_omi_ppm_no2[l,:,:,1] = no2_u2[l,:,:,1]

               
       # average the wrfno2_omi between 13:00 and 14:00 localtime
        wrf_omi_avg_no2 = np.nanmean(wrf_omi_all_no2, axis=3)
        wrf_omi_avg_pres = np.nanmean(wrf_omi_all_pres, axis=3)
        wrf_omi_avg_no2_ppm = np.nanmean(wrf_omi_ppm_no2, axis=3)

        self.wrf_omi_avg['no2'] = wrf_omi_avg_no2[:,:,:]
        self.wrf_omi_avg['pres'] = wrf_omi_avg_pres[:,:,:]
        self.wrf_omi_avg['no2ppm'] = wrf_omi_avg_no2_ppm[:,:,:] 

        print('WRF-Chem:', np.nanmin(wrf_omi_avg_no2), np.nanmax(wrf_omi_avg_no2), np.nanmin(wrf_omi_avg_pres), np.nanmax(wrf_omi_avg_pres))
      

    def utchourfromlatlon(self,lat,lon,year, month, day):
        # convert between UTC time and local time according to lon and lat

        lon2 = lon.values.tolist()
        lat2 = lat.values.tolist()
        timezone_str = tf.timezone_at(lng=lon2, lat=lat2)
        
        if timezone_str == None:
            if lon > -100.0:
                timezone_str = 'America/New_York'
            else:
                timezone_str = 'America/Los_Angeles'
        tz = pytz.timezone(timezone_str)
        d = datetime(year,month, day, 00,00,00)
        uos = tz.utcoffset(d, is_dst=True)
        utchour = uos.seconds/60.0/60.0
        utcday = uos.days
        if utcday < 0:
            utchour = (24-utchour)*-1 # Local - UTC
        return utchour
              

    def readwrfoutput(self, infile):
        # read wrf-chem output, return wrf-chem data directory

        print('--> reading wrf output file of : ', infile   )
        ncfile = Dataset(infile,'r',format = 'NETCDF4_CLASSIC')
        wrfdata = {}
        
        no2data = wrf.getvar(ncfile, 'no2',timeidx=wrf.ALL_TIMES)  # NO2 Mixing ratio, ppmv
        tdata = wrf.getvar(ncfile, 'T',timeidx=wrf.ALL_TIMES)      # K,perturbation potential temperature theta-t0
        pdata = wrf.getvar(ncfile, 'P',timeidx=wrf.ALL_TIMES)      # Pa,perturbation pressure
        pbdata = wrf.getvar(ncfile, 'PB',timeidx=wrf.ALL_TIMES)    # Pa,base state pressure
        phdata = wrf.getvar(ncfile, 'PH',timeidx=wrf.ALL_TIMES)
        phbdata = wrf.getvar(ncfile, 'PHB',timeidx=wrf.ALL_TIMES)
        

        # presure: base state + PB (KSMP)
        pb2data = np.zeros([np.shape(pdata)[0], np.shape(pdata)[1], np.shape(pdata)[2], np.shape(pdata)[3]],dtype=np.float)
        pb2data[:,:,:,:] = pdata[:,:,:,:]+ pbdata[:,:,:,:]

        # convert the perturbation potential temperature (from 300K reference) to temp
        tbdata = np.zeros([np.shape(tdata)[0], np.shape(tdata)[1], np.shape(tdata)[2], np.shape(tdata)[3]],dtype=np.float)
        tbdata[:,:,:,:] =(300.0+tdata[:,:,:,:])*((pb2data[:,:,:,:]/1.0e5)**0.286)
        
        wrfdata = {'no2': no2data, 'pb2':pb2data, 'ph':phdata, 'phb': phbdata, 't2': tbdata}

        ncfile.close()        
        return wrfdata
        
#--- 
class tropomi_process:
    def __init__(self, year, month, day):
        self.omidata_alltrack = []
        self.scatwts_alltrack = []
        self.year = year
        self.month = month
        self.day = day
        
    def avgtropomi(self):
        year = self.year
        month = self.month
        day = self.day

        fm = file_management()
        timeindex = '{:04d}'.format(year)+'{:02d}'.format(month) + '{:02d}'.format(day)
        subdirlist = fm.subdirlist(Basedir_tropomi)

        omidata_alltrack = []
        
        for subdir in subdirlist:
            fflist = fm.filelist(subdir, keyword='____'+timeindex)
            tracknumber = len(fflist)
            for infile in fflist:
                # identify if the swath cover the WRF-Chem domain
                llregion = self.identifyregion(infile)
                if llregion == True:
                    omidata = self.readtropomi(infile)
                    omidata_alltrack.append(omidata)
                else:
                    pass
                
        self.omidata_alltrack = omidata_alltrack
        
    def identifyregion(self, infile):
        # identify if the swath cover the WRF-Chem domain
        # wrf-chem domain
        wrflonlat = wrfcoord
        wrflon = wrflonlat['lon']
        wrflat = wrflonlat['lat']

        # identify if the file is located in the region or not
        ds = Dataset(infile,"r")
        variable = np.squeeze(wrf.to_np(ds.groups['PRODUCT']['latitude'][:,:,:]))
        lat_ind = variable
        variable = np.squeeze(wrf.to_np(ds.groups['PRODUCT']['longitude'][:,:,:]))
        lon_ind = variable
        ds.close()

        # determine the location in wrfchem
        locind = extractwrfcoord(lats=lat_ind, lons=lon_ind)

        wrf_ind_x = int(np.shape(wrflon)[0]) # lat
        wrf_ind_y = int(np.shape(wrflon)[1]) # lon

        llregion = False
        xtomi = locind[0,:]
        ytomi = locind[1,:]

        ind = np.where((xtomi >= 0.0) & (ytomi >= 0.0) & (xtomi < wrf_ind_y) & (ytomi < wrf_ind_x))
        if np.shape(ind)[1] == 0:
            pass
        else:
           llregion = True
        return llregion
                    
    def readtropomi(self,infile):
        # read tropomi swath L2 NO2 data, return OMI NO2 data directory
        omidata = {}
        print('--> reading tropomi datafile of : ', infile)
        ds = Dataset(infile,"r")
        
        variable = ds.groups['PRODUCT']['nitrogendioxide_tropospheric_column']
        omidata['nitrogendioxide_tropospheric_column'] = variable
        
        variable = ds.groups['PRODUCT']['qa_value']
        omidata['qa_value'] = variable        
        
        variable = ds.groups['PRODUCT']['averaging_kernel']
        omidata['averaging_kernel'] = variable
        
        variable = ds.groups['PRODUCT']['air_mass_factor_total']
        omidata['air_mass_factor_total'] = variable
        
        variable = ds.groups['PRODUCT']['air_mass_factor_troposphere']
        omidata['air_mass_factor_troposphere'] = variable
        
        variable = ds.groups['PRODUCT']['latitude']
        omidata['latitude'] = variable
        
        variable = ds.groups['PRODUCT']['longitude']
        omidata['longitude'] = variable
        
        variable = ds.groups['PRODUCT']['tm5_constant_a']
        omidata['tm5_constant_a'] = variable
        
        variable = ds.groups['PRODUCT']['tm5_constant_b']
        omidata['tm5_constant_b'] = variable
        
        variable = ds.groups['PRODUCT']['tm5_tropopause_layer_index']
        omidata['tm5_tropopause_layer_index'] = variable
        
        variable = ds.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']['surface_pressure']
        omidata['surface_pressure'] = variable
        
        variable = ds.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']['cloud_fraction_crb']
        omidata['cloud_fraction_crb'] = variable

        variable = ds.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['nitrogendioxide_stratospheric_column']
        omidata['nitrogendioxide_stratospheric_column'] = variable

        variable = ds.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['air_mass_factor_stratosphere']
        omidata['air_mass_factor_stratosphere'] = variable
        
        variable = ds.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['nitrogendioxide_slant_column_density']
        omidata['nitrogendioxide_slant_column_density'] = variable
       
        variable = ds.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['latitude_bounds']
        omidata['latitude_bounds'] = variable #time x scanline x groud_pixel x corners, 1x3245x450x4

        variable = ds.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['longitude_bounds']
        omidata['longitude_bounds'] = variable #time x scanline x groud_pixel x corners, 1x3245x450x4

        # calculate the preslev
        pleva = omidata['tm5_constant_a']
        plevb = omidata['tm5_constant_b']
        
        spre = np.squeeze(omidata['surface_pressure'], axis=0)
        aks = omidata['averaging_kernel']
                
        FillValue = omidata['surface_pressure']._FillValue
        preslev = np.copy(aks)
        preslev[:,:,:,:]  = np.nan
        aks = None # to save memory
        del aks
        
        spre[np.where(spre == FillValue)] = np.nan

        for l in range(np.shape(preslev)[3]):
            preslev[0,:,:,l] = (pleva[l,0]+spre[:,:]*plevb[l,0] + pleva[l,1]+spre[:,:]*plevb[l,1]) / 2.0 # center of the vertical layer            
        omidata['preslev'] = preslev
        
        return omidata
        ds.close()
               
#---


if __name__ == "__main__":
    # year, month, day
    print(sys.argv[1], sys.argv[2], sys.argv[3])
    main(int(sys.argv[1]), int(sys.argv[2]),int(sys.argv[3]))

