
# This code is written to process monthly averaged map for both wrfchem and TROPOMI
# --- Meng Li, 2019. 5. 9
# --- Contact: meng.li@noaa.gov; meng.li.atm@gmail.com

import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import wrf
from os import listdir
from os.path import isfile


Basedir_wrfoutput = os.environ.get('Basedir_wrfoutput')

'''
===================================================
Main program: wrfoutput, tropomi, and evaluation
===================================================
'''

#=================Preparation Codes==================
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

#===========MAIN PROGRAM STARTS HERE===============
    
# GET THE WRF COORDIATE INFORMATION
wrfcoord = extractwrfcoord()
#---
class output_config:   
    def __init__(self):
        SMALL_SIZE=12
        MEDIUM_SIZE = 16
        BIG_SIZE = 18
        plt.rc('font', size=SMALL_SIZE )
        plt.rc('axes', titlesize=SMALL_SIZE)
        plt.rc('axes', labelsize=MEDIUM_SIZE)
        plt.rc('xtick', labelsize=SMALL_SIZE)
        plt.rc('ytick', labelsize=SMALL_SIZE)
        plt.rc('legend', fontsize=MEDIUM_SIZE)
        plt.rc('figure', titlesize=BIG_SIZE)
        plt.rc('font',**{'family':'sans-serif','sans-serif':['arial']})

        
    def outputnc_2d(self, fn, value, valuename,valueunit):      
        print('--> output 2d data for:', valuename)
        ds = Dataset(fn, 'w', format = 'NETCDF4_CLASSIC')
        ds.createDimension('longitude', np.shape(value)[1])
        ds.createDimension('latitude', np.shape(value)[0])
        dlong = ds.createVariable('longitude', 'f4', ['latitude','longitude'])
        dlat = ds.createVariable('latitude', 'f4', ['latitude','longitude'])
        dsec = ds.createVariable(valuename, 'f4', ['latitude','longitude'])
    					
        #lat, lon = wrf.latlon_coords(value)
        lon = wrfcoord['lon']
        lat = wrfcoord['lat']
        ds.variables['longitude'][:,:] = lon[:,:]
        ds.variables['latitude'][:,:] = lat[:,:]
        ds.longitude = 'Edge of grids, West to East'
        ds.latitude = 'Edge of grids, South to North'
        ds.variables[valuename][:,:] = value[:,:]
        ds.valuename = valueunit
        ds.close()       
    

    def outputnc_3d(self, fn, lon,lat,time, valuename, value, valueunit):       
    	ds = Dataset(fn, 'w', format = 'NETCDF4_CLASSIC')   
    	ds.createDimension('longitude', np.shape(lon)[0])
    	ds.createDimension('latitude', np.shape(lat)[0])
    	ds.createDimension('time', np.shape(time)[0])
    
    	dlong = ds.createVariable('longitude', 'f4', ['longitude'])
    	dlat = ds.createVariable('latitude', 'f4', ['latitude'])
    	dmonth = ds.createVariable('time', 'f4', ['time'])	
    	dsec = ds.createVariable(valuename, 'f4', ['time','latitude','longitude'])
    					
    	ds.variables['longitude'][:] = lon[:]
    	ds.variables['latitude'][:] = lat[:]
    	ds.variables['time'][:] = time[:]
    	ds.longitude = 'Edge of grids, West to East'
    	ds.latitude = 'Edge of grids, South to North'
    	ds.month = 'Time'
    
    	ds.variables[valuename][:,:,:] = value[:,:,:]		
    	ds.valuename = valueunit
    	ds.close()
        
    def plot_2dmap(self, fn, value,valuename, valueunit,  mindata=0.0, maxdata = 0.0):    
        print('--> plotting 2d map for:', valuename       )
        from wrf import (to_np, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords)
        import cartopy.crs as crs        
        from cartopy.feature import NaturalEarthFeature
        from matplotlib.cm import get_cmap
        from cartopy.io.shapereader import Reader
        from cartopy.feature import ShapelyFeature
        
	# get the cartopy mapping object
        cart_proj = get_cartopy(wrfcoord['lon'])
        lats, lons = latlon_coords(wrfcoord['lon'])    
        # create a figure
        fig = plt.figure(figsize=(12,6))    
        # set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)
        # download and add the states and coastlines
        # states = NaturalEarthFeature(category="cultural",scale="50m",
        #                             facecolor="none", name="admin_1_states_provinces_shp")
        states_reader = Reader('/scratch1/BMC/rcm2/mli/tech/ne_50m_admin_1_states_provinces/ne_50m_admin_1_states_provinces.shp')
        states = ShapelyFeature(states_reader.geometries(), crs.PlateCarree())
        ax.add_feature(states,linewidth=0.5, edgecolor="black", facecolor='none')
        #ax.coastlines('50m',linewidth=0.8)
        coast_reader = Reader('/scratch1/BMC/rcm2/mli/tech/ne_50m_coastline/ne_50m_coastline.shp')
        coast = ShapelyFeature(coast_reader.geometries(), crs.PlateCarree())
        ax.add_feature(coast, linewidth=0.8, edgecolor='black', facecolor='none')

        # make the contour outlines and filled contoures for the value
        #plt.contour(lons, lats, value, 10, colors="black", transform=crs.PlateCarree())
        #plt.contourf(lons, lats, value, vmin=mindata, vmax=maxdata, transform=crs.PlateCarree(), cmap='jet')
        if (mindata == 0.0 and maxdata == 0.0):
            mindata = np.min(value)
            maxdata = np.max(value)
            
        plt.pcolormesh(lons, lats, value, vmin = mindata, vmax = maxdata, cmap='jet',transform=crs.PlateCarree() )
        #plt.imshow(lons, lats, value)
        cb = plt.colorbar(ax=ax, shrink=.98)
        cb.ax.tick_params(labelsize=18, length=8)
        # set the map bounds
        ax.set_xlim(cartopy_xlim(wrfcoord['lon']))
        ax.set_ylim(cartopy_ylim(wrfcoord['lat']))
        
        # ad the grid lines
        ax.gridlines(color="black", linestyle = "dotted")
        plt.title(valuename + ', unit: ' +valueunit, fontsize=22, fontweight='bold')
        plt.savefig(fn, dpi=300)
        #plt.show()
        plt.clf()
        plt.close()

    def plot_2dmap_ccolbar(self, fn, value,valuename, valueunit,  mindata=0.0, maxdata = 0.0):
        print('--> plotting 2d map for:', valuename)
        from wrf import (to_np, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords)
        import cartopy.crs as crs
        from cartopy.feature import NaturalEarthFeature
        from matplotlib.cm import get_cmap
        import matplotlib as mpl
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
        from cartopy.io.shapereader import Reader
        from cartopy.feature import ShapelyFeature

        # get the cartopy mapping object
        cart_proj = get_cartopy(wrfcoord['lon'])
        lats, lons = latlon_coords(wrfcoord['lon'])
        # create a figure
        fig = plt.figure(figsize=(12,6))
        # set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)
        # download and add the states and coastlines
        states_reader = Reader('/scratch1/BMC/rcm2/mli/tech/ne_50m_admin_1_states_provinces/ne_50m_admin_1_states_provinces.shp')
        states = ShapelyFeature(states_reader.geometries(), crs.PlateCarree())
        ax.add_feature(states,linewidth=0.5, edgecolor="black", facecolor='none')
        #ax.coastlines('50m',linewidth=0.8)
        coast_reader = Reader('/scratch1/BMC/rcm2/mli/tech/ne_50m_coastline/ne_50m_coastline.shp')
        coast = ShapelyFeature(coast_reader.geometries(), crs.PlateCarree())
        ax.add_feature(coast, linewidth=0.8, edgecolor='black', facecolor='none')

        # make the contour outlines and filled contoures for the value
        #plt.contour(lons, lats, value, 10, colors="black", transform=crs.PlateCarree())
        #plt.contourf(lons, lats, value, vmin=mindata, vmax=maxdata, transform=crs.PlateCarree(), cmap='jet')
        if (mindata == 0.0 and maxdata == 0.0):
            mindata = np.min(value)
            maxdata = np.max(value)
        # define a custom colorbar
        #cmap = plt.cm.RdBu_r
        #cmap = plt.get_cmap('bwr')
        cmap = plt.get_cmap('PuBuGn')
        # extract all colors from the .jet map
        cmaplist = [cmap(i) for i in range(cmap.N)]
        # CREATE THE NEW MAP
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
        # define the bins and normalize
        #bounds = np.linspace(0,20, 21)
        bounds = np.linspace(mindata,maxdata,11)
        norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=cmap.N)

        plt.pcolormesh(lons, lats, value, vmin = mindata, vmax = maxdata, norm=norm,cmap=cmap,transform=crs.PlateCarree() )
        #plt.imshow(lons, lats, value)
        cb = plt.colorbar(ax=ax, shrink=.98,extend='both', ticks=bounds)
        cb.ax.tick_params(labelsize=18, length=8)
        # set the map bounds
        #ax.set_xlim(cartopy_xlim(wrfcoord['lon']))
        #ax.set_ylim(cartopy_ylim(wrfcoord['lat']))
        ax.set_xlim(cartopy_xlim(wrfcoord['lon']))
        ax.set_ylim(cartopy_ylim(wrfcoord['lat']))

        # ad the grid lines
        ax.gridlines(color="black", linestyle = "dotted")
        plt.title(valuename + ', unit: ' +valueunit, fontsize=22, fontweight='bold' )
        plt.savefig(fn, dpi=300)
        #plt.show()
        plt.clf()
        plt.close()

    def plot_scatter(self,fn, title, x, y, mindata=0.0, maxdata=1e17, mb=None, nmb=None, label='$\mathregular{NO_2}$ column',
            xlabel='TROPOMI $\mathregular{NO_2}$ column, $\mathregular{10^{15}}$ molec c$\mathregular{m^{-2}}$',
            ylabel='WRF-Chem $\mathregular{NO_2}$ column, $\mathregular{10^{15}}$ molec c$\mathregular{m^{-2}}$'):
        from scipy import stats
        #fig = plt.figure(figsize=(6,6))
        fig, ax = plt.subplots(figsize=(6,6))
        x = x/1e15
        y = y/1e15
        mindata = mindata/1e15
        maxdata = maxdata/1e15
        slope, intercept, r_value, p_value, stderr=stats.linregress(x, y)
        plt.scatter(x,y, marker = 'o',facecolors='cornflowerblue', edgecolors='b',label=label, s=150) # s=50
        plt.xlim(mindata, maxdata)
        plt.ylim(mindata, maxdata)
        plt.xlabel(xlabel,fontsize=18, weight=500)
        plt.ylabel(ylabel,fontsize=18, weight=500)

        xarr = np.arange(start=0.0, stop=1e17,step=1e15)
        plt.plot(xarr, xarr,'k--',label='_nolegend')

        y2 = xarr*slope + intercept
        plt.plot(xarr, y2, 'r-',label='liner regression')
        div = (maxdata - mindata)/18.0

        if intercept > 0.0:
            eqstr = 'Y='+'{:5.2f}'.format(slope)+'X+'+'{:5.2f}'.format(intercept)+'e15'#'$\mathregular{10^{15}}$'
        else:
            eqstr = 'Y='+'{:5.2f}'.format(slope)+'X-'+'{:5.2f}'.format(intercept*(-1.0))+'e15'
        plt.text(maxdata/2.0, mindata+div*5.0,eqstr,color='r',fontname = 'arial',fontsize=18) # fontweight='bold', fontstyle='italic'
        r_value = r_value * r_value
        rstr = '$\mathregular{R^{2}}$:'+'{:5.2f}'.format(r_value)
        plt.text(maxdata/2.0, mindata+div*3.5,rstr,color='r',fontsize=18)

        plt.text(maxdata*0.2/6.0, maxdata-div*5.0, 'N: ' + str(len(x)),fontname='arial',fontsize=18)
        #plt.text(maxdata/2.0, maxdata-div*2.0, 'slope: '+ '{:5.2f}'.format(slope), **csfont)
        #plt.text(maxdata/2.0, maxdata-div*3.0, 'intercept: '+ '{:5.2f}'.format(intercept/1e15) + 'e15', **csfont)
        if mb != None:
            plt.text(maxdata*0.2/6.0, maxdata-div*6.5, 'MB: ' + '{:5.2f}'.format(mb/1e15)+'e15',fontname='arial', fontsize=18)
        if nmb != None:
            plt.text(maxdata*0.2/6.0, maxdata-div*8.0, 'NMB: ' + '{:5.2f}'.format(nmb*100.0)+'%',fontname='arial', fontsize=18)
 
        plt.legend(loc='upper left',fontsize=16)
        #plt.title(title,fontname='arial', fontsize=16)
        ax.tick_params(length=8, width=1, labelsize=18)
        plt.savefig(fn, bbox_inches = 'tight', dpi=300)
        #plt.show()
        plt.clf()
        plt.close()
        print('--> scatter: r,slope,intercept:', r_value, slope, intercept/1e15)
        
    def plot_minus(self, fn, value,valuename, valueunit,  mindata=-1.0e16, maxdata = 1.0e16):
        print('--> plotting 2d minus map for:', valuename)
        
        from matplotlib import cm
        from wrf import (to_np, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords)
        import cartopy.crs as crs
        from cartopy.feature import NaturalEarthFeature
        from matplotlib.cm import get_cmap
        from cartopy.io.shapereader import Reader
        from cartopy.feature import ShapelyFeature       
 
        # get the cartopy mapping object
        cart_proj = get_cartopy(wrfcoord['lon'])
        lats, lons = latlon_coords(wrfcoord['lon'])    
        # create a figure
        fig = plt.figure(figsize=(12,6))    
        # set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)
        # download and add the states and coastlines
        # states = NaturalEarthFeature(category="cultural",scale="50m",
        #                             facecolor="none", name="admin_1_states_provinces_shp")
        states_reader = Reader('/scratch1/BMC/rcm2/mli/tech/ne_50m_admin_1_states_provinces/ne_50m_admin_1_states_provinces.shp')
        states = ShapelyFeature(states_reader.geometries(), crs.PlateCarree())
        ax.add_feature(states,linewidth=0.5, edgecolor="black",facecolor='none')
        #ax.coastlines('50m',linewidth=0.8)
        coast_reader = Reader('/scratch1/BMC/rcm2/mli/tech/ne_50m_coastline/ne_50m_coastline.shp')
        coast = ShapelyFeature(coast_reader.geometries(), crs.PlateCarree())
        ax.add_feature(coast, linewidth=0.8, edgecolor='black', facecolor='none')

        # make the contour outlines and filled contoures for the value
        #plt.contour(lons, lats, value, 10, colors="black", transform=crs.PlateCarree())
        #plt.contourf(lons, lats, value, vmin=mindata, vmax=maxdata, transform=crs.PlateCarree(), cmap='jet')
        if (mindata == 0.0 and maxdata == 0.0):
            mindata = np.min(value)
            maxdata = np.max(value)
        
        cmap = plt.get_cmap('bwr')
        plt.pcolormesh(lons, lats, value, vmin = mindata, vmax = maxdata, cmap=cmap,transform=crs.PlateCarree() )       
        # color bar
        cb = plt.colorbar(ax=ax, shrink=.98)
        cb.ax.tick_params(labelsize=18, length=8)
        # set the map bounds
        ax.set_xlim(cartopy_xlim(wrfcoord['lon']))
        ax.set_ylim(cartopy_ylim(wrfcoord['lat']))
        
        # ad the grid lines
        ax.gridlines(color="black", linestyle = "dotted")
        plt.title(valuename + ', unit: ' +valueunit, fontsize=22, fontweight='bold')
        plt.savefig(fn, dpi=300)
        #plt.show()
        plt.clf()
        plt.close()
        

#---
if __name__ == '__main__':
        main()
    
 



