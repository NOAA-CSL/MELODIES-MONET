# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
'''
Plot_2D.py
this code is designed for plotting CESM output 
can be used for either finite volume or spectral element (+ regional refinement)
'''

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as cm
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.collections import PolyCollection
import matplotlib
from matplotlib import ticker

class Plot_2D(object):
    
    def __init__(self, var, lons=None, lats=None, lon_range=[-180,180], lat_range=[-90,90],
                 scrip_file="", ax=None, cmap=None, projection=ccrs.PlateCarree(), center_180=False, 
                 grid_line=False, grid_line_lw=1, coast=True, country=True, state=False, 
                 resolution="10m", feature_line_lw=0.5, feature_color="black",
                 lonlat_info=True, lonlat_line=True, lon_interval=None, lat_interval=None,
                 lon_labels=None, lat_labels=None,
                 font_family="STIXGeneral", label_size=15, colorbar=True, 
                 log_scale=False, log_scale_min=None, diff=False, orientation="horizontal", 
                 shrink=0.8, pad=0.12, fraction=0.1, extend='both',
                 colorticks=None, colorlabels=None, pretty_tick=True, nticks=None, 
                 cmax=None, cmin=None, title="", title_size=20, title_bold=False,
                 unit="", unit_size=15, unit_bold=False, unit_italic=True, unit_offset=[0.0,0.0],
                 verbose=False):
        
        # ========================================================================
        # ===== Error check and pass input values to class-accessible values =====
        # ========================================================================
        # variable dimension check
        if (np.ndim(var) > 2) or (np.ndim(var) < 1):
            raise ValueError( '"var" must be 1-D (SE) or 2-D (FV) array' )
        
        # xarray, lon, and lat check
        if type(var) in [ xr.core.dataset.Dataset, xr.core.dataarray.DataArray ]:
            if np.ndim(var) == 2: # FV results
                self.model_type = 'FV'
                if verbose:
                    print( '"var" is a xarray variable, longitude and latitude values ' + \
                           'are being automatically assigned by xarray dimension variables')
                if lons != None:
                    print( 'Warning: "lons" is assigned but not used ' + \
                           'because xarray itself has longitude values')
                if lats != None:
                    print( 'Warning: "lats" is assigned but not used ' + \
                           'because xarray itself has longitude values')                
                self.var = np.copy( var.values )
                self.lon = np.copy( var.lon.values )
                self.lat = np.copy( var.lat.values )
            else: # SE results
                self.model_type = 'SE'
                self.var = np.copy( var.values )
        else:
            self.var = np.copy(var)
            if np.ndim(var) == 2: # FV results
                self.model_type = 'FV'
                if np.shape(lons) == ():
                    raise ValueError( '"lons" must be provided for FV model output' )
                else:
                    self.lon = np.copy( lons )
                if np.shape(lats) == ():
                    raise ValueError( '"lats" must be provided for FV model output' )
                else:
                    self.lat = np.copy( lats )
            else: # SE results
                self.model_type = 'SE'        
        
        # lon_range dimension check
        if len(lon_range) != 2:
            raise ValueError( 'Check lon_range!' + '\n' + \
                              'Current Values:', lon_range)
        else:
            if (lon_range==[-180,180]) & (lat_range==[-90,90]):
                self.regional = False
            else:
                self.regional = True
            
            self.lon_range = lon_range
        
        # lat_range dimension check
        if len(lat_range) != 2:
            raise ValueError( 'Check lat_range!' + '\n' + \
                              'Current Values:', lat_range)
        else:
            self.lat_range = lat_range

        # To shift map by 180 degree in x-axis
        if center_180:
            if projection==ccrs.PlateCarree():
                projection = ccrs.PlateCarree(central_longitude=180)
            else:
                if projection.proj4_params['lon_0'] != 180:
                    raise ValueError( 'central_longitude must be set to 180 if center_180 is True' )
            
            if self.lon_range[1] < self.lon_range[0]:
                self.lon_range[0] -= 180
                self.lon_range[1] += 180

                    
        # Read scrip file in case of SE model output
        if self.model_type == 'SE':
            if type(scrip_file) == xr.core.dataset.Dataset:
                ds_scrip = scrip_file
                if verbose:
                    print( "use xarray dataset for scrip file" )
            else:
                if scrip_file == "":
                    raise ValueError( '"scrip_file" must be specified for SE model output' )
                if type(scrip_file) != str:
                    raise ValueError( '"scrip_file" must be provided as "string"' )
                if verbose:
                    print( "Read SCRIP file:", scrip_file )
                ds_scrip = xr.open_dataset( scrip_file )
                
            self.corner_lon = np.copy( ds_scrip.grid_corner_lon.values )
            self.corner_lat = np.copy( ds_scrip.grid_corner_lat.values )
            self.center_lon = np.copy( ds_scrip.grid_center_lon.values )
            self.center_lat = np.copy( ds_scrip.grid_center_lat.values )

            
        # Color map check
        if cmap == None:
            if diff:
                #self.cmap = cm.seismic
                self.cmap = cm.bwr
            else:
                self.cmap = cm.jet
                self.cmap = Cmap.cmap
        else:
            self.cmap = cmap
        
        # axis check
        if ax == None:
            self.fig = plt.figure( figsize=(8,5) )
            ax = self.fig.add_subplot(1,1,1, projection=projection)
        else:
            self.fig = ax.figure
        self.ax = ax
        
        # nticks check, assign the base value if None
        if nticks == None:
            if log_scale:
                self.nticks = 7
            else:
                self.nticks = 5
        else:
            self.nticks = nticks
            
        
        # Pass input keywords
        self.scrip_file = scrip_file
        self.font_family = font_family
        self.projection = projection
        self.center_180 = center_180
        self.verbose = verbose
        self.grid_line = grid_line
        self.grid_line_lw = grid_line_lw
        self.label_size = label_size
        self.coast = coast
        self.country = country
        self.state = state
        self.resolution = resolution
        self.feature_line_lw = feature_line_lw
        self.feature_color = feature_color
        self.lonlat_info = lonlat_info
        self.lonlat_line = lonlat_line
        self.lon_interval = lon_interval
        self.lat_interval = lat_interval
        self.lon_labels = lon_labels
        self.lat_labels = lat_labels
        self.colorbar = colorbar
        self.log_scale = log_scale
        self.log_scale_min = log_scale_min
        self.orientation = orientation
        self.shrink = shrink
        self.pad = pad
        self.fraction = fraction
        self.extend = extend
        self.colorticks = colorticks
        self.colorlabels = colorlabels
        self.pretty_tick = pretty_tick
        self.cmax = cmax
        self.cmin = cmin
        self.title = title
        self.title_size = title_size
        self.title_bold = title_bold
        self.unit = unit
        self.unit_size = unit_size
        self.unit_bold = unit_bold
        self.unit_italic = unit_italic
        self.unit_offset = unit_offset
        # === END Error check and pass input values to class-accessible values ===
        # ========================================================================
        
        
        # =======================================================================
        # ============================ Initial Setup ============================
        # =======================================================================        
        # If lon_range doesn't match longitude values, 
        # shift longitude values by 180 degree
        if self.model_type == 'FV': # 2D FV model output
            if ( (np.min(self.lon_range) < 0) & (np.max(self.lon) > 180) ):
                if not self.center_180:
                    self.lon[self.lon > 180.] -= 360.
                    if verbose:
                        print( "FV model: Shift longitude values by 180 degree" )
        else: # 1D SE model output
            if ( (np.min(self.lon_range) < 0) & (np.max(self.corner_lon) > 180) ):
                if not self.center_180:
                    self.corner_lon[self.corner_lon > 180.] -= 360.
                    if verbose:
                        print( "SE model: Shift longitude values by 180 degree" )

        # automatically set longitude and latitude intervals        
        if self.lonlat_info:
            if self.lon_interval == None:
                lon_length = lon_range[1] - lon_range[0]
                self.lon_interval = np.around( lon_length / 6. )
                if self.lon_interval < 1:
                    self.lon_interval = 1
            if self.lat_interval == None:
                lat_length = lat_range[1] - lat_range[0]
                self.lat_interval = np.around( lat_length / 6. )
                if self.lat_interval < 1:
                    self.lat_interval = 1

        # set vertices for SE model output
        if self.model_type == 'SE':
            
            self.lons_corners = np.copy( self.corner_lon.reshape( self.corner_lon.shape[0],
                                                                  self.corner_lon.shape[1],1) )
            self.lats_corners = np.copy( self.corner_lat.reshape( self.corner_lat.shape[0],
                                                                  self.corner_lat.shape[1],1) )

            if not self.center_180:
                self.lons_corners[ self.lons_corners > 180. ] -= 360
                self.center_lon[ self.center_lon > 180. ] -= 360
            
            self.lons_corners_add = []
            self.lats_corners_add = []
            self.var_add = []
            # For longitudes -180, 180
            for i, cenlon in enumerate( self.center_lon ):
                lon_maxmin = np.max( self.lons_corners[i,:,:] ) - \
                             np.min( self.lons_corners[i,:,:] )
                if self.center_180:
                    if ( lon_maxmin > 180 ):
                        if np.mean( self.lons_corners[i,:,:] ) <= 180:  
                            inds2 = np.where( self.lons_corners[i,:,:] < 180)[0]
                            tmp_lons_corners = np.copy( self.lons_corners[i,:] )
                            tmp_lons_corners[inds2] = 360.
                            self.lons_corners_add.append( tmp_lons_corners )
                            self.lats_corners_add.append( self.lats_corners[i,:] )

                            inds = np.where( self.lons_corners[i,:,:] > 180 )[0]
                            self.lons_corners[i,inds] = 0.

                            self.var_add.append( self.var[i] )

                        elif np.mean( self.lons_corners[i,:,:] ) > 180:
                            inds2 = np.where( self.lons_corners[i,:,:] > 180)[0]
                            tmp_lons_corners = np.copy( self.lons_corners[i,:] )
                            tmp_lons_corners[inds2] = 0
                            self.lons_corners_add.append( tmp_lons_corners )
                            self.lats_corners_add.append( self.lats_corners[i,:] )                     

                            inds = np.where( self.lons_corners[i,:,:] > 180 )[0]
                            self.lons_corners[i,inds] = 360.

                            self.var_add.append( self.var[i] )
                            
                else:
                    if ( lon_maxmin > 180 ):
                        if np.mean( self.lons_corners[i,:,:] ) <= 0:  
                            inds2 = np.where( self.lons_corners[i,:,:] < 0)[0]
                            tmp_lons_corners = np.copy( self.lons_corners[i,:] )
                            tmp_lons_corners[inds2] = 180.
                            self.lons_corners_add.append( tmp_lons_corners )
                            self.lats_corners_add.append( self.lats_corners[i,:] )

                            inds = np.where( self.lons_corners[i,:,:] > 0 )[0]
                            self.lons_corners[i,inds] = -180.

                            self.var_add.append( self.var[i] )

                        elif np.mean( self.lons_corners[i,:,:] ) > 0:
                            inds2 = np.where( self.lons_corners[i,:,:] > 0)[0]
                            tmp_lons_corners = np.copy( self.lons_corners[i,:] )
                            tmp_lons_corners[inds2] = -180.
                            self.lons_corners_add.append( tmp_lons_corners )
                            self.lats_corners_add.append( self.lats_corners[i,:] )                     

                            inds = np.where( self.lons_corners[i,:,:] < 0 )[0]
                            self.lons_corners[i,inds] = 180.

                            self.var_add.append( self.var[i] )
                            
                            
            self.lons_corners_add = np.array(self.lons_corners_add)        
            self.lats_corners_add = np.array(self.lats_corners_add)
                            
            if self.center_180:
                self.lons_corners[ self.lons_corners > 180. ] -= 360
                self.lons_corners_add[ self.lons_corners_add > 180. ] -= 360
                self.center_lon[ self.center_lon > 180. ] -= 360
                
                for lc1 in np.arange(len(self.lons_corners[:,0,0])):
                    for lc2 in np.arange(len(self.lons_corners[0,:,0])):
                        for lc3 in np.arange(len(self.lons_corners[0,0,:])):
                            if self.lons_corners[lc1,lc2,lc3] >= 0:
                                self.lons_corners[lc1,lc2,lc3] -= 180
                            elif self.lons_corners[lc1,lc2,lc3] < 0:
                                self.lons_corners[lc1,lc2,lc3] += 180
                                
                for lc1 in np.arange(len(self.lons_corners_add[:,0,0])):
                    for lc2 in np.arange(len(self.lons_corners_add[0,:,0])):
                        for lc3 in np.arange(len(self.lons_corners_add[0,0,:])):
                            if self.lons_corners_add[lc1,lc2,lc3] >= 0:
                                self.lons_corners_add[lc1,lc2,lc3] -= 180
                            elif self.lons_corners_add[lc1,lc2,lc3] < 0:
                                self.lons_corners_add[lc1,lc2,lc3] += 180
                
                for lc1 in np.arange(len(self.center_lon)):
                    if self.center_lon[lc1] >= 0:
                        self.center_lon[lc1] -= 180
                    elif self.center_lon[lc1] < 0:
                        self.center_lon[lc1] += 180
                
                

            if self.lons_corners_add != []:
                self.lons_corners = np.concatenate( (self.lons_corners, 
                                                     self.lons_corners_add), axis=0 )
                self.lats_corners = np.concatenate( (self.lats_corners, 
                                                     self.lats_corners_add), axis=0 )
                
            self.var = np.concatenate( (self.var, np.array(self.var_add)), axis=0 )
            self.verts = np.concatenate( ( self.lons_corners, self.lats_corners ), axis=2)
            
                
            
            
        # set plot color properties (FV model output)
        if self.model_type == 'FV':
            self.lon_inds = np.where( ( self.lon >= self.lon_range[0] ) & \
                                      ( self.lon <= self.lon_range[-1] ) )[0]
            self.lat_inds = np.where( ( self.lat >= self.lat_range[0] ) & \
                                      ( self.lat <= self.lat_range[-1] ) )[0]
            self.var_slice = self.var[ self.lat_inds[0]:self.lat_inds[-1]+1,
                                       self.lon_inds[0]:self.lon_inds[-1]+1 ]
        elif self.model_type == 'SE':
            self.ncol_inds = np.where( ( self.center_lon >= self.lon_range[0] ) & \
                                       ( self.center_lon <= self.lon_range[-1] ) & \
                                       ( self.center_lat >= self.lat_range[0] ) & \
                                       ( self.center_lat <= self.lat_range[-1] ) )[0]
            self.corner_lon_slice = self.corner_lon[ self.ncol_inds, : ]
            self.corner_lat_slice = self.corner_lat[ self.ncol_inds, : ]
            self.center_lon_slice = self.center_lon[ self.ncol_inds ]
            self.center_lat_slice = self.center_lat[ self.ncol_inds ]
        
            self.var_slice = self.var[ self.ncol_inds ]

            
        # set colorbar properties
        self.kwd_pretty_tick = {}
        if not self.log_scale:
            if self.cmax == None:
                self.cmax = np.max( self.var_slice )
            else:
                self.kwd_pretty_tick['max_set'] = self.cmax

            if self.cmin == None:
                self.cmin = np.min( self.var_slice )
            else:
                self.kwd_pretty_tick['min_set'] = self.cmin
        
        #if self.colorbar:
        # Automatically set tick values
        if self.pretty_tick:
            if np.shape(self.colorticks) == ():
                if self.log_scale:
                    if self.cmin == None:
                        self.cmin_od = \
                            np.floor(np.log10(np.abs(np.min( \
                                self.var_slice[self.var_slice != 0]))))
                        self.cmin_sign = np.sign(np.min(self.var_slice[self.var_slice != 0])) 
                    else:
                        self.cmin_od = np.floor(np.log10(np.abs(self.cmin))) 
                        self.cmin_sign = np.sign(self.cmin)
                    if self.cmax == None:
                        self.cmax_od = \
                            np.floor(np.log10(np.abs(np.max( \
                                self.var_slice[self.var_slice != 0]))))
                        self.cmax_sign = np.sign( np.max( \
                                self.var_slice[self.var_slice != 0]) )
                    else:
                        self.cmax_od = np.floor(np.log10(np.abs(self.cmax)))
                        self.cmax_sign = np.sign(self.cmax)

                    self.sign_sum = self.cmin_sign + self.cmax_sign
                    if self.sign_sum in [1,2]:
                        if self.cmax_od > 3:
                            nticks_init = self.cmax_od + 1
                        else:
                            nticks_init = self.cmax_od - self.cmin_od + 1

                        if nticks_init > 6:
                            Nticks_list = [5,6,7,4]
                            check_loop = True
                            for nticks in Nticks_list:
                                if self.cmax_od > 3:
                                    tmparray = np.linspace(0, self.cmax_od, 
                                                           np.int(nticks))
                                else:
                                    tmparray = np.linspace(self.cmin_od, self.cmax_od, 
                                                           np.int(nticks))
                                checksum = np.sum( np.abs(tmparray - tmparray.astype('i')) )
                                if np.abs(checksum) < 1e-7:
                                    check_loop = False
                                    break

                            if check_loop:
                                interval = np.ceil( self.cmax_od / 3 ).astype('I')
                                tick_start = self.cmax_od
                                colorticks_h = []
                                nticks = 0
                                while (tick_start > 0):
                                    colorticks_h.append( tick_start )
                                    tick_start -= interval
                                self.colorticks = np.array( [0] + \
                                                           list(10**(np.flip(colorticks_h))) )
                                nticks = len(self.colorticks)
                            else:
                                if (self.cmax_od > 3):
                                    linticks = np.linspace( 0, self.cmax_od, num=np.int(nticks) )
                                    self.colorticks = np.array( [0] + \
                                                             list(10**(linticks[1:])) )
                                else:
                                    self.colorticks = np.logspace( self.cmin_od, self.cmax_od, 
                                                                   num=np.int(nticks) )                                        
                        else:
                            nticks = nticks_init
                            linticks = np.linspace( self.cmin_od, self.cmax_od, 
                                                    np.int(nticks) )                               
                            zero_ind = np.where( linticks == 0 )
                            self.colorticks = 10**( linticks )
                            self.colorticks[zero_ind] = 0


                        if self.cmin_od < 0:
                            self.linthresh = 10**(self.cmin_od)
                        else:
                            self.linthresh = 1e1
                        self.linscale = 1.5

                    elif self.sign_sum == 0:
                        if (self.cmin_od > 0) or (self.cmax_od > 0):
                            if diff:
                                self.max_od = np.max( [self.cmin_od, self.cmax_od] )
                                self.cmax_od = self.max_od
                                self.cmin_od = self.max_od
                            nticks_init = self.cmax_od - self.cmin_od * self.cmin_sign + 1

                            if nticks_init > 6:
                                if diff:
                                    Nticks_list = [4,5,3,6,2]
                                else:
                                    Nticks_list = [9,8,7,6,5,10,11,12,13,4,3,2]
                                check_loop = True
                                for nticks in Nticks_list:
                                    if diff:
                                        tmparray = np.linspace( 0,self.cmax_od )
                                    else:
                                        tmparray = np.linspace(self.cmin_od * self.cmin_sign, 
                                                               self.cmax_od, np.int(nticks))

                                    checksum = np.sum( np.abs(tmparray-tmparray.astype('i')) )
                                    if np.abs(checksum) < 1e-7:
                                        check_loop = False
                                        break

                                if diff:
                                    if check_loop:
                                        interval = np.ceil( self.cmax_od / 3 ).astype('I')
                                        tick_start = self.cmax_od
                                        colorticks_h = []
                                        nticks = 0
                                        while (tick_start > 0):
                                            colorticks_h.append( tick_start )
                                            tick_start -= interval
                                            nticks += 1
                                        self.colorticks = list(-10**(np.array(colorticks_h))) + \
                                                       [0] + list( 10**(np.flip(colorticks_h) ) )
                                        nticks = len(self.colorticks)
                                    else:
                                        linticks_hl = np.linspace( -self.cmin_od, 0, nticks )
                                        linticks_hr = np.linspace( 0, self.cmax_od, nticks )[1:]
                                        self.colorticks = 10**( linticks_hl ) * -1 + \
                                                          10**( linticks_hr ) 
                                        self.colorticks[nticks] = 0
                                        nticks = nticks * 2 - 1
                                else:
                                    linticks = tmparray
                                    zero_ind = np.where( linticks == 0 )
                                    self.colorticks = 10**( np.abs(linticks) ) * np.sign(linticks)
                                    self.colorticks[zero_ind] = 0                                            

                            else:
                                nticks = nticks_init
                                linticks = np.linspace( self.cmin_od * self.cmin_sign, 
                                                    self.cmax_od, np.int(nticks) )
                                zero_ind = np.where( linticks == 0 )
                                self.colorticks = 10**( linticks ) * np.sign(linticks)
                                self.colorticks[zero_ind] = 0

                            self.linthresh = 1e1
                            self.linscale = 1.5

                        else:
                            if diff:
                                self.cmax_od, self.cmin_od = \
                                    np.max( [self.cmin_od, self.cmax_od] ), \
                                    np.max( [self.cmin_od, self.cmax_od] )

                            if self.log_scale_min == None:
                                self.min_order = np.min( [self.cmin_od, self.cmax_od] )
                                if self.min_order > 4:
                                    self.cmin_p = 10                           
                                elif self.min_order > 0:
                                    self.cmin_p = 1
                                elif self.min_order <= 0:
                                    self.cmin_p = 10**(self.min_order) * 1e-4
                            else:
                                self.cmin_p = self.log_scale_min

                            self.colorticks = np.array( [ -10**(self.cmin_od),
                                          -np.sqrt( 10**(self.cmin_od) * self.cmin_p ),
                                                          -self.cmin_p,
                                                          0, 
                                                          self.cmin_p,
                                           np.sqrt( 10**(self.cmax_od) * self.cmin_p ),
                                                          10**( self.cmax_od ) ] )
                            nticks = 7
                            self.linthresh = self.cmin_p
                            self.linscale = 1.5

                    elif self.sign_sum in [-2, -1]:
                        nticks_init = self.cmax_od - self.cmin_od + 1
                        if nticks_init > 6:
                            Nticks_list = [9,8,7,6,5,10,11,12,13,4]
                            for nticks in Nticks_list:
                                tmparray = np.linspace(self.cmin_od, self.cmax_od, 
                                                       np.int(nticks))
                                checksum = np.sum( np.abs(tmparray - tmparray.astype('i')) )
                                if np.abs(checksum) < 1e-7:
                                    break
                        else:
                            nticks = nticks_init
                        self.colorticks = np.logspace( -self.cmin_od, -self.cmax_od, 
                                                       num=np.int(nticks) )
                        self.linthresh = -self.cmax
                        self.linscale = -self.cmax

                    self.nticks = nticks
                else:
                    self.cbprop = get_cbar_prop( [self.var_slice], Ntick_set=self.nticks,
                                                **self.kwd_pretty_tick )
                    self.colorticks = self.cbprop.colorticks
                    self.colorlabels = self.cbprop.colorlabels
                self.cmin = self.colorticks[0]
                self.cmax = self.colorticks[-1]
        else:
            if np.shape(self.colorticks) == ():
                if self.log_scale:
                    self.cmin = np.min(self.var_slice[self.var_slice != 0])
                    self.cmax = np.max(self.var_slice[self.var_slice != 0])

                    self.cmin_od = np.floor( np.log10( np.abs(self.cmin) ) )
                    self.cmax_od = np.floor( np.log10( np.abs(self.cmax) ) )
                    self.cmin_sign = np.sign( self.cmin )
                    self.cmax_sign = np.sign( self.cmax )
                    self.min_order = np.min( [self.cmin_od, self.cmax_od] )
                    if self.min_order > 4:
                        self.cmin_p = 10                           
                    elif self.min_order > 0:
                        self.cmin_p = 1
                    elif self.min_order <= 0:
                        self.cmin_p = 10**(self.min_order) * 1e-4

                    self.colorticks = np.array( [ self.cmin,
                                                 -np.sqrt( -self.cmin * self.cmin_p ),
                                                 -self.cmin_p,
                                                  0,
                                                  self.cmin_p,
                                                  np.sqrt( self.cmax * self.cmin_p ),
                                                  self.cmax ] )
                    self.linthresh = self.cmin_p
                    self.linscale = 1.5

                else:                        
                    self.colorticks = np.linspace( self.cmin, self.cmax, self.nticks )
            else:
                self.linthresh = np.min( np.abs( self.colorticks[self.colorticks != 0] ) )
                self.linscale = 1.5

            if np.shape(self.colorlabels) == ():
                if not self.log_scale:
                    self.colorlabels = np.copy( self.colorticks )
            self.cmin = self.colorticks[0]
            self.cmax = self.colorticks[-1]
            

        # Adjust colorticks in case of maximum value > 1e3
        # make colorlabels for log_scale plot
        if self.log_scale:
            if (self.sign_sum in [1,2]) & (np.max(self.colorticks) > 1e3):
                ind_above_1 = np.where( self.colorticks > 1 )
                self.colorticks = [0] + list(self.colorticks[ind_above_1])
                self.nticks = len(self.colorticks)
            elif (self.sign_sum in [-1,-2]) & (np.min(self.colorticks) < 1e-3):
                ind_below_m1 = np.where( self.colorticks < -1 )
                self.colorticks = list( self.colorticks[ind_below_m1] ) + [0]
                self.nticks = len(self.colorticks)

            self.colorlabels = []
            for ct in self.colorticks:
                if ct == 0:
                    lbtmp = '0'
                else:
                    if ct < 0:
                        tmpind = 5
                        sign = "-"
                    elif ct > 0:
                        tmpind = 4
                        sign = ""

                    if format( np.abs(ct), '.4e' )[tmpind-2:tmpind] == '00':
                        lbtmp = format( ct, '.2e' )[tmpind:]
                        lbtmp = sign + '$\\mathdefault{' + lbtmp.replace('e','10^{')[:-2] + \
                                 str(int(lbtmp[-2:])) + '}}$'
                    else:
                        lbtmp = format( ct, '.4e' )[:tmpind] + format( ct, '.2e' )[tmpind:]
                        lbtmp = '$\\mathdefault{' + lbtmp.replace('e','\\times10^{')[:-2] + \
                                 str(int(lbtmp[-2:])) + '}}$'
                self.colorlabels.append( lbtmp.replace( '+', '' ) )

        # Dictionary declaration for keywords in external function call
        self.kwd_pcolormesh = {}
        self.kwd_polycollection = {}
        self.kwd_title = {}
        self.kwd_unit = {}
        
        # Construct pcolormesh dictionary
        if self.grid_line:
            if self.model_type == 'FV':
                self.kwd_pcolormesh['edgecolors'] = 'black'
                self.kwd_pcolormesh['lw'] = self.grid_line_lw
            elif self.model_type == 'SE':
                self.kwd_polycollection['edgecolor'] = 'black'
                self.kwd_polycollection['lw'] = self.grid_line_lw
        else:
            if self.model_type == 'SE':
                self.kwd_polycollection['edgecolor'] = 'face'
            
        if self.log_scale:
            #if self.model_type == 'FV':
            #    # for compability with some matplotlib version
            try:
                self.kwd_pcolormesh['norm'] = \
                    matplotlib.colors.SymLogNorm( linthresh=self.linthresh, 
                                                  linscale=self.linscale,
                                                  vmin=self.cmin, vmax=self.cmax )
                self.kwd_polycollection['norm'] = \
                    matplotlib.colors.SymLogNorm( linthresh=self.linthresh, 
                                                  linscale=self.linscale,
                                                  vmin=self.cmin, vmax=self.cmax )
            except:                    
                self.kwd_pcolormesh['norm'] = \
                    matplotlib.colors.SymLogNorm( linthresh=self.linthresh,
                                                  linscale=self.linscale,
                                                  vmin=self.cmin, vmax=self.cmax,
                                                  base=10 ) 
                self.kwd_polycollection['norm'] = \
                    matplotlib.colors.SymLogNorm( linthresh=self.linthresh,
                                                  linscale=self.linscale,
                                                  vmin=self.cmin, vmax=self.cmax,
                                                  base=10 ) 
                
            #elif self.model_type == 'SE':
            #    raise ValueError( 'log scale for SE model is currently not supported' )
            
        if self.title != "":
            if self.title_bold:
                self.kwd_title['weight'] = 'semibold'
        
        if self.unit != "":
            if self.unit_bold:
                self.kwd_unit['weight'] = 'semibold'
            if self.unit_italic:
                self.kwd_unit['style'] = 'italic'
            
        # ========================== END Initial Setup ==========================
        # =======================================================================

        # Call 2D plot code
        self.plot()

    # ========================================================================
    # =============================== Plotting ===============================
    # ========================================================================
    def plot(self):
        
        # === Set font family for the whole plot first ===
        plt.rcParams['font.family'] = self.font_family
        
        # === FV model output (2D with longitude and latitude values) ===
        if self.model_type == 'FV':
            if self.center_180:
                self.lon += 180
                self.lon[ self.lon > 180. ] -= 360

            self.im = self.ax.pcolormesh(self.lon, self.lat, self.var, 
                                         cmap=self.cmap, transform=self.projection, 
                                         vmin=self.cmin, vmax=self.cmax,
                                         **self.kwd_pcolormesh )
        elif self.model_type == 'SE':
            self.im = PolyCollection( self.verts, cmap=self.cmap,
                                      **self.kwd_polycollection )
            self.im.set_array( self.var )
            self.im.set_clim( vmin=self.cmin, vmax=self.cmax )
            self.ax.add_collection(self.im)
        
        # === Set longitude & Latitude labels & lines ===
        self.ax.set_xlim(self.lon_range)
        self.ax.set_ylim(self.lat_range)
        self.ax.tick_params(labelsize=self.label_size)
        if self.coast:
            self.ax.coastlines(resolution=self.resolution, 
                               lw=self.feature_line_lw, color=self.feature_color )
        if self.country:
            self.ax.add_feature(cfeature.BORDERS.with_scale(self.resolution), 
                                lw=self.feature_line_lw, edgecolor=self.feature_color )
        if self.state:
            self.ax.add_feature(cfeature.STATES.with_scale(self.resolution), 
                                lw=self.feature_line_lw, edgecolor=self.feature_color )
        if self.lonlat_info:
            if self.lon_labels==None:
                self.lonticklabel = np.arange(self.lon_range[0],self.lon_range[1]+0.1,
                                              self.lon_interval)
            else:
                if self.lon_labels[0] > self.lon_labels[-1]:
                    self.lonticklabel = np.copy( self.lon_labels )
                    if self.center_180:
                        self.lonticklabel += 180
                        self.lonticklabel[ self.lonticklabel > 180. ] -= 360    

                else:
                    self.lonticklabel = self.lon_labels
            
            if self.lat_labels==None:
                self.latticklabel = np.arange(self.lat_range[0],self.lat_range[1]+0.1,
                                              self.lat_interval)
            else:
                self.latticklabel = self.lat_labels


            self.ax.set_xticks(self.lonticklabel,crs=self.ax.projection)
            self.ax.set_yticks(self.latticklabel,crs=self.ax.projection)
            self.ax.set_xlabel('')
            self.ax.set_ylabel('')

            self.lon_formatter = LongitudeFormatter(zero_direction_label=True)
            self.lat_formatter = LatitudeFormatter()
            self.ax.xaxis.set_major_formatter(self.lon_formatter)
            self.ax.yaxis.set_major_formatter(self.lat_formatter)

            
            if self.regional:
                # if self.center_180:
                #     self.lonticklabel += 180
                #     self.lonticklabel[ self.lonticklabel > 180. ] -= 360    
                # self.gl2 = self.ax.gridlines( lw=1.0, color='black', alpha=0.5, linestyle=':' )

                # if self.center_180:
                #     self.lonticklabel += 180
                #     self.lonticklabel[ self.lonticklabel > 180. ] -= 360    
                #self.gl = self.ax.gridlines( lw=1.0, color='black', alpha=0.5, linestyle=':' )
                #self.gl.xlocator = ticker.FixedLocator( self.lonticklabel )                
                
                if self.lonlat_line:
                    for lontick in self.lonticklabel:
                        self.ax.plot( [lontick, lontick], self.ax.get_ylim(), 
                                      lw=1.0, linestyle=':', color='black', alpha=0.5 )
                    for lattick in self.latticklabel:
                        self.ax.plot( self.ax.get_xlim(), [lattick, lattick],
                                      lw=1.0, linestyle=':', color='black', alpha=0.5 )
            
            else:
                if self.lonlat_line:
                    self.ax.grid( lw=1.0, color='black', alpha=0.5, linestyle=':')
  
        # === Set colorbar properties ===
        if self.colorbar:
            self.cb = self.fig.colorbar( self.im, ax=self.ax, orientation=self.orientation,
                                         shrink=self.shrink, pad=self.pad, extend=self.extend,
                                         fraction=self.fraction, ticks=self.colorticks )

            self.cb.ax.tick_params(labelsize=self.label_size-1)
            #if not (self.log_scale & self.pretty_tick):
            self.cb.ax.set_xticklabels( self.colorlabels, size=13 )       

            # === Add a unit if specified ===
            if self.unit != "":
                self.cb.ax.text( 1.02+self.unit_offset[0], 1.0+self.unit_offset[1], self.unit,
                                 size=self.unit_size, ha='left', va='top',
                                 transform=self.cb.ax.transAxes, **self.kwd_unit )
        
        # === Add a title if specified ===
        if self.title != "":
            self.ax.set_title( self.title, fontsize=self.title_size, 
                               fontfamily=self.font_family, **self.kwd_title )
            
    # ============================= END Plotting =============================
    # ========================================================================

    # ===== Defining __call__ method =====
    def __call__(self):
        print( '=== var ===')
        print( np.shape(var) )

        
        
class get_cbar_prop(object):
    
    def __init__(self, arrays, colorticks=None, colorlabels=None, 
                       ranges=None, Nticks_list=None, 
                       max_find_method='ceil', tick_find_method='ceil',
                       min_set=None, max_set=None, Ntick_set=None ):
        
        # set ranges
        if not ranges:
            ranges = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 
                      0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 
                      0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 
                      1, 2, 3, 4, 5, 6, 8, 
                      10, 12, 15, 20, 25, 30, 40, 50, 
                      60, 80, 100, 150, 200, 300, 400, 500, 600, 800, 
                      1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000 ]
        self.ranges = ranges
        
        # Calculate max value
        if max_set == None:
            maxval = np.nanmax( arrays )
            self.maxval = maxval
        
        # Calculate min value
        if min_set == None:
            minval = np.nanmin( arrays )
            self.minval = minval
        
        # Calculate plotmax
        if max_set:
            plotmax = max_set
            self.maxval = plotmax
            
        else:
            if maxval in ranges:
                plotmax = maxval

            elif max_find_method=='ratio':
                tmp = np.log10( maxval / np.array(ranges)  )
                plotmax = ranges[ np.argmin( np.abs(tmp) ) ]
                #print( 'ratio', plotmax )

            elif max_find_method=='floor':
                tmp = maxval - np.array(ranges)
                plotmax = ranges[ np.where( tmp >= 0. )[0][-1] ]
                #print( 'floor', plotmax )

            elif max_find_method=='ceil':
                tmp = maxval - np.array(ranges)
                plotmax = ranges[ np.where( tmp >= 0. )[0][-1]+1 ]
                #print( 'ceil', plotmax )

        self.plotmax = plotmax
        
        # Calculate plotmin
        if min_set != None:
            plotmin = min_set
            self.minval = plotmin
        else:
            plotmin = 0
        self.plotmin = plotmin
        
        
        if Ntick_set == None:
            
            # Set Nticks_list
            if not Nticks_list:
                Nticks_list = [ 5, 6, 4, 3, 2 ]
            self.Nticks_list = Nticks_list

            # Construct colorticks array
            for ntick in Nticks_list:
                tmparray = np.linspace(plotmin, plotmax, ntick)
                if maxval < 3:
                    factor = 10**( -np.floor(np.log10( plotmax )) + 1 )
                    checksum = np.sum( np.abs(tmparray*factor - (tmparray*factor).astype('i')) )
                    #print( ntick, checksum, factor )
                else:
                    checksum = np.sum( np.abs(tmparray - tmparray.astype('i')) )
                    #print( ntick, checksum, factor )

                if np.abs(checksum) < 1e-7:
                    #print( 'ntick', ntick )
                    #print( tmparray )
                    break

        else:
            ntick = Ntick_set
            self.Nticks_list = [ Ntick_set ]
            tmparray = np.linspace(plotmin, plotmax, ntick)
                    
                
        colorticks = tmparray
        checksum = np.sum( np.abs( tmparray - tmparray.astype('i') ) )
        if checksum < 1e-7:
            colorlabels = tmparray.astype('i')
        else:
            colorlabels = colorticks.astype('S')
            for i, cbt in enumerate( colorticks ):
                colorlabels[i] = '{:.8}'.format(cbt)
            
                
        self.Nticks = ntick
        self.colorticks = colorticks
        self.colorlabels = colorlabels.astype('U')

                       
    # Defining __call__ method 
    def __call__(self): 
        
        print( 'ranges', self.ranges )    
        print( 'maxval', self.maxval )
        print( 'minval', self.minval )
        print( 'plotmax', self.plotmax )
        print( 'plotmin', self.plotmin )
        print( 'Nticks_list', self.Nticks_list )
        print( 'Nticks', self.Nticks )
        print( 'colorticks', self.colorticks )
        print( 'colorlabels', self.colorlabels )
