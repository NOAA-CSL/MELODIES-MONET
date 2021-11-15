def read_OMPS_nm(fname):
    import xarray as xr
    import h5py
    from datetime import datetime,timedelta
    
    with h5py.File(fname,'r') as f:
        time = f['GeolocationData']['Time'][:]
        to3 = f['ScienceData']['ColumnAmountO3'][:]
        lat = f['GeolocationData']['Latitude'][:]
        lon = f['GeolocationData']['Longitude'][:]
        aprior = f['AncillaryData']['APrioriLayerO3'][:]
        plevs = f['DimPressureLevel'][:]
        layere = f['ScienceData']['LayerEfficiency'][:]
        flags = (f['ScienceData']['QualityFlags'][:])
        cloud_fraction = f['ScienceData']['RadiativeCloudFraction'][:]
        
    to3[((to3 < 50.0)|(to3 > 700.0))] = -999.
    layere[((layere < 0.)|(layere > 10.))] = 0
    time[((time < -5e9)|(time > 1e10))] = -999.
    ssday = datetime(year=1993,day=1,month=1,hour=0)
    time = np.asarray([ssday + timedelta(seconds = i) for i in time])
   
    to3[(cloud_fraction > .3)] = -999.
    layere[(cloud_fraction > .3),:] = 0
    to3[(flags >= 138 )] = -999.
    
    ds = xr.Dataset(
            {
            'ozone_column': (['x','y'],to3),
            'apriori': (['x','y','z'],aprior),
            'layer_efficiency': (['x','y','z'],layere),
            },
            coords={
                'longitude':(['x','y'],lon),
                'latitude':(['x','y'],lat),
                'time':(['x'],time),
                'pressure':(['z'],plevs),
            },
            attrs={'missing_value':-999}
        )         
    return ds