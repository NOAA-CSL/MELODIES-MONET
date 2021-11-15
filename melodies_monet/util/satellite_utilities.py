# File started by Maggie Bruckner. As of 11/12/2021, the functions in this file are 
# only included as examples of things that might go in this file.

import xesmf as xe
import numpy as np
from datetime import datetime,timedelta
from numba import jit

@jit(nopython=True)
def model_co_column(obs_pressure,ak,model_pressure,model_co,model_z,out_array):
    '''Returns CO column for model with averaging kernel applied. Usage of jit should 
        speed up calculation.
    '''
    n1,n2,nlkern = obs_pressure.shape
    fac = 1.0e4*6.023e23/28.97/9.8/1000.
    # initialize intermediates
    coint = np.zeros(nlkern,dtype='float64')
    zint = np.zeros(nlkern,dtype='float64') # set to zeros in loop
    
    for i in range(n1):
        for j in range(n2):
            # calculating column
            dum0 = obs_pressure[i,j,:]
            dum1 = model_pressure[:,i,j]
            index = np.where((dum0 >= np.min(dum1)) & (dum0 <= np.max(dum1)))[0]
            if len(dum0[~np.isnan(dum0)]) > 0 and len(index) > 0: 
                coint[:] = 0
                zint[:] = 0 
                
                amin = np.argmin(index)
                amax = np.argmax(index)
                coint[amin:amax] = np.interp(dum0[amin:amax],dum1[::-1],model_co[::-1,i,j])
                zint[amin:amax] = np.interp(dum0[amin:amax],dum1[::-1],model_z[::-1,i,j])
                plm1 = 0.5 * (dum0[nlkern-1]+dum0[nlkern-2])
                zlm1 = 0.5 * (zint[nlkern-1]+zint[nlkern-2])
                
                for k in np.arange(nlkern-3,-1,-1):
                    pl=.5*(dum0[k+1]+dum0[k])
                    zl=.5*(zint[k+1]+zint[k])

                    dp=plm1-pl
                    dz=zl-zlm1
                    
                    if dp > 0.0 and dz > 0.0:
                        out_array[i,j]=out_array[i,j]+dp*coint[k]*ak[i,j,k]*fac/dz
                    plm1=pl
                    zlm1=zl
                dp=2.*(dum0[0]-plm1)
                dz=2.*(zint[0]-zlm1)
                if dp > 0.0 and dz > 0.0:    
                    out_array[i,j] = out_array[i,j]+dp*coint[0]*ak[i,j,0]*fac/dz
                
    return out_array*1.e-18

def model_to_swath_time(mod,model_times,sat_times,obs):
    ''' Interpolation of model data to satellite observations in time. Also performs
        calculation of CO column.
    '''
    
    import pandas as pd
    import numpy as np
    nf,nz,nj,nk = mod['co'].shape
    file_tstep = (model_times[1] - model_times[0])
    # select indicies of files needed for interpolation
    files_for_interp = []
    for i in range(nf):
        if model_times[i] <= sat_times[0,0]: 
            if np.abs(sat_times[0,0]-model_times[i]) <= file_tstep: files_for_interp.append(i)
        elif model_times[i] >= sat_times[0,0]:
            if np.abs(sat_times[0,0]-model_times[i]) <= file_tstep: files_for_interp.append(i)
            elif np.abs(model_times[i] - sat_times[-1,-1]) <= file_tstep: files_for_interp.append(i)
        else: pass
    
    # calculate model columns, with time interpolation
    column_out = np.full((nj,nk),0,dtype='float32')
    for i in files_for_interp:
        if len(files_for_interp) > 1:
            tind = ((np.abs(sat_times - model_times[i]) <= file_tstep))
            out_array = np.full_like(column_out[tind],0,dtype='float32')
            tfac = np.abs(sat_times-model_times[i])/file_tstep
            column_out[tind] += tfac*model_co_column(obs.obj['pressure_levels'].values[tind],\
                                                     obs.obj['column_averaging_kernel'].values[tind]\
                                                             ,mod['pressure'].values[i],mod['co'].values[i],\
                                                             mod['zdash'].values[i],out_array)
        
    return column_out