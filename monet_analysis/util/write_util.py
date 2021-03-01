import numpy as np

def write_ncf(dset,output_name, title=''):
    """Short summary.

    Parameters
    ----------
    dset : type
        Description of parameter `dset`.
    output_name : type
        Description of parameter `output_name`.

    Returns
    -------
    type
        Description of returned object.

    """
    import pandas as pd
    print('Writing:', output_name)
    comp = dict(zlib=True,complevel=7)
    encoding = {}
    for i in dset.data_vars.keys():
        mask_and_scale = ['float','float32','']
        if dset[i].dtype == 'float'dset.data_vars.keys()
        encoding[i] = comp
    dset.attrs['title'] = title
    dset.attrs['format'] = 'NetCDF-4'
    dset.attrs['date_created'] = pd.to_datetime('today').strftime('%Y-%m-%d')
    dset.to_netcdf(output_name,encoding=encoding)


def compute_scale_and_offset(mn, mx, n,dtype=np.float32):
    """
    min is the minumum of the values
    max is the maximum of the values
    n is the integer bit length  (ie 32 for np.int32 or 16 for np.int16)
    """
    # stretch/compress data to the available packed range
    scale_factor = (mx - mn) / (2 ** n - 1)
    # translate the range to be symmetric about zero
    add_offset = mn + 2 ** (n - 1) * scale_factor
    return (scale_factor.astype(dtype), add_offset.astype(dtype))

def pack_value(values, scale_factor, offset, dtype):
    return ((values - offset) / scale_factor).astype(dtype)

def get_min_max(da):
    return (da.min().compute(), da.max().compute())

def compress_variable(da):
    da = da.fillna(-1)
    mn,mx = get_min_max(da)
    scale_factor,offset=compute_scale_and_offset(mn,mx,32,dtype=da.dtype)
    da.data = pack_value(da,scale_factor,offset,dtype=np.int32).data
    da.attrs['scale_factor'] = scale_factor.values
    da.attrs['add_offset'] = offset.values
    da.attrs['_FillValue'] = -1
    da.attrs['missing_value'] = -1
    return da
