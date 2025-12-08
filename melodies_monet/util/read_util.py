# SPDX-License-Identifier: Apache-2.0
#
import os
import warnings
import pandas as pd
import xarray as xr
import glob


def read_saved_data(analysis, filenames, method, attr, xr_kws={}):
    """Read previously saved dict containing melodies-monet data (:attr:`paired`, :attr:`models`, or :attr:`obs`)
    from pickle file or netcdf file, populating the :attr:`paired`, :attr:`models`, or :attr:`obs` dict.

    Parameters
    ----------
    analysis : melodies_monet.driver.analysis
        Instance of the analysis class from driver script.
    filenames : str or iterable
        str or list for reading in pkl. For netCDF, must be dict with format {group1:str or iterable of filenames, group2:...}
    method : str
        One of either 'pkl' or 'netcdf'.
    attr : str
        The analysis attribute that will be populated with the saved data. One of either 'paired' or 'models' or 'obs'.
    **kwargs : optional
        Additional keyword arguments for xr.open_dataset()

    Returns
    -------
    None
    """
    from glob import glob
    from .. import tutorial
    
    # Determine where to read files from
    if getattr(analysis,'output_dir_read') is not None:
        read_dir = getattr(analysis,'output_dir_read')
    else:
        read_dir = ''
    
    # expand any wildcards in the filenames
    if method=='pkl':
        if isinstance(filenames,str):
            files = sorted([file for sublist in [glob(os.path.join(read_dir,file)) for file in [filenames]] for file in sublist])
        else:
            files = sorted([file for sublist in [glob(os.path.join(read_dir,file)) for file in filenames] for file in sublist])
        if not files:
            raise FileNotFoundError('No such file: ',filenames)
    elif method=='netcdf':
        if isinstance(filenames,dict): 
            files = {}
            for group in filenames.keys():
                if isinstance(filenames[group],str):
                    files[group] = sorted([file for sublist in [glob(os.path.join(read_dir,file)) for file in [filenames[group]]] for file in sublist])
                else:
                     if filenames[group][0].startswith("example:"):
                        files[group] = sorted([file for sublist in [
                            [tutorial.fetch_example(":".join(s.strip() for s in file.split(":")[1:]))] for file in filenames[group]] for file in sublist])
                     else:
                         files[group] = sorted([file for sublist in [glob(os.path.join(read_dir,file)) for file in filenames[group]] for file in sublist])
                if not files[group]:
                    raise FileNotFoundError('No such file: ', filenames[group])
        else:
            raise TypeError('NetCDF format filenames need to be specified as a dict, with format {group1:str or iterable of filenames, group2:...}')
    
    # Set analysis.read such that it now contains expanded filenames so user has list of read files
    expanded_filenames = getattr(analysis,'read')
    expanded_filenames[attr]['filenames'] = files 
    setattr(analysis, 'read', expanded_filenames)

    # for converting name of attribute to name of class for constructing
    class_names = {'paired':'pair','models':'model','obs':'observation'}

    if method=='pkl':
        if len(files)==1:
            setattr(analysis, attr, read_pkl(files[0]))
        elif len(files)>1:
            for count, file in enumerate(files):
                if count==0:
                    attr_out = read_pkl(file)
                else:
                    attr_append = read_pkl(file)
                    for group in attr_out.keys():
                        attr_out[group].obj = xr.merge([attr_out[group].obj,attr_append[group].obj])
            setattr(analysis, attr,  attr_out)

    elif method=='netcdf':
        xr_dict = {}
        for group in files.keys():
            if isinstance(files[group],str):
                group_files = [files[group]]
            else:
                group_files = files[group]
            xr_dict[group] = read_analysis_ncf(group_files,xr_kws)
        setattr(analysis, attr,  xarray_to_class(class_type=class_names[attr],group_ds=xr_dict))

def read_pkl(filename):
    """Function to read a pickle file containing part of the analysis class (models, obs, paired)

    Parameters
    ----------
    filename : type
        Description of parameter `filename`.
        
    Returns
    -------
    obj : type
        Description of returned object.

    """
    from joblib import load
    
    print('Reading:', filename)
    with open(filename, 'rb') as file:
        obj = load(file)
        
    return obj

def read_analysis_ncf(filenames,xr_kws={}):
    """Function to read netcdf4 files containing an object within an attribute of a part of the
    analysis class (models, obs, paired). For example, a single model/obs pairing or a single model. 
    If the object is saved in multiple files, the function will merge the files.

    Parameters
    ----------
    filenames : str or iterable
        Description of parameter `filename`.
    xr_kws : optional
        Additional keyword arguments for xr.open_dataset()
        
    Returns
    -------
    ds_out : type
        Xarray dataset containing merged files.

    """
    if len(filenames)==1:
        print('Reading:', filenames[0])
        ds_out = xr.open_dataset(filenames[0],**xr_kws)
        
    elif len(filenames)>1:
        for count, file in enumerate(filenames):
            print('Reading:', file)

            if count==0:
                ds_out = xr.open_dataset(file,**xr_kws)
                group_name1 =  ds_out.attrs['group_name']

            else:
                ds_append = xr.open_dataset(file,**xr_kws)
                # Test if all the files have the same group to prevent merge issues
                if group_name1 != ds_append.attrs['group_name']:
                    raise Exception('The group names are not consistent between the netcdf files being read.') 
                else:
                    ds_out = xr.merge([ds_out,ds_append])
            
    return ds_out

def xarray_to_class(class_type,group_ds):
    """Remake dict containing driver class instances from dict of xarray datasets. Dict of xarray datasets must contain 
    global attribute that contains json formatted class attributes.

    Parameters
    ----------
    class_type : str
        One of 'model', 'pair' or 'observation'
    group_ds : dict
        dict containing xarray datasets from read_grouped_ncf.

    Returns
    -------
    class_dict
    """
    import json
    from melodies_monet import driver
    
    class_dict = {}
    for group in group_ds.keys():
        if class_type == 'pair':
            c=driver.pair()
        elif class_type == 'model':
            c=driver.model()
        elif class_type == 'observation':
            c=driver.observation()

        obj_dict = json.loads(group_ds[group].attrs['dict_json'])
        
        for attr in obj_dict.keys():
            setattr(c, attr, obj_dict[attr])
        c.obj = group_ds[group]
        class_dict[group]=c

    return class_dict

def read_aircraft_obs_csv(filename,time_var=None):
    """Function to read .csv formatted aircraft observations.

    Parameters
    ----------
    filename : str 
        Filename of .csv file to be read
    time_var : optional
        The variable in the dataset that should be converted to 
        datetime format, renamed to `time` and set as a dimension.
        
    Returns
    -------
    ds_out : xarray.Dataset
        Xarray dataset containing information from .csv file

    """
    df = pd.read_csv(filename)
    if time_var is not None:
        df.rename(columns={time_var:'time'},inplace=True)
        df['time']  = pd.to_datetime(df['time'])
        
    # Sort the values based on time
    df.sort_values(by='time',inplace=True,ignore_index=True)
        
    df.set_index('time',inplace=True)
    
    return xr.Dataset.from_dataframe(df)


def read_site_excel(data_path, site_data, site_number=None, **kwargs):
    """Load a site from from an MS Excel sheet.
    Currently optimized for CDPHE's VOC canister data.

    Parameters
    ----------
    data_path: str
        Path to the excel file containing the data
    site_data: dict
        dict containing the data of the site.
        Required keys are:
            'coords': {'latitude': float, 'longitude': float}
                Coordinates of the site
            'sheet_name': str
                name of the excel sheet containing the site
        Optional keys:
            skiprows: int
                rows that should be skipped when reading the excel sheet.
                Defaults to 0.
            headers: int | list[int]
                rows with headers for building the data frame.
                Take into account that 'skiprows' takes precedence.
                I.e., if the headers are in rows [15,16], but you
                already typed 'skiprows': 15, you should type
                'headers': [0,1].
                Defaults to 0.
            site_id: str
                ID of the site. Defaults to sheet_name.
            qualifier_name: str
                Name of row containing the qualifiers.
                If None, 'Qualifier' is assumed.
            ignore_qualifiers: None
                If None, only data without qualifiers is plotted.
            na_values: scalar | str | list-like | dict | default None
                Values that are NaN
    site_number: int
        Number of site. If only one site is provided, it should be 0.
        This keyword is set for clearer compilation of multiple sites.

    Returns
    -------
    xr.Dataset
        Dataset containing the information of the site
    """
    params = {
        **{
            "skiprows": None,
            "headers": 0,
            "site_id": site_data["sheet_name"],
            "qualifier_name": "Qualifier",
            "keep_qualifiers": None,
            "na_values": None,
            "timezone": "UTC",
            "analyte": ("Analysis", "Analyte"),
            "time_var": ("Sample", "Date"),
            "unit_var": ("Detection", "Units"),
            "results": ('CAS', 'Result'),
            "repeated_values": None,
            "sampling_start_hour": 6,
            "sampling_length": 3,
        },
        **site_data,
    }
    site_number = 0 if site_number is None else site_number
    if isinstance(params['na_values'], str):
        params['na_values'] = [params['na_values']]
    data = pd.read_excel(
        data_path,
        sheet_name=params["sheet_name"],
        skiprows=params["skiprows"],
        header=params["header"],
        na_values=(params["na_values"] + [r'^\s*$']),
        keep_default_na=True,
    ).dropna(how='all')

    data = _apply_qualifiers(data, params['qualifier_name'], params['keep_qualifiers'])
    data.loc[:, "time_local"] = (
            data[params["time_var"]].dt.normalize()
            + pd.Timedelta(params['sampling_start_hour'], 'h')
    )
    data.loc[:, "time_utc"] = (
            data["time_local"].dt.tz_localize(params["timezone"]).dt.tz_convert(None)
    ).copy()
    variables = data[params["analyte"]].unique()
    compiled_data = xr.Dataset()
    for v in variables:
        tmp_data = data.loc[data[params["analyte"]] == v]
        tmp_ds = xr.Dataset()
        tmp_ds["time"] = (("time",), tmp_data["time_utc"].values)
        tmp_ds["x"] = (("x",), [site_number])
        if params["timezone"] not in ['UTC', 'UCT', 'Etc/UTC', 'Etc/UCT', 'Etc/Universal']:
            tmp_ds["time_local"] = (("time", "x"), tmp_data["time_local"].values[..., None])
            time_local = tmp_ds["time_local"].drop_duplicates(dim='time')
        tmp_ds[v] = (("time", "x"), tmp_data[params["results"]].values[..., None])
        try:
            tmp_ds[v].attrs = {"units": tmp_data[params["unit_var"]].unique()}
        except KeyError:
            tmp_ds[v].attrs = {"units": params["unit_var"]}
        try:
            tmp_ds = tmp_ds.groupby('time').mean()
        except ValueError:
            pass
        tmp_ds["time_local"] = time_local
        compiled_data = xr.merge([compiled_data, tmp_ds])

    compiled_data["longitude"] = (("x",), [params["site_coords"]["longitude"]])
    compiled_data["latitude"] = (("x",), [params["site_coords"]["latitude"]])
    compiled_data["siteid"] = (("x",), [params["site_id"]])
    return compiled_data


def _apply_qualifiers(data, qa_name, keep_qualifiers):
    """Apply a qualifier

    Parameters
    ----------
    data: pd.DataFrame
        Dataframe containing the data
    qa_name: str | tuple
        Name of the column containing the qualifiers
    keep_qualifiers: str | list[str] | 
        Qualifiers to keep
    """
    if keep_qualifiers is None or keep_qualifiers != "no":
        return data[data[*qa_name].isnull()]
    elif keep_qualifiers != "all":
        return data[data[*qa_name].isnull() | data[*qa_name].isin(list(keep_qualifiers))]
    return data


def compile_sites_excel(data_path, site_dict):
    """Compiles all sites in a file
    Currently optimized for CDPHE VOC canister data, but should work for most.

    Parameters
    ----------
    data_path : str
        Path to the excel file containing the data.
    site_dict : dict
        Dictionary containing a key per site (ideally, the site's name)
        and a dict as site_value, to pass to read_cdphe_site

    Returns
    -------
    xr.Dataset
        Dataset with all sites compiled
    """

    compiled_data = xr.Dataset()
    for n, s in enumerate(list(site_dict.keys())):
        site_data = site_dict[s]
        ds = read_site_excel(data_path, site_data, site_number=n)
        compiled_data = xr.merge([compiled_data, ds])
    return compiled_data


def control_reading_excel(data_path, site_type, site_dict):
    """Controls the reading and file preparation process.

    Parameters
    ----------
    path : str | list | glob object
        Path to files that should be opened
    site_type : str
        Type of site/excel to read. Currently, only CDPHE VOC Canisters are implemented
    site_dict : dict

    Returns
    -------
    xr.Dataset
        Dataset with excel compiled
    """
    if site_type != "pt_sfc":
        warnings.warn("site_type is not pt_sfc. Will attempt to read anyway")
    return compile_sites_excel(data_path, site_dict)


def read_pandora(path):
    """Calls tools for reading Pandora Global Network text files.
    It is only a wrapper around the pandonia_global_network.py tool


    Parameters
    ----------
    path: str | list[str]
        String containing the paths

    Returns
    -------
    xr.Dataset
        Formatted dataset. Should work for MELODIES-MONET.
    """
    from .pandonia_global_network import open_mfdataset as read_and_format
    return read_and_format(path)


def read_noaa_gml(filename):
    """Function to read .dat and/or netcdf formatted NOAA GML observations.

    Parameters
    ----------
    filename : str
        Filename of .dat or .nc file to be read

    Returns
    -------
    xarray.Dataset
        Xarray dataset containing information from .dat/.nc file

    """

    #NOAA GML files come in two flavors: .dat and .nc
    if filename.lower().endswith('.nc'): #for netcdf version of files

        #open the netcdf file
        dset = xr.open_dataset(filename)
        
        #get and store the time variable
        time = dset.time
        
        #get the id global attribute from the netCDF file and extract the siteid in lower case format
        #and replicate along the time dimension...note that index [3] in the line of code 
        #below assumes that the id global attribute format won't change.
        siteid=[dset.attrs['id'].split("_")[3].lower()]*time.size

        #get the variable in the file.  Note that their should only be two
        #variables in the file: time and the data variable.
        vnames=list(dset.data_vars.keys())

        #Build the xarray dataset to return...currently this is hardwired to 'O3(PPB)', but could be changed
        #to vnames[0] and properly annotated in the conrol file which would allow for different constituents.
        ds = xr.Dataset(data_vars={'O3(PPB)':dset[vnames[0]], 'siteid': (('time'),siteid, {'descripton':'replicated siteid'})},coords={'time':time})

        #close the netcdf file when done
        dset.close()

    else: #for .dat version of files

        # Ignore preamble
        preamble_end_line = 0
        with open(filename, 'r') as f:
            for i, line in enumerate(f):
                if line.strip() and line.strip().split()[0] == 'STN':
                    preamble_end_line = i
                    break

        # Read the file again, skipping the identified preamble lines
        df = pd.read_csv(
            filename,
            skiprows=preamble_end_line,
            sep=r"\s+",
        )
        df = df.rename(
            columns={"YEAR": "year", "MON": "month", "DAY": "day", "HR": "hour", "STN": "siteid"}
        )
        df['time'] = pd.to_datetime(df[["year", "month", "day", "hour"]])
        siteid = df["siteid"][0]
        df = df[df.columns[~df.columns.isin(["year", "month", "day", "hour"])]]

        # Sort the values based on time
        df = df.sort_values(by='time', ignore_index=True)
        df = df.set_index('time')
        ds = xr.Dataset.from_dataframe(df)

    return ds


def read_noaa_gml_multifile(filenames):
    """Function to read .dat formatted NOAA GML observations.

    Parameters
    ----------
    filenames : str | list[str]
        Filenames of .dat file to be read

    Returns
    -------
    xarray.Dataset
        Xarray dataset containing information from .dat file

    """
    if isinstance(filenames, str):
        files = sorted(glob.glob(filenames))
    else:
        files = []
        for file in filenames:
            files.extend(sorted(glob.glob(file)))

    data = read_noaa_gml(files[0])
    if len(files) > 1:
        for file in len(files):
            data2 = read_noaa_gml(file)
            if data['siteid'][0].values != data['siteid'][0].values:
                raise Exception('Only one station at a time is supported')
            data = xr.merge([data, data2], dim='time')
    return data.sortby('time')
