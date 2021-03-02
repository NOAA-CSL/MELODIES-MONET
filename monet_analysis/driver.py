""" This is the overall control file.  It will drive the entire analysis package"""
import monetio as mio
import monet as m
import os
import xarray as xr
from .util import write_ncf

class pair:
    def __init__(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        self.type = 'pt_sfc'
        self.radius_of_influence = 1e6
        self.obs = None
        self.model = None
        self.filename = None

class observation:
    def __init__(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        self.obs = None
        self.label = None
        self.file = None
        self.obj = None
        self.type = 'pt_src'

    def open_obs(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        try:
            if os.path.isfile(self.file):
                _,extension = os.path.splitext(self.file)
                if extension in ['.nc','.ncf','.netcdf','.nc4']:
                    if len(glob(self.file)) > 1:
                        self.obj = xr.open_dataset()
                    self.obj = xr.open_dataset(self.file)
        except ValueError:
            print('something happened opening file')

    def obs_to_df(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        self.obj = self.obj.monet._df_to_da()

class model:
    def __init__(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        self.model = None
        self.file_str = None
        self.files = None
        self.label = None
        self.obj = None
        self.mapping = None

    def glob_files(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        from numpy import sort
        from glob import glob
        self.files = sort(glob(self.file_str))

    def open_model_files(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        self.glob_files()
        if self.model.lower() == 'cmaq':
            if len(self.files) > 1:
                self.obj = mio.cmaq.open_mfdataset(self.files)
            else:
                self.obj = mio.cmaq.open_dataset(self.files[0])
        elif self.model.lower() == 'wrfchem':
            if len(self.files) > 1:
                self.obj = mio.wrfchem.open_mfdataset(self.files)
            else:
                self.obj = mio.wrfchem.open_dataset(self.files)
        elif self.model.lower() == 'rrfs_cmaq':
            if len(self.files) > 1:
                self.obj = mio.rrfs_cmaq.open_mfdataset(self.files)
            else:
                self.obj = mio.rrfs_cmaq.open_dataset(self.files)
        elif self.model.lower() == 'gsdchem':
            if len(self.files) > 1:
                self.obj = mio.fv3chem.open_mfdataset(self.files)
            else:
                self.obj = mio.fv3chem.open_dataset(self.files)

class analysis:
    def __init__(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        self.control = 'control.yaml'
        self.control_dict = None
        self.models = {}
        self.obs = {}
        self.paired = {}

    def read_control(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        import yaml
        with open(self.control, 'r') as stream:
            self.control_dict = yaml.safe_load(stream)

    def open_models(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        if 'model' in self.control_dict:
            # open each model
            for model in self.control_dict['model']:
                m = model()
                m.model = model
                m.label = self.control_dict['model'][model]['label']
                m.glob_files()
                m.mapping = self.control_dict['model'][model]['mapping']
                self.models[m.label] = m

    def open_obs(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        if 'obs' in self.control_dict:
            for obs in self.control_dict['obs']:
                o = observation()
                o.obs = obs
                o.label = obs
                o.obs_type = self.control_dict['obs'][obs]['obs_type']
                o.open_file(self.control_dict['obs'][obs]['filename'])
                self.obs[o.label] = o

    def pair_data(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        pairs = {}
        for m in self.models:
            # Now we have the models we need to loop through the mapping table for each network and pair the data
            # each paired dataset will be output to a netcdf file with 'model_label_network.nc'
            for obs_to_pair in m.mapping.keys():
                # get the variables to pair from the model data (ie don't pair all data)
                model_obj = m.obj[m.mapping[obs_to_pair].keys()]
                ## TODO:  add in ability for simple addition of variables from

                # simplify the objs object with the correct mapping vairables
                obs = self.obs[obs_to_pair]
                # pair the data
                # if pt_sfc (surface point network or monitor)
                if obs.obs_type.lower() == 'pt_sfc':
                    # convert this to pandas dataframe
                    obs.obs_to_df()
                    # now combine obs with
                    paired_data = model_obj.monet.combine_point(obs.obj,radius_of_influence=1e6,suffix=m.label)
                    # this outputs as a pandas dataframe.  Convert this to xarray obj
                    p = pair()
                    p.obs = obs.label
                    p.model = m.label
                    p.filename = '{}_{}.nc'.format(p.obs,p.model)
                    p.obj = paired_data.monet._df_to_da()
                    write_util.write_ncf(p.obj,p.filename) # write out to file
                # TODO: add other network types / data types where (ie flight, satellite etc)

    ### TODO: Create the plotting driver (most complicated one)
    #def plotting(self):
