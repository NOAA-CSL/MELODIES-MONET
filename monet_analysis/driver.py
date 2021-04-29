""" This is the overall control file.  It will drive the entire analysis package"""
import monetio as mio
import monet as m
import os
import xarray as xr
import pandas as pd
import numpy as np
import datetime

# from util import write_ncf


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
        self.model_vars = None
        self.obs_vars = None
        self.filename = None

    def fix_paired_xarray(self, dset=None):

        # first convert to dataframe
        df = dset.to_dataframe().reset_index(drop=True)

        # now get just the single site index
        dfpsite = df.rename({'siteid': 'x'}, axis=1).drop_duplicates(subset=['x'])
        columns = dfpsite.columns  # all columns
        site_columns = [
            'latitude',
            'longitude',
            'x',
            'site',
            'msa_code',
            'cmsa_name',
            'epa_region',
            'state_name',
            'msa_name',
            'site',
            'utcoffset',
        ]  # only columns for single site identificaiton

        # site only xarray obj (no time dependence)
        dfps = dfpsite.loc[:, columns[columns.isin(site_columns)]].set_index(['x']).to_xarray()  # single column index

        # now pivot df and convert back to xarray using only non site_columns
        site_columns.remove('x')  # need to keep x to merge later
        dfx = df.loc[:, df.columns[~df.columns.isin(site_columns)]].rename({'siteid': 'x'}, axis=1).set_index(['time', 'x']).to_xarray()

        # merge the time depenedent and time independent
        out = xr.merge([dfx, dfps])

        # reset x index and add siteid back to the xarray object
        if ~pd.api.types.is_numeric_dtype(out.x):
            siteid = out.x.values
            out['x'] = range(len(siteid))
            out['siteid'] = (('x'), siteid)

        return out


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
        self.variable_dict = None

    def open_obs(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        from glob import glob
        from numpy import sort

        try:
            if os.path.isfile(self.file):
                _, extension = os.path.splitext(self.file)
                if extension in ['.nc', '.ncf', '.netcdf', '.nc4']:
                    if len(glob(self.file)) > 1:
                        self.obj = xr.open_mfdataset(sort(glob(self.file)))
                    self.obj = xr.open_dataset(self.file)
                elif extension in ['.ict', '.icarrt']:
                    self.obj = mio.icarrt.add_data(self.file)
                self.mask_and_scale()  # mask and scale values from the control values
                if 'epa_region' in self.obj.data_vars:
                    self.obj = self.obj.rename({'epa_region': 'REGION_LABEL'})
        except ValueError:
            print('something happened opening file')

    def mask_and_scale(self):
        """Mask and scale obs to convert units and set detection limits"""
        vars = self.obj.data_vars
        if self.variable_dict is not None:
            for v in vars:
                if v in self.variable_dict:
                    d = self.variable_dict[v]
                    if 'unit_scale' in d:
                        scale = d['unit_scale']
                    else:
                        scale = 1
                    if 'unit_scale_method' in d:
                        if d['unit_scale_method'] == '*':
                            self.obj[v].data *= scale
                        elif d['unit_scale_method'] == '/':
                            self.obj[v].data /= scale
                        elif d['unit_scale_method'] == '+':
                            self.obj[v].data += scale
                        elif d['unit_scale_method'] == '-':
                            self.obj[v].data += -1 * scale
                    if 'obs_limit' in d:
                        self.obj[v].data = self.obj[v].where(self.obj[v] >= d['obs_limit'])

    def obs_to_df(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        self.obj = self.obj.to_dataframe().reset_index().drop(['x', 'y'], axis=1)


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
        self.variable_dict = None
        self.figure_kwargs = None

    def glob_files(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        from numpy import sort
        from glob import glob

        print(self.file_str)
        self.files = sort(glob(self.file_str))

    def open_model_files(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        self.glob_files()
        if 'cmaq' in self.model.lower():
            if len(self.files) > 1:
                self.obj = mio.cmaq.open_mfdataset(self.files)
            else:
                self.obj = mio.cmaq.open_dataset(self.files[0])
        elif 'wrfchem' in self.model.lower():
            if len(self.files) > 1:
                self.obj = mio.wrfchem.open_mfdataset(self.files)
            else:
                self.obj = mio.wrfchem.open_dataset(self.files)
        elif 'rrfs' in self.model.lower():
            if len(self.files) > 1:
                self.obj = mio.rrfs_cmaq.open_mfdataset(self.files)
            else:
                self.obj = mio.rrfs_cmaq.open_dataset(self.files)
        elif 'gsdchem' in self.model.lower():
            if len(self.files) > 1:
                self.obj = mio.fv3chem.open_mfdataset(self.files)
            else:
                self.obj = mio.fv3chem.open_dataset(self.files)
        self.mask_and_scale()

    def mask_and_scale(self):
        """Mask and scale obs to convert units and set detection limits"""
        vars = self.obj.data_vars
        if self.variable_dict is not None:
            for v in vars:
                if v in self.variable_dict:
                    d = self.variable_dict[v]
                    if 'unit_scale' in d:
                        scale = d['unit_scale']
                    else:
                        scale = 1
                    if 'unit_scale_method' in d:
                        if d['unit_scale_method'] == '*':
                            self.obj[v].data *= scale
                        elif d['unit_scale_method'] == '/':
                            self.obj[v].data /= scale
                        elif d['unit_scale_method'] == '+':
                            self.obj[v].data += scale
                        elif d['unit_scale_method'] == '-':
                            self.obj[v].data += -1 * scale


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
        self.start_time = None
        self.end_time = None

    def read_control(self, control=None):
        """Reads the yaml control file.  If not set assumes control file is control.yaml

        Parameters
        ----------
        control : type
            Description of parameter `control`.

        Returns
        -------
        type
            Description of returned object.

        """
        import yaml

        if control is not None:
            self.control = control

        with open(self.control, 'r') as stream:
            self.control_dict = yaml.safe_load(stream)

        # set analysis time
        self.start_time = pd.Timestamp(self.control_dict['analysis']['start_time'])
        self.end_time = pd.Timestamp(self.control_dict['analysis']['end_time'])

    def open_models(self):
        """Opens all models and creates model instances for monet-analysis"""
        if 'model' in self.control_dict:
            # open each model
            for mod in self.control_dict['model']:
                # create a new model instance
                m = model()
                # this is the model type (ie cmaq, rapchem, gsdchem etc)
                m.model = self.control_dict['model'][mod]['mod_type']
                # set the model label in the dictionary and model class intance
                m.label = mod
                # create file string (note this can include hot strings)
                m.file_str = self.control_dict['model'][mod]['files']
                # create mapping
                m.mapping = self.control_dict['model'][mod]['mapping']
                # add variable dict
                print(mod)
                print(self.control_dict['model'][mod])
                if 'variables' in self.control_dict['model'][mod].keys():
                    m.variable_dict = self.control_dict['model'][mod]['variables']
                if 'figure_kwargs' in self.control_dict['model'][mod].keys():
                    m.figure_kwargs = self.control_dict['model'][mod]['figure_kwargs']
                # open the model
                m.open_model_files()
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
                o.file = self.control_dict['obs'][obs]['filename']
                o.open_obs()
                if 'variables' in self.control_dict['obs'][obs]:
                    self.variable_dict = self.control_dict['obs'][obs]['variables']
                self.obs[o.label] = o

    def pair_data(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        pairs = {}
        for model_label in self.models:
            mod = self.models[model_label]
            # Now we have the models we need to loop through the mapping table for each network and pair the data
            # each paired dataset will be output to a netcdf file with 'model_label_network.nc'
            for obs_to_pair in mod.mapping.keys():
                # get the variables to pair from the model data (ie don't pair all data)
                keys = [key for key in mod.mapping[obs_to_pair].keys()]
                obs_vars = [mod.mapping[obs_to_pair][key] for key in keys]

                model_obj = mod.obj[keys]
                ## TODO:  add in ability for simple addition of variables from

                # simplify the objs object with the correct mapping vairables
                obs = self.obs[obs_to_pair]

                # pair the data
                # if pt_sfc (surface point network or monitor)
                if obs.obs_type.lower() == 'pt_sfc':
                    # convert this to pandas dataframe
                    obs.obs_to_df()
                    # now combine obs with
                    paired_data = model_obj.monet.combine_point(obs.obj, radius_of_influence=1e6, suffix=mod.label)
                    # print(paired_data)
                    # this outputs as a pandas dataframe.  Convert this to xarray obj
                    p = pair()
                    p.obs = obs.label
                    p.model = mod.label
                    p.model_obj = mod
                    p.model_vars = keys
                    p.obs_vars = obs_vars
                    p.filename = '{}_{}.nc'.format(p.obs, p.model)
                    p.obj = paired_data.monet._df_to_da()
                    label = "{}_{}".format(p.obs, p.model)
                    self.paired[label] = p
                    p.obj = p.fix_paired_xarray(dset=p.obj)
                    # write_util.write_ncf(p.obj,p.filename) # write out to file
                # TODO: add other network types / data types where (ie flight, satellite etc)

    ### TODO: Create the plotting driver (most complicated one)
    # def plotting(self):
    def plotting(self):
        """This function will cycle through all the plots and control variables needed to make the plots necessary

        Returns
        -------
        type
            Description of returned object.

        """
        from plots import surfplots as splots

        # first get the plotting dictionary from the yaml file
        plot_dict = self.control_dict['plotting']

        # now we are going to loop through each plot_group (note we can have multiple plot groups)
        # a plot group can have
        #     1) a singular plot type
        #     2) multiple paired datasets or model datasets depending on the plot type
        #     3) kwargs for creating the figure ie size and marker (note the default for obs is 'x')
        for grp in plot_dict.keys():
            grp_dict = plot_dict[grp]  # this is the plot group

            pair_labels = grp_dict['data']
            # get the plot type
            plot_type = grp_dict['type']

            current_obsvar = None

            # first get the observational obs labels
            pair1 = self.paired[list(an.paired.keys())[0]]
            obs_vars = pair1.obs_vars

            # loop through obs variables
            for obsvar in obs_vars:
                for p_index, p in enumerate(pair_labels):
                    # find the pair model label that matches the obs var
                    index = p.obs_vars.index(obsvar)
                    modvar = p.model_vars[index]

                    if plot_type.lower() == 'timeseries':
                        pairdf = p.obj.to_dataframe().reset_index(drop=True).dropna(subset=[modvar])
                        if p_index == 0:
                            ax = splots.timeseries(pairdf, column=obsvar, label=p.obs)
                            if p.model_obj.figure_kwargs is not None:
                                ax = splots.timeseries(pairdf, column=modvar, label=p.model, ax=ax, **p.model_obj.figure_kwargs)
                            else:
                                ax = splots.timeseries(pairdf, column=modvar, label=p.model, ax=ax)
                        else:
                            if p.model_obj.figure_kwargs is not None:
                                ax = splots.timeseries(pairdf, column=modvar, label=p.model, ax=ax, **p.model_obj.figure_kwargs)
                            else:
                                ax = splots.timeseries(pairdf, column=modvar, label=p.model, ax=ax)
