""" This is the overall control file.  It will drive the entire analysis package"""
import monetio as mio
import monet as m
import os
import xarray as xr
import pandas as pd
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
        self.filename = None

    def fix_paired_xarray(self,dset=None):

        # first convert to dataframe
        df = dset.to_dataframe().reset_index(drop=True)

        # now get just the single site index
        dfpsite = df.rename({'siteid':'x'},axis=1).drop_duplicates(subset=['x'])
        columns = dfpsite.columns # all columns
        site_columns = ['latitude','longitude','x','site','msa_code','cmsa_name','epa_region','state_name','msa_name','site','utcoffset'] # only columns for single site identificaiton

        # site only xarray obj (no time dependence)
        dfps = dfpsite.loc[:,columns[columns.isin(site_columns)]].set_index(['x']).to_xarray() # single column index

        # now pivot df and convert back to xarray using only non site_columns
        site_columns.remove('x') # need to keep x to merge later
        dfx = df.loc[:,df.columns[~df.columns.isin(site_columns)]].rename({'siteid':'x'},axis=1).set_index(['time','x']).to_xarray()

        # merge the time depenedent and time independent
        out = xr.merge([dfx,dfps])

        # reset x index and add siteid back to the xarray object
        if ~pd.api.types.is_numeric_dtype(out.x):
            siteid = out.x.values
            out['x'] = range(len(siteid))
            out['siteid'] = (('x'),siteid)

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
                _,extension = os.path.splitext(self.file)
                if extension in ['.nc','.ncf','.netcdf','.nc4']:
                    if len(glob(self.file)) > 1:
                        self.obj = xr.open_mfdataset(sort(glob(self.file)))
                    self.obj = xr.open_dataset(self.file)
                elif extension in ['.ict','.icarrt']:
                    self.obj = mio.icarrt.add_data(self.file)
        except ValueError:
            print('something happened opening file')

    def obs_to_df(self):
        """Short summary.

        Returns
        -------
        type
            Description of returned object.

        """
        self.obj = self.obj.to_dataframe().reset_index().drop(['x','y'],axis=1)

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
        self.end_time = pd.Timestamp(self.control_dict['analysis']['start_time'])

    def open_models(self):
        """Opens all models and creates model instances for monet-analysis

        """
        if 'model' in self.control_dict:
            # open each model
            for mod in self.control_dict['model']:
                # create a new model instance
                m = model()
                # this is the model type (ie cmaq, rapchem, gsdchem etc)
                m.model = self.control_dict['model'][mod]['mod_type']
                #set the model label in the dictionary and model class intance
                m.label = mod
                # create file string (note this can include hot strings)
                m.file_str = self.control_dict['model'][mod]['files']
                #create mapping
                m.mapping = self.control_dict['model'][mod]['mapping']
                #open the model
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
                    paired_data = model_obj.monet.combine_point(obs.obj,radius_of_influence=1e6,suffix=mod.label)
                    print(paired_data)
                    # this outputs as a pandas dataframe.  Convert this to xarray obj
                    p = pair()
                    p.obs = obs.label
                    p.model = mod.label
                    p.filename = '{}_{}.nc'.format(p.obs,p.model)
                    p.obj = paired_data.monet._df_to_da()
                    label = "{}_{}".format(p.obs,p.model)
                    self.paired[label] = p
                    p.obj = p.fix_paired_xarray(dset=p.obj)
                    # write_util.write_ncf(p.obj,p.filename) # write out to file
                # TODO: add other network types / data types where (ie flight, satellite etc)

    ### TODO: Create the plotting driver (most complicated one)
    #def plotting(self):
    def plotting(self):
        """This function will cycle through all the plots and control variables needed to make the plots necessary

        Returns
        -------
        type
            Description of returned object.

        """
        # Make the plots for each paired dataset
        for paired_label in self.paired:
            paired = self.paired[paired_label]
            mapping_table = self.models[paired.model].mapping[paired.obs]

            # check what type of observation this is (pt_sfc.... etc)
            obs_type = self.control_dict['obs'][paired.obs]['obs_type'].lower()

            #get the plots kwarg from the obs
            plots_yaml = self.control_dict['obs'][paired.obs]['plots']
            # get the basename of the plots string
            basename = plots_yaml['basename']

            # first check if the plotting type is a point surface observation
            if obs_type == 'pt_sfc':
                # TODO: add new plot types here and below
                known_plot_types = ['taylor','spatial_bias','timeseries','spatial_overlay']
                # now we want to loop through each plot type
                plot_types = [i.lower() for i in self.control_dict['obs'][paired.obs]['plots'].keys()]
                # only loop over plot types that are in the pt_sfc
                good_to_go = plot_types[[i in known_plot_types for i in plot_types]]

                # loop over good_to_go plot types:
                for plot_type in good_to_go:
                    plot_type_dict = plot_types[plot_type]

                    # first do domain plots
                    if plot_type == 'talyor':
                        _taylor_plot(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=False)
                    if plot_type == 'spatial_bias':
                        _spatial_bias(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=False)
                    if plot_type == 'timeseries':
                        _timeseries(paried, plot_dict=plot_type_dict, region=None, epa_regulatory=False)
                    if plot_type == 'spatial_overlay':
                        _spatial_overlay(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=False)
                    # TODO: Add additional plot types here

                    if 'regulatory' in plot_type_dict:
                        if plot_type_dict['regulatory']:
                            if plot_type == 'talyor':
                                _taylor_plot(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=True)
                            if plot_type == 'spatial_bias':
                                _spatial_bias(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=True)
                            if plot_type == 'timeseries':
                                _timeseries(paried, plot_dict=plot_type_dict, region=None, epa_regulatory=True)
                            if plot_type == 'spatial_overlay':
                                _spatial_overlay(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=True)




            if 'epa_regions' in self.control_dict['obs'][paired.obs]['plots']:

            subset = self.control_dict['obs'][paired.obs]['plots']['epa_regions']
            reg = self.control_dict['obs'][paired.obs]['plots']['regulatory']
            df = paired.obj.to_dataframe()
            regions = self.control_dict['obs'][paired.obs]['plots']['epa_regions_names']
            out_name = paired_label + '_' + self.control_dict['obs'][paired.obs]['plots']['plots_basename']
            #For each model species in the mapping table make plots.
            # if pt_sfc (surface point network or monitor)
            if self.control_dict['obs'][paired.obs]['obs_type'].lower() == 'pt_sfc':
                if self.control_dict['obs'][paired.obs]['plots']['taylor_diagram'] == True:
                    scale_ty = self.control_dict['obs'][paired.obs]['plots']['taylor_diagram_scale']
                    import taylor_plots as taylor
                    if subset == True:
                        for region in regions:
                            taylor.make_taylor_plots(df,out_name,subset,self.start_time,self.end_time,reg,scale_ty,mapping_table,region)
                    else:
                        taylor.make_taylor_plots(df,out_name,subset,self.start_time,self.end_time,reg,scale_ty,mapping_table)
    def _taylor_plot(paired, plot_dict=None):
        # TODO: Create wrapper function for the taylor plots
        a = 1

    def _spatial_bias(paired):
        # TODO: create wrapper for spatial bias
        a = 1

    def _timeseries(paired):
        # TODO: create wrapper for timeseries
        a = 1

    def _spatial_overlay(paired):
        # TODO: write wrapper for overlay plots
        a = 1
