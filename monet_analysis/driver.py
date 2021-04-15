""" This is the overall control file.  It will drive the entire analysis package"""
import monetio as mio
import monet as m
import os
import xarray as xr
import pandas as pd
import numpy as np
import datetime
from monet.util.tools import get_relhum
from plots import surfplots as splot
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
        self.end_time = pd.Timestamp(self.control_dict['analysis']['end_time'])

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
            sub_map = {mapping_table[i]: i for i in mapping_table} 
            #I reorder this as all the species in the plots are based on obs species as the key
            
            # check what type of observation this is (pt_sfc.... etc)
            obs_type = self.control_dict['obs'][paired.obs]['obs_type'].lower()

            #get the plots kwarg from the obs
            plots_yaml = self.control_dict['obs'][paired.obs]['plots']
            # get the basename of the plots string
            basename = plots_yaml['plots_basename']

            # first check if the plotting type is a point surface observation
            if obs_type == 'pt_sfc':
                # TODO: add new plot types here and below
                known_plot_types = ['taylor','spatial_bias','timeseries','spatial_overlay']
                # now we want to loop through each plot type
                plot_types = [i.lower() for i in self.control_dict['obs'][paired.obs]['plots']['plot_types'].keys()]
                # only loop over plot types that are in the pt_sfc
                good_to_go = list(set(known_plot_types) & set(plot_types))
                print('Will loop good_to_go')
                # loop over good_to_go plot types:
                for plot_type in good_to_go:
                    print(plot_type)
                    plot_type_dict = self.control_dict['obs'][paired.obs]['plots']['plot_types'][plot_type]
                    #Loop over species here and apply unit corrections for each observational type.
                    species = sub_map.keys()
                    #create the df
                    df = paired.obj.to_dataframe()
                    df_replace = df.replace(-1.0, np.nan)  # Replace all exact -1.0 values with nan, 
                    #need to confirm this with Barry, do -1.0 values always reflect nan's in MONET? 
                    #Or is this new to compressed format of obs? Maybe something to set for each observations?
                    for sp in species:
                        print(sp)
                        df_drop = df_replace.dropna(subset=[sp, sub_map.get(sp)])  # Drops all rows with obs species = NaN
                        #For each observations apply corrections for units here, but eventually should make units universal when bring in obs and model output.
                        if paired.obs == 'airnow':
                            if sp == 'WS':
                                df_drop.loc[:, 'WS'] = df_drop.loc[:, 'WS'] * 0.514  # convert obs knots-->m/s
                                df_drop.query('WS > 0.2', inplace=True
                                         )  # Filter out calm WS obs (< 0.2 m/s), should not be trusted--creates artificially larger postive  model bias
                            elif sp == 'BARPR':
                                df_drop.loc[:, 'PRSFC'] = df_drop.loc[:, 'PRSFC'] * 0.01  # convert model Pascals-->millibars
                            elif sp == 'PRECIP':
                                df_drop.loc[:, 'PRECIP'] = df_drop.loc[:, 'PRECIP'] * 0.1  # convert obs mm-->cm
                            elif sp == 'TEMP':
                                df_drop.loc[:, 'TEMP2'] = df_drop.loc[:, 'TEMP2'] - 273.16  # convert model K-->C
                            elif sp == 'RHUM':
                            # convert model mixing ratio to relative humidity
                                df_drop.loc[:, 'Q2'] = get_relhum(df_drop.loc[:, 'TEMP2'], df_drop.loc[:, 'PRSFC'],
                                                                  df_drop.loc[:, 'Q2'])
                            # df2.rename(index=str,columns={"Q2": "RH_mod"},inplace=True)
                            elif sp == 'CO':
                                df_drop.loc[:, 'CO'] = df_drop.loc[:, 'CO'] * 1000.0  # convert obs ppm-->ppb
#figure out the start dates and output name
                        #If applicable, select start and end time. 
                        if self.start_time != None and self.end_time != None:
                            dfnew = df_drop.loc[self.start_time:self.end_time]
                            startdatename = str(datetime.datetime.strftime(self.start_time, '%Y-%m-%d_%H'))
                            enddatename = str(datetime.datetime.strftime(self.end_time, '%Y-%m-%d_%H'))
                            outname = "{}.{}.{}.{}.{}".format(basename, paired_label, sp, startdatename, enddatename)
                            if sp == 'PM2.5':
                                outname = outname.replace('PM2.5', 'PM2P5')
                        else:
                            dfnew = df_drop
                            outname = "{}.{}.{}".format(basename, paired_label, sp)
                            if sp == 'PM2.5':
                                outname = outname.replace('PM2.5', 'PM2P5')
                        print(outname)
                        # first do domain plots
                        if plot_type == 'taylor':
                            print('entered taylor plt')
                            if self.start_time != None and self.end_time != None:
                                splot.make_taylor_plot(dfnew, sub_map.get(sp), sp , outname, paired.obs, paired.model,
                                             plot_dict=plot_type_dict, region=None, epa_regulatory=False, time_avg=True)
                            else: #If model analysis times are not provided plot each hour of the output rather than an average.
                                splot.make_taylor_plot(dfnew, sub_map.get(sp), sp , outname,
                                             plot_dict=plot_type_dict, region=None, epa_regulatory=False, time_avg=False)
                        if plot_type == 'spatial_bias':
                            splot.make_spatial_bias(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=False)
                        if plot_type == 'timeseries':
                            splot.make_timeseries(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=False)
                        if plot_type == 'spatial_overlay':
                            splot.make_spatial_overlay(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=False)
                    # TODO: Add additional plot types here

                    #if 'regulatory' in plot_type_dict:
                    #    # Converts OZONE, PM10, or PM2.5 dataframe to NAAQS regulatory values
                    #    if jj == 'OZONE' and reg is True:
                    #        df2 = make_8hr_regulatory(df_drop, [jj, sub_map.get(jj)]).rename(
                    #            index=str, columns={jj + '_y': jj, sub_map.get(jj) + '_y': sub_map.get(jj)}
                    #        )
                    #    elif jj == 'PM2.5' and reg is True:
                    #        df2 = make_24hr_regulatory(df_drop, [jj, sub_map.get(jj)]).rename(
                    #              index=str, columns={jj + '_y': jj, sub_map.get(jj) + '_y': sub_map.get(jj)}
                    #          )
                    #    elif jj == 'PM10' and reg is True:
                    #        df2 = make_24hr_regulatory(df_drop, [jj, sub_map.get(jj)]).rename(
                    #            index=str, columns={jj + '_y': jj, sub_map.get(jj) + '_y': sub_map.get(jj)}
                    #        )
                    #    else:
                    #        df2 = df_drop
                        #add reg to output name
                    #    outname = "{}.{}".format(outname, 'reg')
                    #    
                    #    if plot_type_dict['regulatory']:
                    #        if plot_type == 'talyor':
                    #            _taylor_plot(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=True)
                    #        if plot_type == 'spatial_bias':
                    #            _spatial_bias(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=True)
                    #        if plot_type == 'timeseries':
                    #            _timeseries(paried, plot_dict=plot_type_dict, region=None, epa_regulatory=True)
                    #       if plot_type == 'spatial_overlay':
                    #            _spatial_overlay(paired, plot_dict=plot_type_dict, region=None, epa_regulatory=True)




            #if 'epa_regions' in self.control_dict['obs'][paired.obs]['plots']:
            #    a = 1