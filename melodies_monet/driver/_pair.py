# SPDX-License-Identifier: Apache-2.0
#
import xarray as xr
import pandas as pd


class pair:
    """The pair class.

    The pair class pairs model data
    directly with observational data along time and space.
    """

    def __init__(self):
        """Initialize a :class:`pair` object."""
        self.type = "pt_sfc"
        self.radius_of_influence = 1e6
        self.obs = None
        self.model = None
        self.model_vars = None
        self.obs_vars = None
        self.filename = None

    def __repr__(self):
        return (
            f"{type(self).__name__}(\n"
            f"    type={self.type!r},\n"
            f"    radius_of_influence={self.radius_of_influence!r},\n"
            f"    obs={self.obs!r},\n"
            f"    model={self.model!r},\n"
            f"    model_vars={self.model_vars!r},\n"
            f"    obs_vars={self.obs_vars!r},\n"
            f"    filename={self.filename!r},\n"
            ")"
        )

    def fix_paired_xarray(self, dset):
        """Reformat the paired dataset.

        Parameters
        ----------
        dset : xarray.Dataset

        Returns
        -------
        xarray.Dataset
            Reformatted paired dataset.
        """
        # first convert to dataframe
        df = dset.to_dataframe().reset_index(drop=True)

        # now get just the single site index
        dfpsite = df.rename({"siteid": "x"}, axis=1).drop_duplicates(subset=["x"])
        columns = dfpsite.columns  # all columns
        site_columns = [
            "latitude",
            "longitude",
            "x",
            "site",
            "msa_code",
            "cmsa_name",
            "epa_region",
            "state_name",
            "msa_name",
            "site",
            "utcoffset",
        ]  # only columns for single site identificaiton

        # site only xarray obj (no time dependence)
        dfps = (
            dfpsite.loc[:, columns[columns.isin(site_columns)]].set_index(["x"]).to_xarray()
        )  # single column index

        # now pivot df and convert back to xarray using only non site_columns
        site_columns.remove("x")  # need to keep x to merge later
        dfx = (
            df.loc[:, df.columns[~df.columns.isin(site_columns)]]
            .rename({"siteid": "x"}, axis=1)
            .set_index(["time", "x"])
            .to_xarray()
        )

        # merge the time dependent and time independent
        out = xr.merge([dfx, dfps])

        # reset x index and add siteid back to the xarray object
        if ~pd.api.types.is_numeric_dtype(out.x):
            siteid = out.x.values
            out["x"] = range(len(siteid))
            out["siteid"] = (("x"), siteid)

        return out
