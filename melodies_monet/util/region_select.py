# SPDX-License-Identifier: Apache-2.0
#

"""
Masking arbitrary regions with regionmask
"""

import warnings
from functools import lru_cache

import pandas as pd

from melodies_monet.util.tools import get_epa_region_bounds, get_giorgi_region_bounds

try:
    import geopandas as gpd
    import regionmask
    from shapely.geometry import MultiPolygon, Polygon
except ImportError:
    regionmask = None


@lru_cache(None)
def get_regions(url_or_path, **kwargs):
    return regionmask.from_geopandas(gpd.read_file(url_or_path), **kwargs)


def _create_custom_mask(data, mask_info):
    """Creates mask using regionmask.

    Parameters
    ----------
    data : xr.Dataset
        Dataset with the data that the user wishes to mask.
    mask_info : list[list[float, float]] | dict[str, list[list[float, float]]]
        List containing the vertices of the polygon to mask.

    Returns
    -------
    xr.Dataset
        masked data
    """
    if not isinstance(mask_info, (list, dict)):
        raise TypeError(f"mask_info type={type(mask_info)}, not valid for create_custom_mask.")
    if isinstance(mask_info, list):
        if isinstance(mask_info[0][0], (float, int)):
            poly = Polygon(mask_info)
        else:
            polys = []
            for p in mask_info:
                polys.append(Polygon(p))
            poly = MultiPolygon(polys)
        regions = regionmask.Regions([poly])
    else:
        all_regions = []
        abbrevs = []
        for rname, rpolygon in mask_info.items():
            abbrevs.append(rname)
            all_regions.append(rpolygon)
        regions = regionmask.Regions(all_regions, abbrevs=abbrevs, name="custom_regions")
    # Regionmask requires "lat" and "lon"
    region_mask = regions.mask(data.rename({"latitude": "lat", "longitude": "lon"}))
    # But MM requires "latitude" and "longitude"
    if "lat" in region_mask.coords:
        region_mask = region_mask.rename({"lat": "latitude", "lon": "longitude"})
    masked_data = data.where(region_mask.notnull())
    return masked_data, regions.bounds_global.tolist()


def _create_predefined_mask(data, name_regiontype, region=None):
    """Creates mask using regionmask.

    Parameters
    ----------
    data : xr.Dataset
        Dataset with the data that the user wishes to mask.
    name_regiontype : str
        name of regionmask defined regions (e.g., "srex", "giorgi")
    region : str
        region to mask. If it can be parsed to an integer, that is used.
        Otherwise, it is assumed to be an abbrev.

    Returns
    -------
    xr.Dataset, list[float]
        mask for input data, and bounding box
    """
    name_regiontype_split = name_regiontype.split(".")
    all_regions = regionmask.defined_regions
    for r in name_regiontype_split:
        all_regions = getattr(all_regions, r)
    # Regionmask requires "lat" and "lon"
    region_mask = all_regions.mask(data.rename({"latitude": "lat", "longitude": "lon"}))
    # But MM requires "latitude" and "longitude"
    if "lat" in region_mask.coords:
        region_mask = region_mask.rename({"lat": "latitude", "lon": "longitude"})
    try:
        selected_region = data.where(region_mask == int(region))
    except ValueError:
        selected_region = data.where(region_mask.cf == region)
    return selected_region, list(all_regions[region].bounds)


def _create_shapefile_mask(data, mask_path=None, mask_url=None, region_name=None, **kwargs):
    """Creates mask from shapefile using regionmask and geopandas.

    Parameters
    ----------
    data : xr.Dataset
        Dataset with the data that the user wishes to mask.
    mask_path : str | pathobject
        Path to the shapefiles. They can be zipped.
    mask_url : str
        url to download the shapefiles. Requires pooch.
        If both mask_url and mask_path are provided, this will be
        ignored.
    region_name : str
        region to mask. If it can be parsed to an integer, that is used.
        Otherwise, it is assumed to be an abbrev.
    **kwargs
        Arguments to pass to regionmask.from_geopandas

    Returns
    -------
    xr.Dataset, list[float]
        mask for the input data, and bounding box
    """

    if mask_url is not None and mask_path is not None:
        warnings.warn(
            "mask_url and mask_path provided. Only one can be used. "
            "Selecting mask_path and discarding URL."
        )

    if mask_path is not None:
        url_or_path = mask_path
    elif mask_url is not None:
        url_or_path = mask_url
    else:
        raise ValueError("Either mask_path or mask_url have to be provided")

    regions = get_regions(url_or_path, **kwargs)

    # Regionmask requires "lat" and "lon"
    region_mask = regions.mask(data.rename({"latitude": "lat", "longitude": "lon"}))
    # But MM requires "latitude" and "longitude"
    if "lat" in region_mask.coords:
        region_mask = region_mask.rename({"lat": "latitude", "lon": "longitude"})

    try:
        selected_region = data.where(region_mask == int(region_name))
    except ValueError:
        selected_region = data.where(region_mask.cf == region_name)
    return selected_region, list(regions[region_name].bounds)


def control_custom_mask(data, domain_type, domain_info=None, **kwargs):
    """Parses region information to return the right type of data.

    Parameters
    ----------
    data : xr.Dataset
        data to be masked
    domain_type : str
        type of data. Used to decide which function to apply. Should begin with
        "custom".
    domain_info : dict
        Dictionary containing relevant information on domain, like url, name, etc.
    **kwargs:
        Extra kwargs to pass to regionmask

    Returns
    -------
    xr.Dataset, list[float]
        masked Dataset, list of bounds
    """
    if regionmask is None:
        raise ImportError(
            "regionmask is not installed, try alternative functions."
            + " create_autoregion can probably do the trick."
        )
    if domain_info is None:
        raise KeyError("If regionmask is used, domain_info must exist.")
    if "custom" not in domain_type:
        raise ValueError("If regionmask is used, the domain_type should be starting with 'custom'")
    if "polygon" in domain_type:
        masked_data, bounds = _create_custom_mask(data, domain_info["mask_info"])
    elif "defined-region" in domain_type:
        name_regiontype = domain_info["name_regiontype"]
        region = domain_info["region"]
        masked_data, bounds = _create_predefined_mask(data, name_regiontype, region)
    elif "file" in domain_type:
        params = domain_info
        params["mask_path"] = domain_info.get("mask_path", None)
        params["mask_url"] = domain_info.get("mask_url", None)
        params["region_name"] = domain_info.get("region_name", None)
        params["abbrevs"] = domain_info.get("abbrevs", "_from_name")
        masked_data, bounds = _create_shapefile_mask(data, **params, **kwargs)
    else:
        raise ValueError(
            "Could not identify the type of domain. Should be 'polygon',"
            + f" 'defined-region' or 'file'. You asked for {domain_type!r}"
        )
    return masked_data, reorder_bounds_for_plotting(bounds)


def reorder_bounds_for_plotting(bounds):
    """Reorders from lonmin, latmin, lonmax, latmax (regionmask)
    to lonmin, lonmax, latmin, latmax (cartopy)

    Parameters
    ----------
    bounds : list[float]
        list in the format [lonmin, latmin, lonmax, latmax]

    Returns
    -------
    list[float]
        list in the format [lonmin, lonmax, latmin, latmax]
    """
    return [bounds[0], bounds[2], bounds[1], bounds[3]]


def create_autoregion(data, domain_type, domain_name, domain_info=None):
    """Selects a region using predefined boundaries.

    Parameters
    ----------
    data : xr.Dataset | pd.DataFrame
        data to be masked
    domain_type : str
        type of data. Used to decide which function to apply.
        If domain_type == 'auto-region:', domain_info is required.
    domain_name : str
        This is used as the region name, or to read the info.
    domain_info: None | dict[str, tuple[float, float, float, float]]
        if not None, dict containing the domain name and a tuple with
        lonmin, lonmax, latmin, latmax. Only required if domain_type
        is auto-region:

    Returns
    -------
    xr.Dataset | pd.DataFrame
        Data as it was provided to the function
    """
    auto_region_id = domain_type.split(":")[1].lower()
    if auto_region_id == "epa":
        bounds = get_epa_region_bounds(acronym=domain_name)
    elif auto_region_id == "giorgi":
        bounds = get_giorgi_region_bounds(acronym=domain_name)
    elif auto_region_id == "box":
        bounds = domain_info["bounds"]
    else:
        raise ValueError(
            "Currently, auto-region selections without a domain query have only "
            "been implemented for Giorgi and EPA regions. You asked for "
            f"{domain_type!r}. If you need more capabilities, check out the custom: "
            "regions capabilities. Be aware that they require regionmask."
        )
    if isinstance(data, pd.DataFrame):
        data_all = data.loc[
            (data["longitude"] >= bounds[0])
            & (data["longitude"] <= bounds[1])
            & (data["latitude"] >= bounds[2])
            & (data["latitude"] <= bounds[3])
        ]
    else:
        data_all = data.where(
            (data["longitude"] >= bounds[0])
            & (data["longitude"] <= bounds[1])
            & (data["latitude"] >= bounds[2])
            & (data["latitude"] <= bounds[3]),
        )
    return data_all, bounds


def select_region(data, domain_type, domain_name, domain_info=None, **kwargs):
    """Selects a region in whichever format it was provided

    Parameters
    ----------
    data : xr.Dataset | pd.DataFrame
        data to be masked
    domain_type : str
        type of data. Used to decide which function to apply.
    domain_name : str
        This is used as the region name, or to read the info.
    domain_info : dict
        Dict containing the domain_name and other relevant information, e. g.,
        lonlat box, mask_url, mask_file, etc.
    **kwargs:
        extra kwargs to pass to the selector, depending on the type of
        data or region.

    Returns
    -------
    xr.Dataset | pd.DataFrame
        Region selected. The type will be the same as the input data.
    """

    if domain_type == "all":
        return data, [-180, 180, -90, 90]
    if domain_type.startswith("auto-region") or (domain_type == "custom:box"):
        data_masked, bounds = create_autoregion(data, domain_type, domain_name, domain_info)
    elif domain_type.startswith("custom"):
        data_masked, bounds = control_custom_mask(data, domain_type, domain_info, **kwargs)
    else:
        if isinstance(data, pd.DataFrame):
            data_masked = data.query(domain_type + " == " + '"' + domain_name + '"')
            lon = data_masked["longitude"]
            lat = data_masked["latitude"]
            bounds = [lon.min(), lon.max(), lat.min(), lat.max()]
        else:
            data_masked = data.where(data[domain_type] == domain_name)
            lon = data_masked["longitude"].where(data_masked[domain_type]==domain_name)
            lat = data_masked["latitude"].where(data_masked[domain_type]==domain_name)
            bounds = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]
    return data_masked, bounds
