# SPDX-License-Identifier: Apache-2.0
#
"""
This file calculates a number of meteorological variables. For hydrometeorological variables, 
variables including dewpoint and relative humidity are not universally provided by the obs. 
For models, wind speed is usually not directly provided. 

This script calculates all of them in one place. 

Users can calculate variables directly from the observations by using the extra_calc dictionary. Potential 
temperature is the only variable with this supported capability. 

See: control_aircraft_looping_AEROMMA_UFSAQM.yaml

Author: Liam Thompson
"""

try:
    import metpy
    from metpy.units import units
except ImportError:
    metpy = None
    units = None

# if future iterations want to calc dewpoint/relh on the observations, this can be done with other metpy
# libraries.

# calc dewpoint 
def dewpoint(obj, varmap = None, output_key = "dewpoint"):
    """Calculates dewpoint in K
    
    Parameters
    ----------
    obj : xarray dataset
        Model/obs data
    varmap : dictionary
        Dictionary defining the data column names to use
        Please provide variable names for "pres_calc" in Pa and "specific_hum" in kg/kg
    output_key : string
        String to name the new dewpoint ouput
    Returns
    -------
    xarray dataset
        Xarray dataset with applied calculation as a new data array
        
    """
    
    if metpy is None:
        raise ImportError(
            "metpy is required for extra_calc. "
            "Install with: conda install -c conda-forge metpy "
            "or with: pip install metpy"
        )
    
    # grab variable names from the yaml and error if not provided
    try:
        pressure = obj[varmap['pres_calc']]
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'pres_calc' in yaml file under 'dewpoint' under 'extra_calc'") from e
    try:
        specific_hum = obj[varmap['specific_hum']]
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'specific_hum' in yaml file under 'dewpoint' under 'extra_calc'") from e

    dpt = (metpy.calc.dewpoint_from_specific_humidity(
        pressure * units.Pa,
        specific_hum * units("kg/kg")
    )).metpy.convert_units("K").metpy.dequantify()
    #Fancy dequantify method puts back in standard xarray format

    obj[output_key] = dpt
        
    return obj

        
# calc relative humidity
def relh(obj, varmap=None, output_key="rel_hum"):
    """Calculates Relative Humidity in %
    
    Parameters
    ----------
    obj : xarray dataset
        Model/obs data
    varmap : dictionary
        Dictionary defining the data column names to use
        Please provide variable names for "pres_calc" in Pa, "temp_calc" in K, "specific_hum" in kg/kg
    output_key : string
        String to name the new relative humidity ouput
    Returns
    -------
    xarray dataset
        Xarray dataset with applied calculation as a new data array
        
    """
    
    if metpy is None:
        raise ImportError(
            "metpy is required for extra_calc. "
            "Install with: conda install -c conda-forge metpy "
            "or with: pip install metpy"
        )

    # grab variable names from the yaml and error if not provided
    try:
        pressure = obj[varmap['pres_calc']]
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'pres_calc' in yaml file under 'rel_hum' under 'extra_calc'") from e
    try:
        specific_hum = obj[varmap['specific_hum']]
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'specific_hum' in yaml file under 'rel_hum' under 'extra_calc'") from e
    try:
        temperature = obj[varmap["temp_calc"]]
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'temp_calc' in yaml file under 'rel_hum' under 'extra_calc'") from e

    rlh = (metpy.calc.relative_humidity_from_specific_humidity(
        pressure * units.Pa,
        temperature * units.kelvin, 
        specific_hum * units("kg/kg"),
        phase = 'auto' #Needed for cold temperatures in stratosphere
    )).metpy.convert_units("%").metpy.dequantify()
    #Fancy dequantify method puts back in standard xarray format

    # Set any value above 100% to NaN
    # I think this happens because the divisions by 0 
    # (e.g. very tinyyyyyy numbers once you get high enough in the atmos) 
    # are handled  weirdly in the metpy library. 
    rlh = rlh.where(rlh <= 100)

    obj[output_key] = rlh
        
    return obj
        
# calc windspeed
def wspd(obj, varmap = None, output_key = "windspeed"):
    """Calculates windspeed in m/s
    
    Parameters
    ----------
    obj : xarray dataset
        Model/obs data
    varmap : dictionary
        Dictionary defining the data column names to use
        Please provide variable names for "u_comp" in m/s, "v_comp" in m/s
    output_key : string
        String to name the new windspeed ouput
    Returns
    -------
    xarray dataset
        Xarray dataset with applied calculation as a new data array
        
    """
    
    if metpy is None:
        raise ImportError(
            "metpy is required for extra_calc. "
            "Install with: conda install -c conda-forge metpy "
            "or with: pip install metpy"
        )

    # grab variable names from the yaml and error if not provided
    try:
        u = obj[varmap["u_comp"]]
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'u_comp' in yaml file under 'windspeed' under 'extra_calc'") from e
    try:
        v = obj[varmap["v_comp"]]
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'v_comp' in yaml file under 'windspeed' under 'extra_calc'") from e

    wspd = (metpy.calc.wind_speed(
        u * units("m/s"),
        v * units("m/s")
        )).metpy.convert_units("m/s").metpy.dequantify()
    #Fancy dequantify method puts back in standard xarray format

    obj[output_key] = wspd
        
    return obj

# calc wind direction
def wdir(obj, varmap = None, output_key = "winddir"):
    """Calculates wind direction in degrees
    
    Parameters
    ----------
    obj : xarray dataset
        Model/obs data
    varmap : dictionary
        Dictionary defining the data column names to use
        Please provide variable names for "u_comp" in m/s, "v_comp" in m/s
    output_key : string
        String to name the new wind direction ouput
    Returns
    -------
    xarray dataset
        Xarray dataset with applied calculation as a new data array
        
    """

    if metpy is None:
        raise ImportError(
            "metpy is required for extra_calc. "
            "Install with: conda install -c conda-forge metpy "
            "or with: pip install metpy"
        )

    # grab variable names from the yaml and error if not provided
    try:
        u = obj[varmap["u_comp"]].compute()
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'u_comp' in yaml file under 'winddir' under 'extra_calc'") from e
    try:
        v = obj[varmap["v_comp"]].compute()
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'v_comp' in yaml file under 'winddir' under 'extra_calc'") from e
    #Unfortunately, wind_direction does not work with dask in metpy.
    #Maybe find another solution as we move to optimize memory usage, 
    #but for now just compute these.

    wdir = (metpy.calc.wind_direction(
        u * units("m/s"),
        v * units("m/s")
        )).metpy.convert_units("degree").metpy.dequantify()
    #Fancy dequantify method puts back in standard xarray format

    obj[output_key] = wdir
        
    return obj

# calc potential temperature
def ptemp(obj, varmap=None, output_key="ptemp"):
    """Calculates potential temperature in K
    
    Parameters
    ----------
    obj : xarray dataset
        Model/obs data
    varmap : dictionary
        Dictionary defining the data column names to use
        Please provide variable names for "pres_calc" in Pa, "temp_calc" in K
    output_key : string
        String to name the new potential temperature ouput
    Returns
    -------
    xarray dataset
        Xarray dataset with applied calculation as a new data array
        
    """

    if metpy is None:
        raise ImportError(
            "metpy is required for extra_calc. "
            "Install with: conda install -c conda-forge metpy "
            "or with: pip install metpy"
        )

    # grab variable names from the yaml and error if not provided
    try:
        pres = obj[varmap['pres_calc']]
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'pres_calc' in yaml file under 'ptemp_mod' or 'ptemp_obs' under 'extra_calc'") from e
    try:
        temp = obj[varmap["temp_calc"]]
    except (KeyError,TypeError) as e:
        raise Exception("Please specify 'temp_calc' in yaml file under 'ptemp_mod' or 'ptemp_obs' under 'extra_calc'") from e
    
    ptmp = (metpy.calc.potential_temperature(
        pres * units.Pa,
        temp * units("K")
    )).metpy.convert_units("K").metpy.dequantify()
    #Fancy dequantify method puts back in standard xarray format

    obj[output_key] = ptmp
        
    return obj