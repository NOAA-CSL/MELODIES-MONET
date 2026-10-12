# https://github.com/NCAR/MELODIES-MONET/issues/175
# last update mentioned moving this to the tools section. That file seems really long already. 
# could be easier to just have a file dedicated to regulatory specific functions 

# Regulatory calculations 
# moved from surf plots

def make_24hr_regulatory(df, col=None):
    """Calculates 24-hour averages

    Parameters
    ----------
    df : dataframe
        Model/obs pair of hourly data
    col : str
        Column label of observation variable to apply the calculation
    Returns
    -------
    dataframe
        dataframe with applied calculation

    """
    #return calc_24hr_ave(df, col)
    return calc_24hr_ave_v1(df, col)


def calc_24hr_ave_v1(df, col=None):
    df.index = df.time_local
    # select sites with nobs >=18, 75% completeness
    df_24hr_ave = (df.groupby("siteid")[col].resample("D").sum(min_count=18, numeric_only=True)/df.groupby("siteid")[col].resample("D").count()).reset_index().dropna()
    df = df.reset_index(drop=True)
    return df.merge(df_24hr_ave, on=["siteid", "time_local"])


def make_8hr_regulatory(df, col=None):
    """Calculates 8-hour rolling average daily

    Parameters
    ----------
    df : dataframe
        Model/obs pair of hourly data
    col : str
        Column label of observation variable to apply the calculation
    Returns
    -------
    dataframe
        dataframe with applied calculation

    """
    #return calc_8hr_rolling_max(df, col, window=8)
    return calc_8hr_rolling_max_v1(df, col, window=8)


def calc_8hr_rolling_max_v1(df, col=None, window=None):
    df.index = df.time_local
    df_rolling = df.groupby("siteid")[col].rolling(window,min_periods=6,center=True, win_type="boxcar").mean(numeric_only=True).reset_index().dropna()
    # JianHe: select sites with nobs >=18, 75% completeness based on EPA
    df_rolling.index = df_rolling.time_local
    df_rolling_max = df_rolling.groupby("siteid").resample("D").max(min_count=18, numeric_only=True).reset_index().dropna()
    df = df.reset_index(drop=True)
    return df.merge(df_rolling_max, on=["siteid", "time_local"])


def get_utcoffset(lat,lon):
    """get UTC offset in hour based on a point (lat/lon)

    Parameters
    ----------
    lat :
        Latitude (deg; -90. to 90.)
    lon :
        Longitude (deg; -180. to 180.)

    Returns
    -------
    UTC offset in hour

    """
    import datetime
    import pytz
    from timezonefinder import TimezoneFinder

    tf = TimezoneFinder()

    timezone_str = tf.timezone_at(lng=lon, lat=lat)

    if timezone_str is None:
        #print('None timezone: ', lat, lon)
        if lon > -100.0:
            timezone_str = 'America/New_York'
        else:
            timezone_str = 'America/Los_Angeles'

        tz = pytz.timezone(timezone_str)
        d=datetime.utcnow()
        uos = tz.utcoffset(d, is_dst=False)
        utchour = uos.seconds/60.0/60.0
        utcday = uos.days

    elif timezone_str.startswith({'Etc','GMT'}):
        #print('Ocean timezone: ', timezone_str)
        tz = pytz.timezone(timezone_str)
        d=datetime.utcnow()
        uos = tz.utcoffset(d, is_dst=False)
        utchour = uos.seconds/60.0/60.0
        utcday = uos.days

    else:
        #print('Land timezone: ', timezone_str)
        tz = pytz.timezone(timezone_str)
        d=datetime.utcnow()
        uos = tz.utcoffset(d, is_dst=True)
        utchour = uos.seconds/60.0/60.0
        utcday = uos.days

    if utcday < 0:
       utchour = (24-utchour)*(-1) # Local - UTC

    return utchour
    