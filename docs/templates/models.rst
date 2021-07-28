******
Models
******

monetio is capable of opening output from several different models including CMAQ,
HYSPLIT and GEFS-Aerosol.  monetio opens all N-Dimensional output in the :py:class:`~xarray.Dataset`.

CMAQ
----

CMAQ is a 3D photochemical grid model developed at the U.S. EPA to simulate air
composition.  monetio is able to read the output IOAPI output and format it to be
compatible with it's datastream.

As an example, lets open some CMAQ data from the Hawaiian volcanic eruption in 2018.
First we will set the path to the data files


.. code-block:: python

    import monetio as mio

    cmaqfile = monetio.__path__ + '/../data/aqm.t12z.aconc.ncf'

    c = mio.cmaq.open_dataset(cmaqfile)


This will return an :py:class:`~xarray.Dataset`.

.. code:: python

  print(c)
  <xarray.Dataset>
  Dimensions:    (DATE-TIME: 2, VAR: 41, time: 45, x: 80, y: 52, z: 1)
  Coordinates:
      latitude   (y, x) float32 ...
      longitude  (y, x) float32 ...
    * time       (time) datetime64[ns] 2018-05-17T12:00:00 2018-05-17T13:00:00 ...
  Dimensions without coordinates: DATE-TIME, VAR, x, y, z
  Data variables:
      TFLAG      (time, VAR, DATE-TIME) int32 dask.array<shape=(45, 41, 2), chunksize=(45, 41, 2)>
      O3         (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      NO2        (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      NO         (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      NO3        (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      N2O5       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      HNO3       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      HONO       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      PNA        (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      CO         (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      FORM       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ALD2       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      PAN        (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      NTR        (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      XO2N       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      SO2        (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ASO4I      (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ASO4J      (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ANH4I      (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ANH4J      (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ANO3I      (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ANO3J      (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      AORGAI     (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      AORGAJ     (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      AORGPAI    (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      AORGPAJ    (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      AORGBI     (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      AORGBJ     (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      AECI       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      AECJ       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      A25I       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      A25J       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      NUMATKN    (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      NUMACC     (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      SRFATKN    (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      SRFACC     (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      AH2OI      (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      AH2OJ      (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ACLI       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ACLJ       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ANAI       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
      ANAJ       (time, z, y, x) float32 dask.array<shape=(45, 1, 52, 80), chunksize=(45, 1, 52, 80)>
  Attributes:
      IOAPI_VERSION:  $Id: @(#) ioapi library version 3.1 $                    ...
      EXEC_ID:        ????????????????                                         ...
      FTYPE:          1
      CDATE:          2018142
      CTIME:          135716
      WDATE:          2018142
      WTIME:          135716
      SDATE:          2018137
      STIME:          120000
      TSTEP:          10000
      NTHIK:          1
      NCOLS:          80
      NROWS:          52
      NLAYS:          1
      NVARS:          41
      GDTYP:          2
      P_ALP:          19.0
      P_BET:          21.0
      P_GAM:          -157.5
      XCENT:          -157.5
      YCENT:          20.53
      XORIG:          -480000.0
      YORIG:          -312000.0
      XCELL:          12000.0
      YCELL:          12000.0
      VGTYP:          1
      VGTOP:          200.0
      VGLVLS:         [1.       0.089794]
      GDNAM:          AQF_HI
      UPNAM:          OPACONC
      VAR-LIST:       O3              NO2             NO              NO3      ...
      FILEDESC:       Concentration file output                                ...
      HISTORY:

All monetio xarray objects have common coordinate names (latitude and longitude) and dimension names (time, x, y, z).  It retains the
original attributes of the file and variable names.  monetio will precalculate some variables while loading the data in a lazy fashion, i.e. it
will not actually do the computation (not stored in memory) until needed:

.. code:: python

    pm25 = c.PM25

where pm25 is a :py:class:`~xarray.DataArray` as it is a single variable.

Prep-Chem-Sources
-----------------

A reader and writer was built into 


HYSPLIT
-----------------
An example with HYSPLIT run for 2019 eruption of Reventador Volcano.

.. code-block:: python

    import monetio as mio

    hysplitfile = monetio.__path__ + '/../data/cdump.bin

    hxr = mio.hysplit.open_dataset(hysplitfile)


This will return an :py:class:`~xarray.Dataset`.


.. code-block:: python

    print(hxr)

    <xarray.Dataset>
    Dimensions:    (time: 14, x: 47, y: 174, z: 5)
    Coordinates:
      * time       (time) datetime64[ns] 2019-02-25T17:00:00 ... 2019-02-26T06:00:00
      * x          (x) int64 8987 8988 8989 8990 8991 ... 9029 9030 9031 9032 9033
      * y          (y) int64 4332 4333 4334 4335 4336 ... 4501 4502 4503 4504 4505
      * z          (z) int64 2000 3000 4000 5000 6000
        longitude  (y, x) float64 -77.85 -77.84 -77.83 ... -77.41 -77.4 -77.39
        latitude   (y, x) float64 -1.774 -1.774 -1.774 ... -0.04456 -0.04456
    Data variables:
        p006       (time, z, y, x) float32 nan nan nan nan nan ... nan nan nan nan
        par2       (time, z, y, x) float32 nan nan nan nan nan ... nan nan nan nan
        par3       (time, z, y, x) float32 nan nan nan nan nan ... nan nan nan nan
        par4       (time, z, y, x) float32 nan nan nan nan nan ... nan nan nan nan
    Attributes:
        Starting Locations:           [(-0.077, -77.656)]
        Source Date:                  [datetime.datetime(2019, 2, 25, 16, 30)]
        Meteorological Model ID:      ERA5
        Number Start Locations:       2
        Number of Levels:             5
        Level top heights (m):        [2000 3000 4000 5000 6000]
        Number of Species:            4
        Sampling Time:                1:00:00
        sample time hours:            1.0
        Species ID:                   ['p006', 'par2', 'par3', 'par4']
        Concentration Grid:           {'Number Lat Points': 9001, 'Number Lon Poi...'}
        Coordinate time description:  Beginning of sampling time



If there are multiple species in the file

.. code-block:: python

    hxr2  = hysplit.add_species(hxr)


returns an '~xarray.DataArray' which is total concentration.


To combine multiple cdump files, use the following command.



.. code-block:: python


    file1 = (monetio.__path__ + '/../data/cdump.bin, 'S1', 'ERA5e0')
    file2 = (monetio.__path__ + '/../data/cdump2.bin, 'S1', 'ERA5e1')

    d1 = datetime.datetime(2019,2,25,18)
    d2 = datetime.datetime(2019,2,25,20)

    hda = hysplit.combine_dataset([file1,file2], drange=[d1,d2])


returns an '~xarray.DataArray' which is total concentration over all species.
The data-array has dimensions of x,y,z,time, source, ens.
The ens dimension tags what ensemble member the run belongs to.
The source dimension tags which source term the run used.


.. code-block:: python


   print(hda)


   <xarray.DataArray (source: 1, ens: 2, time: 2, z: 3, y: 38, x: 26)>
    array([[[[[[0., ..., 0.],
               ...,
               [0., ..., 0.]],

              ...,
              ...,

              [[0., ..., 0.],
               ...,
               [0., ..., 0.]]]]]], dtype=float32)
    Coordinates:
      * y          (y) int64 4458 4459 4460 4461 4462 ... 4491 4492 4493 4494 4495
      * x          (x) int64 9000 9001 9002 9003 9004 ... 9021 9022 9023 9024 9025
      * time       (time) datetime64[ns] 2019-02-25T18:00:00 2019-02-25T19:00:00
      * z          (z) int64 3000 4000 5000
      * ens        (ens) <U6 'Era5e0' 'Era5e1'
      * source     (source) <U2 'S1'
        latitude   (y, x) float64 -0.5145 -0.5145 -0.5145 ... -0.1445 -0.1445
        longitude  (y, x) float64 -77.72 -77.71 -77.7 ... -77.49 -77.48 -77.47
    Attributes:
        sample time hours:  1.0

To calcluate mass loading

.. code-block:: python

   massload =  hypslit.hysp_massload(hxr, threshold=0, mult=1e10, species=None)

All points with value below or equal to threshold will be returned as 0.
mult is a multiplicative factor applied before the thresholding.
species can be a list of values from the "Species ID" attribute. 
If it is None then all species will be used.

To find top heights

.. code-block:: python

   massload =  hypslit.hysp_hysp_heights(hxr, threshold=0, height_mult=1/1000.0, mult=1e10, mass_load=False)

returns xarray DataArray which gives top height of each level which contains mass loading higher
than the given threshold value. mult is a mutiplicative factor applied before thresholding.
height_mult is a multiplicative factor used to convert heights from meters to some other unit. 
In this example heights are converted to km.
mass_load is a boolean which indicates whether the height should be determined from the mass loading value (True)
or the concentration value (False). 




