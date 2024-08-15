Troubleshooting
===============

Installation problems
---------------------
* Conda installation fails:
    * Often the problem is in the installation of wrf-python. Check that your computer does not have an Apple Silicon CPU (Apple Intel should be fine) and that the Python version is compatible with the wrf-python Conda package.
 
Runtime issues
--------------
* Pairing (:meth:`~melodies_monet.driver.analysis.pair_models`) takes too long:
    * ``analysis.pair_models()`` is one of the most expensive functions in MELODIES MONET, and you might be running out of memory. If you have access to more RAM, try it with that. A Time Chunking functionality is being developed to deal with this issue.
* The plots are not produced. The error message references ``leg.legendHandles``:
    * You are probably using matplotlib 3.9+ with Pandas 1.x. Downgrade matplotlib to 3.8 (upgrading Pandas should also work, but you might run into some incompatibilities for some specific functionalities, especially those related to downloading observational data with MONETIO. Check :doc:`/getting_started/installation`).
* Failure downoloading maps:
    * MELODIES-MONET uses Cartopy.crs for mapping. Cartopy.crs downloads shapefiles automatically to create the maps, and if you are offline or working on a machine with download restrictions, this might fail. Check :doc:`/appendix/machine-specific-install`. This is the situatino for *NOAA HPC Hera*, and the solution discussed there should work.
