# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(
    name='melodies_monet',
    version='0.1',
    url='https://github.com/noaa-csl/melodies_monet',
    license='MIT',
    include_package_data=True,
    author='Barry Baker',
    author_email='barry.baker@noaa.gov',
    maintainer='Barry Baker',
    maintainer_email='barry.baker@noaa.gov',
    packages=find_packages(),
    package_data={'': ['data/*.txt', 'data/*.dat', 'data/*.hdf', 'data/*.ncf', 'data/*.jpg', 'data/*.png', 'data/*.ncf']},
    keywords=['model', 'verification', 'hysplit', 'cmaq', 'atmosphere', 'camx', 'evaluation', 'ufs', 'rrfs'],
    description='MELODIES-MONET Unified verification package',
    install_requires=[
        'pandas',
        'netcdf4',
        'xarray',
        'dask',
        'pyresample',
        'matplotlib',
        'seaborn',
        'cartopy',
        'pydecorate',
        'global_land_mask',
        'pyyaml',
    ],
    extra_requires={'xesmf;platform_system!="Windows"'},
)
