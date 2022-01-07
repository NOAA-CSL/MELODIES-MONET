import os
import sys
import logging
import numpy as np
import pyhdf.SD as hdf


hdftypes = {'int16': hdf.SDC.INT16, 'uint16': hdf.SDC.UINT16,
            'int32': hdf.SDC.INT32, 'uint32': hdf.SDC.UINT32,
            'float32': hdf.SDC.FLOAT32, 'float64': hdf.SDC.FLOAT64}


def hdf_open(filename):
    """
    filename - file name
    return - file id
    """
    try:
        fileid = hdf.SD(filename, hdf.SDC.READ)
        logging.debug('hdfio.hdf_open:' + filename)
    except IOError:
        logging.error('hdfio.hdf_open:' + filename)
        sys.exit(1)
    return fileid


def hdf_create(filename):
    """
    filename - file name
    return - file id
    """
    try:
        fileid = hdf.SD(filename,
            hdf.SDC.WRITE | hdf.SDC.CREATE | hdf.SDC.TRUNC)
        logging.debug('hdfio.hdf_create:' + filename)
    except IOError:
        logging.error('hdfio.hdf_create:' + filename)
        sys.exit(1)
    return fileid


def hdf_close(fileid):
    """
    fileid - file id
    """
    fileid.end()


def hdf_list(fileid):
    """
    fileid - file id
    return datasets
    """
    datasets = sorted(fileid.datasets())
    indices = list()
    for dataset in datasets:
        index = hdf.SD.nametoindex(fileid, dataset)
        indices.append(index)
        logging.debug(
            'hdfio.hdf_list:' + str(index) + ':' + dataset)
    return datasets, indices


def hdf_read(fileid, varname):
    """
    fileid - file id
    varname - variable name
    return data
    """
    logging.debug('hdfio.hdf_read:' + varname)
    varid = fileid.select(varname)
    data = varid[:]
    return data


def hdf_write_coord(fileid, coordname, data):
    """
    fileid - file id
    coordname - coordinate variable name
    data - coordinate array
    """
    logging.debug('hdfio.hdf_write_coord:' + coordname)
    coordid = fileid.create(
        coordname, hdftypes[str(data.dtype)], data.shape)
    dimid = coordid.dim(0)
    dimid.setname(coordname)
    coordid[:] = data
    coordid.endaccess()


def hdf_write_field(fileid, fieldname, coordnames, data, fill=None):
    """
    fileid - file id
    fieldname - field variable name
    coordnames - tuple of coordinate variable names
    data - field array
    fill - fill value
    """
    logging.debug('hdfio.hdf_write_field:' + fieldname)
    fieldid = fileid.create(
        fieldname, hdftypes[str(data.dtype)], data.shape)
    for i in range(len(coordnames)):
        dimid = fieldid.dim(i)
        dimid.setname(coordnames[i])
    fieldid[:] = data
    if fill is not None:
        fieldid.setfillvalue(fill)
    fieldid.endaccess()
