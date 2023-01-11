import os
import logging
from glob import glob


def fill_date_template(template_str, date_str):
    """
    Replace date template parameters with values from date string

    Parameters
        template_str (str): template string
        date_str (str yyyy-mm-m_abbr-dd-ddd): date string

    Returns
        None
    """

    yyyy_str, mm_str, m_abbr_str, dd_str, ddd_str \
        = tuple(date_str.split('-'))

    if 'DDD' in template_str:
        return template_str.replace(
            'YYYY', yyyy_str).replace(
            'DDD', ddd_str)
    else:
        return template_str.replace(
            'YYYY', yyyy_str).replace(
            'MM', mm_str).replace(
            'M_ABBR', m_abbr_str).replace(
            'DD', dd_str)


def find_file(datadir, filestr):
    """
    Parameters
        datadir (str): data directory
        filestr (str): filename regular expression

    Returns
        filename (str): complete path of matching filename in data directory
    """

    pattern = os.path.join(os.path.expandvars(datadir), filestr)
    files = glob(pattern)

    if len(files) == 0:
        raise Exception('no file matches for %s' % pattern)
    if len(files) > 1:
        raise Exception('more than one file match %s' % pattern)

    filename = files[0]
    logging.info(filename)

    return filename
