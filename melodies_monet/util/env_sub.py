import os

def env_sub(file_str):
    """Substitute envirnoment variables prefixed by $
    in file location string.

    Parameters
    -------
    -------
    file_str (str): file string

    Returns
    -------
    file_str (str): modified file string
    """
    for subpath in file_str.split("/"):
        if "$" in subpath:
            envvar = subpath.replace("$", "")
            envval = os.getenv(envvar)
            if envval is None:
                print("environment variable not defined: " + subpath)
                exit(1)
            else:
                file_str = file_str.replace(subpath, envval)

    return file_str
