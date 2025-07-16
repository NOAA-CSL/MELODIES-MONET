"""
TEMPO downloader
Created by NSF - National Center for Atmospheric Research, in
the Atmospheric Chemistry Observations and Modeling Laboratory
"""

try:
    import earthaccess
except ImportError:
    raise ImportError(
        "The earthaccess module is required for accessing TEMPO data. "
        + "Download and install it with pip ( pip install earthaccess ) or "
        + "conda ( conda install -c conda-forge earthaccess )"
    )


# """"""""""""""""""""""
# "    User Changes    "
# """"""""""""""""""""""

name_key = "TEMPO_NO2_l2"
start_time = "2024-08-22 18:00:00"
end_time = "2024-08-22 19:00:00"
path = "."


# """"""""""""""""""""""
# "  End User Changes  "
# """"""""""""""""""""""


def main() -> None:
    earthaccess.login(strategy="interactive")
    results = earthaccess.search_data(short_name=name_key, temporal=(start_time, end_time))
    files = [granule["meta"]["native-id"] for granule in results]
    user_input = input(f"Do you want to download {files}? (yes/no) : ")
    if user_input.lower() in ["yes", "y"]:
        print("Downloading")
        earthaccess.download(results, path)
    else:
        print("Did not answer yes, exiting")


if __name__ == "__main__":
    main()
