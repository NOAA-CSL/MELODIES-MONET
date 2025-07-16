from getpass import getpass
from http.client import NOT_FOUND, UNAUTHORIZED
from pathlib import Path

from requests import Session

url_to_path = lambda url, output_dir: output_dir.joinpath(url.split('/')[-1])
print("Welcome to the ASDC Download Script!\nThis script downloads data from https://asdc.larc.nasa.gov/data/")
with Session() as session:
    # get login
    url = input("Enter the top level URL (you can also enter 'test' to download a small dataset)\n\turl: ")
    if url == "test":
        url = "https://asdc.larc.nasa.gov/data/AJAX/CH2O_1/"
    token = getpass("Enter your token, if you don't have a token, get one from https://urs.earthdata.nasa.gov/\n\ttoken: ")
    if not token:
        print("Token cannot be blank, exiting.")
        exit()
    session.headers = {"Authorization": f"Bearer {token}"}\

    # verify login works
    response = session.get(url)
    if not response.ok:
        if response.status_code == UNAUTHORIZED:
            print(f"Earthdata Login reponded with Unauthorized, did you enter a valid token?")
            exit()
        if response.status_code == NOT_FOUND:
            print(f"The top level URL does not exist, select a URL within https://asdc.larc.nasa.gov/data/")
            exit()
    output_dir = Path('data')

    # get a list of all urls
    pages = [url]
    file_urls = []
    print("Getting file links")
    for i, page in enumerate(pages):
        print(f"Checking {page} for links", end="\r", flush=True)
        response = session.get(page)
        if not response.ok:
            if response.status_code == NOT_FOUND:
                print(f"The follwoing page was not found: {url}")
            else:
                print(f"Recieved {response.reason} status for {page}")
        content = response.content.decode('utf-8')
        if '<table id="indexlist">' not in content:
            print(f"Data table not found for {page}")
            continue
        table_content = content.split('<table id="indexlist">')[-1].split('</table>')[0]
        hrefs = {part.split('"')[0] for i, part in enumerate(table_content.split('href="')) if i}
        for href in hrefs:
            if href.endswith('/'):
                pages.append(page + href)
            else:
                file_urls.append(page + href)
    if not file_urls:
        print("No files found, exiting.")
        exit()
    
    # offer to remove existing data
    output_dir.mkdir(exist_ok=True)
    if output_dir.exists() and len(list(output_dir.iterdir())):
        if input(f"There's already data in {output_dir.absolute()}, \n\tRemove it? [y/n]: ") == "y":
            for path in output_dir.iterdir():
                path.unlink()

    # get a list of new files (ignore already download files if the size is the same)
    print("Getting size")
    total_size = 0
    file_count = len(file_urls)
    new_files = []
    for i, url in enumerate(file_urls):
        print(f"Getting size for file {i+1} of {file_count}", end="\r", flush=True)
        _response = session.head(url)
        if url_to_path(url, output_dir).exists() and _response.headers.get('content-length') != url_to_path(url, output_dir).stat().st_size:
            continue
        total_size += int(_response.headers.get('content-length'))
        new_files.append(url)
    if not new_files:
        print("No new files, exiting.")
        exit()
    if input(f"Found {len(new_files)} files totaling {total_size // 1024**2} MB in {output_dir.absolute()}.\n\tDownload [y/n]: ") == 'n':
        exit()
    
    # downlaod files
    for i, url in enumerate(new_files):
        print(f"Downloading file {i+1} of {len(file_urls)}", end="\r", flush=True)
        _response = session.get(url)
        with url_to_path(url, output_dir).open('wb') as file:
            file.write(_response._content)
    print("\nDownload Complete")

