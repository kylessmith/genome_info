import os
from os.path import join
import os
import requests
import shutil
from tqdm.auto import tqdm
import pickle


def download_file(urlpath: str,
                  name: str):
    """
    Download a file from the data catalog

    Parameters
    ----------
        urlpath : str
            URL path to the file to download
        name : str
            Name of the file to download

    Returns
    -------
        None
    """

    # Determine the path to the data directory
    destdir = os.path.split(os.path.realpath(__file__))[0]

    # make an HTTP request within a context manager
    with requests.get(urlpath, stream=True) as r:
        
        # check header to get content length, in bytes
        total_length = int(r.headers.get("Content-Length"))
        
        # implement progress bar via tqdm
        with tqdm.wrapattr(r.raw, "read", total=total_length, desc="")as raw:
        
            # save the output to a file
            with open(join(destdir, name), 'wb') as output:
                shutil.copyfileobj(raw, output)

    return None


def get_base_file(genome_name: str):
    """
    Get the base file name.

    Parameters
    ----------
        genome_name : str
            The genome name.

    Returns
    -------
        base_file : str
            The base file name.
    """

    urlpath = "https://raw.githubusercontent.com/kylessmith/" + genome_name + "_info/master/base/" + genome_name + ".pickle"
    download_file(urlpath, join("base", genome_name + ".pickle"))

    return None


def download_genome_files(genome_name: str):
    """
    Download the genome files.

    Parameters
    ----------
        genome_name : str
            The genome name.

    Returns
    -------
        None
    """

    from zipfile import ZipFile

    # Read base file
    base_directory = join(os.path.split(os.path.realpath(__file__))[0], "base")
    with open(join(base_directory, genome_name + '.pickle'), 'rb') as openfile:
        pickle_object = pickle.load(openfile)

    # Download files
    data_directory = join(os.path.split(os.path.realpath(__file__))[0], "data")
    for i in range(pickle_object["package_number"]):
        n = str(i + 1)
        filename = join(data_directory, genome_name + "_" + n + ".zip")
        urlpath = "https://raw.githubusercontent.com/kylessmith/" + genome_name + "_info/master/data/" + genome_name + "_" + n + ".zip"
        download_file(urlpath, filename)
        
        # Unzip files
        with ZipFile(filename, 'r') as zObject:
            zObject.extractall(path=join(data_directory, genome_name))
        
        # Clean up
        os.remove(join(data_directory, genome_name + "_" + n + ".zip"))

    # Download external files
    external_directory = join(join(data_directory, genome_name), "external")
    os.mkdir(external_directory)
    for file in pickle_object["external_downloads"]:
        urlpath = pickle_object["external_downloads"][file]
        download_file(urlpath, join(external_directory, file))

    return None


def check_genome(genome_name: str) -> bool:
    """
    Check if the genome is downloaded.

    Parameters
    ----------
        genome_name : str
            The genome name.

    Returns
    -------
        is_downloaded : bool
            Whether the genome is downloaded.
    """

    # Check if base file exists
    current_dir = os.path.split(os.path.realpath(__file__))[0]
    base_file = join(current_dir, "base", genome_name + ".pickle")
    if not os.path.exists(base_file):
        return False
    
    # Check if data directory exists
    data_directory = join(current_dir, "data", genome_name)
    if not os.path.exists(data_directory):
        return False
    
    # Check if external directory exists
    external_directory = join(data_directory, "external")
    if not os.path.exists(external_directory):
        return False
    
    return True


def download_genome(genome_name: str):
    """
    Download the genome.

    Parameters
    ----------
        genome_name : str
            The genome name.

    Returns
    -------
        None
    """

    # Check if genome is downloaded
    if check_genome(genome_name):
        return None
    
    # Download files
    get_base_file(genome_name)
    download_genome_files(genome_name)

    return None
