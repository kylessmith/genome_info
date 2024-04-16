import numpy as np
import pandas as pd
import os
import glob
import json
from os.path import join
from typing import List, Any, Dict
from intervalframe import IntervalFrame

# Local imports
from .readers import *
from .download import *


class InfoReader(object):
    """
    Class for reading info files.
    """

    def __init__(self, genome_name: str) -> None:
        """
        Initialize the class.
        """

        # Check if genome is downloaded
        download_genome(genome_name)

        # Read functions
        self.genome_name = genome_name
        self.readers = {"IntervalFrame":IntervalFrame.read_parquet,
                        "DataFrame":pd.read_parquet,
                        "pfm":read_pfm,
                        "pickle":read_pickle,
                        "multi_IntervalFrame":read_muiltiIntervalFrame}
        
        # Find data directory
        self.current_dir = os.path.split(os.path.realpath(__file__))[0]
        self.data_directory = join(self.current_dir, "data")
        self.base_directory = join(self.current_dir, "base")

        # Read base file
        self.info_files = {}
        self.base_object = read_pickle(join(self.base_directory, genome_name + '.pickle'))
        self.name = self.base_object["name"]
        self.version = self.base_object["version"]

        # Find files
        for key in self.base_object["files"]:
            values = self.base_object["files"][key]
            if values[1] in self.readers:
                found_file = glob.glob(join(self.data_directory, genome_name, "*", values[0]))
                if len(found_file) == 1:
                    self.info_files[key] = found_file[0]
                elif len(found_file) == 0:
                    raise ValueError("No file found for %s" % values[0])
                else:
                    raise ValueError("Multiple file found: %s" % found_file)
            else:
                raise ValueError("Unknown reader type: %s" % values[1])
            
        # Assign keys
        self.keys = list(self.info_files.keys()) + list(self.base_object["external_downloads"].keys())
            
        return None
    

    def __getitem__(self, key: str) -> Any:
        """
        Get an item.
        """

        # Read file
        try:
            read_type = self.base_object["files"][key][1]
            if read_type in self.readers:
                data = self.readers[read_type](self.info_files[key])
            else:
                raise ValueError("Unknown reader type: %s" % key)
        except KeyError:
            if key in self.base_object["external_downloads"]:
                data = self.base_object["external_downloads"][key]
            else:
                raise KeyError("Unknown key: %s" % key)
        
        return data
    
