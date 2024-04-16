import os
import glob
from typing import List, Any, Dict
from intervalframe import IntervalFrame


def read_pfm(directory_name: str) -> List[tuple]:
    """
    Read a PFM files from a directory.

    Parameters
    ----------
        directory_name : str
            The directory name.

    Returns
    -------
        pfm : List[tuple]
            The list of PFM tuples.
    """

    import MOODS.parsers
    import MOODS.tools

    # Read matrices
    bg = MOODS.tools.flat_bg(4)
    pseudocount = 0.0001
    filenames = glob.glob(os.path.join(directory_name,"*.pfm"))
    matrices = [MOODS.parsers.pfm_to_log_odds(f, bg, pseudocount) for f in filenames]
    pfm_names = [os.path.splitext(os.path.basename(f))[0].split(".")[0] for f in filenames]

    return pfm_names, matrices


def read_pickle(pickle_name: str) -> Any:
    """
    Read a pickle file.

    Parameters
    ----------
        pickle_name : str
            The file name.

    Returns
    -------
        pickle_object : Any
            Object read from pickle.
    """

    import pickle

    # Read pickle
    with open(pickle_name, "rb") as pickle_file:
        pickle_object = pickle.load(pickle_file)

    return pickle_object


def read_muiltiIntervalFrame(directory_name: str) -> Dict[str, IntervalFrame]:
    """
    Read a multiIntervalFrame from a directory.

    Parameters
    ----------
        directory_name : str
            The directory name.

    Returns
    -------
        multiIntervalFrame : Dict[IntervalFrame]
            The multiIntervalFrame.
    """

    # Read files
    results = {}
    filenames = glob.glob(os.path.join(directory_name,"*.parquet"))
    for f in filenames:
        name = os.path.splitext(os.path.basename(f))[0].split(".parquet")[0]
        results[name] = IntervalFrame.read_parquet(f)

    return results

