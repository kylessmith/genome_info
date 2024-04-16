import numpy as np
from intervalframe import IntervalFrame
from ailist import LabeledIntervalArray
import pandas as pd


def add_upstream(intervals: IntervalFrame,
                 n: int,
                 selection: np.ndarray | None = None) -> IntervalFrame:
    """
    Add upstream to intervals

    Parameters
    ----------
        intervals : IntervalFrame
            Intervals to add upstream to
        n : int
            Number of bases to add upstream
        selection : np.ndarray
            Selection of intervals to add upstream to

    Returns
    -------
        iframe : IntervalFrame
            Intervals with upstream added
    """

    # Extract arrays
    starts = intervals.index.starts
    ends = intervals.index.ends
    labels = intervals.index.labels

    # Add upstream 
    if selection is None:
        starts = starts - n
    else:
        starts[selection] = starts[selection] - n

    # Correct zeros
    starts[starts < 0] = 0

    # Construct IntervalFrame
    new_intervals = LabeledIntervalArray()
    new_intervals.add(starts, ends, labels)
    iframe = IntervalFrame(new_intervals, intervals.df)

    return iframe


def add_downstream(intervals: IntervalFrame,
                   n: int,
                   selection: np.ndarray | None = None) -> IntervalFrame:
    """
    Add downstream to intervals

    Parameters
    ----------
        intervals : IntervalFrame
            Intervals to add downstream to
        n : int
            Number of bases to add downstream
        selection : np.ndarray
            Selection of intervals to add downstream to
    
    Returns
    -------
        iframe : IntervalFrame
            Intervals with downstream added
    """

    # Extract arrays
    starts = intervals.index.starts
    ends = intervals.index.ends
    labels = intervals.index.labels

    # Add downstream 
    if selection is None:
        ends = ends + n
    else:
        ends[selection] = ends[selection] + n

    # Construct IntervalFrame
    new_intervals = LabeledIntervalArray()
    new_intervals.add(starts, ends, labels)
    iframe = IntervalFrame(new_intervals, intervals.df)

    return iframe


def adjust_bounds(iframe: IntervalFrame,
                  upstream: int,
                  downstream: int,
                  filter_duplicates: bool = True) -> IntervalFrame:
    """
    Adjust bounds of intervals

    Parameters
    ----------
        iframe : IntervalFrame
            Intervals to adjust
        upstream : int
            Number of bases to add upstream
        downstream : int
            Number of bases to add downstream
        filter_duplicates : bool
            Filter duplicates

    Returns
    -------
        iframe : IntervalFrame
            Intervals with bounds adjusted
    """

    # Filter duplicates
    if filter_duplicates:
        selection1 = np.logical_and(iframe.df.loc[:,"Strand"].values=="+",
                                    ~pd.Index(iframe.df.loc[:,"gene_name"].values).duplicated(keep="first"))
        selection2 = np.logical_and(iframe.df.loc[:,"Strand"].values=="-",
                                    ~pd.Index(iframe.df.loc[:,"gene_name"].values).duplicated(keep="last"))
        is_unique = np.logical_or(selection1, selection2)
        iframe = iframe.iloc[is_unique,:]

    if upstream != 0:
        selection = iframe.df.loc[:,"Strand"].values=="+"
        iframe = add_upstream(iframe, upstream, selection)
        selection = iframe.df.loc[:,"Strand"].values=="-"
        iframe = add_upstream(iframe, upstream, selection)

    if downstream != 0:
        selection = iframe.df.loc[:,"Strand"].values=="+"
        iframe = add_downstream(iframe, downstream, selection)
        selection = iframe.df.loc[:,"Strand"].values=="-"
        iframe = add_downstream(iframe, downstream, selection)

    return iframe