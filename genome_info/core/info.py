import pandas as pd
import numpy as np
import os
from os.path import join
from typing import List, Any, Dict
from ailist import LabeledIntervalArray
from intervalframe import IntervalFrame
from tqdm import tqdm

# Local imports
from ..genomes.genomes import InfoReader
from ..kmers.kmer_reader import read_kmers, read_sequence, read_sequence_code, read_sequence_code_intervals, gc_percent
from .utilities import adjust_bounds


def get_include():
	"""
	Get file directory if C headers

	Parameters
	----------
		None

	Returns
	-------
		location : str
			Directory to header files
	"""

	# Grab file location
	location = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]

	return location


class GenomeInfo(object):
    """
    Class for reading genome info files.
    """

    def __init__(self, genome_name: str) -> None:
        """
        Initialize the class.
        """

        # Read info file
        self.info = InfoReader(genome_name)
        self.name = self.info.name
        self.version = self.info.version
        self.keys = self.info.keys
        self.seq_file = join(self.info.data_directory, self.info.genome_name, "external", self.info.genome_name + ".2bit")
        self.pfm_scanner = None

        return None
    
    def __repr__(self):
        """
        Print the info file.
        """

        repr_string = "Genome info file for %s version %s\n" % (self.name, self.version)

        return repr_string
    

    def __getitem__(self, key):
        """
        Get an item from the info file.
        """

        try:
            value = self.info[key]
        except KeyError:
            value = self.info.base_object["properties"][key]

        return value

    
    def sequence(self,
                chromosome: str,
                start: int,
                end: int) -> str:
        """
        Get sequence

        Parameters
        ----------
            chromosome : str
                Chromosome name
            start : int
                Start position
            end : int
                End position

        Returns
        -------
            sequence : str
        """

        # Get sequence
        sequence = read_sequence(self.seq_file, chromosome, start, end)

        return sequence
    

    def sequence_code(self,
                        chromosome: str,
                        start: int,
                        end: int) -> str:
        """
        Get sequence

        Parameters
        ----------
            chromosome : str
                Chromosome name
            start : int
                Start position
            end : int
                End position

        Returns
        -------
            sequence : str
        """

        # Get sequence
        sequence = read_sequence_code(self.seq_file, chromosome, start, end)

        return sequence
    

    def sequence_code_intervals(self,
                                intervals: LabeledIntervalArray,
                                max_length: int = 1000) -> np.ndarray:
        """
        Get sequence

        Parameters
        ----------
            chromosome : str
                Chromosome name
            start : int
                Start position
            end : int
                End position

        Returns
        -------
            sequence : str
        """

        # Get sequence
        sequence = read_sequence_code_intervals(self.seq_file, intervals, max_length)

        return sequence

    
    def interval_kmers(self,
                        intervals: LabeledIntervalArray,
                        k: int = 2,
                        last_n: int = 0) -> Dict[str,int]:
        """
        Get kmers from intervals

        Parameters
        ----------
            intervals : LabeledIntervalArray
                Intervals
            k : int
                Kmer length
            last_n : int
                Last n bases to query

        Returns
        -------
            kmers : Dict[str,int]
                Kmer counts
        """

        # Calculate kmers
        kmers = read_kmers(self.seq_file, intervals, k, last_n)

        return kmers


    def count_kmers(self,
                    chrom: str,
                    start: int,
                    end: int,
                    k: int = 2) -> Dict[str,int]:
        """
        Count kmers in interval

        Parameters
        ----------
            chrom : str
                Chromosome
            start : int
                Start position
            end : int
                End position
            k : int
                Kmer length

        Returns
        -------
            kmers : Dict[str,int]
                Kmer counts
        """

        # Calculate kmers
        intervals = LabeledIntervalArray()
        intervals.add(start, end, chrom)
        kmers = read_kmers(self.seq_file, intervals, k)

        return kmers


    def load_pfm_scanner(self,
                         pvalue: float = 5e-05,
                         pseudocounts = 0.0001) -> None:
        """
        Load PFM scanner

        Parameters
        ----------
            None

        Returns
        -------
            None
        """

        import MOODS.parsers
        import MOODS.tools
        import MOODS.scan

        # Get PFM scanner
        pfms = self["pfm"]
        n_motifs = len(pfms[0])
        bg = MOODS.tools.flat_bg(4)

        # Create matrices
        matrices = [None] * 2 * n_motifs
        thresholds = [None] * 2 * n_motifs
        for i, motif in enumerate(pfms[1]):
            #matrices[i] = MOODS.parsers.pfm_to_log_odds(motif, bg, pseudocounts)
            matrices[i] = motif
            matrices[i+n_motifs] = MOODS.tools.reverse_complement(matrices[i])

            thresholds[i] = MOODS.tools.threshold_from_p(motif, bg, pvalue)
            thresholds[i+n_motifs] = thresholds[i]
        
        # Create scanner
        self.pfm_scanner = MOODS.scan.Scanner(7)
        self.pfm_scanner.set_motifs(matrices = matrices,
                                    bg = bg,
                                    thresholds = thresholds)
        self.pfm_names = pfms[0] + pfms[0]

        return None
    

    def pfm_scan(self,
                 sequence: str):
        """
        PFM scan sequence

        Parameters
        ----------
            sequence : str
                Sequence to scan

        Returns
        -------
            results : Dict[float]
                Dictionary of results
        """

        # Load scanner
        if self.pfm_scanner is None:
            self.load_pfm_scanner()

        # Scan sequence
        results = self.pfm_scanner.scan(sequence)

        # Process results
        pfm_results = {}
        for i, rs in enumerate(results):
            for r in rs:
                try:
                    pfm_results[self.pfm_names[i]]
                    pfm_results[self.pfm_names[i]] = max(pfm_results[self.pfm_names[i]], r.score)
                except KeyError:
                    pfm_results[self.pfm_names[i]] = r.score

        return pfm_results
    

    def _motif_match(self,
                    intervals: LabeledIntervalArray) -> IntervalFrame:
        """
        Match motifs to intervals

        Parameters
        ----------
            intervals : LabeledIntervalArray
                Intervals to match

        Returns
        -------
            motif_matches : IntervalFrame
                Motif matches
        """

        # Load scanner
        if self.pfm_scanner is None:
            self.load_pfm_scanner()

        # Match motifs
        names = pd.unique(np.array(self.pfm_names))
        n_motifs = len(names)
        motif_matches = pd.DataFrame(np.zeros((len(intervals), n_motifs),
                                              dtype = np.uint8),
                                     columns = names)

        for i, interval in enumerate(intervals):
            sequence = self.sequence(interval.label, interval.start, interval.end)
            results = self.pfm_scan(sequence)
            for motif in results:
                motif_matches.loc[i, motif] = 1

        # Convert to IntervalFrame
        iframe = IntervalFrame(intervals = intervals, df = motif_matches)
        
        return iframe
    

    def motif_match(self,
                    intervals: LabeledIntervalArray,
                    n_jobs: int = 1) -> IntervalFrame:
        """
        Match motifs to intervals

        Parameters
        ----------
            intervals : LabeledIntervalArray
                Intervals to match
            n_jobs : int
                Number of jobs

        Returns
        -------
            motif_matches : IntervalFrame
                Motif matches
        """

        # Set scanner
        self.pfm_scanner = None

        if n_jobs == 1:
            iframe = self._motif_match(intervals)
        else:
            import math
            from joblib import Parallel, delayed

            def process_interval_chunk(intervals, g):
                m = self._motif_match(intervals)
                
                return m
            
            def chunkify(data, n_jobs):
                """Split data into chunks based on the number of jobs."""
                chunk_size = math.ceil(len(data) / n_jobs)
                for i in range(0, len(data), chunk_size):
                    yield data[i:i + chunk_size]

            interval_chunks = chunkify(intervals, n_jobs)
            results = Parallel(n_jobs=n_jobs)(
                delayed(process_interval_chunk)(chunk, self) for chunk in interval_chunks
            )

            r = results[0]
            iframe = r.concat(results[1:])

        return iframe


    def get_intervals(self,
                        key: str,
                        upstream: int = 0,
                        downstream: int = 0,
                        filter_column: str = None,
                        filter_selection: str = None,
                        filter_duplicates: bool = True) -> IntervalFrame:
        """
        Get an item from the info file.

        Parameters
        ----------
            key : str
                Key to get
            upstream : int
                Upstream of intervals
            downstream : int
                Downstream of intervals
            filter_column : str
                Column to filter on
            filter_selection : str
                Selection to filter on
            filter_duplicates : bool
                Flag to filter duplicates

        Returns
        -------
            value : IntervalFrame
                Intervals
        """
        
        # Read intervals
        value = self.info[key]
        if not isinstance(value, IntervalFrame):
            return None
        
        # Filter intervals
        if filter_column is not None:
            if filter_selection is not None:
                chosen = value.df.loc[:,filter_column].values == filter_selection
                value = value.iloc[chosen,:]

        # Add upstream and downstream
        value = adjust_bounds(value,
                                upstream,
                                downstream,
                                filter_duplicates)

        return value


    def calculate_bias(self,
                        intervals: LabeledIntervalArray,
                        include_blacklist: bool = True,
                        include_repeat: bool = True,
                        include_gc: bool = True,
                        include_mappability: bool = True) -> IntervalFrame:
        """
        Calculate bias per interval
        
        Parameters
        ----------
            intervals : LabeledIntervalArray
                Labeled intervals
            include_blacklist : bool
                Flag to include blacklist
            include_repeat : bool
                Flag to include repeat
            include_gc : bool
                Flag to include gc
            include_mappability : bool
                Flag to include mappability

        Returns
        -------
            bias_record : IntervalFrame
                Bias for given intervals
        """


        # Initialize bias records
        bias_record = IntervalFrame(intervals=intervals)

        # Calculate blacklist
        if include_blacklist:
            blacklist = self["blacklist"]
            bias_record.df.loc[:,"blacklist"] = intervals.percent_coverage(blacklist.index)

        # Calculate repeat
        if include_repeat:
            repeat = self["repeat"]
            bias_record.df.loc[:,"repeat"] = intervals.percent_coverage(repeat.index)

        # Calculate gc
        if include_gc:
            bias_record.df.loc[:,"gc"] = gc_percent(self.seq_file, intervals)

        # Calculate mappability
        if include_mappability:
            mappability = self["mappability"]
            bias_record.df.loc[:,"mappability"] = intervals.percent_coverage(mappability.index)

        return bias_record


    def calculate_bin_bias(self,
                            bin_size: int = 100000,
                            include_blacklist: bool = True,
                            include_repeat: bool = True,
                            include_gc: bool = True,
                            include_mappability: bool = True) -> IntervalFrame:
        """
        Calculate bias per bin
        
        Parameters
        ----------
            bin_size : int
                Size of bins
            include_blacklist : bool
                Flag to include blacklist
            include_repeat : bool
                Flag to include repeat
            include_gc : bool
                Flag to include gc
            include_mappability : bool
                Flag to include mappability

        Returns
        -------
            bias_record : IntervalFrame
                Bias for given bin size
        """

        # Check previous calculation
        key = "bin_bias_%d" % bin_size
        if key in self.keys:
            return self[key]

        # Initialize bins
        intervals = LabeledIntervalArray.create_bin(self["chrom_sizes"], bin_size=bin_size)

        # Initialize bias records
        bias_record = self.calculate_bias(intervals,
                                            include_blacklist,
                                            include_repeat,
                                            include_gc,
                                            include_mappability)

        return bias_record
    

    def convert_gene_names(self,
                           names: np.ndarray) -> np.ndarray:
        """
        Convert gene names

        Parameters
        ----------
            names : np.ndarray
                Gene names
        
        Returns
        -------
            converted_names : np.ndarray
                Converted gene names
        """

        # Read gene names
        gene_names = self["gene_names"]

        # Convert gene names
        converted_names = np.array([gene_names[name] for name in names])

        return converted_names
    

    def annotate_regions(self,
                        intervals: LabeledIntervalArray,
                        promoter_upstream: int = 1000,
                        promoter_downstream: int = 1000) -> np.ndarray:
        """
        Annotate regions

        Parameters
        ----------
            intervals : LabeledIntervalArray
                Labeled intervals
            promoter_upstream : int
                Upstream of promoter
            promoter_downstream : int
                Downstream of promoter

        Returns
        -------
            anno_record : np.ndarray
                Annotation for given intervals
        """


        # Initialize bias records
        anno_record = IntervalFrame(intervals=intervals)

        # Calculate exons
        regions = self["exon"].merge(1)
        anno_record.df.loc[:,"exon"] = intervals.percent_coverage(regions.index)

        # Calculate introns
        regions = self["intron"].merge(1)
        anno_record.df.loc[:,"intron"] = intervals.percent_coverage(regions.index)

        # Calculate promoters
        regions = self.get_intervals(key = "tss",
                                     upstream = promoter_upstream,
                                     downstream = promoter_downstream,
                                     filter_column = "gene_type",
                                     filter_selection = "protein_coding",
                                     filter_duplicates = True).merge(1)
        anno_record.df.loc[:,"promoter"] = intervals.percent_coverage(regions.index)

        # Determine calls
        anno_call = np.repeat("intergenic", anno_record.df.shape[0])
        anno_call[anno_record.df.loc[:,"intron"].values > 0.5] = "intron"
        anno_call[anno_record.df.loc[:,"exon"].values > 0.5] = "exon"
        anno_call[anno_record.df.loc[:,"promoter"].values > 0.5] = "promoter"

        return anno_call
    

    def filter_blacklist(self,
                         intervals: IntervalFrame,
                         blacklist_cutoff: float = 0.1) -> IntervalFrame:
        """
        Filter intervals by blacklist

        Parameters
        ----------
            intervals : IntervalFrame
                Intervals to filter
            blacklist_cutoff : float
                Cutoff for blacklist

        Returns
        -------
            filtered_intervals : IntervalFrame
                Filtered intervals
        """


        # Calculate blacklist
        blacklist = self["blacklist"]
        blacklist_overlap = intervals.index.percent_coverage(blacklist.index)
        filtered_intervals = intervals.iloc[blacklist_overlap < blacklist_cutoff,:]

        return filtered_intervals