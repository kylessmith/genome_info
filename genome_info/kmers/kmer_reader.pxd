import numpy as np
cimport numpy as np
np.import_array()
cimport cython
from libc.stdint cimport uint32_t, int32_t, int64_t, uint16_t
from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t
from ailist.array_query_core cimport pointer_to_numpy_array


cdef extern from "interval_kmer.c":
    # C is include here so that it doesn't need to be compiled externally
    pass

cdef extern from "2bit.c":
    # C is include here so that it doesn't need to be compiled externally
    pass

cdef extern from "2bit.h":
    # C is include here so that it doesn't need to be compiled externally
    pass

cdef extern from "interval_kmer.h":
	
    ctypedef struct kmer_t:
        char *name
        int count

    ctypedef struct kmer_count_t:
        int max_kmers
        int n_kmers
        kmer_t *kmers
        void *kmer_lookup

    ctypedef struct base_freq_t:
        float *A
        float *T
        float *G
        float *C
        int n_intervals
        int n_bases
        int up
        int down

    ctypedef struct tribase_freq_t:
        float *AAA
        float *AAT
        float *AAG
        float *AAC
        float *ATA
        float *ATT
        float *ATG
        float *ATC
        float *AGA
        float *AGT
        float *AGG
        float *AGC
        float *ACA
        float *ACT
        float *ACG
        float *ACC
        float *TAA
        float *TAT
        float *TAG
        float *TAC
        float *TTA
        float *TTT
        float *TTG
        float *TTC
        float *TGA
        float *TGT
        float *TGG
        float *TGC
        float *TCA
        float *TCT
        float *TCG
        float *TCC
        float *GAA
        float *GAT
        float *GAG
        float *GAC
        float *GTA
        float *GTT
        float *GTG
        float *GTC
        float *GGA
        float *GGT
        float *GGG
        float *GGC
        float *GCA
        float *GCT
        float *GCG
        float *GCC
        float *CAA
        float *CAT
        float *CAG
        float *CAC
        float *CTA
        float *CTT
        float *CTG
        float *CTC
        float *CGA
        float *CGT
        float *CGG
        float *CGC
        float *CCA
        float *CCT
        float *CCG
        float *CCC
        int n_intervals
        int n_bases
        int up
        int down

    ctypedef struct interval_base_freq_t:
        base_freq_t *start
        base_freq_t *end

    ctypedef struct interval_tribase_freq_t:
        tribase_freq_t *start
        tribase_freq_t *end

    #-------------------------------------------------------------------------------------
    # interval_kmer.c
    #=====================================================================================

    void interval_base_freq_destroy(interval_base_freq_t *ibf) nogil
    interval_base_freq_t *read_interval_base_freq(labeled_aiarray_t *laia, char *fname, int n_bases) nogil

    void interval_tribase_freq_destroy(interval_tribase_freq_t *itbf) nogil
    interval_tribase_freq_t *read_interval_tribase_freq(labeled_aiarray_t *laia, char *fname, int n_bases) nogil

    void kmer_count_destroy(kmer_count_t *kc) nogil
    int fetch_kmer(kmer_count_t *kc, char *seq) nogil
    kmer_count_t *interval_kmer_count(labeled_aiarray_t *laia, char *fname, int kmer, int last_n) nogil
    char *fetch_sequence(char *fname, char *name, int start, int end) nogil
    void fetch_sequence_code(char *fname, char *name, int start, int end, int *seq_code) nogil
    void gc_content(labeled_aiarray_t *laia, char *fname, float gc[]) nogil


cdef kmer_count_t *_read_kmers(char *fname, labeled_aiarray_t *laia, int k, int last_n)
cdef bytes _fetch_sequence(char *fname, char *name, int start, int end)
cdef void _fetch_sequence_code(char *fname, char *name, int start, int end, int[::1] seq_code)
cdef void _gc_percent(char *fname, labeled_aiarray_t *laia, float[::1] gc)