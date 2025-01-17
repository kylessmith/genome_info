#cython: embedsignature=True
#cython: profile=False

import pandas as pd
import numpy as np
cimport numpy as np
np.import_array()

from libc.stdint cimport uint32_t, int32_t, int64_t, uint16_t
from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t
from ailist.array_query_core cimport pointer_to_numpy_array
from ..genomes.genomes import InfoReader


cdef kmer_count_t *_read_kmers(char *fname, labeled_aiarray_t *laia, int k, int last_n):
    cdef kmer_count_t *kc = interval_kmer_count(laia, fname, k, last_n)

    return kc


def read_kmers(str seq_fn, LabeledIntervalArray laia, int k, int last_n = 0):

    cdef str twobit_name = seq_fn
    cdef bytes fname = twobit_name.encode()
    cdef kmer_count_t *kc = _read_kmers(fname, laia.laia, k, last_n)

    result = {}
    cdef bytes kmer_seq
    cdef int count
    cdef uint32_t t
    for t in range(kc.n_kmers):
        kmer_seq = kc.kmers[t].name
        count = fetch_kmer(kc, kmer_seq)

        result[kmer_seq.decode()] = count

    # Delete
    kmer_count_destroy(kc)

    return result


cdef bytes _fetch_sequence(char *fname, char *name, int start, int end):
    cdef bytes seq = fetch_sequence(fname, name, start, end)

    return seq

def read_sequence(str seq_fn, str chrom, int start, int end):

    cdef str twobit_name = seq_fn
    cdef bytes fname = twobit_name.encode()
    cdef bytes name = chrom.encode()

    cdef bytes seq = _fetch_sequence(fname, name, start, end)
    cdef str py_seq = seq.decode()

    return py_seq


cdef void _fetch_sequence_code(char *fname, char *name, int start, int end, int[::1] seq_code):
    fetch_sequence_code(fname, name, start, end, &seq_code[0])

    return 


def read_sequence_code(str seq_fn, str chrom, int start, int end):
    
    cdef str twobit_name = seq_fn
    cdef bytes fname = twobit_name.encode()
    cdef bytes name = chrom.encode()

    cdef np.ndarray seq_arr = np.zeros(end - start, dtype=np.int32)
    cdef int[::1] seq_code = seq_arr
    _fetch_sequence_code(fname, name, start, end, seq_code)

    return seq_arr


def read_sequence_code_intervals(str seq_fn, LabeledIntervalArray laia, int max_len = 1000):
    
    cdef str twobit_name = seq_fn
    cdef bytes fname = twobit_name.encode()

    cdef np.ndarray seq_arr = np.zeros((laia.size, max_len), dtype=np.int8)
    seq_arr[:, :] = -1

    cdef int i
    for i in range(laia.size):
        interval = laia[i]
        arr = read_sequence_code(seq_fn, interval.label, interval.start, interval.end)
        arr = arr.astype(np.int8)
        seq_arr[i, :len(arr)] = arr

    return seq_arr


cdef void _gc_percent(char *fname, labeled_aiarray_t *laia, float[::1] gc):
    gc_content(laia, fname, &gc[0])

    return


def gc_percent(str seq_fn, LabeledIntervalArray laia, str genome_version = "hg38"):
    cdef str twobit_name = seq_fn
    cdef bytes fname = twobit_name.encode()

    cdef np.ndarray gc = np.zeros(laia.size, dtype=np.single)
    cdef float[::1] gc_mem = gc

    _gc_percent(fname, laia.laia, gc_mem)

    return gc


def read_bounds_base_freq(str seq_fn, LabeledIntervalArray laia, int n_bases):

    cdef str twobit_name = seq_fn
    cdef bytes fname = twobit_name.encode()

    cdef np.ndarray start_freq = np.zeros((4, n_bases), dtype=float)
    cdef np.ndarray end_freq = np.zeros((4, n_bases), dtype=float)
    cdef interval_base_freq_t *ibf = read_interval_base_freq(laia.laia, fname, n_bases)
    cdef int i
    for i in range(n_bases):
        start_freq[0, i] = ibf.start.A[i]
        start_freq[1, i] = ibf.start.T[i]
        start_freq[2, i] = ibf.start.G[i]
        start_freq[3, i] = ibf.start.C[i]
        end_freq[0, i] = ibf.end.A[i]
        end_freq[1, i] = ibf.end.T[i]
        end_freq[2, i] = ibf.end.G[i]
        end_freq[3, i] = ibf.end.C[i]

    interval_base_freq_destroy(ibf)

    return start_freq, end_freq


def read_bounds_tribase_freq(str seq_fn, LabeledIntervalArray laia, int n_bases):
    
    cdef str twobit_name = seq_fn
    cdef bytes fname = twobit_name.encode()

    cdef np.ndarray start_freq = np.zeros((64, n_bases), dtype=float)
    cdef np.ndarray end_freq = np.zeros((64, n_bases), dtype=float)

    cdef interval_tribase_freq_t *ibf = read_interval_tribase_freq(laia.laia, fname, n_bases)
    cdef int i
    for i in range(n_bases):
        start_freq[0, i] = ibf.start.AAA[i]
        start_freq[1, i] = ibf.start.AAC[i]
        start_freq[2, i] = ibf.start.AAG[i]
        start_freq[3, i] = ibf.start.AAT[i]
        start_freq[4, i] = ibf.start.ACA[i]
        start_freq[5, i] = ibf.start.ACC[i]
        start_freq[6, i] = ibf.start.ACG[i]
        start_freq[7, i] = ibf.start.ACT[i]
        start_freq[8, i] = ibf.start.AGA[i]
        start_freq[9, i] = ibf.start.AGC[i]
        start_freq[10, i] = ibf.start.AGG[i]
        start_freq[11, i] = ibf.start.AGT[i]
        start_freq[12, i] = ibf.start.ATA[i]
        start_freq[13, i] = ibf.start.ATC[i]
        start_freq[14, i] = ibf.start.ATG[i]
        start_freq[15, i] = ibf.start.ATT[i]
        start_freq[16, i] = ibf.start.CAA[i]
        start_freq[17, i] = ibf.start.CAC[i]
        start_freq[18, i] = ibf.start.CAG[i]
        start_freq[19, i] = ibf.start.CAT[i]
        start_freq[20, i] = ibf.start.CCA[i]
        start_freq[21, i] = ibf.start.CCC[i]
        start_freq[22, i] = ibf.start.CCG[i]
        start_freq[23, i] = ibf.start.CCT[i]
        start_freq[24, i] = ibf.start.CGA[i]
        start_freq[25, i] = ibf.start.CGC[i]
        start_freq[26, i] = ibf.start.CGG[i]
        start_freq[27, i] = ibf.start.CGT[i]
        start_freq[28, i] = ibf.start.CTA[i]
        start_freq[29, i] = ibf.start.CTC[i]
        start_freq[30, i] = ibf.start.CTG[i]
        start_freq[31, i] = ibf.start.CTT[i]
        start_freq[32, i] = ibf.start.GAA[i]
        start_freq[33, i] = ibf.start.GAC[i]
        start_freq[34, i] = ibf.start.GAG[i]
        start_freq[35, i] = ibf.start.GAT[i]
        start_freq[36, i] = ibf.start.GCA[i]
        start_freq[37, i] = ibf.start.GCC[i]
        start_freq[38, i] = ibf.start.GCG[i]
        start_freq[39, i] = ibf.start.GCT[i]
        start_freq[40, i] = ibf.start.GGA[i]
        start_freq[41, i] = ibf.start.GGC[i]
        start_freq[42, i] = ibf.start.GGG[i]
        start_freq[43, i] = ibf.start.GGT[i]
        start_freq[44, i] = ibf.start.GTA[i]
        start_freq[45, i] = ibf.start.GTC[i]
        start_freq[46, i] = ibf.start.GTG[i]
        start_freq[47, i] = ibf.start.GTT[i]
        start_freq[48, i] = ibf.start.TAA[i]
        start_freq[49, i] = ibf.start.TAC[i]
        start_freq[50, i] = ibf.start.TAG[i]
        start_freq[51, i] = ibf.start.TAT[i]
        start_freq[52, i] = ibf.start.TCA[i]
        start_freq[53, i] = ibf.start.TCC[i]
        start_freq[54, i] = ibf.start.TCG[i]
        start_freq[55, i] = ibf.start.TCT[i]
        start_freq[56, i] = ibf.start.TGA[i]
        start_freq[57, i] = ibf.start.TGC[i]
        start_freq[58, i] = ibf.start.TGG[i]
        start_freq[59, i] = ibf.start.TGT[i]
        start_freq[60, i] = ibf.start.TTA[i]
        start_freq[61, i] = ibf.start.TTC[i]
        start_freq[62, i] = ibf.start.TTG[i]
        start_freq[63, i] = ibf.start.TTT[i]
        end_freq[0, i] = ibf.end.AAA[i]
        end_freq[1, i] = ibf.end.AAC[i]
        end_freq[2, i] = ibf.end.AAG[i]
        end_freq[3, i] = ibf.end.AAT[i]
        end_freq[4, i] = ibf.end.ACA[i]
        end_freq[5, i] = ibf.end.ACC[i]
        end_freq[6, i] = ibf.end.ACG[i]
        end_freq[7, i] = ibf.end.ACT[i]
        end_freq[8, i] = ibf.end.AGA[i]
        end_freq[9, i] = ibf.end.AGC[i]
        end_freq[10, i] = ibf.end.AGG[i]
        end_freq[11, i] = ibf.end.AGT[i]
        end_freq[12, i] = ibf.end.ATA[i]
        end_freq[13, i] = ibf.end.ATC[i]
        end_freq[14, i] = ibf.end.ATG[i]
        end_freq[15, i] = ibf.end.ATT[i]
        end_freq[16, i] = ibf.end.CAA[i]
        end_freq[17, i] = ibf.end.CAC[i]
        end_freq[18, i] = ibf.end.CAG[i]
        end_freq[19, i] = ibf.end.CAT[i]
        end_freq[20, i] = ibf.end.CCA[i]
        end_freq[21, i] = ibf.end.CCC[i]
        end_freq[22, i] = ibf.end.CCG[i]
        end_freq[23, i] = ibf.end.CCT[i]
        end_freq[24, i] = ibf.end.CGA[i]
        end_freq[25, i] = ibf.end.CGC[i]
        end_freq[26, i] = ibf.end.CGG[i]
        end_freq[27, i] = ibf.end.CGT[i]
        end_freq[28, i] = ibf.end.CTA[i]
        end_freq[29, i] = ibf.end.CTC[i]
        end_freq[30, i] = ibf.end.CTG[i]
        end_freq[31, i] = ibf.end.CTT[i]
        end_freq[32, i] = ibf.end.GAA[i]
        end_freq[33, i] = ibf.end.GAC[i]
        end_freq[34, i] = ibf.end.GAG[i]
        end_freq[35, i] = ibf.end.GAT[i]
        end_freq[36, i] = ibf.end.GCA[i]
        end_freq[37, i] = ibf.end.GCC[i]
        end_freq[38, i] = ibf.end.GCG[i]
        end_freq[39, i] = ibf.end.GCT[i]
        end_freq[40, i] = ibf.end.GGA[i]
        end_freq[41, i] = ibf.end.GGC[i]
        end_freq[42, i] = ibf.end.GGG[i]
        end_freq[43, i] = ibf.end.GGT[i]
        end_freq[44, i] = ibf.end.GTA[i]
        end_freq[45, i] = ibf.end.GTC[i]
        end_freq[46, i] = ibf.end.GTG[i]
        end_freq[47, i] = ibf.end.GTT[i]
        end_freq[48, i] = ibf.end.TAA[i]
        end_freq[49, i] = ibf.end.TAC[i]
        end_freq[50, i] = ibf.end.TAG[i]
        end_freq[51, i] = ibf.end.TAT[i]
        end_freq[52, i] = ibf.end.TCA[i]
        end_freq[53, i] = ibf.end.TCC[i]
        end_freq[54, i] = ibf.end.TCG[i]
        end_freq[55, i] = ibf.end.TCT[i]
        end_freq[56, i] = ibf.end.TGA[i]
        end_freq[57, i] = ibf.end.TGC[i]
        end_freq[58, i] = ibf.end.TGG[i]
        end_freq[59, i] = ibf.end.TGT[i]
        end_freq[60, i] = ibf.end.TTA[i]
        end_freq[61, i] = ibf.end.TTC[i]
        end_freq[62, i] = ibf.end.TTG[i]
        end_freq[63, i] = ibf.end.TTT[i]

    interval_tribase_freq_destroy(ibf)

    tribases = ["AAA", "AAC",
                "AAG", "AAT",
                "ACA", "ACC",
                "ACG", "ACT",
                "AGA", "AGC",
                "AGG", "AGT",
                "ATA", "ATC",
                "ATG", "ATT",
                "CAA", "CAC",
                "CAG", "CAT",
                "CCA", "CCC",
                "CCG", "CCT",
                "CGA", "CGC",
                "CGG", "CGT",
                "CTA", "CTC",
                "CTG", "CTT",
                "GAA", "GAC",
                "GAG", "GAT",
                "GCA", "GCC",
                "GCG", "GCT",
                "GGA", "GGC",
                "GGG", "GGT",
                "GTA", "GTC",
                "GTG", "GTT",
                "TAA", "TAC",
                "TAG", "TAT",
                "TCA", "TCC",
                "TCG", "TCT",
                "TGA", "TGC",
                "TGG", "TGT",
                "TTA", "TTC",
                "TTG", "TTT"]
    start_results = pd.DataFrame(start_freq, index=tribases)
    end_results = pd.DataFrame(end_freq, index=tribases)

    return start_results, end_results