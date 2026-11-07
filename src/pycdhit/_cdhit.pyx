# coding: utf-8
# cython: language_level=3
"""Bindings to CD-HIT, an accelerated sequence clustering program.

References:
    - Li, W., & Godzik, A. (2006).
      Cd-hit: a fast program for clustering and comparing large sets of
      protein or nucleotide sequences. *Bioinformatics (Oxford, England),
      22(13)*, 1658–1659. :doi:`10.1093/bioinformatics/btl158` Epub 2006
      May 26. :pmid:`16731699`.
    - Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012).
      CD-HIT: accelerated for clustering the next-generation sequencing data.
      *Bioinformatics (Oxford, England), 28(23)*, 3150–3152.
      :doi:`10.1093/bioinformatics/bts565` :pmid:`23060610`.

"""

# --- C imports --------------------------------------------------------------

from libc.string cimport strncmp

from cpython.unicode cimport PyUnicode_AsUTF8AndSize

cimport scoring_matrices
from scoring_matrices cimport ScoringMatrix

cimport cdhit
from cdhit.options cimport Options
from cdhit.seq cimport SequenceDB, Sequence

cdef Options DEFAULT_OPTIONS = Options()

# --- Python imports ---------------------------------------------------------

import array
import dataclasses
import os
import typing

__version__ = PROJECT_VERSION

# --- Functions --------------------------------------------------------------

@dataclasses.dataclass(frozen=True)
class Clustering:
    clusters: array.array
    representatives: array.array

    @property
    def n_clusters(self):
        """`int`: The number of clusters found in the clustering.
        """
        return len(self.representatives)

def cluster(
    sequences: typing.Iterable[str],
    *,
    double identity_threshold = DEFAULT_OPTIONS.cluster_thd,
    int cpus = 0,
    gap_open: int = 11,
    gap_extend: int = 1,
    ScoringMatrix scoring_matrix not None = ScoringMatrix.from_name("BLOSUM62")
) -> Clustering:
    """Cluster protein sequences with the CD-HIT algorithm.

    Arguments:
        sequences (iterable of `str`, `bytes` or buffer-like object): The
            sequences to cluster with CD-HIT.

    Keyword Arguments:
        identity_threshold (`float`): The identity threshold with which to
            cluster the sequences (``-c`` parameter from CD-HIT).
        cpus (`int`): The number of parallel threads to use to clustering.
            Pass 0 to use all available system threads, as given
            by `os.cpu_count`.
        gap_open (`int`): The penalty for opening a new gap in the pairwise
            alignment.
        gap_extend (`int`): The penalty for extending an existing gap in the
            pairwise alignment.
        scoring_matrix (`~scoring_matrices.ScoringMatrix`): The scoring 
            matrix to use for pairwise alignments. Defaults to BLOSUM62.

    """
    cdef ssize_t     length
    cdef size_t      i
    cdef const char* s
    cdef str         sequence
    cdef Sequence*   seq      = NULL
    cdef SequenceDB  db       = SequenceDB()
    cdef Options     options  = Options()

    options.max_memory = 0  # don't attempt limiting memory
    options.threads = 1     # use a single thread (for now)
    options.input = "n/a"   # bogus filename so CD-HIT doesn't complain
    options.output = "n/a"  # bogus filename so CD-HIT doesn't complain

    # Set parallelism from parameters
    if cpus < 0:
        raise ValueError(f"`cpus` must be positive or null, got {cpus!r}")
    elif cpus == 0:
        _cpus = os.cpu_count() or 1
        options.threads = max(_cpus, 1)
    else:
        options.threads = cpus

    # Update clustering parameters
    if identity_threshold >= 1 or identity_threshold <= 0:
        raise ValueError(f"`identity_threshold` must be a float between 0 and 1, got {identity_threshold!r}")
    options.cluster_thd = identity_threshold

    # Update Smith-Waterman parameters
    if gap_open <= 0:
        raise ValueError(f"`gap_open` must be a positive number, got {gap_open!r}")
    if gap_extend <= 0:
        raise ValueError(f"`gap_extend` must be a positive number, got {gap_extend!r}")
    options.mat.set_gap(-gap_open, -gap_extend)
    if scoring_matrix is not None:
        if not scoring_matrix.is_integer():
            raise ValueError(f"`scoring_matrix` requires an integer scoring matrix")
        elif strncmp(scoring_matrix.alphabet_ptr(), cdhit.aa, cdhit.MAX_AA) != 0:
            raise ValueError(f"`scoring_matrix` has an invalid alphabet")
        for i in range(cdhit.MAX_AA):
            for j in range(cdhit.MAX_AA):
                ai = chr(cdhit.aa[i])
                aj = chr(cdhit.aa[j])
                options.mat.matrix[i][j] = <int> scoring_matrix[ai, aj]

    # Use CD-HIT validation code to ensure options are consistent
    try:
        options.Validate()
    except Exception as exc:
        raise ValueError("Failed to validate parameters") from exc

    # Initialize the degeneration rates
    options.InitNAA(cdhit.MAX_UAA)
    options.NAAN = options.NAAN_array[options.NAA]
    db.NAAN = options.NAAN_array[options.NAA]

    # Append sequence data to a new sequence database (seq_db.Read(...))
    for i, sequence in enumerate(sequences):
        s = PyUnicode_AsUTF8AndSize(sequence, &length)
        seq = new Sequence()
        seq[0] = s
        seq.tot_length = length
        seq.index = db.sequences.size()
        db.sequences.Append(seq)

    # Run clustering algorithm with GIL relaxed
    with nogil:
        db.SortDivide(options)
        db.DoClustering(options)

    # Extract flat clusters and representatives
    clusters = array.array('Q')
    for k in range(db.sequences.size()):
        clusters.append(db.sequences[k].cluster_id)
    representatives = array.array('Q')
    for k in range(db.rep_seqs.size()):
        representatives.append(db.rep_seqs[k])
    return Clustering(
        clusters=clusters,
        representatives=representatives,
    )



# --- Init code ----------------------------------------------------------------

