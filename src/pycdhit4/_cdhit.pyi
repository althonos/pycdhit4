import array
import dataclasses
import typing

from scoring_matrices import ScoringMatrix

__version__: str

@dataclasses.dataclass(frozen=True)
class Clustering:
    clusters: typing.Sequence[int]
    representatives: typing.Sequence[int]

    @property
    def n_clusters(self) -> int: ...

def cluster(
    sequences: typing.Iterable[str],
    *,
    identity_threshold: float = 0.9,
    cpus: int = 0,
    gap_open: int = 11,
    gap_extend: int = 1,
    scoring_matrix: ScoringMatrix = ...,
) -> Clustering: ...
