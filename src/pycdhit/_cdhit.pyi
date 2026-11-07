import array
import dataclasses
import typing

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
    cpus: typing.Optional[int] = None,
) -> Clustering: ...
