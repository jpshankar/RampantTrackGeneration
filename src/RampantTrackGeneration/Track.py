from dataclasses import dataclass

from .edges.data import EdgeVertexInfo

from uuid import uuid4

from voronout.Point import Point

@dataclass(frozen=True)
class Track:
    nodes: dict[uuid4, Point]
    stops: dict[EdgeVertexInfo, tuple[Point]]
    edges: tuple[EdgeVertexInfo]