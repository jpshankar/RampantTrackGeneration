from dataclasses import dataclass

from .data import NodeInfo

from .edges.data import EdgeVertexInfo

from uuid import uuid4

from voronout.Point import Point

@dataclass(frozen=True)
class Track:
    nodes: dict[uuid4, Point]
    stops: dict[uuid4, tuple[Point]]
    edges: dict[uuid4, EdgeVertexInfo]

    startNodeId: uuid4
    destinationNodeId: uuid4
    
    nodeInfo: dict[uuid4, NodeInfo]