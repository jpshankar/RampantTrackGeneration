from dataclasses import dataclass

from uuid import uuid4
from voronout.Point import Point

@dataclass(frozen=True)
class GraphNode:
    nodeId: uuid4
    nodePoint: Point