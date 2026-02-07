from dataclasses import dataclass

from .EdgeVertexInfo import EdgeVertexInfo

from voronout.Point import Point

@dataclass(frozen=True)
class EdgeIntersectionInfo:
    intersectionPoint: Point
    intersectedEdge: EdgeVertexInfo