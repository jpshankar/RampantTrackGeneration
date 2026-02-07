from dataclasses import dataclass

from .EdgeVertexInfo import EdgeVertexInfo

@dataclass(frozen=True)
class EdgesMakingAngle:
    edge0: EdgeVertexInfo
    edge1: EdgeVertexInfo