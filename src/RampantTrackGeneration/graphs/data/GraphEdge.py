from dataclasses import dataclass

from ...edges.data.EdgeVertexInfo import EdgeVertexInfo

from uuid import uuid4

@dataclass(frozen=True)
class GraphEdge:
    edgeId: uuid4
    edgeVertices: EdgeVertexInfo
    
    edgeLength: float