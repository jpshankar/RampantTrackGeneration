from ..edges.data import EdgeVertexInfo

from networkx import Graph
from networkx.algorithms import approximation as nxApproximation

from uuid import uuid4

from voronout.Point import Point

class Ops:
    @staticmethod
    def graphVertex(graph: Graph, vertexId: uuid4) -> Point:
        if graph.has_node(n = vertexId):
            vertex = graph.nodes[vertexId]
            return Point(x = vertex["x"], y = vertex["y"])
        else:
            return None
        
    @staticmethod
    def graphEdgeLength(graph: Graph, edge: EdgeVertexInfo) -> float:
        return graph.edges[edge.vertex0Id, edge.vertex1Id]["edgeLength"]
        
    # Scaling p1 and p2 back down to Voronout's 0 -> 1 scale, to maintain congruity.
    @staticmethod
    def scaledGraphPointDistance(graph: Graph, p1: Point, p2: Point):
        graphWidthScalar = 1 / graph.graph["width"]
        graphHeightScalar = 1 / graph.graph["height"]

        scaledP1 = p1.scale(widthScalar = graphWidthScalar, heightScalar = graphHeightScalar)
        scaledP2 = p2.scale(widthScalar = graphWidthScalar, heightScalar = graphHeightScalar)

        return Point.distance(p1 = scaledP1, p2 = scaledP2)
    
    @staticmethod
    def addVertexToGraph(graph: Graph, vertexId: uuid4, vertexX: float, vertexY: float):
        if not graph.has_node(n = vertexId):
            graph.add_node(node_for_adding = vertexId, x = vertexX, y = vertexY)

    @staticmethod
    def addConnectionToGraph(graph: Graph, connection: EdgeVertexInfo, vertexDistance: float):
        vertex0Id = connection.vertex0Id
        vertex1Id = connection.vertex1Id

        if not graph.has_edge(u = vertex0Id, v = vertex1Id):
            graph.add_edge(u_of_edge = vertex0Id, v_of_edge = vertex1Id, edgeLength = vertexDistance)

    @staticmethod
    def removeEdgeAndCleanUpNodes(edgeGraph: Graph, existingEdge: EdgeVertexInfo):
        existingEdgeVertex0Id = existingEdge.vertex0Id
        existingEdgeVertex1Id = existingEdge.vertex1Id

        edgeGraph.remove_edge(u = existingEdgeVertex0Id, v = existingEdgeVertex1Id)

        if not bool(tuple(edgeGraph.neighbors(n = existingEdgeVertex0Id))):
            edgeGraph.remove_node(n = existingEdgeVertex0Id)

        if not bool(tuple(edgeGraph.neighbors(n = existingEdgeVertex1Id))):
            edgeGraph.remove_node(n = existingEdgeVertex1Id)

    @staticmethod
    def removeEdgeAndReturnDisconnected(edgeGraph: Graph, existingEdge: EdgeVertexInfo) -> tuple[tuple[EdgeVertexInfo]]:
        existingEdgeVertex0Id = existingEdge.vertex0Id
        existingEdgeVertex1Id = existingEdge.vertex1Id

        beforeExistingConnectivities = nxApproximation.all_pairs_node_connectivity(G = edgeGraph)

        vertex0ConnectivityBeforeDeletion = beforeExistingConnectivities[existingEdgeVertex0Id]
        vertex1ConnectivityBeforeDeletion = beforeExistingConnectivities[existingEdgeVertex1Id]

        Ops.removeEdgeAndCleanUpNodes(edgeGraph = edgeGraph, existingEdge = existingEdge)

        afterExistingConnectivities = nxApproximation.all_pairs_node_connectivity(G = edgeGraph)

        afterExistingConnectivitiesVertex0 = afterExistingConnectivities[existingEdgeVertex0Id] if existingEdgeVertex0Id in afterExistingConnectivities else {}
        afterExistingConnectivitiesVertex1 = afterExistingConnectivities[existingEdgeVertex1Id] if existingEdgeVertex1Id in afterExistingConnectivities else {}

        disconnectedByDeletionVertex0 = tuple((maybeDisconnectedVertex for (maybeDisconnectedVertex, numConnectionsLeft) in afterExistingConnectivitiesVertex0.items() if numConnectionsLeft == 0 and vertex0ConnectivityBeforeDeletion[maybeDisconnectedVertex] > 0))
        disconnectedByDeletionVertex1 = tuple((maybeDisconnectedVertex for (maybeDisconnectedVertex, numConnectionsLeft) in afterExistingConnectivitiesVertex1.items() if numConnectionsLeft == 0 and vertex1ConnectivityBeforeDeletion[maybeDisconnectedVertex] > 0))

        return (disconnectedByDeletionVertex0, disconnectedByDeletionVertex1)

