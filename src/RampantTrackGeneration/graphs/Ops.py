from .data import GraphEdge, GraphNode
from ..edges.data import EdgeVertexInfo

import rustworkx
from rustworkx import PyGraph as Graph

from uuid import uuid4

from voronout.Point import Point

class Ops:
    @staticmethod
    def vertexIdToGraphNodeInd(graph: Graph, vertexId: uuid4) -> int:
        filterOnVertexId = lambda gn: gn.nodeId == vertexId
        filteredNodes = graph.filter_nodes(filter_function = filterOnVertexId)

        return filteredNodes[0] if filteredNodes else None
        
    @staticmethod
    def vertexIdToGraphNode(graph: Graph, vertexId: uuid4) -> GraphNode:
        graphNodeInd = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = vertexId)
        return graph.get_node_data(node = graphNodeInd) if graphNodeInd != None else None
        
    @staticmethod
    def graphVertex(graph: Graph, vertexId: uuid4) -> Point:
        graphNode = Ops.vertexIdToGraphNode(graph = graph, vertexId = vertexId)
        return graphNode.nodePoint if graphNode else None
    
    @staticmethod
    def graphEdge(graph: Graph, edgeVertexInfo: EdgeVertexInfo) -> GraphEdge:
        vertex0Ind = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = edgeVertexInfo.vertex0Id)
        vertex1Ind = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = edgeVertexInfo.vertex1Id)

        vertex0IndValid = vertex0Ind == 0 or bool(vertex0Ind)
        vertex1IndValid = vertex1Ind == 0 or bool(vertex1Ind)

        if vertex0IndValid and vertex1IndValid:
            graphEdgeIndices = graph.edge_indices_from_endpoints(node_a = vertex0Ind, node_b = vertex1Ind)
            if graphEdgeIndices:
                graphEdgeInd = graphEdgeIndices[0]
                return graph.get_edge_data_by_index(edge_index = graphEdgeInd)
            else:
                return None
        else:
            return None
    
    @staticmethod
    def graphVertexNeighbors(graph: Graph, vertexId: uuid4) -> tuple[uuid4]:
        vertexNodeInd = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = vertexId)
        if vertexNodeInd != None:
            vertexNeighborInds = graph.neighbors(node = vertexNodeInd)
            return tuple((graph.get_node_data(node = vertexNeighborInd).nodeId for vertexNeighborInd in vertexNeighborInds))
        else:
            return tuple()
        
    @staticmethod
    def graphEdgeLength(graph: Graph, edge: EdgeVertexInfo) -> float:
        graphEdge = Ops.graphEdge(graph = graph, edgeVertexInfo = edge)
        return graphEdge.edgeLength if graphEdge else None
        
    # Scaling p1 and p2 back down to Voronout's 0 -> 1 scale, to maintain congruity.
    @staticmethod
    def scaledGraphPointDistance(graph: Graph, p1: Point, p2: Point):
        graphWidthScalar = 1 / graph.attrs["width"]
        graphHeightScalar = 1 / graph.attrs["height"]

        scaledP1 = p1.scale(widthScalar = graphWidthScalar, heightScalar = graphHeightScalar)
        scaledP2 = p2.scale(widthScalar = graphWidthScalar, heightScalar = graphHeightScalar)

        return Point.distance(p1 = scaledP1, p2 = scaledP2)
    
    @staticmethod
    def addVertexToGraph(graph: Graph, vertexNode: GraphNode):
        if not Ops.vertexIdToGraphNode(graph = graph, vertexId = vertexNode.nodeId):
            graph.add_node(obj = vertexNode)

    @staticmethod
    def addConnectionToGraph(graph: Graph, connectionEdge: GraphEdge):
        if not Ops.graphEdge(graph = graph, edgeVertexInfo = connectionEdge.edgeVertices):
            vertex0Ind = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = connectionEdge.edgeVertices.vertex0Id)
            vertex1Ind = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = connectionEdge.edgeVertices.vertex1Id)

            graph.add_edge(node_a = vertex0Ind, node_b = vertex1Ind, edge = connectionEdge)
    
    @staticmethod
    def cleanUpNodeIfPossible(graph: Graph, nodeInd: int):
        if not tuple(graph.neighbors(node = nodeInd)):
            graph.remove_node(node = nodeInd)

    @staticmethod
    def removeEdgeAndCleanUpNodes(graph: Graph, existingEdge: EdgeVertexInfo):
        existingEdgeVertex0Id = existingEdge.vertex0Id
        existingEdgeVertex1Id = existingEdge.vertex1Id

        vertex0Ind = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = existingEdgeVertex0Id)
        vertex1Ind = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = existingEdgeVertex1Id)

        graph.remove_edge(node_a = vertex0Ind, node_b = vertex1Ind)
        
        Ops.cleanUpNodeIfPossible(graph = graph, nodeInd = vertex0Ind)
        Ops.cleanUpNodeIfPossible(graph = graph, nodeInd = vertex1Ind)

    @staticmethod
    def getCurrentEdgeNodeConnectivities(graph: Graph, edgeNode0Ind: int, edgeNode1Ind: int) -> tuple[tuple[uuid4]]:
        currentEdgeNodeConnectivities = rustworkx.all_pairs_all_simple_paths(graph)

        vertex0ConnectivityBeforeDeletionInds = tuple(currentEdgeNodeConnectivities[edgeNode0Ind])
        vertex0ConnectivityBeforeDeletionIds = tuple((graph.get_node_data(node = vertex0ConnectivityBeforeDeletionInd).nodeId for vertex0ConnectivityBeforeDeletionInd in vertex0ConnectivityBeforeDeletionInds))

        vertex1ConnectivityBeforeDeletionInds = tuple(currentEdgeNodeConnectivities[edgeNode1Ind])
        vertex1ConnectivityBeforeDeletionIds = tuple((graph.get_node_data(node = vertex1ConnectivityBeforeDeletionInd).nodeId for vertex1ConnectivityBeforeDeletionInd in vertex1ConnectivityBeforeDeletionInds))

        return (vertex0ConnectivityBeforeDeletionIds, vertex1ConnectivityBeforeDeletionIds)
    
    @staticmethod
    def getDisconnectedByEdgeRemoval(graph: Graph, edgeVertex0Id: uuid4, edgeVertex1Id: uuid4, edgeNode0BeforeConnectivities: tuple[uuid4], edgeNode1BeforeConnectivities: tuple[uuid4]) -> tuple[tuple[uuid4]]:
        afterVertex0NodeInd = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = edgeVertex0Id)
        afterVertex1NodeInd = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = edgeVertex1Id)

        afterExistingConnectivities = rustworkx.all_pairs_all_simple_paths(graph)

        vertex0ConnectivityAfterDeletionInds = tuple(afterExistingConnectivities[afterVertex0NodeInd]) if afterVertex0NodeInd and afterVertex0NodeInd in afterExistingConnectivities else tuple()
        vertex0ConnectivityAfterDeletionIds = tuple((graph.get_node_data(vertex0ConnectivityAfterDeletionInd).nodeId for vertex0ConnectivityAfterDeletionInd in vertex0ConnectivityAfterDeletionInds))
        
        vertex1ConnectivityAfterDeletionInds = tuple(afterExistingConnectivities[afterVertex1NodeInd]) if afterVertex1NodeInd and afterVertex1NodeInd in afterExistingConnectivities else tuple()
        vertex1ConnectivityAfterDeletionIds = tuple((graph.get_node_data(vertex1ConnectivityAfterDeletionInd).nodeId for vertex1ConnectivityAfterDeletionInd in vertex1ConnectivityAfterDeletionInds))
    
        disconnectedByDeletionVertex0 = tuple((vertex0ConnectivityBeforeDeletionId for vertex0ConnectivityBeforeDeletionId in edgeNode0BeforeConnectivities if vertex0ConnectivityBeforeDeletionId not in vertex0ConnectivityAfterDeletionIds and vertex0ConnectivityBeforeDeletionId != edgeVertex1Id))
        disconnectedByDeletionVertex1 = tuple((vertex1ConnectivityBeforeDeletionId for vertex1ConnectivityBeforeDeletionId in edgeNode1BeforeConnectivities if vertex1ConnectivityBeforeDeletionId not in vertex1ConnectivityAfterDeletionIds and vertex1ConnectivityBeforeDeletionId != edgeVertex0Id))

        return (disconnectedByDeletionVertex0, disconnectedByDeletionVertex1)

    @staticmethod
    def removeEdgeAndReturnDisconnected(graph: Graph, existingEdge: EdgeVertexInfo) -> tuple[tuple[uuid4]]:
        existingEdgeVertex0Id = existingEdge.vertex0Id
        existingEdgeVertex1Id = existingEdge.vertex1Id

        beforeVertex0NodeInd = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = existingEdgeVertex0Id)
        beforeVertex1NodeInd = Ops.vertexIdToGraphNodeInd(graph = graph, vertexId = existingEdgeVertex1Id)

        (vertex0ConnectivityBeforeDeletionIds, vertex1ConnectivityBeforeDeletionIds) = Ops.getCurrentEdgeNodeConnectivities(graph = graph, edgeNode0Ind = beforeVertex0NodeInd, edgeNode1Ind = beforeVertex1NodeInd)

        Ops.removeEdgeAndCleanUpNodes(graph = graph, existingEdge = existingEdge)
        return Ops.getDisconnectedByEdgeRemoval(graph = graph, edgeVertex0Id = existingEdgeVertex0Id, edgeVertex1Id = existingEdgeVertex1Id, edgeNode0BeforeConnectivities = vertex0ConnectivityBeforeDeletionIds, edgeNode1BeforeConnectivities = vertex1ConnectivityBeforeDeletionIds)