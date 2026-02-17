from .edges.data import EdgesMakingAngle, EdgeVertexInfo

from .edges.Ops import Ops as EdgesOps

from .graphs.data import GraphEdge, GraphNode
from .graphs.Ops import Ops as GraphOps

from math import floor

from rustworkx import articulation_points, PyGraph as Graph

from voronout.Point import Point
from voronout.VoronoiDiagram import VoronoiDiagram

from .Track import Track

import numpy
import random

from uuid import uuid4

class TrackGenerator:
    @staticmethod
    def _edgeAngleViability(angle: float, minAcceptableAngle: float) -> bool:
        # Neither angle nor its complement should be < minAcceptableAngle.
        return angle >= minAcceptableAngle and (180 - angle) >= minAcceptableAngle
    
    @staticmethod
    def _calculateEdgeAnglesWithGraphVertexEdges(graph: Graph, edge: EdgeVertexInfo, vertexId: uuid4) -> dict[EdgesMakingAngle, float]:
        vertexNeighborIds = GraphOps.graphVertexNeighbors(graph = graph, vertexId = vertexId)
        vertexEdgesMakingAngles = []

        for vertexNeighborId in vertexNeighborIds:
            otherEdge = EdgeVertexInfo(vertex0Id = vertexId, vertex1Id = vertexNeighborId)
            edgesMakingAngle = EdgesMakingAngle(edge0 = edge, edge1 = otherEdge)
            vertexEdgesMakingAngles.append(edgesMakingAngle)

        return EdgesOps.calculateEdgeAnglesWithExtantVertexEdges(edgeGraph = graph, edgesMakingAngles = tuple(vertexEdgesMakingAngles))

    @staticmethod
    def _findPossibleViableConnection(graph: Graph, possibleConnections: tuple[EdgeVertexInfo], minAcceptableAngle: float) -> EdgeVertexInfo:
        possibleViableConnections = []
        for possibleConnection in possibleConnections:
            zeroConnections = TrackGenerator._calculateEdgeAnglesWithGraphVertexEdges(graph = graph, edge = possibleConnection, vertexId = possibleConnection.vertex0Id)
            oneConnections = TrackGenerator._calculateEdgeAnglesWithGraphVertexEdges(graph = graph, edge = possibleConnection, vertexId = possibleConnection.vertex1Id)

            allZeroConnectionsValid = all((TrackGenerator._edgeAngleViability(angle = angle, minAcceptableAngle = minAcceptableAngle) for angle in zeroConnections.values()))
            allOneConnectionsValid = all((TrackGenerator._edgeAngleViability(angle = angle, minAcceptableAngle = minAcceptableAngle) for angle in oneConnections.values()))

            if allZeroConnectionsValid and allOneConnectionsValid:
                possibleViableConnections.append(possibleConnection)

        return random.choice(possibleViableConnections) if possibleViableConnections else None
    
    @staticmethod
    def _updateVertexPotentialConnections(existingConnectionsGraph: Graph, potentialConnectionsGraph: Graph, vertexToUpdateId: uuid4): 
        verticesConnectedTo = GraphOps.graphVertexNeighbors(graph = existingConnectionsGraph, vertexId = vertexToUpdateId)
        
        potentialConnectingVertexNodeIds = tuple((maybePotentialConnectingVertex.nodeId for maybePotentialConnectingVertex in existingConnectionsGraph.nodes()))
        potentialConnectingVertexIds = tuple((potentialConnectingVertexId for potentialConnectingVertexId in potentialConnectingVertexNodeIds if potentialConnectingVertexId not in verticesConnectedTo and potentialConnectingVertexId != vertexToUpdateId))

        vertexToUpdateNode = GraphOps.vertexIdToGraphNode(graph = existingConnectionsGraph, vertexId = vertexToUpdateId)
        GraphOps.addVertexToGraph(graph = potentialConnectionsGraph, vertexNode = vertexToUpdateNode)
        
        for potentialConnectingVertexId in potentialConnectingVertexIds:
            potentialConnectionVertex = GraphOps.vertexIdToGraphNode(graph = existingConnectionsGraph, vertexId = potentialConnectingVertexId)

            GraphOps.addVertexToGraph(graph = potentialConnectionsGraph, vertexNode = potentialConnectionVertex)

            potentialConnectionVertexInfo = EdgeVertexInfo(vertex0Id = vertexToUpdateId, vertex1Id = potentialConnectingVertexId)
            
            potentialConnectionP1 = vertexToUpdateNode.nodePoint
            potentialConnectionP2 = potentialConnectionVertex.nodePoint

            potentialConnectionDistance = GraphOps.scaledGraphPointDistance(graph = potentialConnectionsGraph, p1 = potentialConnectionP1, p2 = potentialConnectionP2)

            potentialConnection = GraphEdge(edgeId = uuid4(), edgeVertices = potentialConnectionVertexInfo, edgeLength = potentialConnectionDistance)

            GraphOps.addConnectionToGraph(graph = potentialConnectionsGraph, connectionEdge = potentialConnection)

    @staticmethod
    def _pruneIntersectionEdgesFromExistingConnectionsGraph(existingConnectionsGraph: Graph, intersectionEdges: tuple[EdgeVertexInfo]):
        nodesThatWouldDisconnect = set((existingConnectionsGraph.get_node_data(node = articulationPointInd) for articulationPointInd in articulation_points(graph = existingConnectionsGraph)))

        safeToRemoveEdges = tuple((intersectionEdge for intersectionEdge in intersectionEdges if intersectionEdge.vertex0Id not in nodesThatWouldDisconnect and intersectionEdge.vertex1Id not in nodesThatWouldDisconnect))

        for safeToRemoveEdge in safeToRemoveEdges:
            safeToRemoveGraphEdge = GraphOps.graphEdge(graph = existingConnectionsGraph, edgeVertexInfo = safeToRemoveEdge)

            safeToRemoveEdgeVertex0Id = safeToRemoveEdge.vertex0Id
            safeToRemoveEdgeVertex1Id = safeToRemoveEdge.vertex1Id

            safeToRemoveVertex0Ind = GraphOps.vertexIdToGraphNodeInd(graph = existingConnectionsGraph, vertexId = safeToRemoveEdgeVertex0Id)
            safeToRemoveVertex1Ind = GraphOps.vertexIdToGraphNodeInd(graph = existingConnectionsGraph, vertexId = safeToRemoveEdgeVertex1Id)

            (vertex0ConnectivityBeforeDeletionIds, vertex1ConnectivityBeforeDeletionIds) = GraphOps.getCurrentEdgeNodeConnectivities(graph = existingConnectionsGraph, edgeNode0Ind = safeToRemoveVertex0Ind, edgeNode1Ind = safeToRemoveVertex1Ind)
            
            existingConnectionsGraph.remove_edge(node_a = safeToRemoveVertex0Ind, node_b = safeToRemoveVertex1Ind)

            (disconnectedVertex0, disconnectedVertex1) = GraphOps.getDisconnectedByEdgeRemoval(graph = existingConnectionsGraph, edgeVertex0Id = safeToRemoveEdgeVertex0Id, edgeVertex1Id = safeToRemoveEdgeVertex1Id, edgeNode0BeforeConnectivities = vertex0ConnectivityBeforeDeletionIds, edgeNode1BeforeConnectivities = vertex1ConnectivityBeforeDeletionIds)

            if len(disconnectedVertex0) > 0 or len(disconnectedVertex1) > 0:
                GraphOps.addConnectionToGraph(graph = existingConnectionsGraph, connectionEdge = safeToRemoveGraphEdge)
            else:
                GraphOps.cleanUpNodeIfPossible(graph = existingConnectionsGraph, nodeInd = safeToRemoveVertex0Ind)
                GraphOps.cleanUpNodeIfPossible(graph = existingConnectionsGraph, nodeInd = safeToRemoveVertex1Ind)

    @staticmethod
    def _handleLonelyExistingConnections(existingConnectionsGraph: Graph, lonelyConnectionMinLengthQuantile: float):
        existingConnectionEdgeLengths = tuple((existingConnectionEdge.edgeLength for existingConnectionEdge in existingConnectionsGraph.edges()))
        existingConnectionEdgeMinQuantile = numpy.quantile(a = existingConnectionEdgeLengths, q = lonelyConnectionMinLengthQuantile)

        maybeHandleableExistingConnections = tuple((existingConnectionEdge for existingConnectionEdge in existingConnectionsGraph.edges() if existingConnectionEdge.edgeLength < existingConnectionEdgeMinQuantile))

        for maybeHandleableExistingConnection in maybeHandleableExistingConnections:
            existingConnectionEdgeInfo = maybeHandleableExistingConnection.edgeVertices

            existingEdgeVertex0Id = existingConnectionEdgeInfo.vertex0Id
            existingEdgeVertex1Id = existingConnectionEdgeInfo.vertex1Id

            vertex0Ind = GraphOps.vertexIdToGraphNodeInd(graph = existingConnectionsGraph, vertexId = existingEdgeVertex0Id)
            vertex1Ind = GraphOps.vertexIdToGraphNodeInd(graph = existingConnectionsGraph, vertexId = existingEdgeVertex1Id)

            zeroNeighbors = tuple(existingConnectionsGraph.neighbors(vertex0Ind))
            oneNeighbors = tuple(existingConnectionsGraph.neighbors(vertex1Ind))

            if not len(zeroNeighbors) > 1 or not len(oneNeighbors) > 1:
                GraphOps.removeEdgeAndCleanUpNodes(graph = existingConnectionsGraph, existingEdge = existingConnectionEdgeInfo)

    # https://math.stackexchange.com/questions/134112/find-a-point-on-a-line-segment-located-at-a-distance-d-from-one-endpoint
    @staticmethod
    def _generateStopsOnConnection(connectionGraph: Graph, connection: EdgeVertexInfo, connectionLengthVertexPadding: float, connectionLengthNodeBuffer: float, minEdgeLengthForStopsGeneration: float) -> tuple[Point]:
        connectionLength = GraphOps.graphEdgeLength(graph = connectionGraph, edge = connection)
        
        # Only generate stops if length is at least larger than ((connectionLengthVertexPadding + connectionLengthNodeBuffer) * 100) percent of edges.
        if (connectionLength >= minEdgeLengthForStopsGeneration):
            connectionVertexZeroId = connection.vertex0Id
            connectionVertexZero = GraphOps.graphVertex(graph = connectionGraph, vertexId = connectionVertexZeroId)

            connectionVertexOneId = connection.vertex1Id
            connectionVertexOne = GraphOps.graphVertex(graph = connectionGraph, vertexId = connectionVertexOneId)

            zeroIsFirstVertex = connectionVertexZero.x < connectionVertexOne.x

            firstVertex = connectionVertexZero if zeroIsFirstVertex else connectionVertexOne
            secondVertex = connectionVertexZero if not zeroIsFirstVertex else connectionVertexOne

            connectionXRangeMinInterval = (connectionLengthVertexPadding + connectionLengthNodeBuffer) * connectionLength
            connectionXRangeMaxInterval = connectionLength - connectionXRangeMinInterval

            connectionNodes = []

            while connectionXRangeMinInterval < connectionXRangeMaxInterval:
                d = random.uniform(connectionXRangeMinInterval, connectionXRangeMaxInterval)        
                td = d / connectionLength

                nextX = firstVertex.x + ((secondVertex.x - firstVertex.x) * td)
                nextY = firstVertex.y + ((secondVertex.y - firstVertex.y) * td)
                
                nextPoint = Point(x = nextX, y = nextY)
                connectionNodes.append(nextPoint)

                connectionXRangeMinInterval = d + (connectionLengthNodeBuffer * connectionLength)
            
            return tuple(connectionNodes)
        else:
            return tuple()
                
    @staticmethod
    def generateTrack(
        diagramWidth: int,
        diagramHeight: int,
        numWalkersOnTrack: int,
        numDestinationsOnTrack: int,
        diagramEdgePercentageToProcess: float, 
        newConnectionAngleMinQuantile: float, 
        lonelyConnectionMinLengthQuantile: float, 
        connectionLengthVertexPadding: float, 
        connectionLengthNodeBuffer: float
    ) -> Track:
        numVoronoiDiagramRegions = (numWalkersOnTrack * numDestinationsOnTrack) * 2

        diagramRegionSites = tuple((Point(x = random.random(), y = random.random()) for _ in range(numVoronoiDiagramRegions)))

        voronoiDiagram = VoronoiDiagram(basePoints = diagramRegionSites, planeWidth = diagramWidth, planeHeight = diagramHeight)

        voronoiDiagramEdgesToProcess = {}
    
        for voronoiDiagramEdge in voronoiDiagram.diagramEdges:
            edgeId = voronoiDiagramEdge.edgeId
            edge = EdgeVertexInfo(vertex0Id = voronoiDiagramEdge.vertex0Id, vertex1Id = voronoiDiagramEdge.vertex1Id)

            if edgeId not in voronoiDiagramEdgesToProcess:
                voronoiDiagramEdgesToProcess[edgeId] = edge

        voronoiVertices = voronoiDiagram.vertices

        # existingConnections models the current state of the diagram -> Track transformation.
        existingConnections = Graph(attrs={"width": diagramWidth, "height": diagramHeight})

        # potentialConnections models all possible reconnections.
        potentialConnections = Graph(attrs={"width": diagramWidth, "height": diagramHeight})

        for existingConnectionNode in voronoiVertices.keys():
            edgesWithVertex = { voronoiDiagramEdgeId: voronoiDiagramEdge for (voronoiDiagramEdgeId, voronoiDiagramEdge) in voronoiDiagramEdgesToProcess.items() if voronoiDiagramEdge.vertex0Id == existingConnectionNode or voronoiDiagramEdge.vertex1Id == existingConnectionNode}
            
            for (edgeId, edge) in edgesWithVertex.items():                
                pointEdgeVertex0Id = edge.vertex0Id
                pointEdgeVertex0 = voronoiVertices[pointEdgeVertex0Id]

                pointEdgeVertex1Id = edge.vertex1Id
                pointEdgeVertex1 = voronoiVertices[pointEdgeVertex1Id]

                vertex0Node = GraphNode(nodeId = pointEdgeVertex0Id, nodePoint = pointEdgeVertex0)
                
                GraphOps.addVertexToGraph(graph = existingConnections, vertexNode = vertex0Node)

                vertex1Node = GraphNode(nodeId = pointEdgeVertex1Id, nodePoint = pointEdgeVertex1)
                GraphOps.addVertexToGraph(graph = existingConnections, vertexNode = vertex1Node)

                pointEdgeVertexDistance = GraphOps.scaledGraphPointDistance(graph = existingConnections, p1 = pointEdgeVertex0, p2 = pointEdgeVertex1)
                pointEdge = GraphEdge(edgeId = edgeId, edgeVertices = edge, edgeLength = pointEdgeVertexDistance)

                GraphOps.addConnectionToGraph(graph = existingConnections, connectionEdge = pointEdge)

        existingEdgeVertexAngles = {}

        for existingConnection in existingConnections.edges():
            existingEdge = existingConnection.edgeVertices

            vertex0EdgeAngles = TrackGenerator._calculateEdgeAnglesWithGraphVertexEdges(graph = existingConnections, edge = existingEdge, vertexId = existingEdge.vertex0Id)
            for vertex0AngleEdges in vertex0EdgeAngles:
                if vertex0AngleEdges not in existingEdgeVertexAngles:
                    existingEdgeVertexAngles[vertex0AngleEdges] = vertex0EdgeAngles[vertex0AngleEdges]

            vertex1EdgeAngles = TrackGenerator._calculateEdgeAnglesWithGraphVertexEdges(graph = existingConnections, edge = existingEdge, vertexId = existingEdge.vertex1Id)
            for vertex1AngleEdges in vertex1EdgeAngles:
                if vertex1AngleEdges not in existingEdgeVertexAngles:
                    existingEdgeVertexAngles[vertex1AngleEdges] = vertex1EdgeAngles[vertex1AngleEdges]

        initialDiagramEdgeAngles = tuple(existingEdgeVertexAngles.values())
        initialDiagramMinAcceptableAngle = numpy.quantile(a = initialDiagramEdgeAngles, q = newConnectionAngleMinQuantile)

        for existingConnectionNode in existingConnections.nodes():
            TrackGenerator._updateVertexPotentialConnections(existingConnectionsGraph = existingConnections, potentialConnectionsGraph = potentialConnections, vertexToUpdateId = existingConnectionNode.nodeId)

        diagramEdgesSorted = sorted(tuple(voronoiDiagramEdgesToProcess.values()), key = lambda vde: GraphOps.graphEdgeLength(graph = existingConnections, edge = vde))
        numDiagramEdgesToProcess = floor(len(diagramEdgesSorted) * diagramEdgePercentageToProcess)

        diagramEdgesToProcess = diagramEdgesSorted[:numDiagramEdgesToProcess]

        intersectionEdges = []

        for diagramEdgeToProcess in diagramEdgesToProcess:
            existingConnectionsEdge = GraphOps.graphEdge(graph = existingConnections, edgeVertexInfo = diagramEdgeToProcess)

            if existingConnectionsEdge:
                edgeToProcessLength = existingConnectionsEdge.edgeLength      
                (disconnectedByDeletionVertex0, disconnectedByDeletionVertex1) = GraphOps.removeEdgeAndReturnDisconnected(graph = existingConnections, existingEdge = diagramEdgeToProcess)
                
                possibleConnections = []

                for vertexDisconnected0 in disconnectedByDeletionVertex0:
                    for vertexDisconnected1 in disconnectedByDeletionVertex1:
                        potentialConnection = EdgeVertexInfo(vertex0Id = vertexDisconnected0, vertex1Id = vertexDisconnected1)
                        # Do not add zero-length (vertex0 == vertex1) edges or edges that already exist.
                        if vertexDisconnected0 != vertexDisconnected1 and potentialConnection not in diagramEdgesToProcess:                            
                            connectionAlreadyAdded = GraphOps.graphEdge(graph = existingConnections, edgeVertexInfo = potentialConnection)

                            if not connectionAlreadyAdded:
                                potentialConnectionEdge = GraphOps.graphEdge(graph = potentialConnections, edgeVertexInfo = potentialConnection)
                                # Ensure that an edge that isn't currently a potentialConnection will fail the length check.
                                potentialConnectionEdgeLength = potentialConnectionEdge.edgeLength if potentialConnectionEdge else 0.0

                                # The reconnection should be longer than the deleted edge.
                                if potentialConnectionEdgeLength > edgeToProcessLength:
                                    possibleConnections.append(potentialConnection)

                maybeViableConnection = TrackGenerator._findPossibleViableConnection(graph = existingConnections, possibleConnections = tuple(possibleConnections), minAcceptableAngle = initialDiagramMinAcceptableAngle)

                if maybeViableConnection:
                    maybeViableGraphEdge = GraphOps.graphEdge(graph = potentialConnections, edgeVertexInfo = maybeViableConnection)

                    maybeViableEdgeId = maybeViableGraphEdge.edgeId
                    maybeViableVertexLength = maybeViableGraphEdge.edgeLength

                    maybeViableVertex0Id = maybeViableConnection.vertex0Id
                    maybeViableVertex1Id = maybeViableConnection.vertex1Id

                    GraphOps.removeEdgeAndCleanUpNodes(graph = potentialConnections, existingEdge = maybeViableConnection)
                    
                    potentialIntersections = tuple((existingConnectionEdge.edgeVertices for existingConnectionEdge in existingConnections.edges() if (existingConnectionEdge.edgeVertices.vertex0Id != maybeViableVertex0Id and existingConnectionEdge.edgeVertices.vertex0Id != maybeViableVertex1Id) and (existingConnectionEdge.edgeVertices.vertex1Id != maybeViableVertex0Id and existingConnectionEdge.edgeVertices.vertex1Id != maybeViableVertex1Id)))
                    edgeIntersectionsInfo = EdgesOps.calculateEdgesIntersectionInfo(edgeGraph = existingConnections, intersectingEdge = maybeViableConnection, otherEdges = potentialIntersections)
                                            
                    if edgeIntersectionsInfo:
                        # We expect this to be sorted along the v0 -> v1 path.
                        precedingVertexId = maybeViableVertex0Id
                        precedingVertex = GraphOps.graphVertex(graph = existingConnections, vertexId = maybeViableVertex0Id)

                        for edgeIntersectionInfo in edgeIntersectionsInfo:
                            edgeIntersectionPoint = edgeIntersectionInfo.intersectionPoint
                            distanceToIntersection = GraphOps.scaledGraphPointDistance(graph = existingConnections, p1 = precedingVertex, p2 = edgeIntersectionPoint)

                            edgeIntersectionPointId = uuid4()
                            edgeIntersectionNode = GraphNode(nodeId = edgeIntersectionPointId, nodePoint = edgeIntersectionPoint)

                            precedingIntersectionVertices = EdgeVertexInfo(vertex0Id = precedingVertexId, vertex1Id = edgeIntersectionPointId)
                            precedingIntersectionEdge = GraphEdge(edgeId = uuid4(), edgeVertices = precedingIntersectionVertices, edgeLength = distanceToIntersection)

                            GraphOps.addVertexToGraph(graph = existingConnections, vertexNode = edgeIntersectionNode)
                            GraphOps.addConnectionToGraph(graph = existingConnections, connectionEdge = precedingIntersectionEdge)
 
                            intersectionEdges.append(precedingIntersectionVertices)
                            edgeIntersected = edgeIntersectionInfo.intersectedEdge

                            edgeIntersectedVertex0Id = edgeIntersected.vertex0Id
                            edgeIntersectedVertex1Id = edgeIntersected.vertex1Id
                                                        
                            edgeIntersectedVertex0 = GraphOps.graphVertex(graph = existingConnections, vertexId = edgeIntersectedVertex0Id)
                            edgeIntersectedVertex1 = GraphOps.graphVertex(graph = existingConnections, vertexId = edgeIntersectedVertex1Id)

                            iv0DistanceToIntersection = GraphOps.scaledGraphPointDistance(graph = existingConnections, p1 = edgeIntersectedVertex0, p2 = edgeIntersectionPoint)
                            iv1DistanceToIntersection = GraphOps.scaledGraphPointDistance(graph = existingConnections, p1 = edgeIntersectionPoint, p2 = edgeIntersectedVertex1)

                            edgeIntersectedVertex0Ind = GraphOps.vertexIdToGraphNodeInd(graph = existingConnections, vertexId = edgeIntersectedVertex0Id)
                            edgeIntersectedVertex1Ind = GraphOps.vertexIdToGraphNodeInd(graph = existingConnections, vertexId = edgeIntersectedVertex1Id)

                            # No need to remove " orphaned " nodes here, since they'll subsequently be re-added with the handling of intersection edges.
                            existingConnections.remove_edge(node_a = edgeIntersectedVertex0Ind, node_b = edgeIntersectedVertex1Ind)
                            
                            if edgeIntersected in intersectionEdges:
                                intersectionEdges.remove(edgeIntersected)

                            intersectionEdge0 = EdgeVertexInfo(vertex0Id = edgeIntersectedVertex0Id, vertex1Id = edgeIntersectionPointId)
                            intersectionEdge0Graph = GraphEdge(edgeId = uuid4(), edgeVertices = intersectionEdge0, edgeLength = iv0DistanceToIntersection)

                            GraphOps.addConnectionToGraph(graph = existingConnections, connectionEdge = intersectionEdge0Graph)                            
                            intersectionEdges.append(intersectionEdge0)

                            intersectionEdge1 = EdgeVertexInfo(vertex0Id = edgeIntersectionPointId, vertex1Id = edgeIntersectedVertex1Id)
                            intersectionEdge1Graph = GraphEdge(edgeId = uuid4(), edgeVertices = intersectionEdge1, edgeLength = iv1DistanceToIntersection)

                            GraphOps.addConnectionToGraph(graph = existingConnections, connectionEdge = intersectionEdge1Graph)
                            intersectionEdges.append(intersectionEdge1)

                            TrackGenerator._updateVertexPotentialConnections(existingConnectionsGraph = existingConnections, potentialConnectionsGraph = potentialConnections, vertexToUpdateId = edgeIntersectionPointId)

                            precedingVertexId = edgeIntersectionPointId
                            precedingVertex = edgeIntersectionPoint

                        # Finish up with (last intersection point) -> vertex1.
                        maybeViableVertex1 = GraphOps.graphVertex(graph = existingConnections, vertexId = maybeViableVertex1Id)
                        distanceToEnd = GraphOps.scaledGraphPointDistance(graph = existingConnections, p1 = precedingVertex, p2 = maybeViableVertex1)

                        toEndEdge = EdgeVertexInfo(vertex0Id = precedingVertexId, vertex1Id = maybeViableVertex1Id)
                        toEndGraphEdge = GraphEdge(edgeId = uuid4(), edgeVertices = toEndEdge, edgeLength = distanceToEnd)

                        GraphOps.addConnectionToGraph(graph = existingConnections, connectionEdge = toEndGraphEdge)

                        intersectionEdges.append(toEndEdge)
                    else:
                        viableVertices = EdgeVertexInfo(vertex0Id = maybeViableVertex0Id, vertex1Id = maybeViableVertex1Id)
                        viableGraphEdge = GraphEdge(edgeId = maybeViableEdgeId, edgeVertices = viableVertices, edgeLength = maybeViableVertexLength)

                        GraphOps.addConnectionToGraph(graph = existingConnections, connectionEdge = viableGraphEdge)

        TrackGenerator._pruneIntersectionEdgesFromExistingConnectionsGraph(existingConnectionsGraph = existingConnections, intersectionEdges = intersectionEdges)
        TrackGenerator._handleLonelyExistingConnections(existingConnectionsGraph = existingConnections, lonelyConnectionMinLengthQuantile = lonelyConnectionMinLengthQuantile)

        existingConnectionEdges = existingConnections.edges()

        edgeLengths = tuple((existingConnectionEdge.edgeLength for existingConnectionEdge in existingConnectionEdges))
        minEdgeStopLengthQuantile = connectionLengthVertexPadding + connectionLengthNodeBuffer

        minEdgeStopLength = numpy.quantile(a = edgeLengths, q = minEdgeStopLengthQuantile)
        
        stops = { existingConnectionEdge.edgeId: TrackGenerator._generateStopsOnConnection(connectionGraph = existingConnections, connection = existingConnectionEdge.edgeVertices, connectionLengthVertexPadding = connectionLengthVertexPadding, connectionLengthNodeBuffer = connectionLengthNodeBuffer, minEdgeLengthForStopsGeneration = minEdgeStopLength) for existingConnectionEdge in existingConnectionEdges}

        nodes = { existingConnectionNode.nodeId: GraphOps.graphVertex(graph = existingConnections, vertexId = existingConnectionNode.nodeId) for existingConnectionNode in existingConnections.nodes()}

        edges = { existingConnectionsEdge.edgeId: existingConnectionsEdge.edgeVertices for existingConnectionsEdge in existingConnectionEdges }

        return Track(nodes = nodes, stops = stops, edges = edges)
