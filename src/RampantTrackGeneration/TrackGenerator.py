from .edges.data import EdgesMakingAngle, EdgeVertexInfo

from .edges.Ops import Ops as EdgesOps
from .graphs.Ops import Ops as GraphOps

from math import floor

from networkx import Graph, all_node_cuts as networkXAllNodeCuts

from voronout.Point import Point
from voronout.VoronoiDiagram import VoronoiDiagram

from .Track import Track

import numpy
import random

from uuid import uuid4

class TrackGenerator:
    @staticmethod
    def _edgeAngleViability(angle: float, minAcceptableAngle: float) -> bool:
        # neither angle nor its complement should be < minAcceptableAngle
        return angle >= minAcceptableAngle and (180 - angle) >= minAcceptableAngle
    
    @staticmethod
    def _calculateEdgeAnglesWithGraphVertexEdges(graph: Graph, edge: EdgeVertexInfo, vertexId: uuid4) -> dict[EdgesMakingAngle, float]:
        vertexNeighborIds = graph.neighbors(n = vertexId)
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
        verticesConnectedTo = existingConnectionsGraph.neighbors(n = vertexToUpdateId)
        potentialConnectingVertexIds = tuple((potentialConnectingVertexId for potentialConnectingVertexId in existingConnectionsGraph.nodes if potentialConnectingVertexId not in verticesConnectedTo and potentialConnectingVertexId != vertexToUpdateId))

        vertexToUpdate = existingConnectionsGraph.nodes[vertexToUpdateId]
        GraphOps.addVertexToGraph(graph = potentialConnectionsGraph, vertexId = vertexToUpdateId, vertexX = vertexToUpdate["x"], vertexY = vertexToUpdate["y"])
        
        for potentialConnectingVertexId in potentialConnectingVertexIds:
            potentialConnectionVertex = existingConnectionsGraph.nodes[potentialConnectingVertexId]

            GraphOps.addVertexToGraph(graph = existingConnectionsGraph, vertexId = potentialConnectingVertexId, vertexX = potentialConnectionVertex["x"], vertexY = potentialConnectionVertex["y"])

            potentialConnection = EdgeVertexInfo(vertex0Id = vertexToUpdateId, vertex1Id = potentialConnectingVertexId)
            
            potentialConnectionP1 = Point(x = vertexToUpdate["x"], y = vertexToUpdate["y"])
            potentialConnectionP2 = Point(x = potentialConnectionVertex["x"], y = potentialConnectionVertex["y"])

            potentialConnectionDistance = GraphOps.scaledGraphPointDistance(graph = potentialConnectionsGraph, p1 = potentialConnectionP1, p2 = potentialConnectionP2)

            GraphOps.addConnectionToGraph(graph = potentialConnectionsGraph, connection = potentialConnection, vertexDistance = potentialConnectionDistance)

    @staticmethod
    def _pruneIntersectionEdgesFromExistingConnectionsGraph(existingConnectionsGraph: Graph, intersectionEdges: tuple[EdgeVertexInfo]):
        existingConnectionsNodeCutSets = networkXAllNodeCuts(G = existingConnectionsGraph)
        nodesThatWouldDisconnect = set()

        for existingConnectionNodeCutSet in existingConnectionsNodeCutSets:
            for existingConnectionNodeCut in existingConnectionNodeCutSet:
                nodesThatWouldDisconnect.add(existingConnectionNodeCut)

        safeToRemoveEdges = tuple((intersectionEdge for intersectionEdge in intersectionEdges if intersectionEdge.vertex0Id not in nodesThatWouldDisconnect and intersectionEdge.vertex1Id not in nodesThatWouldDisconnect))

        for safeToRemoveEdge in safeToRemoveEdges:
            safeToRemoveEdgeLength = GraphOps.graphEdgeLength(graph = existingConnectionsGraph, edge = safeToRemoveEdge)
            (disconnectedVertex0, disconnectedVertex1) = GraphOps.removeEdgeAndReturnDisconnected(edgeGraph = existingConnectionsGraph, existingEdge = safeToRemoveEdge)
            if len(disconnectedVertex0) > 0 or len(disconnectedVertex1) > 0:
                # Need to disconnect/reconnect like this because refreshing existingConnectionsNodeCutSets after every removal is computationally too expensive
                existingConnectionsGraph.add_edge(u_of_edge = safeToRemoveEdge.vertex0Id, v_of_edge = safeToRemoveEdge.vertex1Id, edgeLength = safeToRemoveEdgeLength)

    @staticmethod
    def _handleLonelyExistingConnections(existingConnectionsGraph: Graph, lonelyConnectionMinLengthQuantile: float):
        existingConnectionsEdgeLengthMapping = {}
        for (ev0Id, ev1Id) in existingConnectionsGraph.edges:
            existingConnectionEdgeInfo = EdgeVertexInfo(vertex0Id = ev0Id, vertex1Id = ev1Id)
            existingConnectionsEdgeLengthMapping[existingConnectionEdgeInfo] = GraphOps.graphEdgeLength(graph = existingConnectionsGraph, edge = existingConnectionEdgeInfo)

        existingConnectionEdgeMinQuantile = numpy.quantile(tuple(existingConnectionsEdgeLengthMapping.values()), lonelyConnectionMinLengthQuantile)

        for (existingConnectionEdgeInfo, existingConnectionEdgeLength) in existingConnectionsEdgeLengthMapping.items():
            existingEdgeVertex0Id = existingConnectionEdgeInfo.vertex0Id
            existingEdgeVertex1Id = existingConnectionEdgeInfo.vertex1Id

            zeroNeighbors = tuple(existingConnectionsGraph.neighbors(existingEdgeVertex0Id))
            oneNeighbors = tuple(existingConnectionsGraph.neighbors(existingEdgeVertex1Id))

            if not len(zeroNeighbors) > 1 or not len(oneNeighbors) > 1:
                if existingConnectionEdgeLength < existingConnectionEdgeMinQuantile:
                    GraphOps.removeEdgeAndCleanUpNodes(edgeGraph = existingConnectionsGraph, existingEdge = existingConnectionEdgeInfo)

    # https://math.stackexchange.com/questions/134112/find-a-point-on-a-line-segment-located-at-a-distance-d-from-one-endpoint
    @staticmethod
    def _generateStopsOnConnection(connectionGraph: Graph, connection: EdgeVertexInfo, connectionLengthVertexPadding: float, connectionLengthNodeBuffer: float, minEdgeLengthForStopsGeneration: float) -> tuple[Point]:
        connectionLength = GraphOps.graphEdgeLength(graph = connectionGraph, edge = connection)
        
        # Only meaningfully generate stops if length is at least larger than ((connectionLengthVertexPadding + connectionLengthNodeBuffer) * 100) percent of edges.
        if (connectionLength >= minEdgeLengthForStopsGeneration):
            connectionVertexZeroId = connection.vertex0Id
            connectionVertexZero = connectionGraph.nodes[connectionVertexZeroId]

            connectionVertexOneId = connection.vertex1Id
            connectionVertexOne = connectionGraph.nodes[connectionVertexOneId]

            zeroIsFirstVertex = connectionVertexZero["x"] < connectionVertexOne["x"]

            firstVertex = connectionVertexZero if zeroIsFirstVertex else connectionVertexOne
            secondVertex = connectionVertexZero if not zeroIsFirstVertex else connectionVertexOne

            connectionXRangeMinInterval = (connectionLengthVertexPadding + connectionLengthNodeBuffer) * connectionLength
            connectionXRangeMaxInterval = connectionLength - connectionXRangeMinInterval

            connectionNodes = []

            while connectionXRangeMinInterval < connectionXRangeMaxInterval:
                d = random.uniform(connectionXRangeMinInterval, connectionXRangeMaxInterval)        
                td = d / connectionLength

                nextX = firstVertex["x"] + ((secondVertex["x"] - firstVertex["x"]) * td)
                nextY = firstVertex["y"] + ((secondVertex["y"] - firstVertex["y"]) * td)
                
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

        voronoiDiagramEdgesToProcess = []
        for voronoiDiagramEdge in voronoiDiagram.diagramEdges:
            edge = EdgeVertexInfo(vertex0Id = voronoiDiagramEdge.vertex0Id, vertex1Id = voronoiDiagramEdge.vertex1Id)
            if edge not in voronoiDiagramEdgesToProcess:
                voronoiDiagramEdgesToProcess.append(edge)

        voronoiDiagramEdges = tuple(voronoiDiagramEdgesToProcess)

        voronoiVertices = voronoiDiagram.vertices

        existingConnections = Graph(width = diagramWidth, height = diagramHeight)
        potentialConnections = Graph(width = diagramWidth, height = diagramHeight)

        for vertexId in voronoiVertices.keys():
            edgesWithVertex = set((maybeEdgeWithPoint for maybeEdgeWithPoint in voronoiDiagramEdges if maybeEdgeWithPoint.vertex0Id == vertexId or maybeEdgeWithPoint.vertex1Id == vertexId))
            
            for edgeWithPoint in edgesWithVertex:
                pointEdgeVertex0Id = edgeWithPoint.vertex0Id
                pointEdgeVertex0 = voronoiVertices[pointEdgeVertex0Id]

                pointEdgeVertex1Id = edgeWithPoint.vertex1Id
                pointEdgeVertex1 = voronoiVertices[pointEdgeVertex1Id]

                GraphOps.addVertexToGraph(graph = existingConnections, vertexId = pointEdgeVertex0Id, vertexX = pointEdgeVertex0.x, vertexY = pointEdgeVertex0.y)
                GraphOps.addVertexToGraph(graph = existingConnections, vertexId = pointEdgeVertex1Id, vertexX = pointEdgeVertex1.x, vertexY = pointEdgeVertex1.y)

                pointEdgeVertexDistance = GraphOps.scaledGraphPointDistance(graph = existingConnections, p1 = pointEdgeVertex0, p2 = pointEdgeVertex1)
                GraphOps.addConnectionToGraph(graph = existingConnections, connection = edgeWithPoint, vertexDistance = pointEdgeVertexDistance)

        existingEdgeVertexAngles = {}

        for (vertex0Id, vertex1Id) in existingConnections.edges:
            existingEdge = EdgeVertexInfo(vertex0Id = vertex0Id, vertex1Id = vertex1Id)

            vertex0EdgeAngles = TrackGenerator._calculateEdgeAnglesWithGraphVertexEdges(graph = existingConnections, edge = existingEdge, vertexId = vertex0Id)
            for vertex0AngleEdges in vertex0EdgeAngles:
                if vertex0AngleEdges not in existingEdgeVertexAngles:
                    existingEdgeVertexAngles[vertex0AngleEdges] = vertex0EdgeAngles[vertex0AngleEdges]

            vertex1EdgeAngles = TrackGenerator._calculateEdgeAnglesWithGraphVertexEdges(graph = existingConnections, edge = existingEdge, vertexId = vertex1Id)
            for vertex1AngleEdges in vertex1EdgeAngles:
                if vertex1AngleEdges not in existingEdgeVertexAngles:
                    existingEdgeVertexAngles[vertex1AngleEdges] = vertex1EdgeAngles[vertex1AngleEdges]

        initialDiagramEdgeAngles = tuple(existingEdgeVertexAngles.values())
        initialDiagramMinAcceptableAngle = numpy.quantile(a = initialDiagramEdgeAngles, q = newConnectionAngleMinQuantile)

        for vertexId in existingConnections.nodes:
            TrackGenerator._updateVertexPotentialConnections(existingConnectionsGraph = existingConnections, potentialConnectionsGraph = potentialConnections, vertexToUpdateId = vertexId)

        diagramEdgesSorted = sorted(voronoiDiagramEdges, key = lambda vde: GraphOps.graphEdgeLength(graph = existingConnections, edge = vde))
        numDiagramEdgesToProcess = floor(len(diagramEdgesSorted) * diagramEdgePercentageToProcess)

        diagramEdgesToProcess = diagramEdgesSorted[:numDiagramEdgesToProcess]

        intersectionEdges = []

        for diagramEdgeToProcess in diagramEdgesToProcess:
            edgeVertex0Id = diagramEdgeToProcess.vertex0Id
            edgeVertex1Id = diagramEdgeToProcess.vertex1Id

            if existingConnections.has_edge(u = edgeVertex0Id, v = edgeVertex1Id):
                edgeToProcessLength = GraphOps.graphEdgeLength(graph = existingConnections, edge = diagramEdgeToProcess)

                (disconnectedByDeletionVertex0, disconnectedByDeletionVertex1) = GraphOps.removeEdgeAndReturnDisconnected(edgeGraph = existingConnections, existingEdge = diagramEdgeToProcess)
                possibleConnections = []

                for vertexDisconnected0 in disconnectedByDeletionVertex0:
                    for vertexDisconnected1 in disconnectedByDeletionVertex1:
                        connectionAlreadyAdded = existingConnections.has_edge(u = vertexDisconnected0, v = vertexDisconnected1)

                        if not connectionAlreadyAdded:
                            potentialConnection = EdgeVertexInfo(vertex0Id = vertexDisconnected0, vertex1Id = vertexDisconnected1)
                            
                            if GraphOps.graphEdgeLength(graph = potentialConnections, edge = potentialConnection) > edgeToProcessLength:
                                possibleConnections.append(potentialConnection)

                maybeViableConnection = TrackGenerator._findPossibleViableConnection(graph = existingConnections, possibleConnections = tuple(possibleConnections), minAcceptableAngle = initialDiagramMinAcceptableAngle)

                if maybeViableConnection:
                    maybeViableVertexLength = GraphOps.graphEdgeLength(graph = potentialConnections, edge = maybeViableConnection) 

                    maybeViableVertex0Id = maybeViableConnection.vertex0Id
                    maybeViableVertex1Id = maybeViableConnection.vertex1Id

                    GraphOps.removeEdgeAndCleanUpNodes(edgeGraph = potentialConnections, existingEdge = maybeViableConnection)
                    
                    # calculate intersections
                    potentialIntersections = tuple((EdgeVertexInfo(vertex0Id = vertex0Id, vertex1Id = vertex1Id) for (vertex0Id, vertex1Id) in existingConnections.edges if (vertex0Id != maybeViableVertex0Id and vertex0Id != maybeViableVertex1Id) and (vertex1Id != maybeViableVertex0Id and vertex1Id != maybeViableVertex1Id)))
                    edgeIntersectionsInfo = EdgesOps.calculateEdgesIntersectionInfo(edgeGraph = existingConnections, intersectingEdge = maybeViableConnection, otherEdges = potentialIntersections)
                                            
                    if edgeIntersectionsInfo:
                        # We expect this to be sorted along the v0 -> v1 path.
                        precedingVertexId = maybeViableVertex0Id
                        precedingVertex = GraphOps.graphVertex(graph = existingConnections, vertexId = maybeViableVertex0Id)

                        for edgeIntersectionInfo in edgeIntersectionsInfo:
                            edgeIntersectionPoint = edgeIntersectionInfo.intersectionPoint
                            distanceToIntersection = GraphOps.scaledGraphPointDistance(graph = existingConnections, p1 = precedingVertex, p2 = edgeIntersectionPoint)

                            edgeIntersectionPointId = uuid4()

                            existingConnections.add_node(node_for_adding = edgeIntersectionPointId, x = edgeIntersectionPoint.x, y = edgeIntersectionPoint.y)

                            existingConnections.add_edge(u_of_edge = precedingVertexId, v_of_edge = edgeIntersectionPointId, edgeLength = distanceToIntersection)
                            intersectionEdges.append(EdgeVertexInfo(vertex0Id = precedingVertexId, vertex1Id = edgeIntersectionPointId))

                            edgeIntersected = edgeIntersectionInfo.intersectedEdge

                            edgeIntersectedVertex0Id = edgeIntersected.vertex0Id
                            edgeIntersectedVertex1Id = edgeIntersected.vertex1Id
                                                        
                            edgeIntersectedVertex0 = GraphOps.graphVertex(graph = existingConnections, vertexId = edgeIntersectedVertex0Id)
                            edgeIntersectedVertex1 = GraphOps.graphVertex(graph = existingConnections, vertexId = edgeIntersectedVertex1Id)

                            iv0DistanceToIntersection = GraphOps.scaledGraphPointDistance(graph = existingConnections, p1 = edgeIntersectedVertex0, p2 = edgeIntersectionPoint)
                            iv1DistanceToIntersection = GraphOps.scaledGraphPointDistance(graph = existingConnections, p1 = edgeIntersectionPoint, p2 = edgeIntersectedVertex1)

                            # No need to remove " orphaned " nodes here, since they'll subsequently be re-added with the handling of intersection edges.
                            existingConnections.remove_edge(u = edgeIntersectedVertex0Id, v = edgeIntersectedVertex1Id)
                            
                            if edgeIntersected in intersectionEdges:
                                intersectionEdges.remove(edgeIntersected)

                            existingConnections.add_edge(u_of_edge = edgeIntersectedVertex0Id, v_of_edge = edgeIntersectionPointId, edgeLength = iv0DistanceToIntersection)
                            intersectionEdges.append(EdgeVertexInfo(vertex0Id = edgeIntersectedVertex0Id, vertex1Id = edgeIntersectionPointId))

                            existingConnections.add_edge(u_of_edge = edgeIntersectionPointId, v_of_edge = edgeIntersectedVertex1Id, edgeLength = iv1DistanceToIntersection)
                            intersectionEdges.append(EdgeVertexInfo(vertex0Id = edgeIntersectionPointId, vertex1Id = edgeIntersectedVertex1Id))

                            TrackGenerator._updateVertexPotentialConnections(existingConnectionsGraph = existingConnections, potentialConnectionsGraph = potentialConnections, vertexToUpdateId = edgeIntersectionPointId)

                            precedingVertexId = edgeIntersectionPointId
                            precedingVertex = edgeIntersectionPoint

                        # finish up with precedingVertex (last intersection point) -> vertex1
                        maybeViableVertex1 = GraphOps.graphVertex(graph = existingConnections, vertexId = maybeViableVertex1Id)

                        distanceToEnd = GraphOps.scaledGraphPointDistance(graph = existingConnections, p1 = precedingVertex, p2 = maybeViableVertex1)

                        existingConnections.add_edge(u_of_edge = precedingVertexId, v_of_edge = maybeViableVertex1Id, edgeLength = distanceToEnd)
                        intersectionEdges.append(EdgeVertexInfo(vertex0Id = precedingVertexId, vertex1Id = maybeViableVertex1Id))
                    else:
                        existingConnections.add_edge(u_of_edge = maybeViableVertex0Id, v_of_edge = maybeViableVertex1Id, edgeLength = maybeViableVertexLength)

        TrackGenerator._pruneIntersectionEdgesFromExistingConnectionsGraph(existingConnectionsGraph = existingConnections, intersectionEdges = intersectionEdges)
        TrackGenerator._handleLonelyExistingConnections(existingConnectionsGraph = existingConnections, lonelyConnectionMinLengthQuantile = lonelyConnectionMinLengthQuantile)

        edges = {uuid4(): EdgeVertexInfo(vertex0Id = vertex0Id, vertex1Id = vertex1Id) for (vertex0Id, vertex1Id) in existingConnections.edges}

        edgeLengths = tuple((GraphOps.graphEdgeLength(graph = existingConnections, edge = edge) for edge in edges.values()))
        minEdgeStopLengthQuantile = connectionLengthVertexPadding + connectionLengthNodeBuffer

        minEdgeStopLength = numpy.quantile(a = edgeLengths, q = minEdgeStopLengthQuantile)
        
        stops = { edgeId: TrackGenerator._generateStopsOnConnection(connectionGraph = existingConnections, connection = edge, connectionLengthVertexPadding = connectionLengthVertexPadding, connectionLengthNodeBuffer = connectionLengthNodeBuffer, minEdgeLengthForStopsGeneration = minEdgeStopLength) for (edgeId, edge) in edges.items()}

        nodes = { nodeId: GraphOps.graphVertex(graph = existingConnections, vertexId = nodeId) for nodeId in existingConnections.nodes}

        return Track(nodes = nodes, stops = stops, edges = edges)
