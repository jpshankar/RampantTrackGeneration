from .data import EdgesMakingAngle, EdgeIntersectionInfo, EdgeVertexInfo

from ..graphs.Ops import Ops as GraphOps

from math import atan, degrees

from rustworkx import PyGraph as Graph

from voronout.Point import Point

class Ops:
    @staticmethod
    def _edgeSlope(edgePoint1: Point, edgePoint2: Point) -> float:
        if edgePoint2.x == edgePoint1.x:
            return 0.0
        else:
            return (edgePoint2.y - edgePoint1.y)/(edgePoint2.x - edgePoint1.x)
        
    # https://www.mathstopia.net/coordinate-geometry/angle-two-lines
    @staticmethod
    def _angleBetweenEdges(edgeGraph: Graph, edge0: EdgeVertexInfo, edge1: EdgeVertexInfo) -> float:
        edge0Slope = Ops._edgeSlope(edgePoint1 = GraphOps.graphVertex(graph = edgeGraph, vertexId = edge0.vertex0Id), edgePoint2 = GraphOps.graphVertex(graph = edgeGraph, vertexId = edge0.vertex1Id))
        edge1Slope = Ops._edgeSlope(edgePoint1 = GraphOps.graphVertex(graph = edgeGraph, vertexId = edge1.vertex0Id), edgePoint2 = GraphOps.graphVertex(graph = edgeGraph, vertexId = edge1.vertex1Id))

        slopeNumerator = edge0Slope - edge1Slope
        slopeDenominator = 1 + (edge0Slope * edge1Slope)

        return degrees(atan(abs(slopeNumerator / slopeDenominator)))
    
    @staticmethod
    def _calculateEdgeAngles(edgeGraph: Graph, angleEdges: EdgesMakingAngle) -> float:
        return Ops._angleBetweenEdges(edgeGraph = edgeGraph, edge0 = angleEdges.edge0, edge1 = angleEdges.edge1)

    def calculateEdgeAnglesWithExtantVertexEdges(edgeGraph: Graph, edgesMakingAngles: tuple[EdgesMakingAngle]) -> dict[EdgesMakingAngle, float]:
        return { edgesMakingAngle: Ops._calculateEdgeAngles(edgeGraph = edgeGraph, angleEdges = edgesMakingAngle) for edgesMakingAngle in edgesMakingAngles }
    
    # https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line
    @staticmethod
    def _calculateEdgeIntersectionInfo(edgeGraph: Graph, intersectingEdge: EdgeVertexInfo, intersectedEdge: EdgeVertexInfo) -> EdgeIntersectionInfo:
        intersectingVertex0 = GraphOps.graphVertex(graph = edgeGraph, vertexId = intersectingEdge.vertex0Id)
        intersectingVertex1 = GraphOps.graphVertex(graph = edgeGraph, vertexId = intersectingEdge.vertex1Id)

        intersectedVertex0 = GraphOps.graphVertex(graph = edgeGraph, vertexId = intersectedEdge.vertex0Id)
        intersectedVertex1 = GraphOps.graphVertex(graph = edgeGraph, vertexId = intersectedEdge.vertex1Id)

        x1dx3 = intersectingVertex0.x - intersectedVertex0.x
        y3dy4 = intersectedVertex0.y - intersectedVertex1.y

        y1dy3 = intersectingVertex0.y - intersectedVertex0.y
        x3dx4 = intersectedVertex0.x - intersectedVertex1.x

        x1dx2 = intersectingVertex0.x - intersectingVertex1.x
        y1dy2 = intersectingVertex0.y - intersectingVertex1.y

        denominator = (x1dx2 * y3dy4) - (y1dy2 * x3dx4)

        if denominator != 0.0:
            tNumerator = (x1dx3 * y3dy4) - (y1dy3 * x3dx4)
            t = tNumerator / denominator

            uNumerator = (x1dx2 * y1dy3) - (y1dy2 * x1dx3)
            u = -(uNumerator / denominator)

            if 0 <= t <= 1 and 0 <= u <= 1:
                intersectionX = intersectingVertex0.x + (t * -x1dx2)
                intersectionY = intersectingVertex0.y + (t * -y1dy2)

                intersectionPoint = Point(x = intersectionX, y = intersectionY)
                return EdgeIntersectionInfo(intersectionPoint = intersectionPoint, intersectedEdge = intersectedEdge)
            else:
                return None
        else:
            return None
    
    @staticmethod
    def calculateEdgesIntersectionInfo(edgeGraph: Graph, intersectingEdge: EdgeVertexInfo, otherEdges: tuple[EdgeVertexInfo]) -> tuple[EdgeIntersectionInfo]:
        edgeIntersectionData = []
        for otherEdge in otherEdges:
            edgeIntersectionInfo = Ops._calculateEdgeIntersectionInfo(edgeGraph = edgeGraph, intersectingEdge = intersectingEdge, intersectedEdge = otherEdge)
            if edgeIntersectionInfo:
                edgeIntersectionData.append(edgeIntersectionInfo)

        intersectingVertex0Id = intersectingEdge.vertex0Id
        intersectingVertex0 = GraphOps.graphVertex(graph = edgeGraph, vertexId = intersectingVertex0Id)

        edgeIntersectionData.sort(key = lambda eii: GraphOps.scaledGraphPointDistance(graph = edgeGraph, p1 = intersectingVertex0, p2 = eii.intersectionPoint))

        return edgeIntersectionData