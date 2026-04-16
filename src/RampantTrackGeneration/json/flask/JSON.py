from flask.json.provider import DefaultJSONProvider

from ... import Track

from dataclasses import asdict

class FlaskRampantTrackGenerationJSONProvider(DefaultJSONProvider):
    def dumps(self, obj, **kw):
        if isinstance(obj, Track):
            nodesObj = {}
            for (nodeId, node) in obj.nodes.items():
                nodeObj = asdict(node)
                nodeObj["nodeInfo"] = obj.nodeInfo[nodeId]

                nodesObj[str(nodeId)] = nodeObj

            edgesObj = {}
            for (edgeId, edge) in obj.edges.items():
                edgeObj = asdict(edge)

                edgeObj["vertex0Id"] = str(edgeObj["vertex0Id"])
                edgeObj["vertex1Id"] = str(edgeObj["vertex1Id"])

                edgeObj["edgeInfo"] = obj.edgeInfo[edgeId]

                edgesObj[str(edgeId)] = edgeObj

            return {
                'nodes': nodesObj,
                'edges': edgesObj,
                'startNodeId': str(obj.startNodeId),
                'destinationNodeId': str(obj.destinationNodeId)
            }
        else:
            return super().dumps(obj, **kw)