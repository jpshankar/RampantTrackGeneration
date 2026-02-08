from flask.json.provider import DefaultJSONProvider

from ... import Track

from uuid import UUID

from voronout.Point import Point

class FlaskRampantTrackGenerationJSONProvider(DefaultJSONProvider):
    def _handlePointDict(self, pointDict: dict[UUID, Point]) -> dict[str, dict[str, float]]:
        return {str(key): self.loads(repr(value)) for (key, value) in pointDict.items()}
        

    def dumps(self, obj, **kw):
        if isinstance(obj, Track):
            stopObj = {str(stopEdgeId): self.loads(f'[{",".join((repr(lineStop) for lineStop in edgeStops))}]') for (stopEdgeId, edgeStops) in obj.stops.items()}
            edgesObj = {str(edgeId): self.loads(repr(edge)) for (edgeId, edge) in obj.edges.items()}
            
            return {
                'nodes': self._handlePointDict(obj.nodes),
                'stops': stopObj,
                'edges': edgesObj
            }
        else:
            return super().dumps(obj, **kw)