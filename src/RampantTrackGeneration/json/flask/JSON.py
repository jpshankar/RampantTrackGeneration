from flask.json.provider import DefaultJSONProvider

from ... import Track

from uuid import UUID

from voronout.Point import Point

class FlaskRampantTrackGenerationJSONProvider(DefaultJSONProvider):
    def _handlePointDict(self, pointDict: dict[UUID, Point]) -> dict[str, dict[str, float]]:
        return {str(key): self.loads(repr(value)) for (key, value) in pointDict.items()}
        

    def dumps(self, obj, **kw):
        if isinstance(obj, Track):
            stopObj = {self.loads(repr(stopLine)): list((self.loads(repr(lineStop)) for lineStop in lineStops)) for (stopLine, lineStops) in obj.stops.items()}
            
            return {
                'nodes': self._handlePointDict(obj.nodes),
                'stops': stopObj,
                'edges': list(stopObj.values())
            }
        else:
            return super().dumps(obj, **kw)