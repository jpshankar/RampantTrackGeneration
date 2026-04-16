from . import StopInfo

from dataclasses import dataclass
from uuid import uuid4

@dataclass(frozen=True)
class EdgeInfo:
    fuelCost: float
    edgeLengthProportion: float

    edgeStopInfo: tuple[StopInfo]

    edgeStopsFromVertex0: tuple[uuid4]
    edgeStopsFromVertex1: tuple[uuid4]

    def __repr__(self):
        edgeStopsFromVertex0Repr = [str(edgeStopFromVertex0) for edgeStopFromVertex0 in self.edgeStopsFromVertex0]
        edgeStopsFromVertex1Repr = [str(edgeStopFromVertex1) for edgeStopFromVertex1 in self.edgeStopsFromVertex1]

        return f'{{"fuelCost": {self.fuelCost}, "edgeLengthProportion": {self.edgeLengthProportion}, "stopInfo": {list(self.edgeStopInfo)}, "edgeStopsFromVertex0": {edgeStopsFromVertex0Repr}, "edgeStopsFromVertex1": {edgeStopsFromVertex1Repr}}}'