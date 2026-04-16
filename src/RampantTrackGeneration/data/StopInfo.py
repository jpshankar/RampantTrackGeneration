from dataclasses import dataclass
from uuid import uuid4

from voronout.Point import Point

@dataclass(frozen=True)
class StopInfo:
    stopPoint: Point
    stopId: uuid4

    fuelAvailable: float

    def __repr__(self):
        return f'{{"stopPoint": {self.stopPoint}, "stopId": "{self.stopId}", "fuelAvailable": {self.fuelAvailable}}}'