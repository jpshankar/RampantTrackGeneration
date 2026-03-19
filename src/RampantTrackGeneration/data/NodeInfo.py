from dataclasses import dataclass

@dataclass(frozen=True)
class NodeInfo:
    distanceToDestination: float
    numNeighbors: int

    def __repr__(self):
        return f'{{"distanceToDestination": {self.distanceToDestination}, "numNeighbors": {self.numNeighbors}}}'