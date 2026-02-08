from dataclasses import dataclass
from uuid import uuid4

@dataclass(frozen=True)
class EdgeVertexInfo:
    vertex0Id: uuid4
    vertex1Id: uuid4

    def __eq__(self, other) -> bool:
        return (self.vertex0Id == other.vertex0Id and self.vertex1Id == other.vertex1Id) or (self.vertex0Id == other.vertex1Id and self.vertex1Id == other.vertex0Id)

    def __repr__(self):
        return f'{{"vertex0Id": "{str(self.vertex0Id)}", "vertex1Id": "{str(self.vertex1Id)}"}}'