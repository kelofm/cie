# --- STL Imports ---
import math


offset: tuple[float,float] = (1.0, 2.0)
domainRadius: float = 2.0
dirichletRadius: float = 1.0
resolution: int = 50
angleIncrement: float = 2 * math.pi / resolution


with open("domain_hull.csv", "w") as domainHullFile:
    with open("dirichlet_constraints.csv", "w") as dirichletFile:
        for iSegment in range(resolution):
            angleBegin: float   = iSegment * angleIncrement
            angleEnd: float     = angleBegin + angleIncrement

            # Domain's boundary polygon that later has to be meshed.
            begin: tuple[float,float] = (
                offset[0] + domainRadius * math.cos(angleBegin),
                offset[1] + domainRadius * math.sin(angleBegin))
            end: tuple[float,float] = (
                offset[0] + domainRadius * math.cos(angleEnd),
                offset[1] + domainRadius * math.sin(angleEnd))
            domainHullFile.write(f"{begin[0]:.4E},{begin[1]:.4E},{end[0]:.4E},{end[1]:.4E}\n")

            # Dirichlet conditions on a curve.
            valueBegin: float   = math.sin(angleBegin)
            valueEnd: float     = math.sin(angleEnd)
            begin: tuple[float,float] = (
                offset[0] + dirichletRadius * math.cos(angleBegin),
                offset[1] + dirichletRadius * math.sin(angleBegin))
            end: tuple[float,float] = (
                offset[0] + dirichletRadius * math.cos(angleEnd),
                offset[1] + dirichletRadius * math.sin(angleEnd))
            dirichletFile.write(f"{begin[0]:.4E},{begin[1]:.4E},{end[0]:.4E},{end[1]:.4E},{valueBegin:.4E},{valueEnd:.4E}\n")
