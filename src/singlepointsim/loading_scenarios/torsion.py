#!/usr/bin/env python3

import numpy.typing as npt
from .base import Scenario


class TorsionScenario(Scenario):
    def get_loading_directions(self) -> tuple[int, int, int]:
        print(
            f"Torsion in the {self.loading_direction_i}-{self.loading_direction_j} plane\n"
        )
        return (
            self.loading_direction_i - 1,
            self.loading_direction_j - 1,
            5 - self.loading_direction_i - self.loading_direction_j,
        )

    def update_dfgrd(self, i: int, j: int, k: int) -> None:
        delta_Djj: npt.NDArray = -self.stress[j] / self.ddsdde[j][j]
        self.dfgrd1[j][j] /= 1 - delta_Djj

    def get_stress_tester(self, stress: npt.NDArray, i: int, j: int, k: int) -> float:
        return abs(stress[j] / stress[i + j + 1])

    def perform_loading(self, i: int, j: int, k: int) -> None:
        m: int = i + j + 1
        self.dfgrd1[i][j] = self.dfgrd0[i][j] + self.velocities[0] * self.dtime
        delta_Dij: npt.NDArray = (
            self.velocities[0] * self.dtime / 2.0 / self.dfgrd1[i][i]
        )
        delta_Djj: npt.NDArray = (
            -1 * (self.stress[j] + self.ddsdde[j][m] * delta_Dij) / self.ddsdde[j][j]
        )

        self.dfgrd1[j][j] = self.dfgrd0[j][j] / (1.0 - delta_Djj)
