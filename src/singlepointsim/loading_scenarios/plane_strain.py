#!/usr/bin/env python3

import numpy.typing as npt
from .base import Scenario


class PlaneStrainScenario(Scenario):
    def get_loading_directions(self) -> tuple[int, int, int]:
        print(
            f"Plane-strain in the {self.loading_direction_i}-{self.loading_direction_j} plane\n"
        )
        return (
            self.loading_direction_i - 1,
            self.loading_direction_j - 1,
            5 - self.loading_direction_i - self.loading_direction_j,
        )

    def update_dfgrd(self, i: int, j: int, k: int) -> None:
        delta_Dii: npt.NDArray = -self.stress[i] / self.ddsdde[i][i]
        self.dfgrd1[i][i] /= 1 - delta_Dii

    def get_stress_tester(self, stress: npt.NDArray, i: int, j: int, k: int) -> float:
        return abs(stress[i] / stress[j])

    def perform_loading(self, i: int, j: int, k: int) -> None:
        self.dfgrd1[j][j] = self.dfgrd0[j][j] + self.velocities[0] * self.dtime
        delta_Djj: npt.NDArray = (self.dfgrd1[j][j] - self.dfgrd0[j][j]) / self.dfgrd1[
            j
        ][j]
        delta_Dii: npt.NDArray = (
            -1 * (self.stress[i] + self.ddsdde[i][j] * delta_Djj) / self.ddsdde[i][i]
        )
        self.dfgrd1[i][i] = self.dfgrd0[i][i] / (1 - delta_Dii)
