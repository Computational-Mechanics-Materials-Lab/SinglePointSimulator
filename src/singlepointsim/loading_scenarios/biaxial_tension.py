#!/usr/bin/env python3

import numpy as np
from .base import Scenario

class BiaxialTensionScenario(Scenario):

    def get_loading_directions(self) -> tuple:
        print(f"Biaxial Tension in the {self.loading_direction_i}-{self.loading_directionj} plane.\n")
        return (
            self.loading_direction_i - 1,
            self.loading_direction_j - 1,
            5 - self.loading_direction_i - self.loading_direction_j,
        )

    def update_dfgrd(self, i: int, j: int, k: int) -> None:
        delta_Dkk = -self.stress[k] / self.ddsdde[k][k]
        self.dfgrd1[k][k] = self.dfgrd1[k][k] / (1.0 - delta_Dkk)

    def get_stress_tester(self, stress: np.ndarray, i: int, j: int, k: int) -> float:
        return 2.0 * abs(stress[k]) / (abs(stress[i]) + abs(stress[j]))

    def perform_loading(self, i: int, j: int, k: int) -> None:
        self.dfgrd1[i][i] = self.dfgrd0[i][i] + self.velocities[0] * self.dtime
        self.dfgrd1[j][j] = self.dfgrd0[j][j] + self.velocities[1] * self.dtime
        delta_Dii = (self.dfgrd1[i][i] - self.dfgrd0[i][i]) / self.dfgrd1[i][i]
        delta_Djj = (self.dfgrd1[j][j] - self.dfgrd0[j][j]) / self.dfgrd1[j][j]
        delta_Dkk = (
            self.stress[k]
            - self.ddsdde[i][i] * delta_Dii
            - self.ddsdde[j][j] * delta_Djj
        ) / self.ddsdde[k][k]

        self.dfgrd1[k][k] = self.dfgrd0[k][k] / (1 - delta_Dkk)
