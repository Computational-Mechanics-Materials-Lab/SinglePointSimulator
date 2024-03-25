#!/usr/bin/env python3

import numpy as np
from .base import Scenario

class TensionCompressionScenario(Scenario):
    def get_loading_directions(self) -> tuple:
        if self.loading_direction_i != self.loading_direction_j:
            # write exception handling later
            raise Exception(
                "WARNING! Loading directions not equivalent in input file.\n"
            )
        print(
            f"Simultaneous Tension and Compression in the {self.loading_direction_i}-{self.loading_direction_j} plane\n"
        )
        return self.loading_direction_j % 3, (self.loading_direction_j + 1) % 3, self.loading_direction_j - 1

    def update_dfgrd(self, i: int, j: int, k: int) -> None:
        delta_Dii = -(self.stress[i] + self.stress[j]) / (
            self.ddsdde[i][i]
            + self.ddsdde[i][j]
            + self.ddsdde[j][i]
            + self.ddsdde[j][j]
        )
        self.dfgrd1[i][i] /= 1 - delta_Dii
        self.dfgrd1[j][j] /= 1 - delta_Dii

    def get_stress_tester(self, stress: np.ndarray, i: int, j: int, k: int) -> float:
        return abs(stress[i] + stress[j]) / abs(stress[k])

    def perform_loading(self, i: int, j: int, k: int) -> None:
        self.dfgrd1[k][k] = self.dfgrd0[k][k] + self.velocities[0] * self.dtime
        delta_Dkk = (self.dfgrd1[k][k] - self.dfgrd0[k][k]) / self.dfgrd1[k][k]
        delta_Dii = (
            -1
            * (
                self.stress[i]
                + self.stress[j]
                + (self.ddsdde[i][k] + self.ddsdde[j][k]) * delta_Dkk
            )
            / (
                self.ddsdde[i][i]
                + self.ddsdde[i][j]
                + self.ddsdde[j][i]
                + self.ddsdde[j][j]
            )
        )

        self.dfgrd1[i][i] = self.dfgrd0[i][i] / (1.0 - delta_Dii)
        self.dfgrd1[j][j] = self.dfgrd0[j][j] / (1.0 - delta_Dii)


