#!/usr/bin/env python3

from base import Scenario

class TorsionScenario(Scenario):

    def get_loading_directions(self) -> tuple:
        print(f"Torsion in the {self.loading_i}-{self.loading_j} plane\n")
        return (
            self.loading_i - 1,
            self.loading_j - 1,
            5 - self.loading_i - self.loading_j,
        )

    def update_dfgrd(self, i: int, j: int, k: int) -> None:
        delta_Djj = -self.stress[j] / self.ddsdde[j][j]
        self.dfgrd1[j][j] /= 1 - delta_Djj

    def get_stress_tester(self, stress: np.ndarray, i: int, j: int, k: int) -> float:
        return abs(stress[j] / stress[i + j + 1])

    def perform_loading(self, i: int, j: int, k: int) -> None:
        m = i + j + 1
        self.dfgrd1[i][j] = self.dfgrd0[i][j] + self.v_1 * self.dtime
        delta_Dij = self.v_1 * self.dtime / 2.0 / self.dfgrd1[i][i]
        delta_Djj = (
            -1 * (self.stress[j] + self.ddsdde[j][m] * delta_Dij) / self.ddsdde[j][j]
        )

        self.dfgrd1[j][j] = self.dfgrd0[j][j] / (1.0 - delta_Djj)
