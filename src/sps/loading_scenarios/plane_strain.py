#!/usr/bin/env python3

from base import Scenario

class PlaneStrainScenario(Scenario):

    def get_loading_directions(self) -> tuple:
        print(f"Plane-strain in the {self.loading_i}-{self.loading_j} plane\n")
        return (
            self.loading_i - 1,
            self.loading_j - 1,
            5 - self.loading_i - self.loading_j,
        )

    def update_dfgrd(self, i: int, j: int, k: int) -> None:
        delta_Dii = -self.stress[i] / self.ddsdde[i][i]
        self.dfgrd1[i][i] /= 1 - delta_Dii

    def get_stress_tester(self, stress: np.ndarray, i: int, j: int, k: int) -> float:
        return abs(stress[i] / stress[j])

    def perform_loading(self, i: int, j: int, k: int) -> None:
        self.dfgrd1[j][j] = self.dfgrd0[j][j] + self.v_1 * self.dtime
        delta_Djj = (self.dfgrd1[j][j] - self.dfgrd0[j][j]) / self.dfgrd1[j][j]
        delta_Dii = (
            -1 * (self.stress[i] + self.ddsdde[i][j] * delta_Djj) / self.ddsdde[i][i]
        )
        self.dfgrd1[i][i] = self.dfgrd0[i][i] / (1 - delta_Dii)
