#!/usr/bin/env python3

from base import Scenario

class CompressionTorsionScenario(Scenario):

    def get_loading_directions(self) -> tuple:
        print(
            f"Simultaneous Compression and Torsion in the {self.loading_i}-{self.loading_j} plane\n"
        )
        return (
            self.loading_i - 1,
            5 - self.loading_i - self.loading_j,
            self.loading_j - 1,
        )

    def update_dfgrd(self, i: int, j: int, k: int) -> None:
        det = (
            self.ddsdde[i][i] * self.ddsdde[j][j]
            - self.ddsdde[i][j] * self.ddsdde[j][i]
        )
        b1 = -self.stress[i]
        b2 = -self.stress[j]
        delta_Dii = (b1 * self.ddsdde[j][j] - b2 * self.ddsdde[i][j]) / det
        delta_Djj = (self.ddsdde[i][i] * b2 - self.ddsdde[j][i] * b1) / det
        self.dfgrd1[i][i] /= 1 - delta_Dii
        self.dfgrd1[j][j] /= 1 - delta_Djj

    def get_stress_tester(self, stress: np.ndarray, i: int, j: int, k: int) -> float:
        # if (abs(stress[k]) + abs(stress[i + j + 1])):
        #    raise Exception
        return (abs(stress[i]) + abs(stress[j])) / (
            abs(stress[k]) + abs(stress[i + j + 1])
        )

    def perform_loading(self, i: int, j: int, k: int) -> None:
        m = i + k + 2
        self.dfgrd1[k][k] = self.dfgrd0[k][k] + self.v_1 * self.dtime
        self.dfgrd1[i][k] = self.dfgrd0[i][k] + self.v_2 * self.dtime

        delta_Dkk = (self.dfgrd1[k][k] - self.dfgrd0[k][k]) / self.dfgrd1[k][k]
        delta_Dik = self.v_2 * self.dtime / 2.0 / self.dfgrd1[i][i]

        det = (
            self.ddsdde[i][i] * self.ddsdde[j][j]
            - self.ddsdde[i][j] * self.ddsdde[j][i]
        )
        b1 = (
            -1 * self.stress[i]
            - self.ddsdde[i][k] * delta_Dkk
            - self.ddsdde[i][m] * delta_Dik
        )
        b2 = (
            -1 * self.stress[j]
            - self.ddsdde[j][k] * delta_Dkk
            - self.ddsdde[j][m] * delta_Dik
        )
        delta_Dii = (b1 * self.ddsdde[j][j] - b2 * self.ddsdde[i][j]) / det
        delta_Djj = (self.ddsdde[i][i] * b2 - self.ddsdde[j][i] * b1) / det
        self.dfgrd1[i][i] = self.dfgrd0[i][i] / (1 - delta_Dii)
        self.dfgrd1[j][j] = self.dfgrd0[j][j] / (1 - delta_Djj)
