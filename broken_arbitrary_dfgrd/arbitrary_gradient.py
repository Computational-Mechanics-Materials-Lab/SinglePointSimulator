#!/usr/bin/env python3

from .base import Scenario
import numpy as np
import numpy.typing as _npt
from scipy.linalg import fractional_matrix_power

class ArbitraryGradientScenario(Scenario):

    def get_loading_directions(self) -> tuple[int, int, int]:
        print("Arbitrary Loading along the diagonal\n")
        return 0, 1, 2

    def update_dfgrd(self, i: int, j: int, k: int) -> None:
        delta_Dii: _npt.NDArray = -(self.stress[i] + self.stress[j]) / (
            self.ddsdde[i][i]
            + self.ddsdde[i][j]
            + self.ddsdde[j][i]
            + self.ddsdde[j][j]
        )
        self.dfgrd1[i][i] /= 1 - delta_Dii
        self.dfgrd1[j][j] /= 1 - delta_Dii

    def get_stress_tester(self, stress: np.ndarray, i: int, j: int, k: int) -> float:
        return abs(stress[i] + stress[j]) / abs(stress[k])

    def perform_arbitrary_gradient_loading(input_dfgrd: np.ndarray) -> None:
        self.dfgrd1 = fractional_matrix_power(input_dfgrd, self.time[0] / self.time_max)

    def run_simulation(self):
        self.reset_variables()
        output_obj = OutputVar()
        ############################
        # Start new time step
        out_of_time = False
        #################################
        # Set up loading directions
        i, j, k = self.get_loading_directions()

        ##########################
        # Call umat
        sse = spd = scd = rpl = drpldt = cmname = 0
        pnewdt = celent = noel = npt = layer = kspt = kstep = kinc = 0

        self.stress, self.statev, self.ddsdde = self.umat(
            self.stress,
            self.statev,
            self.ddsdde,
            sse,
            spd,
            scd,
            rpl,
            self.ddsddt,
            self.drplde,
            drpldt,
            self.strain,
            self.dstrain,
            self.time,
            self.dtime,
            self.temp,
            self.dtemp,
            self.predef,
            self.dpred,
            cmname,
            self.NDI,
            self.NSHR,
            self.props,
            self.coords,
            self.drot,
            pnewdt,
            celent,
            self.dfgrd0,
            self.dfgrd1,
            noel,
            npt,
            layer,
            kspt,
            kstep,
            kinc,
            ntens=self.NTENS,
            nstatv=self.nstatv,
            nprops=self.nprops,
        )

        self.statev_ref[: self.nstatv] = self.statev[: self.nstatv]

        while not out_of_time:
            #######################
            # Increment one time step
            self.time += self.dtime
            # loading
            self.perform_arbitrary_gradient_loading(input_dfgrd)

            ########################
            # Start new convergence iteration.
            cont = False
            Kinc = 0
            test = np.inf
            while (test > self.TOLERANCE) and (Kinc <= self.MAX_ITR):  # Label 500
                Kinc += 1

                # Calculate F_dot * time
                # Note: this is seemingly faster in modern Python than relying on numpy for small dims - JK
                array1 = np.array(
                    [
                        [self.dfgrd1[m][n] - self.dfgrd0[m][n] for n in range(3)]
                        for m in range(3)
                    ]
                )
                array3 = np.array(
                    [
                        [
                            (self.dfgrd1[m][n] + self.dfgrd0[m][n]) / 2.0
                            for n in range(3)
                        ]
                        for m in range(3)
                    ]
                )

                # Multiply F_dot * F_inverse * dtime
                array2 = np.linalg.inv(np.array(array3))
                array3 = np.matmul(array1, array2)
                # look into putting these into a memoized function - I suspect repetitive input - JK
                D_dt = (array3 + np.transpose(array3)) / 2
                W_dt = (array3 + np.transpose(array3)) / 2

                # Store D_dt in dstrain
                self.dstrain[:3] = np.diag(D_dt)
                self.dstrain[3:5] = 2 * np.diag(D_dt, k=1)  # superdiagonal of D_dt
                self.dstrain[5] = 2 * D_dt[1][2]

                # Convert spin to drot[i,j] array for the UMAT
                self.drot = self.spin_to_matrix(W_dt)

                # Call umat to get stress
                self.stress, self.statev, self.ddsdde = self.umat(
                    self.stress,
                    self.statev,
                    self.ddsdde,
                    sse,
                    spd,
                    scd,
                    rpl,
                    self.ddsddt,
                    self.drplde,
                    drpldt,
                    self.strain,
                    self.dstrain,
                    self.time,
                    self.dtime,
                    self.temp,
                    self.dtemp,
                    self.predef,
                    self.dpred,
                    cmname,
                    self.NDI,
                    self.NSHR,
                    self.props,
                    self.coords,
                    self.drot,
                    pnewdt,
                    celent,
                    self.dfgrd0,
                    self.dfgrd1,
                    noel,
                    npt,
                    layer,
                    kspt,
                    kstep,
                    kinc,
                )

                #######
                # Check to see if another iteration is needed.
                ######
                test = self.get_stress_tester(self.stress, i, j, k)

                output_obj.vals.append(test)

                if test > self.TOLERANCE:
                    self.statev[:] = self.statev_ref
                    self.update_dfgrd(i, j, k)
                if Kinc > self.MAX_ITR:
                    self.time -= self.dtime
                    self.dtime /= 2.0
                    cont = True

            if cont:
                continue  # restarts outer while loop
            # Finished Increment!
            # Calc effective delta strain and add it to E_eff
            sm = np.sum(np.square(D_dt))
            dE_eff = np.sqrt(2.0 * sm / 3.0)
            self.E_eff += dE_eff

            # Update strain gradient
            self.strain += self.dstrain

            # Update deformation gradient
            self.dfgrd0[:] = self.dfgrd1

            # Update statev_ref
            self.statev_ref[: self.nstatv] = self.statev[: self.nstatv]

            # save outputs
            output_obj.stress.append(self.von_mises_stress())
            output_obj.strain.append(self.strain.copy())
            output_obj.time.append(self.time[1])
            output_obj.Eeff.append(self.E_eff)
            output_obj.all_dfgrd.append(self.dfgrd1.copy())

            ##########
            # Finish if out of time
            ##########
            if self.time[0] != self.dtime:
                self.dtime = self.dtime * 1.5
            if self.dtime > self.dtime_max:
                self.dtime = self.dtime_max
            if (self.time_max - self.time[0]) < self.dtime:
                self.dtime = self.time_max - self.time[0]
                out_of_time = True

        self.test_matrices(output_obj.stress, "Stress")
        self.test_matrices(output_obj.strain, "Strain")
        self.test_matrices(output_obj.all_dfgrd, "Deformation Gradient")
        return output_obj
