#!/usr/bin/env python3

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from math import sin, cos, sqrt

@dataclass(slots=True)
class SPSOutputs:
    stress: list[float] = field(default_factory=list)
    strain: list[float] = field(default_factory=list)
    time: list[float] = field(default_factory=list)
    Eeff: list[float] = field(default_factory=list)
    all_dfgrd: list[float] = field(default_factory=list)
    vals: list[float] = field(default_factory=list)

class Scenario:
    __slots__ = (
        "dtime",
        "time_max",
        "dtime_max",
        "props",
        "nstatv",
        "displacements",
        "loading_direction_i",
        "loading_direction_j",
        "temp",
        "dtemp",
        "velocities",
        "dfgrd",

        "umat",
        "props_ref",
        "dtime_ref",
        "statev",
        "statev_ref",
        "ddsdde",
        "ddsddt",
        "dpred",
        "drplde",
        "dstrain",
        "predef",
        "strain",
        "stress",
        "coords",

        "time",
        "E_eff",
        "dfgrd0",
        "dfgrd1",
        "drot",

        "NTENS",
        "MAX_ITR",
        "TOLERANCE",
        "NSHR",
        "NDI",
        )

    def __init__(self, step, umat):
        self.dtime = step.dtime
        self.time_max = step.time_max
        self.dtime_max = step.dtime_max
        self.props = step.props
        self.nprops = step.nprops
        self.nstatv = step.nstatv
        self.temp = step.temp
        self.dtemp = step.dtemp

        # Potentially None VVV
        self.displacements = step.displacements
        self.loading_direction_i = step.loading_direction_i
        self.loading_direction_j = step.loading_direction_j
        self.velocities = step.velocities
        self.dfgrd = step.dfgrd
        # Potentially None ^^^

        self.umat = umat

        # TODO Parameterize/Configure?
        self.NTENS = 6
        self.MAX_ITR = 4
        self.TOLERANCE = 0.001
        self.NSHR = 3
        self.NDI = 3

        self.props_ref = np.copy(self.props)
        self.dtime_ref = np.copy(self.dtime)

        self.statev = np.zeros(self.nstatv)
        self.statev_ref = np.zeros(self.nstatv)

        self.ddsdde = np.zeros((self.NTENS, self.NTENS))
        self.ddsddt = np.zeros(self.NTENS)
        self.dpred = np.zeros(1)
        self.drplde = np.zeros(self.NTENS)
        self.dstrain = np.zeros(self.NTENS)
        self.predef = np.zeros(1)
        self.strain = np.zeros(self.NTENS)
        self.stress = np.zeros(self.NTENS)
        self.coords = np.zeros(3)

        # Initialize time and deformation gradients
        self.time = np.zeros(2)
        self.E_eff = 0.0
        self.dfgrd0 = np.asfortranarray(np.identity(3))
        self.dfgrd1 = np.asfortranarray(np.identity(3))
        self.drot = np.asfortranarray(np.identity(3))


    def reset_variables(self) -> None:
        """resets all variables at the beginning of the simulation (in case of multiple function calls)
        without reassigning the variables in memory"""
        self.props[:] = self.props_ref
        self.ddsdde.fill(0)
        self.ddsddt.fill(0)
        self.dpred.fill(0)
        self.drplde.fill(0)
        self.dstrain.fill(0)
        self.predef.fill(0)
        self.strain.fill(0)
        self.stress.fill(0)
        self.coords.fill(0)
        self.time.fill(0)
        self.E_eff = 0.0
        self.dfgrd0 = np.asfortranarray(np.identity(3))
        self.dfgrd1 = np.asfortranarray(np.identity(3))
        self.drot = np.asfortranarray(np.identity(3))
        self.statev.fill(0)
        self.statev_ref.fill(0)
        self.dtime = self.dtime_ref

    def von_mises_stress(self) -> float:
        return np.sqrt(
            0.5
            * (
                np.power((self.stress[0] - self.stress[1]), 2)
                + np.power((self.stress[1] - self.stress[2]), 2)
                + np.power((self.stress[2] - self.stress[0]), 2)
            )
            + 3 * (np.sum(np.square(self.stress[3:6])))
        )

    def test_matrices(self, matrix: np.ndarray, var_name: str) -> None:
        # np.isreal will flag true any NaN or inf values, so additional conditions required
        if not np.all(np.isreal(matrix)) or (
            np.any(np.isnan(matrix)) or np.any(np.isinf(matrix))
        ):
            raise Exception(
                f"ERROR: {var_name} matrix contains NaN, infinity, or a non-real number"
            )

    def spin_to_matrix(self, a):
        """
        Converts spin tensor to a rotation matrix.
        """
        p1 = a[2][1]
        p2 = a[0][2]
        p3 = a[1][0]
        ang = sqrt(p1 * p1 + p2 * p2 + p3 * p3)

        s = sin(ang)
        c = cos(ang)

        # Normalize vector
        if ang < 1e-300:
            p1 = 0
            p2 = 0
            p3 = 1.0
        else:
            p1 = p1 / ang
            p2 = p2 / ang
            p3 = p3 / ang

        b = np.zeros((3, 3))
        b[0][0] = c + (1.0 - c) * p1**2
        b[0][1] = (1.0 - c) * p1 * p2 - s * p3
        b[0][2] = (1.0 - c) * p1 * p3 + s * p2
        b[1][0] = (1.0 - c) * p2 * p1 + s * p3
        b[1][1] = c + (1.0 - c) * p2**2
        b[1][2] = (1.0 - c) * p2 * p3 - s * p1
        b[2][0] = (1.0 - c) * p3 * p1 - s * p2
        b[2][1] = (1.0 - c) * p3 * p2 + s * p1
        b[2][2] = c + (1.0 - c) * p3 * p3

        return b

    def run_simulation(self):
        self.reset_variables()
        output_obj = SPSOutputs()
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
            self.perform_loading(i, j, k)

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

        names = list(output_obj.__slots__)
        names.remove("vals")

        final_result = pd.DataFrame({
            "time": output_obj.time,
            "stress": output_obj.stress,
            "strain": output_obj.strain,
            "Eeff": output_obj.Eeff,
            "all_dfgrd": output_obj.all_dfgrd,
            })
        return final_result
