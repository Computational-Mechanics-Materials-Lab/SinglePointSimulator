#!/usr/bin/env python3
from __future__ import annotations
import numpy as np
import numpy.typing as _npt
from math import sin, cos, sqrt
from typing import Callable, Self, TYPE_CHECKING
from copy import copy
from abc import abstractmethod

from ..sps_outputs import SPSOutputs

if TYPE_CHECKING:
    from ..sps_step import SPSStep


class Scenario:
    __slots__ = (
        "dtime",
        "time_max",
        "dtime_max",
        "props",
        "nstatv",
        "loading_direction_i",
        "loading_direction_j",
        "temp",
        "dtemp",
        "velocities",
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

    def __init__(self: Self, step: SPSStep, umat: Callable) -> None:
        self.dtime: float = step.dtime
        self.time_max: float = step.time_max
        self.dtime_max: float = step.dtime_max
        self.props: _npt.NDArray = np.array(step.props)
        self.nprops: int = step.nprops
        self.nstatv: int = step.nstatv
        self.temp: float = step.temp
        self.dtemp: float = step.dtemp

        self.loading_direction_i: int = step.loading_direction_i
        self.loading_direction_j: int = step.loading_direction_j
        self.velocities: list[float] = step.velocities

        self.umat: Callable = umat

        # TODO Parameterize/Configure?
        self.NTENS: int = 6
        self.MAX_ITR: int = 4
        self.TOLERANCE: float = 0.001
        self.NSHR: int = 3
        self.NDI: int = 3

        self.props_ref: _npt.NDArray = np.copy(self.props)
        self.dtime_ref: float = copy(self.dtime)

        self.statev: _npt.NDArray = np.zeros(self.nstatv)
        self.statev_ref: _npt.NDArray = np.zeros(self.nstatv)

        self.ddsdde: _npt.NDArray = np.zeros((self.NTENS, self.NTENS))
        self.ddsddt: _npt.NDArray = np.zeros(self.NTENS)
        self.dpred: _npt.NDArray = np.zeros(1)
        self.drplde: _npt.NDArray = np.zeros(self.NTENS)
        self.dstrain: _npt.NDArray = np.zeros(self.NTENS)
        self.predef: _npt.NDArray = np.zeros(1)
        self.strain: _npt.NDArray = np.zeros(self.NTENS)
        self.stress: _npt.NDArray = np.zeros(self.NTENS)
        self.coords: _npt.NDArray = np.zeros(3)

        # Initialize time and deformation gradients
        self.time: _npt.NDArray = np.zeros(2)
        self.E_eff: float = 0.0
        self.dfgrd0: _npt.NDArray = np.asfortranarray(np.identity(3))
        self.dfgrd1: _npt.NDArray = np.asfortranarray(np.identity(3))
        self.drot: _npt.NDArray = np.asfortranarray(np.identity(3))

    def reset_variables(self: Self) -> None:
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

    def von_mises_stress(self: Self) -> float:
        return np.sqrt(
            0.5
            * (
                np.power((self.stress[0] - self.stress[1]), 2)
                + np.power((self.stress[1] - self.stress[2]), 2)
                + np.power((self.stress[2] - self.stress[0]), 2)
            )
            + 3 * (np.sum(np.square(self.stress[3:6])))
        )

    def test_matrices(self: Self, matrix: _npt.NDArray, var_name: str) -> None:
        # np.isreal will flag true any NaN or inf values, so additional conditions required
        if not np.all(np.isreal(matrix)) or (
            np.any(np.isnan(matrix)) or np.any(np.isinf(matrix))
        ):
            raise Exception(
                f"ERROR: {var_name} matrix contains NaN, infinity, or a non-real number"
            )

    @abstractmethod
    def get_loading_directions(self: Self) -> tuple[int, int, int]: ...

    @abstractmethod
    def get_stress_tester(
        self: Self, stress: _npt.NDArray, i: int, j: int, k: int
    ) -> float: ...

    @abstractmethod
    def perform_loading(self: Self, i: int, j: int, k: int) -> None: ...

    @abstractmethod
    def update_dfgrd(self: Self, i: int, j: int, k: int) -> None: ...

    def spin_to_matrix(self: Self, a: _npt.NDArray) -> _npt.NDArray:
        """
        Converts spin tensor to a rotation matrix.
        """
        p1: float = a[2][1]
        p2: float = a[0][2]
        p3: float = a[1][0]
        ang: float = sqrt(p1 * p1 + p2 * p2 + p3 * p3)

        s: float = sin(ang)
        c: float = cos(ang)

        # Normalize vector
        if ang < 1e-300:
            p1 = 0
            p2 = 0
            p3 = 1.0
        else:
            p1 = p1 / ang
            p2 = p2 / ang
            p3 = p3 / ang

        b: _npt.NDArray = np.zeros((3, 3))
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

    def run_simulation(self: Self) -> SPSOutputs:
        self.reset_variables()
        output_obj: SPSOutputs = SPSOutputs()
        ############################
        # Start new time step
        out_of_time: bool = False
        #################################
        # Set up loading directions
        i: int
        j: int
        k: int
        i, j, k = self.get_loading_directions()

        ##########################
        # Call umat
        sse: float = 0.0
        spd: float = 0.0
        scd: float = 0.0
        rpl: float = 0.0
        drpldt: float = 0.0
        cmname: float = 0.0
        pnewdt: float = 0.0
        celent: float = 0.0

        noel: int = 0
        npt: int = 0
        layer: int = 0
        kspt: int = 0
        kstep: int = 0
        kinc: int = 0

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
            cont: bool = False
            Kinc: int = 0
            test: float = np.inf
            while (test > self.TOLERANCE) and (Kinc <= self.MAX_ITR):  # Label 500
                Kinc += 1

                # Calculate F_dot * time
                dfgrd_diff: _npt.NDArray = self.dfgrd1 - self.dfgrd0
                dfgrd_mean: _npt.NDArray = np.mean((self.dfgrd1, self.dfgrd0), axis=0)

                # Multiply F_dot * F_inverse * dtime
                dfgrd_mean_inv: _npt.NDArray = np.linalg.inv(dfgrd_mean)
                dfgrd_diff_times_mean_inv: _npt.NDArray = np.matmul(
                    dfgrd_diff, dfgrd_mean_inv
                )

                D_dt = (dfgrd_diff_times_mean_inv + dfgrd_diff_times_mean_inv.T) / 2
                W_dt = (dfgrd_diff_times_mean_inv + dfgrd_diff_times_mean_inv.T) / 2

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
            sm: _npt.NDArray = np.sum(np.square(D_dt))
            dE_eff: float = np.sqrt(2.0 * sm / 3.0)
            self.E_eff += dE_eff

            # Update strain gradient
            self.strain += self.dstrain

            # Update deformation gradient
            self.dfgrd0[:] = self.dfgrd1

            # Update statev_ref
            self.statev_ref[: self.nstatv] = self.statev[: self.nstatv]

            # save outputs
            foo: float = self.von_mises_stress()
            print("Von Mises Stress Type:", type(foo))
            output_obj.stress.append(foo)
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
