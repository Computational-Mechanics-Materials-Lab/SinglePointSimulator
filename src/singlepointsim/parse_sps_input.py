#/usr/bin/env python3

import argparse
import pathlib
from dataclasses import dataclass
from typing import Final

import numpy as np

import json
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from .loading_scenarios import *
from .loading_scenarios.base import Scenario

LOADING_SCENARIO_MAPPING: Final[dict[str, Scenario]] = {
        "tension":             TensionCompressionScenario,
        "compression":         TensionCompressionScenario,

        "torsion":             TorsionScenario,

        "plane strain":        PlaneStrainScenario,
        "plane_strain":        PlaneStrainScenario,
        "plane":               PlaneStrainScenario,
        "strain":              PlaneStrainScenario,
        "planestrain":         PlaneStrainScenario,
        "plane-strain":        PlaneStrainScenario,

        "torsion compression": CompressionTorsionScenario,
        "compression torsion": CompressionTorsionScenario,
        "torsion_compression": CompressionTorsionScenario,
        "compression_torsion": CompressionTorsionScenario,
        "torsioncompression":  CompressionTorsionScenario,
        "compressiontorsion":  CompressionTorsionScenario,
        "compression-torsion": CompressionTorsionScenario,
        "torsion-compression": CompressionTorsionScenario,

        "biaxial tension":     BiaxialTensionScenario,
        "biaxialtension":      BiaxialTensionScenario,
        "biaxial_tension":     BiaxialTensionScenario,
        "biaxial-tension":     BiaxialTensionScenario,

        "arbitrary":           ArbitraryGradientScenario,
        "abitrary gradient":   ArbitraryGradientScenario,
        "arbitrary-gradient":  ArbitraryGradientScenario,
        "arbitrary_gradient":  ArbitraryGradientScenario,
        "arbitrarygradient":   ArbitraryGradientScenario,
}


@dataclass(slots=True)
class SPSStep:
    loading_scenario: str | Scenario
    dtime: float
    time_max: float
    dtime_max: float
    props: list[int | float]
    nstatv: int
    nprops: int = 0
    displacements: list[float] | float | None = None
    loading_direction_i: int | None = None
    loading_direction_j: int | None = None
    temp: float = 290.0
    dtemp: float = 0.0
    velocities: list[float] | None = None
    dfgrd: list[float] | None = None
    randomize_dfgrd: bool = False

    def __post_init__(self):

        self.nprops = len(self.props)

        if self.displacements is not None:
            if not isinstance(self.displacements, list):
                self.displacements = [self.displacements]
            self.velocities = [d / self.time_max for d in self.displacements]

        try:
            self.loading_scenario = LOADING_SCENARIO_MAPPING[self.loading_scenario.lower().strip()]
        except KeyError as e:
            raise Exception(f"{self.loading_scenario} is not a valid loading scenario") from e

        if self.dfgrd is not None:
            if self.randomize_dfgrd:
                self.dfgrd = self.generate_random_gradient()
            else:
                self.dfgrd = self.generate_gradient()


    def generate_gradient(self):
        mask = np.asfortranarray(np.identity(3), dtype=int)
        mask[0][2] = 1
        final_dfgrd = np.asfortranarray(np.zeros((3, 3)))
        mask_indices = np.array(np.where(mask == 1))
        num_to_change = np.size(mask_indices[0])
        current_num = np.array(range(1, num_to_change + 1))
        values_len = np.size(self.dfgrd)
        values_indices = (1 // (values_len ** (num_to_change - current_num))) % values_len

        for i in range(num_to_change):
            final_dfgrd[tuple(mask_indices[:, i])] = self.dfgrd[values_indicies[i]]
        return final_dfgrd


    def generate_random_gradient(self):
        mask = np.asfortranarray(np.identity(3), dtype=int)
        mask[0][2] = 1
        final_dfgrd = np.asfortranarray(np.zeros((3, 3)))
        min_val = min(self.dfgrd)
        max_val = max(self.dfgrd)

        # TODO Parametrize?
        dfgrd_tol = 0.005
        upper_bound = 1 - dfgrd_tol
        lower_bound = 1 +dfgrd_tol

        mask_indices = np.where(mask != 0)

        num_non_zero = np.size(mask_indices[0])

        probs = np.array([upper_bound-min_val, max_val-lower_bound])
        probs /= probs.sum()

        rand_arr = np.array([np.random.choice([np.random.uniform(min_val, upper_bound), np.random.uniform(lower_bound, max_val)], p=probs) for _ in range(num_non_zero)])
        final_dfgrd[mask_indices] = rand_arr
        return final_dfgrd


@dataclass(slots=True)
class SPSInput:
    props: list[int | float]
    nstatv: int
    steps: list[SPSStep] | list[dict]
    nprops: int = 0

    def __post_init__(self):
        self.nprops = len(self.props)
        temp_steps = []
        for s in self.steps:
            s["props"] = self.props
            s["nstatv"] = self.nstatv
            temp_steps.append(SPSStep(**s))

        self.steps = temp_steps


def parse_sps_input():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", nargs="*")
    args = parser.parse_args()
    if len(vars(args)) != 1:
        raise Exception("Expected one .toml or .json input file")
    if not args.input_file:
        raise Exception("No .toml or .json input file was provided")
    input_file = pathlib.Path(args.input_file[0])

    try:
        with open(input_file, "rb") as fp:
            input_dict = tomllib.load(fp)
            umat = pathlib.Path(input_dict.pop("umat"))
            results_directory = pathlib.Path(input_dict.pop("results_directory", ".")
            return umat, results_directory, SPSInput(**input_dict)

    except tomllib.TOMLDecodeError:
        try:
            with open(file_path, "r") as fp:
                input_dict = json.load(fp)
                umat = pathlib.Path(input_dict.pop("umat"))
            results_directory = pathlib.Path(input_dict.pop("results_directory", ".")
                return umat, results_directory, SPSInput(**input_dict)

        except json.decoder.JSONDeocdeError:
            raise Exception(f"{input_file} is not in .toml or .json format")
