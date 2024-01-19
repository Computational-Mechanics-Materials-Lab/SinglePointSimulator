#/usr/bin/env python3

import argparse
import pathlib
from dataclasses import dataclass
from typing import Final

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
        "compression_torsion": CompressionTorsoinScenario,
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
    loading_direction_i: int
    loading_direction_j: int
    velocity_1: float
    props: list[int | float]
    nstatv: int
    temp: float = 290.0
    dtemp: float = 0.0
    velocity_2: float | None = None
    nprops: int = -1

    def __post_init__(self):
        self.nprops = len(self.props)
        if self.velocity_2 is None:
            self.velocity_2 = self.velocity_1

        try:
            self.loading_scenario = LOADING_SCENARIO_MAPPING[self.loading_scenario.lower()]
        except KeyError as e:
            raise Exception(f"{self.loading_scenario} is not a valid loading scenario") from e


@dataclass(slots=True)
class SPSInput:
    props: list[int | float]
    nstatv: int
    steps: list[SPSStep] | list[dict]

    def __post_init__(self):
        self.nprops = len(self.props)
        temp_steps = []
        for s in self.steps:
            s["props"] = self.props
            s["nstatv"] = self.nstatv
            temp_steps.append(SPSStep(**s)]

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
            return umat, SPSInput(**input_dict)

    except tomllib.TOMLDecodeError:
        try:
            with open(file_path, "r") as fp:
                input_dict = json.load(fp)
                umat = pathlib.Path(input_dict.pop("umat"))
                return umat, SPSInput(**input_dict)

        except json.decoder.JSONDeocdeError:
            raise Exception(f"{input_file} is not in .toml or .json format")
