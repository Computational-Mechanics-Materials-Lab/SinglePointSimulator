from __future__ import annotations
from dataclasses import dataclass, field
from typing import Self, Final, Type
import pandas as pd
import numpy.typing as npt
from thick import Thick
from .loading_scenarios.base import Scenario
from .loading_scenarios import (
    TensionCompressionScenario,
    TorsionScenario,
    PlaneStrainScenario,
    CompressionTorsionScenario,
    BiaxialTensionScenario,
)


@dataclass(slots=True)
class SPSStep:
    loading_scenario: str | Type[Scenario]
    dtime: float
    time_max: float
    dtime_max: float
    props: list[int | float]
    nstatv: int
    displacements: list[float] | float
    loading_direction_i: int
    loading_direction_j: int
    nprops: int = 0
    temp: float = 290.0
    dtemp: float = 0.0
    velocities: list[float] = field(default_factory=list)

    def __post_init__(self: Self) -> None:
        loading_scenario_mapping: Final[dict[str, Type[Scenario]]] = Thick.make_thin(
            {
                (
                    "tension",
                    "compression",
                ): TensionCompressionScenario,
                ("torsion",): TorsionScenario,
                (
                    "plane strain",
                    "plane_strain",
                    "plane",
                    "strain",
                    "planestrain",
                    "plane-strain",
                ): PlaneStrainScenario,
                (
                    "torsion compression",
                    "compression torsion",
                    "torsion_compression",
                    "compression_torsion",
                    "torsioncompression",
                    "compressiontorsion",
                    "compression-torsion",
                    "torsion-compression",
                ): CompressionTorsionScenario,
                (
                    "biaxial tension",
                    "biaxialtension",
                    "biaxial_tension",
                    "biaxial-tension",
                ): BiaxialTensionScenario,
                # ("arbitrary", "abitrary gradient", "arbitrary-gradient", "arbitrary_gradient", "arbitrarygradient",): ArbitraryGradientScenario,
            }
        )

        self.nprops = len(self.props)

        if self.displacements is not None:
            if isinstance(self.displacements, float):
                self.displacements = [self.displacements]
            assert isinstance(self.displacements, list)
            self.velocities = [d / self.time_max for d in self.displacements]

        try:
            assert isinstance(self.loading_scenario, str)
            self.loading_scenario = loading_scenario_mapping[
                self.loading_scenario.lower().strip()
            ]
            assert isinstance(self.loading_scenario, Scenario)
        except KeyError as e:
            raise Exception(
                f"{self.loading_scenario} is not a valid loading scenario"
            ) from e


@dataclass(slots=True)
class SPSInput:
    props: list[int | float]
    nstatv: int
    steps: list[SPSStep] | list[dict]
    nprops: int = 0

    def __post_init__(self: Self) -> None:
        self.nprops = len(self.props)
        temp_steps: list[SPSStep] = []
        for s in self.steps:
            assert isinstance(s, dict)
            s["props"] = self.props
            s["nstatv"] = self.nstatv
            temp_steps.append(SPSStep(**s))

        self.steps = temp_steps
        assert isinstance(self.steps, list)


@dataclass(slots=True)
class SPSOutputs:
    stress: list[float] = field(default_factory=list)
    strain: list[npt.NDArray] = field(default_factory=list)
    time: list[float] = field(default_factory=list)
    Eeff: list[float] = field(default_factory=list)
    all_dfgrd: list[npt.NDArray] = field(default_factory=list)
    vals: list[float] = field(default_factory=list)

    def to_df(self: Self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "time": self.time,
                "stress": self.stress,
                "strain": self.strain,
                "Eeff": self.Eeff,
                "all_dfgrd": self.all_dfgrd,
            }
        )
