from __future__ import annotations
from dataclasses import dataclass
from typing import Self
from .sps_step import SPSStep


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
