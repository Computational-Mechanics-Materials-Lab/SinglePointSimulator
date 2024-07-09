from __future__ import annotations
from .parse_sps_input import parse_sps_input
from .loading_scenarios.base import Scenario
from .compile_umat import compile_umat
from .sps_dataclasses import SPSInput, SPSStep, SPSOutputs
from types import ModuleType
from typing import Callable, Type

import pandas as pd
import pathlib


def run_sps() -> None:
    umat: pathlib.Path
    results_dir: pathlib.Path
    input_data: SPSInput
    umat, results_dir, input_data = parse_sps_input()

    if not results_dir.exists():
        results_dir.mkdir(parents=True)
    compiled_umat: ModuleType = compile_umat(umat)
    umat_func: Callable = compiled_umat.umat

    i: int
    step: SPSStep
    num_steps: int = len(input_data.steps)
    for i, step in enumerate(input_data.steps, start=1):
        # TODO input_dfgrd
        sim: Type[Scenario] = step.loading_scenario(step, umat_func)
        outputs: SPSOutputs = sim.run_simulation()
        outputs_df: pd.DataFrame = outputs.to_df()

        outputs_df.to_csv(results_dir / f"{umat.stem}_step{i:0{num_steps}}.csv")
