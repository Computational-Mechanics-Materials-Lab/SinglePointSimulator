from .parse_sps_input import parse_sps_input
from .compile_umat import compile_umat

import pandas as pd
import pathlib

def run_sps():
    umat, results_dir, input_data = parse_sps_input()
    if not results_dir.exists():
        results_dir.mkdir(parents=True)
    compiled_umat = compile_umat(umat)
    umat_func = compiled_umat.umat

    for i, step in enumerate(input_data.steps, start=1):
        # TODO input_dfgrd
        sim = step.loading_scenario(step, umat_func)
        outputs = sim.run_simulation()
        outputs.to_csv(results_dir / f"{umat.stem}_{i}.csv")
