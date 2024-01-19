from .parse_sps_input import parse_sps_input
from .compile_umat import compile_umat

import pandas as pd
import pathlib

def main():
    umat, input_data = parse_sps_input()
    compiled_umat = compile_umat(umat)

    for i, step in enumerate(input_data.steps):
        # TODO input_dfgrd
        sim = step.loading_scenario(step, compiled_umat)
        outputs = sim.run_simulation()
        pd.DataFrame.from_dict({"Time": outputs.time, "Stress": outputs.stress, "Strain": [s[2] for s in outputs.strain]}).to_csv(pathlib.Path("results", f"sps_results_{i}.csv"))


if __name__ == "__main__":
    main()
