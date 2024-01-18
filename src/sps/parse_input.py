#/usr/bin/env python3

import argparse
from dataclasses import dataclass

import json
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

@dataclass(slots=True)
class SPSInput:
    dtime: float
    time_max: float
    dtime_max: float
    loading_i: int
    loading_j: int
    nprops: int = -1
    props: list[int | float]
    nstatv: int
    v_1: int
    temp: float = 290.0
    dtemp: float = 0.0
    v_2: int | None = None

    def __post_init__(self):
        if self.v2 is None:
            self.v2 = self.v1
        self.nprops = len(self.props)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", nargs="*")
    args = parser.parse_args()
    if len(args) != 1:
        raise Exception("Expected one .toml or .json input file")
    input_file = pathlib.Path(args[0])

    try:
        with open(input_file, "rb") as fp:
            input_dict = tomllib.load(fp)
            return SPSInput(**input_dict)

    except tomllib.TOMLDecodeError:
        try:
            with open(file_path, "r") as fp:
                input_dict = json.load(fp)
                return SPSInput(**input_dict)

        except json.decoder.JSONDeocdeError:
            raise Exception(f"{input_file} is not in .toml or .json format")
