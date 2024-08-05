# /usr/bin/env python3

import argparse
import pathlib
from typing import BinaryIO, TextIO, Any


import json

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from .sps_input import SPSInput


class MissingUMATException(Exception):
    pass


def parse_sps_input() -> tuple[pathlib.Path, pathlib.Path, SPSInput]:
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("umat_file", nargs="?")
    args: argparse.Namespace = parser.parse_args()
    if len(vars(args)) > 2:
        raise Exception(
            "Expected one .toml or .json input file and at most one umat file!"
        )
    if not args.input_file:
        raise Exception("No .toml or .json input file was provided")
    input_file: pathlib.Path = pathlib.Path(args.input_file)

    umat_file: pathlib.Path | None = None
    if args.umat_file is not None:
        umat_file = pathlib.Path(args.umat_file)

    input_dict: dict[str, Any]
    try:
        bfp: BinaryIO
        with open(input_file, "rb") as bfp:
            input_dict = tomllib.load(bfp)

    except tomllib.TOMLDecodeError:
        try:
            tfp: TextIO
            with open(input_file, "r") as tfp:
                input_dict = json.load(tfp)
        except json.decoder.JSONDecodeError:
            raise Exception(f"{input_file} is not in .toml or .json format")

    try:
        return dict_to_ret(input_dict, umat_file)
    except MissingUMATException:
        raise MissingUMATException(
            f"UMAT not supplied as either a command line argument or a key-value pair in {input_file}!"
        )


def dict_to_ret(
    input_dict: dict[str, Any], umat_file: pathlib.Path | None
) -> tuple[pathlib.Path, pathlib.Path, SPSInput]:
    umat_path: str | None = input_dict.pop("umat", None)
    if umat_file is None:
        if umat_path is None:
            raise MissingUMATException
        umat_file = pathlib.Path(umat_path)
    results_directory: pathlib.Path = pathlib.Path(
        input_dict.pop("results_directory", ".")
    )
    return umat_file, results_directory, SPSInput(**input_dict)
