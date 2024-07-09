#!/usr/bin/env python3
import pathlib
import importlib
import subprocess
import os
from typing import TextIO
from types import ModuleType
import numpy.f2py
from numpy.f2py.crackfortran import COMMON_FREE_EXTENSIONS, COMMON_FIXED_EXTENSIONS
from .patch_meson import patch_meson


def compile_umat(
    fortran_path: pathlib.Path, regenerate_sig: bool = False
) -> ModuleType:
    module_name: str
    signature_path: pathlib.Path
    module_name, signature_path = get_fortran_metadata(fortran_path)

    try:
        return importlib.import_module(module_name)

    except ModuleNotFoundError:
        pass

    patch_meson()

    if (not signature_path.exists()) or regenerate_sig:
        create_pyf(fortran_path, signature_path, module_name)

    compile_to_so(fortran_path, signature_path, module_name)
    return importlib.import_module(module_name)


def get_fortran_metadata(fortran_path: pathlib.Path) -> tuple[str, pathlib.Path]:
    module_name: str = fortran_path.stem
    signature_path: pathlib.Path = pathlib.Path(f"{module_name}.pyf")
    return module_name, signature_path


def create_pyf(
    fortran_path: pathlib.Path, signature_path: pathlib.Path, module_name: str
) -> None:
    numpy.f2py.run_main(
        [
            "-h",
            str(signature_path),
            str(fortran_path),
            "-m",
            module_name,
            "--overwrite-signature",
            "--verbose",
        ]
    )

    sig_file: TextIO
    with open(signature_path, "r", encoding="utf-8") as sig_file:
        sig_contents: str = sig_file.read()

    sig_contents = sig_contents.replace(" :: stress\n", ", intent(in,out) :: stress\n")
    sig_contents = sig_contents.replace(" :: statev\n", ", intent(in,out) :: statev\n")
    sig_contents = sig_contents.replace(" :: ddsdde\n", ", intent(in,out) :: ddsdde\n")
    sig_contents = sig_contents.replace(" :: time\n", ", intent(in) :: time\n")
    sig_contents = sig_contents.replace(" :: dtime\n", ", intent(in) :: dtime\n")
    sig_contents = sig_contents.replace(" :: dfgrd0\n", ", intent(in) :: dfgrd0\n")
    sig_contents = sig_contents.replace(" :: dfgrd1\n", ", intent(in) :: dfgrd1\n")

    with open(signature_path, "w") as sig_file:
        sig_file.write(sig_contents)


def compile_to_so(
    fortran_path: pathlib.Path, signature_path: pathlib.Path, module_name: str
) -> None:
    #

    old_FC_env = os.environ.get("FC")
    old_FFLAGS_env = os.environ.get("FFLAGS")

    # TODO let the user add/remove these
    new_fflags_list = ["-fcray-pointer", "-cpp", "-fPIC", "-O3", "-funroll-loops"]

    if fortran_path.suffix in COMMON_FIXED_EXTENSIONS:
        new_fflags_list.append("-ffixed-form")

    elif fortran_path.suffix in COMMON_FREE_EXTENSIONS:
        new_fflags_list.append("-ffree-form")

    # TODO let the user select this
    os.environ["FC"] = "gfortran"
    # os.environ["FC"] = "mpifort"
    # os.environ["FC"] = "mpif77"
    os.environ["FFLAGS"] = " ".join(new_fflags_list)

    try:
        ret: subprocess.CompletedProcess = subprocess.run(
            [
                "python",
                "-m",
                "numpy.f2py",
                "-c",
                signature_path,
                fortran_path,
                "--backend",
                "meson",
            ],
            capture_output=True,
        )

        if ret.returncode != 0:
            raise RuntimeError(ret.stderr)

    finally:
        if old_FC_env is None:
            del os.environ["FC"]
        else:
            os.environ["FC"] = old_FC_env

        if old_FFLAGS_env is None:
            del os.environ["FFLAGS"]
        else:
            os.environ["FFLAGS"] = old_FFLAGS_env
