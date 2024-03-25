#!/usr/bin/env python3
import pathlib
import importlib
import subprocess
import shutil


def compile_umat(fortran_path: pathlib.Path) -> None:
    module_name = fortran_path.stem
    signature_path = pathlib.Path(f"{module_name}.pyf")

    try:
        return importlib.import_module(module_name)

    except ModuleNotFoundError:

        # Handle Signature File
        if not signature_path.exists():
            result = subprocess.run([
                "python",
                "-m",
                "numpy.f2py",
                "-h",
                signature_path,
                fortran_path,
                "-m",
                module_name,
                "--overwrite-signature",
                ])

            if result.returncode != 0:
                raise Exception(f"Error generating signature file for {fortran_path} with f2py")

            else:
                with open(signature_path, "r") as sig_file:
                    sig_contents = sig_file.read()

                sig_contents = sig_contents.replace(" :: stress\n", ", intent(in,out) :: stress\n")
                sig_contents = sig_contents.replace(" :: statev\n", ", intent(in,out) :: statev\n")
                sig_contents = sig_contents.replace(" :: ddsdde\n", ", intent(in,out) :: ddsdde\n")
                sig_contents = sig_contents.replace(" :: time\n", ", intent(in) :: time\n")
                sig_contents = sig_contents.replace(" :: dtime\n", ", intent(in) :: dtime\n")
                sig_contents = sig_contents.replace(" :: dfgrd0\n", ", intent(in) :: dfgrd0\n")
                sig_contents = sig_contents.replace(" :: dfgrd1\n", ", intent(in) :: dfgrd1\n")

                with open(signature_path, "w") as sig_file:
                    sig_file.write(sig_contents)

        mpif90 = shutil.which("mpif90")
        mpif77 = shutil.which("mpif77")

        if mpif90 is None:
            raise Exception("Executable for mpif90 could not be found")
        if mpif77 is None:
            raise Exception("Executable for mpif77 could not be found")

        # TODO intel vs gnu vs ...
        fortran_compiler_args = "-ffree-form -ffree-line-length-512 -fcray-pointer -cpp"

        subprocess.run([
            "python",
            "-m",
            "numpy.f2py",
            "-c",
            signature_path,
            fortran_path,
            "-m",
            module_name,
            f"--opt={fortran_compiler_args}",
            f"--f90exec={mpif90}",
            f"--f77exec={mpif77}",
            ])

        return importlib.import_module(module_name)
