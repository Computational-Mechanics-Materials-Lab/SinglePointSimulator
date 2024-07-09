import pathlib
import mesonbuild
import sys
from typing import TextIO


def patch_meson() -> None:
    meson_root_path: pathlib.Path = pathlib.Path(mesonbuild.__file__).parent
    meson_compilers_dir_path: pathlib.Path = meson_root_path / "compilers"
    meson_fortran_path: pathlib.Path = meson_compilers_dir_path / "fortran.py"
    meson_compilers_path: pathlib.Path = meson_compilers_dir_path / "compilers.py"

    meson_patched_f77: str = "#################\n# PATCHED FOR F77\n#################\n"

    # Open the fortran.py file, check for a flag that this program has already modified it. We'll just exit if it has, if not, we'll edit this copy and save it back
    fortranpy_fh: TextIO
    with open(meson_fortran_path, "r") as fortranpy_fh:
        fortranpy_content: str = fortranpy_fh.read()

    if fortranpy_content.startswith(meson_patched_f77):
        return

    user_understands: str = input(
        'WARNING: This program requires changes to the Meson build system which are not yet implemented in the current version of Meson! Meson does not correctly handle compilation of F77/Fixed Form Fortran.\nThis function will patch the current version of Meson in your Python Environment to address this.\nIf you want undo this, you will have to uninstall and then reinstall meson with pip.\nThis may not work if meson is not installed in a virtual environment!\n\nPlease type "I UNDERSTAND" in all CAPS to confirm this change: '
    )
    if user_understands.strip() != "I UNDERSTAND":
        sys.exit(1)

    # Fix the issues in fortran.py

    # Go ahead and re-added that patch note
    fortranpy_content = meson_patched_f77 + fortranpy_content

    # Fix the Sanity Check
    old_sanity_check: str = """        source_name = 'sanitycheckf.f90'\n        code = 'program main; print *, "Fortran compilation is working."; end program\\n'"""

    new_sanity_check: str = """        source_name = 'sanitycheckf.f'\n        code = '       PROGRAM MAIN\\n         PRINT *, "Fortran compilation is working."\\n       END PROGRAM\\n'"""

    fortranpy_content = fortranpy_content.replace(old_sanity_check, new_sanity_check)

    # Fix the library searches
    old_library_check: str = "'stop; end program'"
    new_library_check: str = "'       PROGRAM MAIN\\n       END PROGRAM'"
    fortranpy_content = fortranpy_content.replace(old_library_check, new_library_check)

    with open(meson_fortran_path, "w") as fortranpy_fh:
        fortranpy_fh.write(fortranpy_content)

    # Fix the compiler.py file
    compilerpy_fh: TextIO
    with open(meson_compilers_path, "r") as compilerpy_fh:
        compilerpy_content: str = compilerpy_fh.read()

    # Fix the fortran flags
    old_fortran_flags: str = (
        "'fortran': ('f90', 'f95', 'f03', 'f08', 'f', 'for', 'ftn', 'fpp'),"
    )
    new_fortran_flags: str = (
        "'fortran': ('f', 'for', 'ftn', 'fpp', 'f90', 'f95', 'f03', 'f08'),"
    )
    compilerpy_content = compilerpy_content.replace(
        old_fortran_flags, new_fortran_flags
    )

    with open(meson_compilers_path, "w") as compilerpy_fh:
        compilerpy_fh.write(compilerpy_content)

    print("Meson has been patched to handle F77 correctly!")
