import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py


class build_py(_build_py):
    """Custom build step: run the Fortran Makefile before packaging."""

    def run(self):
        # Path to src/moria/fortran_compile
        here = Path(__file__).parent
        fc_dir = here / "src" / "moria" / "fortran_compile"

        # Run the Makefile to build .xOg executables
        # Adjust F77/flags inside Makefile.txt as needed; we don't touch Fortran.
        subprocess.run(
            ["make", "-f", "Makefile.txt"],
            cwd=fc_dir,
            check=True,
        )

        # Then run the normal build
        super().run()


setup(cmdclass={"build_py": build_py})