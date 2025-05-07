import os
import sys
import subprocess
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        # Ensure CMake is installed
        try:
            subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError("CMake is required to build the C++ extensions")
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        import pybind11  # build-system ensures this is available

        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cfg = "Debug" if self.debug else "Release"
        build_temp = os.path.join(self.build_temp, ext.name)
        os.makedirs(build_temp, exist_ok=True)

        # CMake configure arguments
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-Dpybind11_DIR={pybind11.get_cmake_dir()}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
        ]
        # Prefer Ninja if it’s installed
        try:
            subprocess.check_output(["ninja", "--version"])
            cmake_args += ["-G", "Ninja"]
        except OSError:
            pass

        # Run CMake configure and build
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_temp)
        subprocess.check_call(
            ["cmake", "--build", ".", "--config", cfg, "--parallel", str(os.cpu_count())],
            cwd=build_temp,
        )

setup(
    name="VegasAfterglow",
    version="0.1.3",
    author="Yihan Wang",
    author_email="yihan.wang@example.com",
    description="MCMC tools for astrophysics",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    package_dir={"VegasAfterglow": "pymodule"},                     # your Python folder
    packages=find_packages(where="pymodule"),
    ext_modules=[CMakeExtension("VegasAfterglow.VegasAfterglowC")],
    cmdclass={"build_ext": CMakeBuild},
    install_requires=[
        "numpy>=1.20",
        "scipy>=1.6",
        "pandas>=1.2",
        "emcee>=3.0",
        "corner>=2.2.1",
        "tqdm>=4.0",
    ],
    python_requires=">=3.8",
    include_package_data=True,
    zip_safe=False,
)