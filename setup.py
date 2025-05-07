import os
import sys
import subprocess
import platform
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        # Make sure CMake is available
        try:
            subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError("CMake is required to build the C++ extension")
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        import pybind11

        # Where to put the final .so / .pyd
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cfg = "Debug" if self.debug else "Release"

        build_temp = os.path.join(self.build_temp, ext.name)
        os.makedirs(build_temp, exist_ok=True)

        # Base CMake args
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-Dpybind11_DIR={pybind11.get_cmake_dir()}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
        ]

        # === ARCHITECTURE SELECTION ===
        host_arch = platform.machine().lower()
        if sys.platform == "darwin":
            # On macOS, tell CMake which arch to build for
            cmake_args.append(f"-DCMAKE_OSX_ARCHITECTURES={host_arch}")
        elif sys.platform.startswith("win"):
            # On Windows with Visual Studio generators, pass -A
            # CMake accepts "ARM64" or "x64"
            win_arch = "ARM64" if "arm" in host_arch else "x64"
            cmake_args += ["-A", win_arch]
        elif sys.platform.startswith("linux"):
            # On Linux aarch64 hosts, ensure we compile for armv8-a
            # On x86_64 hosts this is a no-op
            cmake_args += [
                f"-DCMAKE_C_FLAGS=-march={host_arch}",
                f"-DCMAKE_CXX_FLAGS=-march={host_arch}"
            ]
        # Configure
        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args,
            cwd=build_temp,
        )
        # Build in parallel
        build_cmd = ["cmake", "--build", ".", "--config", cfg, "--parallel", str(os.cpu_count())]
        subprocess.check_call(build_cmd, cwd=build_temp)

setup(
    name="VegasAfterglow",
    version="0.1.3",
    author="Yihan Wang",
    author_email="yihan.astro@gmail.com",
    description="MCMC tools for astrophysics",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    packages=["VegasAfterglow"],
    package_dir={"VegasAfterglow": "pymodule"},
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