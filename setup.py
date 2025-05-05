from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import platform
import pybind11

def find_sources():
    cpp_sources = ["pybind/pybind.cpp", "pybind/mcmc.cpp"]
    c_sources = []

    for root, _, files in os.walk("src"):
        for fn in files:
            if fn.endswith(".cpp"):
                cpp_sources.append(os.path.join(root, fn))

    zlib_dir = os.path.join("external", "zlib")
    for fn in os.listdir(zlib_dir):
        if fn.endswith(".c"):
            c_sources.append(os.path.join(zlib_dir, fn))

    return cpp_sources + c_sources, c_sources

extra_compile_args = []
extra_link_args = []

sources, c_sources = find_sources()
system = platform.system()
archflags = os.environ.get("ARCHFLAGS", "")

if system == "Linux":
    compile_args_cpp = ["-std=c++20", "-O3", "-march=native", "-flto", "-w", "-DNDEBUG", "-fPIC", "-ffast-math"]
    compile_args_c = ["-O3", "-flto", "-w", "-DNDEBUG", "-fPIC"]

elif system == "Darwin":
    compile_args_cpp = ["-std=c++20", "-O3", "-flto", "-w", "-DNDEBUG", "-fPIC", "-ffast-math"]
    compile_args_c = ["-O3", "-flto", "-w", "-DNDEBUG", "-fPIC"]
    if "-arch arm64" not in archflags or "-arch x86_64" not in archflags:
        compile_args_cpp.append("-march=native")
        compile_args_c.append("-march=native")
    extra_link_args = ["-undefined", "dynamic_lookup"]

elif system == "Windows":
    compile_args_cpp = ["/std:c++20", "/O2", "/DNDEBUG", "/fp:fast"]
    compile_args_c = ["/O2", "/DNDEBUG"]
    extra_link_args = []
else:
    compile_args_cpp = []
    compile_args_c = []
    extra_link_args = []
class CustomBuildExt(build_ext):
    def build_extensions(self):
        for ext in self.extensions:
            objects = []
            new_sources = []
            for src in ext.sources:
                is_c_file = src.endswith(".c")
                flags = compile_args_c if is_c_file else compile_args_cpp
                flags = flags.copy()  # Important: don't mutate global list
                if is_c_file and self.compiler.compiler_type == "msvc":
                    flags.append("/TC")  # <-- Force C mode for MSVC

                obj = self.compiler.compile(
                    [src],
                    output_dir=self.build_temp,
                    include_dirs=ext.include_dirs,
                    extra_postargs=flags,
                    depends=ext.depends,
                )
                objects.extend(obj)
                if not is_c_file:
                    new_sources.append(src)

            ext.extra_objects = objects
            ext.sources = new_sources  # Keep only C++ sources for linking
        build_ext.build_extensions(self)

ext_modules = [
    Extension(
        "VegasAfterglow.VegasAfterglowC",
        sources=sources,
        include_dirs=[
            pybind11.get_include(),
            "include",
            "external",
            "external/zlib"
        ],
        extra_link_args=extra_link_args,
        language="c++",
    )
]

setup(
    name="VegasAfterglow",
    version="0.1.0",
    description="MCMC tools for astrophysics",
    author="Yihan Wang, Connery Chen & Bing Zhang",
    author_email="yihan.astro@gmail.com",
    packages=["VegasAfterglow"],
    package_dir={"VegasAfterglow": "python"},
    ext_modules=ext_modules,
    cmdclass={"build_ext": CustomBuildExt},
    install_requires=[
        "numpy>=1.19",
        "pandas>=1.1",
        "emcee>=3.0",
        "pybind11>=2.6.0",
        "corner>=2.2.1",
        "tqdm>=4.0"
    ],
    python_requires=">=3.8",
)