from setuptools import setup, Extension
import os
import platform
import pybind11

def find_sources():
    sources = ["pybind/pybind.cpp", "pybind/mcmc.cpp"]
    for root, _, files in os.walk("src"):
        for fn in files:
            if fn.endswith(".cpp"):
                sources.append(os.path.join(root, fn))
    return sources

system = platform.system()
archflags = os.environ.get("ARCHFLAGS", "")

# Base flags
extra_compile_args = []
extra_link_args = []

# Platform-specific settings
if system == "Linux":
    extra_compile_args = ["-std=c++20", "-O3", "-flto", "-w", "-DNDEBUG", "-fPIC", "-ffast-math"]
    extra_link_args = []

elif system == "Darwin":
    extra_compile_args = ["-std=c++20", "-O3", "-flto", "-w", "-DNDEBUG", "-fPIC", "-ffast-math"]
    extra_link_args = ["-undefined", "dynamic_lookup"]
    # Don't use -march=native for universal builds
    #if "-arch arm64" not in archflags or "-arch x86_64" not in archflags:
    #    extra_compile_args.append("-march=native")

elif system == "Windows":
    extra_compile_args = ["/std:c++20", "/O2", "/DNDEBUG", "/fp:fast"]
    extra_link_args = []

ext_modules = [
    Extension(
        "VegasAfterglow.VegasAfterglowC",
        sources=find_sources(),
        include_dirs=[
            pybind11.get_include(),
            "include",
            "external"
        ],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        language="c++",
    )
]

setup(ext_modules=ext_modules)