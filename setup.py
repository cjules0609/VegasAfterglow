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

    zlib_dir = os.path.join("external", "zlib")
    for fn in os.listdir(zlib_dir):
        if fn.endswith(".c"):
            sources.append(os.path.join(zlib_dir, fn))

    return sources

system = platform.system()
archflags = os.environ.get("ARCHFLAGS", "")

# Base flags
extra_compile_args = []
extra_link_args = []

# Platform-specific settings
if system == "Linux":
    extra_compile_args = ["-std=c++20", "-O3", "-march=native", "-flto", "-w", "-DNDEBUG", "-fPIC", "-ffast-math"]
    extra_link_args = []

elif system == "Darwin":
    extra_compile_args = ["-std=c++20", "-O3", "-flto", "-w", "-DNDEBUG", "-fPIC", "-ffast-math"]
    extra_link_args = ["-undefined", "dynamic_lookup"]
    # Don't use -march=native for universal builds
    if "-arch arm64" not in archflags or "-arch x86_64" not in archflags:
        extra_compile_args.append("-march=native")

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
            "external",
            "external/zlib"
        ],
        extra_compile_args=extra_compile_args,
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