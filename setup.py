# setup.py (at root)
from setuptools import setup, find_packages, Extension
import os, platform

def find_sources():
    sources = ["pybind/pybind.cpp", "pybind/mcmc.cpp"]
    for root, _, files in os.walk("src"):
        for fn in files:
            if fn.endswith(".cpp"):
                sources.append(os.path.join(root, fn))
    return sources

extra_compile_args = ["-std=c++20", "-O3", "-march=native", "-flto", "-w", "-DNDEBUG", "-fPIC", "-ffast-math"]
extra_link_args    = ["-lz"]
if platform.system() == "Darwin":
    extra_link_args += ["-undefined", "dynamic_lookup"]
elif platform.system() == "Linux":
    extra_compile_args.append("-fPIC")

ext_modules = [
    Extension(
        "VegasAfterglow.VegasAfterglowC",
        find_sources(),
        include_dirs=["include", "external"],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        language="c++",
    )
]

setup(
    name="VegasAfterglow",
    version="0.1.0",
    description="MCMC tools for astrophysics",
    author="Your Name",
    author_email="you@example.com",

    # <-- here’s the trick:
    packages=["VegasAfterglow"],
    package_dir={"VegasAfterglow": "python"},

    ext_modules=ext_modules,
    install_requires=[
        "numpy>=1.19",
        "pandas>=1.1",
        "emcee>=3.0",
        "pybind11>=2.6.0",
    ],
    python_requires=">=3.7",
)