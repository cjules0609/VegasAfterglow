from setuptools import setup, Extension, find_packages
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
    extra_compile_args = ["-std=c++17", "-O3", "-flto", "-w", "-DNDEBUG", "-fPIC", "-ffast-math", "-pipe"]
    extra_link_args = []

elif system == "Darwin":
    extra_compile_args = ["-std=c++17", "-O3", "-flto", "-w", "-DNDEBUG", "-fPIC", "-ffast-math", "-pipe"]
    extra_link_args = ["-undefined", "dynamic_lookup"]
    # Don't use -march=native for universal builds
    #if "-arch arm64" not in archflags or "-arch x86_64" not in archflags:
    #    extra_compile_args.append("-march=native")

elif system == "Windows":
    extra_compile_args = ["/std:c++17", "/O2", "/DNDEBUG", "/fp:fast", "/MP", "/GL"]
    extra_link_args = ["/LTCG"]

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

package_data = {
    "VegasAfterglow": ["*.py"],
}

setup(
    name="VegasAfterglow",
    version="0.1.2",
    description="MCMC tools for astrophysics",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Yihan Wang, Connery Chen, Bing Zhang",

    license="BSD-3-Clause",                    
    data_files=[
        ("", ["LICENSE"]),
        ("assets", [os.path.join("assets", f) for f in os.listdir("assets") if not f.startswith(".")]),
    ],

    python_requires=">=3.8",
    install_requires=[
        "emcee>=3.0",
        "pybind11>=2.6.0",
        "corner>=2.2.1",
        "tqdm>=4.0",
        "numpy>=1.20",
        "scipy>=1.6",
        "pandas>=1.2"
    ],
    extras_require={
        "dev": ["ninja", "pytest", "black"],
    },

    packages=["VegasAfterglow"],
    package_dir={"VegasAfterglow": "pymodule"},
    package_data=package_data,
    include_package_data=True,

    ext_modules=ext_modules,
    zip_safe=False,
    )