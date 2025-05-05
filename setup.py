from setuptools import setup, Extension
import os
import platform
import pybind11
import multiprocessing

def find_sources():
    sources = ["pybind/pybind.cpp", "pybind/mcmc.cpp"]
    for root, _, files in os.walk("src"):
        for fn in files:
            if fn.endswith(".cpp"):
                sources.append(os.path.join(root, fn))
    return sources

system = platform.system()
archflags = os.environ.get("ARCHFLAGS", "")

# Determine number of parallel build jobs
cpu_count = multiprocessing.cpu_count()
parallel_jobs = max(1, cpu_count - 1)  # Use N-1 cores for compilation

# Base flags
extra_compile_args = []
extra_link_args = []

# Platform-specific settings
if system == "Linux":
    extra_compile_args = [
        "-std=c++20", 
        "-O3", 
        "-flto=auto", 
        "-w", 
        "-DNDEBUG", 
        "-fPIC", 
        "-ffast-math",
        # The following flags can improve compilation speed
        "-pipe",
        f"-j{parallel_jobs}"
    ]
    extra_link_args = ["-flto=auto"]

elif system == "Darwin":
    extra_compile_args = [
        "-std=c++20", 
        "-O3", 
        "-flto=thin", 
        "-w", 
        "-DNDEBUG", 
        "-fPIC", 
        "-ffast-math",
        # The following flags can improve compilation speed
        "-pipe"
    ]
    extra_link_args = ["-undefined", "dynamic_lookup", "-flto=thin"]
    # Apple's clang doesn't support -j flag directly like this

elif system == "Windows":
    extra_compile_args = [
        "/std:c++20", 
        "/O2", 
        "/DNDEBUG", 
        "/fp:fast",
        "/MP",  # Parallel compilation
        "/GL"   # Whole program optimization
    ]
    extra_link_args = ["/LTCG"]  # Link-time code generation

# Set environment variables for faster builds
os.environ["CFLAGS"] = "-O3"
os.environ["CXXFLAGS"] = "-O3"

# Enable ccache if available
if system != "Windows" and os.system("which ccache > /dev/null 2>&1") == 0:
    os.environ["CC"] = "ccache gcc"
    os.environ["CXX"] = "ccache g++"

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

# Pass the parallel flag to setuptools if supported
setup(
    ext_modules=ext_modules,
    options={'build_ext': {'parallel': parallel_jobs}}
)