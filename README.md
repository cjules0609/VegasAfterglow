# VegasAfterglow (under construction, stay tuned!)

<img align="left" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/logo.svg" alt="VegasAfterglow Logo" width="250"/>

[![C++ Version](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17) 
[![PyPI version](https://img.shields.io/pypi/v/VegasAfterglow.svg)](https://pypi.org/project/VegasAfterglow/) 
[![Build Status](https://github.com/YihanWangAstro/VegasAfterglow/actions/workflows/PyPI-build.yml/badge.svg)](https://github.com/YihanWangAstro/VegasAfterglow/actions/workflows/PyPI-build.yml) 
[![License](https://img.shields.io/badge/License-BSD--3--Clause-green.svg)](LICENSE) 
[![Platform](https://img.shields.io/badge/Platform-Linux%20|%20macOS%20|%20Windows-lightgrey.svg)]() 
[![Python Version](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/) 

**VegasAfterglow** is a high-performance C++ framework designed for the comprehensive modeling of gamma-ray burst (GRB) afterglows. It achieves exceptional computational efficiency, enabling the generation of multi-wavelength light curves in milliseconds and facilitating robust Markov Chain Monte Carlo (MCMC) parameter inference in seconds to minutes. The framework incorporates advanced models for shock dynamics (both forward and reverse shocks), diverse radiation mechanisms (synchrotron with self-absorption, and inverse Compton scattering with Klein-Nishina corrections), and complex structured jet configurations. A user-friendly Python wrapper is provided to streamline integration into scientific data analysis workflows.

<br clear="left"/>

---

## Table of Contents

- [Features](#features)
- [Performance Highlights](#performance-highlights)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

---

## Features

<h3 align="center">Shock Dynamics & Evolution</h3>

<div style="display: flex; align-items: flex-start; gap: 20px;">
  <img src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/shock_dynamics.svg" width="400"/>
  
  <div>
    <ul>
      <li><strong>Forward and Reverse Shock Modeling:</strong> Simulates both shocks with arbitrary magnetization levels.</li>
      <li><strong>Relativistic and Non-Relativistic Regimes:</strong> Accurately models shock evolution across all velocity regimes.</li>
      <li><strong>Ambient Medium:</strong> Supports uniform Interstellar Medium (ISM), stellar wind environments, and user-defined density profiles.</li>
      <li><strong>Energy and Mass Injection:</strong> Supports user-defined profiles for continuous energy and/or mass injection into the blast wave.</li>
    </ul>
  </div>
</div>

<h3 align="center">Jet Structure & Geometry</h3>

<div style="display: flex; align-items: flex-start; gap: 20px;">
  <img src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/jet_geometry.svg" width="400"/>
  
  <div>
    <ul>
      <li><strong>Structured Jet Profiles:</strong> Allows user-defined angular profiles for energy distribution, initial Lorentz factor, and magnetization.</li>
      <li><strong>Arbitrary Viewing Angles:</strong> Supports off-axis observers at any viewing angle relative to the jet axis.</li>
      <li><strong>Jet Spreading:</strong> Includes lateral expansion dynamics for realistic jet evolution (experimental).</li>
      <li><strong>Non-Axisymmetric Jets:</strong> Capable of modeling complex, non-axisymmetric jet structures.</li>
    </ul>
  </div>
</div>

<h3 align="center">Radiation Mechanisms</h3>

<div style="display: flex; align-items: flex-start; gap: 20px;">
  <img src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/radiation_mechanisms.svg" width="400"/>

  <div>
    <ul>
      <li><strong>Synchrotron Radiation:</strong> Calculates synchrotron emission from shocked electrons.</li>
      <li><strong>Synchrotron Self-Absorption (SSA):</strong> Includes SSA effects, crucial at low frequencies.</li>
      <li><strong>Inverse Compton (IC) Scattering:</strong> Models IC processes, including:
        <ul>
          <li>Synchrotron Self-Compton (SSC) from both forward and reverse shocks.</li>
          <li>Pairwise IC between forward and reverse shock electron and photon populations (experimental).</li>
          <li>Includes Klein-Nishina corrections for accurate high-energy IC cooling and emission.</li>
        </ul>
      </li>
    </ul>
  </div>
</div>

---

## Performance Highlights

<img align="left" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/convergence_plot.png" width="300"/>

VegasAfterglow delivers exceptional computational performance through deep optimization of its core algorithms:

<div style="display: flex; align-items: flex-start; gap: 20px;">
  <div>
    <ul>
      <li><strong>Ultra-fast Model Evaluation:</strong> Generates a 30-point single-frequency light curve (forward shock & synchrotron only) in approximately 0.6 milliseconds on an Apple M2 chip with a single core.</li>
      <li><strong>Rapid MCMC Exploration:</strong> Enables parameter estimation with 10,000 MCMC steps in record time on an 8-core Apple M2 chip:
        <ul>
          <li>~10 seconds for on-axis structured jet scenarios</li>
          <li>~30 seconds for off-axis modeling scenarios</li>
        </ul>
      </li>
    </ul>
  </div>
</div>

This level of performance is achieved through optimized algorithm implementation and efficient memory access patterns, facilitating comprehensive Bayesian inference on standard laptop hardware in seconds to minutes rather than hours or days. The accelerated convergence speed enables rapid iteration through different physical models and makes VegasAfterglow suitable for both detailed analysis of individual GRB events and large-scale population studies.

<br clear="left"/>

---

## Installation

VegasAfterglow is available as a Python package with C++ source code also provided for direct use.

### Python Installation

To install VegasAfterglow using pip:

```bash
pip install VegasAfterglow
```

This is the recommended method for most users. VegasAfterglow requires Python 3.8 or higher.

<details>
<summary><b>Alternative: Install from Source</b> <i>(click to expand/collapse)</i></summary>
<br>

For cases where pip installation is not viable or when the development version is required:

1. Clone this repository:
```bash
git clone https://github.com/YihanWangAstro/VegasAfterglow.git
```

2. Navigate to the directory and install the Python package:
```bash
cd VegasAfterglow
pip install .
```

Standard development environments typically include the necessary prerequisites (C++17 compatible compiler). For build-related issues, refer to the prerequisites section in [C++ Installation](#c-installation).
</details>

### C++ Installation

For advanced users who need to compile and use the C++ library directly:

<details>
<summary><b>Instructions for C++ Installation</b> <i>(click to expand/collapse)</i></summary>
<br>

1. Clone the repository (if not previously done):
```bash
git clone https://github.com/YihanWangAstro/VegasAfterglow.git
cd VegasAfterglow
```

2. Compile the static library:
```bash
make lib
```

3. (Optional) Compile and run tests:
```bash
make tests
```

Upon successful compilation, you can create custom C++ problem generators using the VegasAfterglow interfaces. For implementation details, refer to the [Creating Custom Problem Generators with C++](#creating-custom-problem-generators-with-c) section or examine the example problem generators in `tests/demo/`.

<details>
<summary><b>Build Prerequisites</b> <i>(click to expand for dependency information)</i></summary>
<br>

The following development tools are required:

- **C++17 compatible compiler**:
  - **Linux**: GCC 7+ or Clang 5+
  - **macOS**: Apple Clang 10+ (with Xcode 10+) or GCC 7+ (via Homebrew)
  - **Windows**: MSVC 19.14+ (Visual Studio 2017 15.7+) or MinGW-w64 with GCC 7+
  
- **Build tools**:
  - Make (GNU Make 4.0+ recommended)
</details>
</details>

---

## Usage

We provide an example MCMC notebook (`script/mcmc.ipynb`) for fitting afterglow light curves and spectra to user-provided data. To avoid conflicts when updating the repository in the future, make a copy of the example notebook in the same directory and work with the copy instead of the original.

The notebook can be run using either Jupyter Notebook or VSCode with the Jupyter extension. Remember to keep your copy in the same directory as the original to ensure all data paths work correctly.

### MCMC Parameter Fitting with VegasAfterglow

This section shows how to use the MCMC module to explore parameter space and determine posterior distributions for GRB afterglow models.

<details>
<summary><b>1. Preparing Data and Configuring the Model</b> <i>(click to expand/collapse)</i></summary>
<br>

```python
from VegasAfterglow import ObsData, Setups, Fitter, ParamDef, Scale
```

VegasAfterglow provides flexible options for loading observational data through the `ObsData` class. You can add light curves (specific flux vs. time) and spectra (specific flux vs. frequency) in multiple ways.

```python
# Create an instance to store observational data
data = ObsData()

# Method 1: Add data directly from lists or numpy arrays

# For light curves
t_data = [1e3, 2e3, 5e3, 1e4, 2e4]  # Time in seconds
flux_data = [1e-26, 8e-27, 5e-27, 3e-27, 2e-27]  # Specific flux in erg/cm²/s/Hz
flux_err = [1e-28, 8e-28, 5e-28, 3e-28, 2e-28]  # Specific flux error in erg/cm²/s/Hz
data.add_light_curve(nu_cgs=4.84e14, t_cgs=t_data, Fnu_cgs=flux_data, Fnu_err=flux_err)

# For spectra
nu_data = [...]  # Frequencies in Hz
spectrum_data = [...] # Specific flux values in erg/cm²/s/Hz
spectrum_err = [...]   # Specific flux errors in erg/cm²/s/Hz
data.add_spectrum(t_cgs=3000, nu_cgs=nu_data, Fnu_cgs=spectrum_data, Fnu_err=spectrum_err)
```

```python
# Method 2: Load from CSV files
import pandas as pd

# Define your bands and files
bands = [2.4e17, 4.84e14]  # Example: X-ray, optical R-band
lc_files = ["data/ep.csv", "data/r.csv"]

# Load light curves from files
for nu, fname in zip(bands, lc_files):
    df = pd.read_csv(fname)
    data.add_light_curve(nu_cgs=nu, t_cgs=df["t"], Fnu_cgs=df["Fv_obs"], Fnu_err=df["Fv_err"])

times = [3000,6000] # Example: time in seconds
spec_files = ["data/spec_1.csv", "data/spec_2.csv"]

# Load spectra from files
for t, fname in zip(times, spec_files):
    df = pd.read_csv(fname)
    data.add_spectrum(t_cgs=t, nu_cgs=df["nu"], Fnu_cgs=df["Fv_obs"], Fnu_err=df["Fv_err"])
```

> **Note:** The `ObsData` interface is designed to be flexible. You can mix and match different data sources, and add multiple light curves at different frequencies as well as multiple spectra at different times.

The `Setups` class defines the global properties and environment for your model. These settings remain fixed during the MCMC process.

```python
cfg = Setups()

# Source properties
cfg.lumi_dist = 3.364e28    # Luminosity distance [cm]  
cfg.z = 1.58               # Redshift

# Physical model configuration
cfg.medium = "wind"        # Ambient medium: "wind", "ISM" (Interstellar Medium) or "user" (user-defined)
cfg.jet = "powerlaw"       # Jet structure: "powerlaw", "gaussian", "tophat" or "user" (user-defined)

# Optional: Advanced grid settings. Default (24, 24, 24) is optimized to converge for most cases.
# cfg.phi_num = 24         # Number of grid points in phi direction
# cfg.theta_num = 24       # Number of grid points in theta direction
# cfg.t_num = 24           # Number of time grid points
```

These settings affect how the model is calculated but are not varied during the MCMC process.
</details>

<details>
<summary><b>2. Defining Parameters and Running MCMC</b> <i>(click to expand/collapse)</i></summary>
<br>

The `ParamDef` class is used to define the parameters for MCMC exploration. Each parameter requires a name, initial value, prior range, and sampling scale:

```python
mc_params = [
    ParamDef("E_iso",    1e52,  1e50,  1e54,  Scale.LOG),       # Isotropic energy [erg]
    ParamDef("Gamma0",     30,     5,  1000,  Scale.LOG),       # Lorentz factor at the core
    ParamDef("theta_c",   0.2,   0.0,   0.5,  Scale.LINEAR),    # Core half-opening angle [rad]
    ParamDef("theta_v",   0.,  None,  None,   Scale.FIXED),     # Viewing angle [rad]
    ParamDef("p",         2.5,     2,     3,  Scale.LINEAR),    # Shocked electron power law index
    ParamDef("eps_e",     0.1,  1e-2,   0.5,  Scale.LOG),       # Electron energy fraction
    ParamDef("eps_B",    1e-2,  1e-4,   0.5,  Scale.LOG),       # Magnetic field energy fraction
    ParamDef("A_star",   0.01,  1e-3,     1,  Scale.LOG),       # Wind parameter
    ParamDef("xi",        0.5,  1e-3,     1,  Scale.LOG),       # Electron acceleration fraction
]
```

**Scale Types:**
- `Scale.LOG`: Sample in logarithmic space (log10) - ideal for parameters spanning multiple orders of magnitude
- `Scale.LINEAR`: Sample in linear space - appropriate for parameters with narrower ranges
- `Scale.FIXED`: Keep parameter fixed at the initial value - use for parameters you don't want to vary

**Parameter Choices:**
The parameters you include depend on your model configuration:
- For "wind" medium: use `A_star` parameter 
- For "ISM" medium: use `n_ism` parameter instead
- Different jet structures may require different parameters

Initialize the `Fitter` class with your data and configuration, then run the MCMC process:

```python
# Create the fitter object
fitter = Fitter(data, cfg)

# Run the MCMC fitting
result = fitter.fit(
    param_defs=mc_params,          # Parameter definitions
    resolution=(24, 24, 24),       # Grid resolution (phi, theta, time)
    total_steps=10000,             # Total number of MCMC steps
    burn_frac=0.3,                 # Fraction of steps to discard as burn-in
    thin=1                         # Thinning factor
)
```

The `result` object contains:
- `samples`: The MCMC chain samples (posterior distribution)
- `labels`: Parameter names
- `best_params`: Maximum likelihood parameter values
</details>

<details>
<summary><b>3. Analyzing Results and Generating Predictions</b> <i>(click to expand/collapse)</i></summary>
<br>

Examine the posterior distribution to understand parameter constraints:

```python
# Print best-fit parameters (maximum likelihood)
print("Best-fit parameters:")
for name, val in zip(result.labels, result.best_params):
    print(f"  {name}: {val:.4f}")

# Compute median and credible intervals
flat_chain = result.samples.reshape(-1, result.samples.shape[-1])
medians = np.median(flat_chain, axis=0)
lower = np.percentile(flat_chain, 16, axis=0)
upper = np.percentile(flat_chain, 84, axis=0)

print("\nParameter constraints (median and 68% credible intervals):")
for i, name in enumerate(result.labels):
    print(f"  {name}: {medians[i]:.4f} (+{upper[i]-medians[i]:.4f}, -{medians[i]-lower[i]:.4f})")
```

Use samples from the posterior to generate model predictions with uncertainties:

```python
# Define time and frequency ranges for predictions
t_out = np.logspace(2, 9, 150)
bands = [2.4e17, 4.84e14] 

# Generate light curves with the best-fit model
lc_best = fitter.light_curves(result.best_params, t_out, bands)

nu_out = np.logspace(6, 20, 150)
times = [3000]
# Generate model spectra at the specified times using the best-fit parameters
spec_best = fitter.spectra(result.best_params, nu_out, times)

# Now you can plot the best-fit model and the uncertainty envelope
```

Corner plots are essential for visualizing parameter correlations and posterior distributions:

```python
import corner

def plot_corner(flat_chain, labels, filename="corner_plot.png"):
    fig = corner.corner(
        flat_chain,
        labels=labels,
        quantiles=[0.16, 0.5, 0.84],  # For median and ±1σ
        show_titles=True,
        title_kwargs={"fontsize": 14},
        label_kwargs={"fontsize": 14},
        truths=np.median(flat_chain, axis=0),  # Show median values
        truth_color='red',
        bins=30,
        fill_contours=True,
        levels=[0.16, 0.5, 0.68],  # 1σ and 2σ contours
        color='k'
    )
    fig.savefig(filename, dpi=300, bbox_inches='tight')

# Create the corner plot
flat_chain = result.samples.reshape(-1, result.samples.shape[-1])
plot_corner(flat_chain, result.labels)
```

Trace plots help verify MCMC convergence:

```python
def plot_trace(chain, labels, filename="trace_plot.png"):
    nsteps, nwalkers, ndim = chain.shape
    fig, axes = plt.subplots(ndim, figsize=(10, 2.5 * ndim), sharex=True)

    for i in range(ndim):
        for j in range(nwalkers):
            axes[i].plot(chain[:, j, i], alpha=0.5, lw=0.5)
        axes[i].set_ylabel(labels[i])
        
    axes[-1].set_xlabel("Step")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)

# Create the trace plot
plot_trace(result.samples, result.labels)
```

**Best Practices for Posterior Exploration:**
1. **Prior Ranges**: Set physically meaningful prior ranges based on theoretical constraints
2. **Convergence Testing**: Check convergence using trace plots and autocorrelation metrics
3. **Parameter Correlations**: Use corner plots to identify degeneracies and correlations
4. **Model Comparison**: Compare different physical models (e.g., wind vs. ISM)
5. **Physical Interpretation**: Connect parameter constraints with physical processes in GRB afterglows
</details>

---

## Advanced Usage: C++ Interface

<details>
<summary><b>Creating Custom Problem Generators with C++</b> <i>(for advanced users)</i></summary>
<br>

After compiling the library, you can create custom applications that use VegasAfterglow's core functionality:

<details>
<summary><b>Working with the C++ API</b> <i>(click to expand/collapse)</i></summary>
<br>

To use VegasAfterglow directly from C++, you will typically perform the following steps:

### 1. Include necessary headers

```cpp
#include "afterglow.h"              // Afterglow models
```

### 2. Define your problem configuration

```cpp
// Example configuration code will be added in future documentation
```

### 3. Compute radiation and create light curves/spectra

```cpp
// Example light curve calculation code will be added in future documentation
```

### 4. Building Custom Applications

Compile your C++ application, linking against the VegasAfterglow library you built:
```bash
g++ -std=c++17 -I/path/to/VegasAfterglow/include -L/path/to/VegasAfterglow/lib -o my_program my_program.cpp -lvegasafterglow
```
Replace `/path/to/VegasAfterglow/` with the actual path to your VegasAfterglow installation.

### 5. Example Problem Generators

The repository includes several example problem generators in the `tests/demo/` directory that demonstrate different use cases. These are the best resources for understanding direct C++ API usage.
</details>
</details>

---

## Contributing

If you encounter any issues, have questions about the code, or want to request new features:

1. **GitHub Issues** - The most straightforward and fastest way to get help:
   - Open an issue at [https://github.com/YihanWangAstro/VegasAfterglow/issues](https://github.com/YihanWangAstro/VegasAfterglow/issues)
   - You can report bugs, suggest features, or ask questions
   - This allows other users to see the problem/solution as well
   - Can be done anonymously if preferred

2. **Pull Requests** - If you've implemented a fix or feature:
   - Fork the repository
   - Create a branch for your changes
   - Submit a pull request with your changes

We value all contributions and aim to respond to issues promptly.

---

## License

VegasAfterglow is released under the **BSD-3-Clause License**.

---

## Citation

If you use VegasAfterglow in your research, please cite the relevant paper(s):


If you use specific modules or features that are described in other publications, please cite those as well according to standard academic practice.


