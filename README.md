# VegasAfterglow (under construction, stay tuned!)

<img align="left" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/logo.svg" alt="VegasAfterglow Logo" width="250"/>

[![C++ Version](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17) 
[![PyPI version](https://img.shields.io/pypi/v/VegasAfterglow.svg)](https://pypi.org/project/VegasAfterglow/) 
[![Build Status](https://github.com/YihanWangAstro/VegasAfterglow/actions/workflows/PyPI-build.yml/badge.svg)](https://github.com/YihanWangAstro/VegasAfterglow/actions/workflows/PyPI-build.yml) 
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE) 
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

VegasAfterglow offers a comprehensive suite of features for detailed GRB afterglow modeling:

**Shock Dynamics & Evolution:**
*   **Forward and Reverse Shock Modeling:** Simulates both shock components with arbitrary magnetization levels.
*   **Relativistic and Non-Relativistic Regimes:** Accurately models shock evolution across all velocity regimes.
*   **Ambient Medium:** Supports various ambient medium profiles including uniform Interstellar Medium (ISM), stellar wind environments, and user-defined density profiles.
*   **Energy and Mass Injection:** Supports user-defined profiles for continuous energy and/or mass injection into the blast wave.

**Jet Structure & Geometry:**
*   **Structured Jet Profiles:** Allows user-defined angular profiles for energy distribution, initial Lorentz factor, and magnetization.
*   **Jet Spreading:** Includes lateral expansion dynamics for realistic jet evolution.
*   **Non-Axisymmetric Jets:** Capable of modeling complex, non-axisymmetric jet structures.

**Radiation Mechanisms:**
*   **Synchrotron Radiation:** Calculates synchrotron emission from shocked electrons.
*   **Synchrotron Self-Absorption (SSA):** Includes SSA effects, crucial at low frequencies.
*   **Inverse Compton (IC) Scattering:** Models IC processes, including:
    *   Synchrotron Self-Compton (SSC) from both forward and reverse shocks.
    *   Pairwise IC between forward and reverse shock electron and photon populations.
    *   Includes Klein-Nishina corrections for accurate high-energy IC cooling and emission.

---

## Performance Highlights

VegasAfterglow delivers exceptional computational performance through deep optimization of its core algorithms:

*   **Ultra-fast Model Evaluation:** Generates a 30-point single-frequency light curve (forward shock & synchrotron only) in approximately **0.6 milliseconds** on an Apple M2 chip.
*   **Rapid MCMC Exploration:** Enables comprehensive parameter estimation (e.g., 10,000 MCMC steps for 8 parameters against 20 data points) in:
    *   Approximately **10 seconds** for on-axis structured jet scenarios.
    *   Approximately **30 seconds** for more complex off-axis scenarios.
*   **Optimized for Interactive Analysis:** Facilitates comprehensive Bayesian inference on standard laptop hardware in seconds to minutes, rather than hours or days, allowing for rapid iteration through different physical models and scenarios.

This level of performance is achieved through meticulous algorithm design, extensive use of vectorization, and careful memory management, making VegasAfterglow highly suitable for both detailed analysis of individual GRB events and large-scale population studies.

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

We provide an example of using MCMC to fit afterglow light curves and spectra to user-provided data. You can run it using either:

<details>
<summary><b>Option 1: Run with Jupyter Notebook</b> <i>(click to expand/collapse)</i></summary>
<br>

1. Install Jupyter Notebook:
```bash
pip install jupyter notebook
```

2. Launch Jupyter Notebook:
```bash
jupyter notebook
```

3. In your browser, open `mcmc.ipynb` inside the `script/` directory
</details>

<details>
<summary><b>Option 2: Run with VSCode + Jupyter Extension</b> <i>(click to expand/collapse)</i></summary>
<br>

1. Install [Visual Studio Code](https://code.visualstudio.com/) and the **Jupyter extension**:
   - Open VSCode
   - Go to the **Extensions** panel (or press `Cmd+Shift+X` on macOS, `Ctrl+Shift+X` on Windows)
   - Search for **"Jupyter"** and click **Install**

2. Open the VegasAfterglow folder in VSCode and navigate to `mcmc.ipynb` in the `script/` directory
</details>

### MCMC Parameter Fitting with VegasAfterglow

This section guides you through using the MCMC (Markov Chain Monte Carlo) module in VegasAfterglow to explore parameter space and determine posterior distributions for GRB afterglow models. Rather than just finding a single best-fit solution, MCMC allows you to quantify parameter uncertainties and understand correlations between different physical parameters.

<details open>
<summary><b>Overview</b> <i>(click to expand/collapse)</i></summary>
<br>

The MCMC module follows these key steps:
1. Create an `ObsData` object to hold your observational data
2. Configure model settings through the `Setups` class
3. Define parameters and their priors for the MCMC process
4. Run the MCMC sampler to explore the posterior distribution
5. Analyze and visualize the results

```python
from VegasAfterglow import ObsData, Setups, Fitter, ParamDef, Scale
```
</details>

<details>
<summary><b>1. Preparing Your Data</b> <i>(click to expand/collapse)</i></summary>
<br>

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
</details>

<details>
<summary><b>2. Configuring the Model</b> <i>(click to expand/collapse)</i></summary>
<br>

The `Setups` class defines the global properties and environment for your model. These settings remain fixed during the MCMC process.

```python
cfg = Setups()

# Source properties
cfg.lumi_dist = 3.364e28    # Luminosity distance [cm]  
cfg.z = 1.58               # Redshift

# Physical model configuration
cfg.medium = "wind"        # Ambient medium: "wind", "ISM" (Interstellar Medium) or "user" (user-defined)
cfg.jet = "powerlaw"       # Jet structure: "powerlaw", "gaussian", "tophat" or "user" (user-defined)

# Optional: Advanced grid settings. 
# cfg.phi_num = 24         # Number of grid points in phi direction
# cfg.theta_num = 24       # Number of grid points in theta direction
# cfg.t_num = 24           # Number of time grid points
```

**Why Configure These Properties?**
- **Source properties:** These parameters define the observer's relation to the source and are typically known from independent measurements
- **Physical model configuration:** These define the fundamental model choices that aren't fitted but instead represent different physical scenarios
- **Grid settings:** Control the numerical precision of the calculations (advanced users). Default is (24, 24, 24). VegasAfterglow is optimized to converge with this grid resolution for most cases.

These settings affect how the model is calculated but are not varied during the MCMC process, allowing you to focus on exploring the most relevant physical parameters.
</details>

<details>
<summary><b>3. Defining Parameters</b> <i>(click to expand/collapse)</i></summary>
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
</details>

<details>
<summary><b>4. Running the MCMC Fitting</b> <i>(click to expand/collapse)</i></summary>
<br>

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
<summary><b>5. Exploring the Posterior Distribution</b> <i>(click to expand/collapse)</i></summary>
<br>

Rather than focusing only on the best-fit parameters, examine the full posterior distribution:

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
</details>

<details>
<summary><b>6. Generating Model Predictions</b> <i>(click to expand/collapse)</i></summary>
<br>

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
</details>

<details>
<summary><b>7. Visualizing Results</b> <i>(click to expand/collapse)</i></summary>
<br>

#### Creating Corner Plots

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

#### Creating Trace Plots

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
</details>

<details>
<summary><b>8. Best Practices for Posterior Exploration</b> <i>(click to expand/collapse)</i></summary>
<br>

1. **Prior Ranges**: Set physically meaningful prior ranges based on theoretical constraints
2. **Convergence Testing**: Check convergence using trace plots and autocorrelation metrics
3. **Parameter Correlations**: Use corner plots to identify degeneracies and correlations
4. **Model Comparison**: Compare different physical models (e.g., wind vs. ISM) using Bayesian evidence
5. **Physical Interpretation**: Connect parameter constraints with physical processes in GRB afterglows
</details>

### Creating Custom Problem Generators with C++

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

VegasAfterglow is released under the **MIT License**.

---

## Citation

If you use VegasAfterglow in your research, please cite the relevant paper(s):


If you use specific modules or features that are described in other publications, please cite those as well according to standard academic practice.


