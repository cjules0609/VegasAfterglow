# VegasAfterglow

**VegasAfterglow** is a **high-performance** C++ framework for modeling **gamma-ray burst (GRB) afterglows**. It supports both **relativistic and non-relativistic** regimes and provides a **flexible, user-configurable** approach for computing light curves and spectra. The framework includes sophisticated shock dynamics, radiation mechanisms, and structured jet models. A Python wrapper is included to streamline configuration and usability.

<div align="center">
<img src="https://img.shields.io/badge/C%2B%2B-20-blue.svg" alt="C++ Version">
<img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
<img src="https://img.shields.io/badge/Platform-Linux%20|%20macOS%20|%20Windows-lightgrey.svg" alt="Platform">
</div>

## Features

- **Forward and Reverse Shock Modeling**  
  - Arbitrary magnetization for both shocks.
  - Works in both relativistic and non-relativistic regimes.

- **Structured Jet with User-Defined Profiles**  
  - Supports custom **energy distribution**, **Lorentz factor**, and **magnetization** profiles.
  - **Jet Spreading** is included for realistic dynamics.
  - **Non-Axisymmetric Jets** allow complex jet structures.

- **Energy Injection Mechanisms**  
  - Allows user-defined energy injection profiles.

- **Synchrotron Radiation with Self-Absorption**  
  - Includes **synchrotron self-absorption (SSA)**.

- **Inverse Compton Scattering (IC)**  
  - Supports **forward SSC**, **reverse SSC**, and **pairwise IC** between forward and reverse shock electrons and photons.
  - Includes **Klein-Nishina corrections** for IC cooling.

## Prerequisites

VegasAfterglow requires the following to build and run:
<details>
<summary>System requirements (expand for details)</summary>

> These are typically provided by your operating system. You can skip this section unless your OS is out of date.

- **C++20 compatible compiler**:
  - **Linux**: GCC 10+ or Clang 10+
  - **macOS**: Apple Clang 12+ (with Xcode 12+) or GCC 10+ (via Homebrew)
  - **Windows**: MSVC 19.27+ (Visual Studio 2019 16.7+) or MinGW-w64 with GCC 10+
  
- **Build tools**:
  - Make (GNU Make 4.0+ recommended) [if you want to compile & run the C++ code]
</details>

## Installation

VegasAfterglow is available as a Python package with C++ source code also provided for direct use.

### Python Installation

<table>
<tr>
<td width="50%" valign="top">

#### Option 1: Install from PyPI (Recommended)

The simplest way to install VegasAfterglow is from PyPI:

```bash
pip install vegasafterglow
```

</td>
<td width="50%" valign="top">

#### Option 2: Install from Source

1. Clone this repository:
```bash
git clone https://github.com/YihanWangAstro/VegasAfterglow.git
```

2. *(Optional)* Create and activate a virtual environment:
```bash
python -m venv vegasafterglow
source vegasafterglow/bin/activate      # On Windows: vegasafterglow\Scripts\activate
```

3. Navigate to the file directory and install the Python package
```bash
cd VegasAfterglow
pip install -e .
```

</td>
</tr>
</table>

### C++ Installation

<details>
<summary>For advanced users who want to compile and use the C++ library directly</summary>

1. Clone the repository (if you haven't already):
```bash
git clone https://github.com/YihanWangAstro/VegasAfterglow.git
cd VegasAfterglow
```

2. Compile the static library:
```bash
make lib
```
Then you can write your own C++ problem generator and use the provided VegasAfterglow interfaces to do afterglow modeling by linking the VegasAfterglow library. (See more details in the [C++ API] section, or just take a look at the example problem generator under `tests/demo/` to see how to use the VegasAfterglow interfaces.)

3. (Optional) Compile and run tests:
```bash
make tests
```
</details>

## Usage

We provide an example of using MCMC to fit afterglow light curves and spectra to user-provided data. You can run it using either:

<table>
<tr>
<td width="50%" valign="top">

### Option 1: Run with Jupyter Notebook

1. Install Jupyter Notebook
```bash
pip install jupyter notebook
```

2. Launch Jupyter Notebook:
```bash
jupyter notebook
```

3. In your browser, open `mcmc.ipynb` inside the `script/` directory

</td>
<td width="50%" valign="top">

### Option 2: Run with VSCode + Jupyter Extension

1. Install [Visual Studio Code](https://code.visualstudio.com/) and the **Jupyter extension**:
  - Open VSCode
  - Go to the **Extensions** panel (or press `Cmd+Shift+X` on MacOS, `Ctrl+Shift+X` on Windows)
  - Search for **"Jupyter"** and click **Install**

2. Open the VegasAfterglow folder in VSCode:
  - Go to **File > Open Folder** and select the cloned `VegasAfterglow` repository.

3. Open `mcmc.ipynb` inside the `script/` directory.

</td>
</tr>
</table>

---

## MCMC Parameter Fitting with VegasAfterglow

This section guides you through using the MCMC (Markov Chain Monte Carlo) module in VegasAfterglow to explore parameter space and determine posterior distributions for GRB afterglow models. Rather than just finding a single best-fit solution, MCMC allows you to quantify parameter uncertainties and understand correlations between different physical parameters.

### Performance Highlights

<div style="background-color: #f0f8ff; padding: 15px; border-radius: 5px; border-left: 5px solid #4682b4;">

VegasAfterglow delivers exceptional computational performance through deep optimization of its core algorithms:

- **Ultra-fast model evaluation**: Generates a 30-point single-frequency light curve (forward shock & synchrotron only) in just **0.6ms** on an Apple M2 chip
- **Rapid MCMC exploration**: Complete 10,000-step parameter estimation with 8 parameters against 20 data points in:
  - **10 seconds** for on-axis GRB cases
  - **30 seconds** for more complex off-axis cases
- **Optimized for interactive analysis**: Perform comprehensive Bayesian inference in seconds rather than hours or days, enabling rapid iteration through different physical scenarios

This extreme performance comes from careful algorithm design, vectorization, and memory optimization, making VegasAfterglow suitable for both individual event analysis and large population studies.
</div>

### Overview

The MCMC module follows these key steps:
1. Create an `ObsData` object to hold your observational data
2. Configure model settings through the `Setups` class
3. Define parameters and their priors for the MCMC process
4. Run the MCMC sampler to explore the posterior distribution
5. Analyze and visualize the results

```python
from VegasAfterglow import ObsData, Setups, Fitter, ParamDef, Scale
```

### 1. Preparing Your Data

VegasAfterglow provides flexible options for loading observational data through the `ObsData` class. You can add light curves (flux vs. time) and spectra (flux vs. frequency) in multiple ways.

<details open>
<summary><b>Using the ObsData Interface:</b></summary>

```python
# Create an instance to store observational data
data = ObsData()

# Method 1: Add data directly from arrays
import numpy as np

# For light curves
t_data = np.array([1e3, 2e3, 5e3, 1e4, 2e4])  # Time in seconds
flux_data = np.array([1e-26, 8e-27, 5e-27, 3e-27, 2e-27])  # Specific flux in erg/cm²/s/Hz
flux_err = np.array([1e-28, 8e-28, 5e-28, 3e-28, 2e-28])  # Specific flux error
data.add_light_curve(nu_cgs=4.84e14, t_cgs=t_data, Fnu_cgs=flux_data, Fnu_err=flux_err)

# For spectra
nu_data = np.logspace(9, 18, 10)  # Frequencies in Hz
spectrum_data = np.array([...])  # Specific flux values in erg/cm²/s/Hz
spectrum_err = np.array([...])   # Specific flux errors in erg/cm²/s/Hz
data.add_spectrum(t_cgs=3000, nu_cgs=nu_data, Fnu_cgs=spectrum_data, Fnu_err=spectrum_err)

# Method 2: Load from CSV files
import pandas as pd

# Define your bands and files
bands = [2.4e17, 4.84e14]  # Example: X-ray, optical R-band
lc_files = ["data/ep.csv", "data/r.csv"]

# Load light curves from files
for nu, fname in zip(bands, lc_files):
    df = pd.read_csv(fname)
    data.add_light_curve(nu_cgs=nu, t_cgs=df["t"], Fnu_cgs=df["Fv_obs"], Fnu_err=df["Fv_err"])
```

> **Note:** The `ObsData` interface is designed to be flexible. You can mix and match different data sources, and add multiple light curves at different frequencies as well as multiple spectra at different times.
</details>

### 2. Configuring the Model

The `Setups` class defines the global properties and environment for your model. These settings remain fixed during the MCMC process.

<details open>
<summary><b>Setting up the model configuration:</b></summary>

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

<div style="background-color: #f0fff0; padding: 10px; border-radius: 5px; border-left: 5px solid #2e8b57;">

**Why Configure These Properties?**
- **Source properties:** These parameters define the observer's relation to the source and are typically known from independent measurements
- **Physical model configuration:** These define the fundamental model choices that aren't fitted but instead represent different physical scenarios
- **Grid settings:** Control the numerical precision of the calculations (advanced users). Default is (24, 24, 24). VegasAfterglow is optimized to converge with this grid resolution for most cases.
</div>

These settings affect how the model is calculated but are not varied during the MCMC process, allowing you to focus on exploring the most relevant physical parameters.
</details>

### 3. Defining Parameters

The `ParamDef` class is used to define the parameters for MCMC exploration. Each parameter requires a name, initial value, prior range, and sampling scale:

<details open>
<summary><b>Setting up MCMC parameters:</b></summary>

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

<div style="background-color: #fff5ee; padding: 10px; border-radius: 5px; border-left: 5px solid #ff7f50;">

**Scale Types:**
- `Scale.LOG`: Sample in logarithmic space (log10) - ideal for parameters spanning multiple orders of magnitude
- `Scale.LINEAR`: Sample in linear space - appropriate for parameters with narrower ranges
- `Scale.FIXED`: Keep parameter fixed at the initial value - use for parameters you don't want to vary
</div>

<div style="background-color: #f5f5f5; padding: 10px; border-radius: 5px; border-left: 5px solid #708090;">

**Parameter Choices:**
The parameters you include depend on your model configuration:
- For "wind" medium: use `A_star` parameter 
- For "ISM" medium: use `n_ism` parameter instead
- Different jet structures may require different parameters
</div>
</details>

### 4. Running the MCMC Fitting

Initialize the `Fitter` class with your data and configuration, then run the MCMC process:

```python
# Create the fitter object
fitter = Fitter(data, cfg)

# Run the MCMC fitting
result = fitter.fit(
    param_defs=mc_params,          # Parameter definitions
    resolution=(24, 24, 24),       # Grid resolution (phi, theta, time)
    total_steps=10000,            # Total number of MCMC steps
    burn_frac=0.3,                 # Fraction of steps to discard as burn-in
    thin=1                         # Thinning factor
)
```

<div style="background-color: #fffaf0; padding: 10px; border-radius: 5px; border-left: 5px solid #ffa500;">

The `result` object contains:
- `samples`: The MCMC chain samples (posterior distribution)
- `labels`: Parameter names
- `best_params`: Maximum likelihood parameter values
</div>

### 5. Exploring the Posterior Distribution

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

### 6. Generating Model Predictions

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

### 7. Visualizing Results

<details>
<summary><b>Creating Corner Plots</b></summary>

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
</details>

<details>
<summary><b>Creating Trace Plots</b></summary>

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

### 8. Tips for Effective Posterior Exploration

<div style="background-color: #f0f0ff; padding: 10px; border-radius: 5px; border-left: 5px solid #6a5acd;">

1. **Prior Ranges**: Set physically meaningful prior ranges based on theoretical constraints
2. **Convergence Testing**: Check convergence using trace plots and autocorrelation metrics
3. **Parameter Correlations**: Use corner plots to identify degeneracies and correlations
4. **Model Comparison**: Compare different physical models (e.g., wind vs. ISM) using Bayesian evidence
5. **Physical Interpretation**: Connect parameter constraints with physical processes in GRB afterglows
</div>

## Directory Structure

```
VegasAfterglow/
│── external/            # External dependencies
│── include/             # Header files
│── pybind/              # Python wrapper files
│── python/              # Python class definitions
│── script/              # MCMC script
│── src/                 # Implementation files
│── tests/               # Unit tests & validation cases
│── docs/                # Documentation
│── Makefile             # Build system
│── pyproject.toml       # Python metadata and build system configuration
│── README.md            # Project documentation
│── setup.py             # Python setup script
```

## Contributing

Contributions are welcome! Feel free to open an issue or submit a pull request.

## License

VegasAfterglow is released under the **MIT License**.


