# VegasAfterglow

<img align="left" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/logo.svg" alt="VegasAfterglow Logo" width="270"/>

[![C++ Version](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![PyPI version](https://img.shields.io/pypi/v/VegasAfterglow.svg)](https://pypi.org/project/VegasAfterglow/)
[![Build Status](https://github.com/YihanWangAstro/VegasAfterglow/actions/workflows/PyPI-build.yml/badge.svg)](https://github.com/YihanWangAstro/VegasAfterglow/actions/workflows/PyPI-build.yml)
[![License](https://img.shields.io/badge/License-BSD--3--Clause-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Linux%20|%20macOS%20|%20Windows-lightgrey.svg)]()
[![Python Version](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Documentation](https://img.shields.io/badge/Documentation-Online-brightgreen.svg)](https://yihanwangastro.github.io/VegasAfterglow/docs/index.html)

**VegasAfterglow** is a high-performance C++ framework with a user-friendly Python interface designed for the comprehensive modeling of gamma-ray burst (GRB) afterglows. It achieves exceptional computational efficiency, enabling the generation of multi-wavelength light curves in milliseconds and facilitating robust Markov Chain Monte Carlo (MCMC) parameter inference in seconds to minutes. The framework incorporates advanced models for shock dynamics (both forward and reverse shocks), diverse radiation mechanisms (synchrotron with self-absorption, and inverse Compton scattering with Klein-Nishina corrections), and complex structured jet configurations.

<br clear="left"/>

---

## Table of Contents

- [VegasAfterglow](#vegasafterglow)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Performance Highlights](#performance-highlights)
  - [Installation](#installation)
    - [Python Installation](#python-installation)
    - [C++ Installation](#c-installation)
  - [Usage](#usage)
    - [Quick Start](#quick-start)
    - [Light Curve \& Spectrum Calculation](#light-curve--spectrum-calculation)
    - [MCMC Parameter Fitting](#mcmc-parameter-fitting)
  - [Documentation](#documentation)
  - [Contributing](#contributing)
  - [License](#license)
  - [Acknowledgments \& Citation](#acknowledgments--citation)

---

## Features

<h3 align="center">Shock Dynamics</h3>

<img align="right" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/shock_dynamics.svg" width="450"/>

- **Forward and Reverse Shock Modeling:** Simulates both shocks via shock crossing dynamics with arbitrary magnetization levels and shell thicknesses.
- **Relativistic and Non-Relativistic Regimes:** Accurately models shock evolution across all velocity regimes.
- **Adiabatic and Radiative Blast Wave:** Supports smooth transition between adiabatic and radiative blast wave.
- **Ambient Medium:** Supports uniform Interstellar Medium (ISM), stellar wind environments, and user-defined density profiles.
- **Energy and Mass Injection:** Supports user-defined profiles for continuous energy and/or mass injection into the blast wave.

<br clear="right"/>

<h3 align="center">Jet Structure & Geometry</h3>

<img align="right" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/jet_geometry.svg" width="450"/>

- **Structured Jet Profiles:** Allows user-defined angular profiles for energy distribution, initial Lorentz factor, and magnetization.
- **Arbitrary Viewing Angles:** Supports off-axis observers at any viewing angle relative to the jet axis.
- **Jet Spreading:** Includes lateral expansion dynamics for realistic jet evolution (experimental).
- **Non-Axisymmetric Jets:** Capable of modeling complex, non-axisymmetric jet structures.

<br clear="right"/>

<h3 align="center">Radiation Mechanisms</h3>

<img align="right" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/radiation_mechanisms.svg" width="450"/>

- **Synchrotron Radiation:** Calculates synchrotron emission from shocked electrons.
- **Synchrotron Self-Absorption (SSA):** Includes SSA effects, crucial at low frequencies.
- **Inverse Compton (IC) Scattering:** Models IC processes, including:
  - Synchrotron Self-Compton (SSC) from both forward and reverse shocks.
  - Pairwise IC between forward and reverse shock electron and photon populations (experimental).
  - Includes Klein-Nishina corrections for accurate synchrotron and IC emission.

<br clear="right"/>

---

## Performance Highlights

<img align="right" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/convergence_plot.png" width="400"/>

VegasAfterglow delivers exceptional computational performance through deep optimization of its core algorithms:

- **Ultra-fast Light Curve Computation:** Generates a 100-point single-frequency light curve (forward shock & synchrotron only) from a structured jet viewed off-axis in approximately 1 milliseconds on an Apple M2 chip with a single core.

- **Rapid MCMC Exploration:** Enables parameter estimation with 10,000 MCMC steps for 8 parameters on 20 data points across multi-wavelength light curves and spectra on an 8-core Apple M2 chip in:
  - ~20 seconds for on-axis structured jet scenarios
  
This level of performance is achieved through optimized algorithm implementation and efficient memory access patterns, facilitating comprehensive Bayesian inference on standard laptop hardware in seconds to minutes rather than hours or days. The accelerated convergence speed enables rapid iteration through different physical models and makes VegasAfterglow suitable for both detailed analysis of individual GRB events and large-scale population studies.

<br clear="right"/>

---

## Installation

VegasAfterglow is available as a Python package with C++ source code also provided for direct use.

### Python Installation

To install VegasAfterglow using pip:

```bash
pip install VegasAfterglow
```

This is the recommended method for most users. VegasAfterglow requires Python 3.7 or higher.

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

Standard development environments typically include the necessary prerequisites (C++20 compatible compiler). For build-related issues, refer to the prerequisites section in [C++ Installation](#c-installation).
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

- **C++20 compatible compiler**:
  - **Linux**: GCC 10+ or Clang 13+
  - **macOS**: Apple Clang 13+ (with Xcode 13+) or GCC 10+ (via Homebrew)
  - **Windows**: MSVC 19.29+ (Visual Studio 2019 16.10+) or MinGW-w64 with GCC 10+
  
- **Build tools**:
  - Make (GNU Make 4.0+ recommended)

</details>
</details>

---

## Usage

### Quick Start

We provide basic example scripts (`script/quick.ipynb` and `script/mcmc.ipynb`) that demonstrate how to set up and run afterglow simulations. This section shows how to calculate light curves and spectra for a simple GRB afterglow model without the need for observational data and perform MCMC parameter fitting with observational data. The notebook can be run using either Jupyter Notebook or VSCode with the Jupyter extension.

To avoid conflicts when updating the repository in the future, make a copy of the example notebook in the same directory and work with the copy instead of the original.

### Light Curve & Spectrum Calculation

The example below walks through the main components needed to model a GRB afterglow, from setting up the physical parameters to producing light curves and spectra.

<details>
<summary><b>Model Setup</b> <i>(click to expand/collapse)</i></summary>
<br>

First, let's set up the physical components of our afterglow model, including the environment, jet, observer, and radiation parameters:

```python
import numpy as np
import matplotlib.pyplot as plt
from VegasAfterglow import ISM, TophatJet, Observer, Radiation, Model

# 1. Define the circumburst environment (constant density ISM)
medium = ISM(n_ism=1) #in cgs unit

# 2. Configure the jet structure (top-hat with opening angle, energy, and Lorentz factor)
jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300) #in cgs unit

# 3. Set observer parameters (distance, redshift, viewing angle)
obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0) #in cgs unit

# 4. Define radiation microphysics parameters
rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)

# 5. Combine all components into a complete afterglow model
model = Model(jet=jet, medium=medium, observer=obs, forward_rad=rad)
```

</details>

<details>
<summary><b>Light Curve Calculation</b> <i>(click to expand/collapse)</i></summary>
<br>

Now, let's compute and plot multi-wavelength light curves to see how the afterglow evolves over time:

```python
# 1. Create logarithmic time array from 10² to 10⁸ seconds (100s to ~3yrs)
times = np.logspace(2, 8, 100)  

# 2. Define observing frequencies (radio, optical, X-ray bands in Hz)
bands = np.array([1e9, 1e14, 1e17])  

# 3. Calculate the afterglow emission at each time and frequency
results = model.specific_flux(times, bands)

# 4. Visualize the multi-wavelength light curves
plt.figure(figsize=(4.8, 3.6),dpi=200)

# 5. Plot each frequency band 
for i, nu in enumerate(bands):
    exp = int(np.floor(np.log10(nu)))
    base = nu / 10**exp
    plt.loglog(times, results['syn'][i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz')

def add_note(plt):
    plt.annotate('jet break',xy=(3e4, 1e-26), xytext=(3e3, 5e-28), arrowprops=dict(arrowstyle='->'))
    plt.annotate(r'$\nu_m=\nu_a$',xy=(8e5, 2e-25), xytext=(7.5e4, 5e-24), arrowprops=dict(arrowstyle='->'))
    plt.annotate(r'$\nu=\nu_a$',xy=(4e6, 4e-25), xytext=(7.5e5, 5e-24), arrowprops=dict(arrowstyle='->'))

add_note(plt)   
plt.xlabel('Time (s)')
plt.ylabel('Flux Density (erg/cm²/s/Hz)')
plt.legend()
plt.title('Light Curves')
plt.savefig('quick-lc.png',dpi=300,bbox_inches='tight')
```

<div align="center">
<img src="assets/quick-lc.png" alt="Afterglow Light Curves" width="600"/>

Running the light curve script will produce this figure showing the afterglow evolution across different frequencies.
</div>
</details>

<details>
<summary><b>Spectrum Analysis</b> <i>(click to expand/collapse)</i></summary>
<br>

We can also examine how the broadband spectrum evolves at different times after the burst:

```python
# 1. Define broad frequency range (10⁵ to 10²² Hz) 
frequencies = np.logspace(5, 22, 100)  

# 2. Select specific time epochs for spectral snapshots 
epochs = np.array([1e2, 1e3, 1e4, 1e5 ,1e6, 1e7, 1e8])

# 3. Calculate spectra at each epoch
results = model.specific_flux(epochs, frequencies)

# 4. Plot broadband spectra at each epoch
plt.figure(figsize=(4.8, 3.6),dpi=200)
colors = plt.cm.viridis(np.linspace(0,1,len(epochs)))

for i, t in enumerate(epochs):
    exp = int(np.floor(np.log10(t)))
    base = t / 10**exp
    plt.loglog(frequencies, results['syn'][:,i], color=colors[i], label=fr'${base:.1f} \times 10^{{{exp}}}$ s')

# 5. Add vertical lines marking the bands from the light curve plot
for i, band in enumerate(bands):
    exp = int(np.floor(np.log10(band)))
    base = band / 10**exp
    plt.axvline(band,ls='--',color='C'+str(i))

plt.xlabel('frequency (Hz)')
plt.ylabel('flux density (erg/cm²/s/Hz)')
plt.legend(ncol=2)
plt.title('Synchrotron Spectra')
plt.savefig('quick-spec.png',dpi=300,bbox_inches='tight')
```

<div align="center">
<img src="assets/quick-spec.png" alt="Broadband Spectra" width="600"/>

The spectral analysis code will generate this visualization showing spectra at different times, with vertical lines indicating the frequencies calculated in the light curve example.
</div>
</details>

These examples demonstrate the core functionality of VegasAfterglow for modeling GRB afterglows. The code is designed to be highly efficient, allowing for rapid exploration of parameter space and comparison with observational data.

### MCMC Parameter Fitting

We provide some example data files in the `data` folder. Remember to keep your copy in the same directory as the original to ensure all data paths work correctly.

<details>
<summary><b>1. Preparing Data and Configuring the Model</b> <i>(click to expand/collapse)</i></summary>
<br>

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import corner
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

data = ObsData()
# Define your bands and files
bands = [2.4e17, 4.84e14, 1.4e14]  # Example: X-ray, optical R-band
lc_files = ["data/ep.csv", "data/r.csv", "data/vt-r.csv"]

# Load light curves from files
for nu, fname in zip(bands, lc_files):
    df = pd.read_csv(fname)
    data.add_light_curve(nu_cgs=nu, t_cgs=df["t"], Fnu_cgs=df["Fv_obs"], Fnu_err=df["Fv_err"])

times = [3000] # Example: time in seconds
spec_files = ["data/ep-spec.csv"]

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
cfg.medium = "wind"        # Ambient medium: "wind", "ism" (Interstellar Medium) or "user" (user-defined)
cfg.jet = "powerlaw"       # Jet structure: "powerlaw", "gaussian", "tophat" or "user" (user-defined)
```

These settings affect how the model is calculated but are not varied during the MCMC process.
</details>

<details>
<summary><b>2. Defining Parameters and Running MCMC</b> <i>(click to expand/collapse)</i></summary>
<br>

The `ParamDef` class is used to define the parameters for MCMC exploration. Each parameter requires a name, prior range, and sampling scale:

```python
mc_params = [
    ParamDef("E_iso",      1e50,  1e54,  Scale.LOG),       # Isotropic energy [erg]
    ParamDef("Gamma0",        5,  1000,  Scale.LOG),       # Lorentz factor at the core
    ParamDef("theta_c",     0.0,   0.5,  Scale.LINEAR),    # Core half-opening angle [rad]
    ParamDef("theta_v",     0.0,   0.0,  Scale.FIXED),     # Viewing angle [rad]
    ParamDef("p",             2,     3,  Scale.LINEAR),    # Shocked electron power law index
    ParamDef("eps_e",      1e-2,   0.5,  Scale.LOG),       # Electron energy fraction
    ParamDef("eps_B",      1e-4,   0.5,  Scale.LOG),       # Magnetic field energy fraction
    ParamDef("A_star",     1e-3,     1,  Scale.LOG),       # Wind parameter
    ParamDef("xi",         1e-3,     1,  Scale.LOG),       # Electron acceleration fraction
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

Check the best-fit parameters and their uncertainties:

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

Use the best-fit parameters to generate model predictions

```python
# Define time and frequency ranges for predictions
t_out = np.logspace(2, 9, 150)
bands = [2.4e17, 4.84e14, 1.4e14] 

# Generate light curves with the best-fit model
lc_best = fitter.light_curves(result.best_params, t_out, bands)

nu_out = np.logspace(6, 20, 150)
times = [3000]
# Generate model spectra at the specified times using the best-fit parameters
spec_best = fitter.spectra(result.best_params, nu_out, times)
```

Now you can plot the best-fit model:

```python
def draw_bestfit(t, lc_fit, nu, spec_fit):
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.5, 7.5))
    
    # Plot light curves
    shifts = [1, 1, 200]
    colors = ['blue', 'orange', 'green']
    
    for i in range(len(lc_files)):
        df = pd.read_csv(lc_files[i])
        ax1.errorbar(df["t"], df["Fv_obs"] * shifts[i], df["Fv_err"] * shifts[i], 
                    fmt='o', color=colors[i], label=lc_files[i])
        ax1.plot(t, np.array(lc_fit[i]) * shifts[i], color=colors[i], lw=1)

    # Plot spectra
    for i in range(len(spec_files)):
        df = pd.read_csv(spec_files[i])
        ax2.errorbar(df["nu"], df["Fv_obs"] * shifts[i], df["Fv_err"] * shifts[i], 
                    fmt='o', color=colors[i], label=spec_files[i])
        ax2.plot(nu, np.array(spec_fit[0]) * shifts[i], color=colors[i], lw=1)

    # Configure axes
    for ax, xlabel, ylabel in [(ax1, 't [s]', r'$F_\nu$ [erg/cm$^2$/s/Hz]'),
                              (ax2, r'$\nu$ [Hz]', r'$F_\nu$ [erg/cm$^2$/s/Hz]')]:
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.legend()

    plt.tight_layout()

draw_bestfit(t_out, lc_best, nu_out, spec_best)
```

Corner plots are essential for visualizing parameter correlations and posterior distributions:

```python
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
        smooth=1,
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

---

## Documentation

Comprehensive documentation is available at **[Documentation](https://yihanwangastro.github.io/VegasAfterglow/docs/index.html)** including:

- **Installation Guide**: Detailed instructions for setting up VegasAfterglow
- **Examples**: Practical examples showing common use cases
- **Python API Reference**: Complete documentation of the Python interface
- **C++ API Reference**: Detailed documentation of C++ classes and functions
- **Contributing Guide**: Information for developers who wish to contribute

The documentation is regularly updated with the latest features and improvements (NOT officially released).

---

## Contributing

If you encounter any issues, have questions about the code, or want to request new features:

1. **GitHub Issues** - The most straightforward and fastest way to get help:
   - Open an issue at [Issues](https://github.com/YihanWangAstro/VegasAfterglow/issues)
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

The BSD 3-Clause License is a permissive open source license that allows you to:

- Freely use, modify, and distribute the software in source and binary forms
- Use the software for commercial purposes
- Integrate the software into proprietary applications

**Requirements:**

- You must include the original copyright notice and the license text
- You cannot use the names of the authors or contributors to endorse derived products
- The license provides no warranty or liability protection

For the full license text, see the [LICENSE](LICENSE) file in the repository.

---

## Acknowledgments & Citation

We would like to thank the contributors who helped improve VegasAfterglow. **Special thanks to Weihua Lei, Shaoyu Fu, Iris Yin, Cuiyuan Dai and Binbin Zhang** for their invaluable work as beta testers, providing feedback and helping with bug fixes during development. We also thank the broader community for their suggestions and support.

If you use VegasAfterglow in your research, please cite the relevant paper(s):

If you use specific modules or features that are described in other publications, please cite those as well according to standard academic practice.
