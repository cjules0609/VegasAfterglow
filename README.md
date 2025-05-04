# VegasAfterglow

**VegasAfterglow** is a **high-performance** C++ framework for modeling **gamma-ray burst (GRB) afterglows**. It supports both **relativistic and non-relativistic** regimes and provides a **flexible, user-configurable** approach for computing light curves and spectra. The framework includes sophisticated shock dynamics, radiation mechanisms, and structured jet models. A Python wrapper is included to streamline configuration and usability.

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

## Installation

VegasAfterglow C++ source code is precompiled and wrapped in Python for immediate usage! 

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

## Usage

VegasAfterglow includes a Jupyter Notebook script to fit afterglow light curves and spectra to user-provided data using MCMC. You can run it using either:

- **Option 1:** Jupyter Notebook (standalone)
- **Option 2:** Visual Studio Code (VSCode) with the Jupyter extension

Choose one of the two methods below:

---

### ✅ Option 1: Run with **Jupyter Notebook**

1. Install Jupyter Notebook
```bash
pip install jupyter notebook
```

2. Launch Jupyter Notebook:
```bash
jupyter notebook
```

3. In your browser, open `mcmc.ipynb` inside the `scripts/` directory

### ✅ Option 2: Run with **VSCode + Jupyter Extension**

1. Install [Visual Studio Code](https://code.visualstudio.com/) and the **Jupyter extension**:
  - Open VSCode
  - Go to the **Extensions** panel (or press `Cmd+Shift+X` on MacOS, `Ctrl+Shift+X` on Windows)
  - Search for **"Jupyter"** and click **Install**

2. Open the VegasAfterglow folder in VSCode:
  - Go to **File > Open Folder** and select the cloned `VegasAfterglow` repository.

3. Open `mcmc.ipynb` inside the `scripts/` directory.

---

Within `mcmc.ipynb`, follow the inline comments for instructions on running the MCMC script.

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
