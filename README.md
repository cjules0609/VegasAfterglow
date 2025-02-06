# VegasAfterglow

**VegasAfterglow** is a **high-performance** C++ framework for modeling **gamma-ray burst (GRB) afterglows**. It supports both **relativistic and non-relativistic** regimes and provides a **flexible, user-configurable** approach for computing light curves and spectra. The framework includes sophisticated shock dynamics, radiation mechanisms, and structured jet models.

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

VegasAfterglow requires **C++17** and **Make**. To compile:

```sh
git clone https://github.com/YihanWangAstro/VegasAfterglow.git
cd VegasAfterglow
make -j$(nproc)  # Use all available cores
```

## Usage

### Simulation Pipeline

1. **Shock Dynamics**  
   - Compute the forward and reverse shock evolution.
   
2. **Electron Distribution**  
   - Solve for the electron energy distribution.

3. **Photon Emission**  
   - Compute synchrotron and inverse Compton radiation.

4. **Observation**  
   - Generate light curves and spectra for a given observer.

### Example Simulation

```cpp
size_t r_num = 32, theta_num = 32, phi_num = 32;
double n_ism = 1 / con::cm3, eps_e = 0.1, eps_B = 1e-3, p = 2.3;
double E_iso = 1e53 * con::erg, Gamma0 = 300, theta_c = 0.1, theta_v = 0.2;

// Define Medium & Jet
auto medium = createISM(n_ism);
auto jet = GaussianJet(theta_c, E_iso, Gamma0);

// Create computational grid
Array t_obs = logspace(1e3 * con::sec, 1e7 * con::sec, 50);
Coord coord = adaptiveGrid(medium, jet, inject::none, t_obs, 0.6, phi_num, theta_num, r_num);

// Generate shocks
Shock f_shock = genForwardShock(coord, medium, jet, inject::none, eps_e, eps_B);

// Compute electron distribution
auto syn_e = genSynElectrons(f_shock, p);

// Compute emitted photons
auto syn_ph = genSynPhotons(f_shock, syn_e);

// Setup observer
Observer obs(coord);
double z = 0.009, lumi_dist = 1.23e26 * con::cm;
obs.observe(f_shock, theta_v, lumi_dist, z);

// Compute flux
double nu_obs = eVtoHz(1 * con::keV);
Array F_nu = obs.specificFlux(t_obs, nu_obs, syn_ph);

// Output results

output(t_obs, "t_obs", con::sec);
output(F_nu, "light_curve", con::erg / con::cm / con::cm / con::sec);

```

## Directory Structure

```
VegasAfterglow/
│── include/             # Header files
│── src/                 # Implementation files
│── tests/               # Unit tests & validation cases
│── boost/               # Header-only Boost.Odeint (manually included)
│── examples/            # Example simulations
│── docs/                # Documentation
│── Makefile             # Build system
│── README.md            # Project documentation
```

## Compilation

The project uses a **Makefile** for building. To compile:

```sh
make -j$(nproc)
```

To clean the build:

```sh
make clean
```

## Contributing

Contributions are welcome! Feel free to open an issue or submit a pull request.

## License

VegasAfterglow is released under the **MIT License**.
