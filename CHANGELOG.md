# Changelog

All notable changes to VegasAfterglow will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.0.0] - 2025-09-06

### Major Features & Enhancements

#### **Advanced MCMC Framework**
- **More model availability**: All built-in jet and medium models now supported in MCMC fitting, including reverse shock, inverse Compton, magnetar injection and etc. See documentations for detials.

#### **Adaptive Mesh Generation**
- **Dynamic Grid Optimization**: New adaptive angular grid generation based on jet properties and viewing angles. Grid points distributed according to Doppler boosting factors for optimal efficiency
- **Performance Gains**: ~5x faster convergence for reverse shocks.

#### **New Jet Models**
- **StepPowerLawJet**: Uniform core with sharp transition to power-law wings for realistic jet structures
- **Enhanced TwoComponentJet**: Separate narrow and wide components with independent energy and Lorentz factor profiles
- **Improved PowerLawJet**: Split power-law indices for energy (`k_e`) and Lorentz factor (`k_g`) angular dependence

### **Performance & Computational Improvements**

#### **Shock Physics Enhancements**
- **Variable Naming Standardization**: Major refactoring for clarity (`EAT_fwd` → `t_obs_fwd`, `Gamma_rel` → `Gamma_th`)
- **Reverse Shock Optimization**: Major code refactoring for reverse shock dynamics. A unified model for shock crossing and post-crossing evolution.
- **Better Crossing Dynamics**: Track the internal energy evolution during shock crossing for improved accuracy with more accurate shock heating and adiabatic cooling. **Note: The enhanced adiabatic cooling and detailed shock heating treatment leads to even weaker reverse shock emission compared to previous versions.**

#### **Numerical & Memory Optimizations**
- **Memory Efficiency**: Reduced memory footprint through optimized array operations and memory access patterns
- **Grid Pre-computation**: Enhanced caching strategies for frequently used calculations

###  **Enhanced Python Interface**

#### **Data Management Improvements**
- **MultiBandData Redesign**: Unified handling of light curves and spectra with flexible weighting
- **Series Calculations**: New methods for calculating flux at specific time-frequency pairs
- **Memory-Efficient Storage**: Optimized data structures for large multi-wavelength datasets

###  **Development & Build System**

#### **Code Quality Enhancements**
- **Enhanced Template System**: Improved compile-time type checking and memory management
- **Better Documentation**: Comprehensive API documentation with detailed parameter descriptions
- **Updated Examples**: New MCMC tutorials with real data fitting demonstrations

### **API Changes & Migration**

#### **Parameter Interface Updates**
- **PowerLawJet**: Split single `k` parameter into separate `k_e` (energy) and `k_g` (Lorentz factor) indices for more flexible modeling
- **TwoComponentJet**: Parameter names changed from `theta_n, E_iso_n, Gamma0_n` to `theta_c, E_iso_c, Gamma0_c` for consistency
- **StepPowerLawJet**: New jet model with parameters `(theta_c, E_iso_c, Gamma0_c, E_iso_w, Gamma0_w, k_e, k_g)`
- **Wind Medium**: Extended with optional `n_ism` and `n_0` parameters for stratified medium modeling: `Wind(A_star, n_ism=0, n_0=inf)`
- **Medium Class**: Simplified from `Medium(rho, mass)` to `Medium(rho)` - removed separate mass parameter
- **Model Methods**: 
  - Removed `specific_flux_sorted_series()` method
  - Added `specific_flux_series_with_expo(t, nu, expo_time, num_points=10)` for exposure time averaging
  - Changed `details(t_obs)` to `details(t_min, t_max)` interface
- **Resolution Parameters**: Default resolution changed from `(0.3, 3.0, 5.0)` to `(0.3, 1, 10)` for optimal performance/accuracy balance
- **MCMC Parameters**: Major restructuring with new parameters for all jet types, medium configurations, and reverse shock physics
  - Added: `theta_v, n_ism, n0, A_star, k_e, k_g, duration, E_iso_w, Gamma0_w, theta_w, L0, t0, q`
  - Added reverse shock parameters: `p_r, eps_e_r, eps_B_r, xi_e_r`
  - Removed: `k_jet` parameter (replaced by `k_e, k_g`)
- **MCMC Data Handling**: Enhanced `MultiBandData` with optional `weights` parameter for both light curves and spectra
- **MCMC Model Interface**: Replaced separate `light_curves()` and `spectra()` methods with unified `specific_flux()` method

---

**Migration Notes**: This release includes some API changes. Most existing code will work with minimal modifications. See the documentation for detailed migration guidance.

**⚠️ Important Physics Changes**: With the improved reverse shock dynamics, we now find that the reverse shock emission is even weaker than what we reported in our previous code paper, which itself was already weaker than the analytical scalings. This further reduction is due to enhanced adiabatic cooling and a more detailed treatment of shock heating. Users fitting reverse shock data may therefore need to re-evaluate their models and parameters.


## [v0.2.8] - 2025-08-15

### Improved

- **Inverse Compton Performance**: Major performance optimization (~10x speedup) for inverse Compton scattering calculations
- **Photon Interface**: Changed parameter name from `P_nu_max` to `I_nu_max` for better consistency in photon interface
- **Spectrum Smoothing**: Enhanced spectrum smoothing at `nu_max` for better numerical stability

### Documentation

- Updated documentation for detailed simulation quantities evolution
- Enhanced examples and API documentation

## [v0.2.7] - 2025-07-31

### Added

- **Internal Quantities Evolution Interface**: New Python interface to check the evolution of internal simulation quantities under various reference frames
- Comprehensive documentation for detailed simulation quantities evolution

### Fixed

- **Performance Issue**: Fixed jet edge detection function that was mistakenly set to 90 degrees, making the computational domain unnecessarily large regardless of the jet profile

### Documentation

- Updated documentation for detailed simulation quantities evolution
- Enhanced examples showing how to track shock dynamics and microphysical parameters

## [v0.2.6] - 2025-07-19

### Added

- **Built-in Two-Component Jet**: Added support for built-in two-component jet configurations

## [v0.2.5] - 2025-07-19

### Improved

- **Jet Edge Detection**: Improved jet edge detection algorithm for user-defined jet profiles

### Documentation

- Updated README with latest features
- Enhanced documentation with additional examples
- Updated API documentation

## [v0.2.4] - 2025-06-28

### Added

- **Magnetar Spin-down Documentation**: Added comprehensive documentation for magnetar spin-down energy injection

### Documentation

- Updated documentation for magnetar energy injection features
- Enhanced README with additional usage examples
- Improved API reference documentation

## [v0.2.3] - 2025-06-27

### Changed

- **API Breaking Change**: Modified Python-level user-defined medium/jet unit system for consistency
- **Parameter Naming**: Changed jet duration parameter name from `T0` (confusing to observer frame) to `duration`

### Improved

- **Data Handling**: Sorted series flux density for better data organization

## [v0.2.2] - 2025-06-23

### Changed

- **Default Resolution**: Changed default resolution settings when unspecified for better performance

### Improved

- **Reverse Shock**: Enhanced reverse shock smoothing algorithms

## [v0.2.1] - 2025-06-22

### Improved

- **Reverse Shock Modeling**:
  - Enhanced reverse shock smoothing
  - Improved post-crossing reverse shock calibration
  - Better single electron peak power calibration

## [v0.1.9] - 2025-06-19

### Changed

- **Python Support**: Removed support for Python 3.7 (minimum requirement now Python 3.8+)

### Fixed

- **Reverse Shock**: Significant corrections to reverse shock calculations

### Documentation

- Updated documentation and examples
- Enhanced API reference

## [v0.1.8] - 2025-06-09

### Improved

- **Reverse Shock**: Refined reverse shock calculations and algorithms
- **Code Quality**: General code cleanup and optimization

## [v0.1.7] - 2025-05-24

### Fixed

- **macOS Compatibility**: Fixed macOS-specific bugs
- **Deep Newtonian Regime**: Refinement for deep Newtonian regime calculations
  - Shocked electron calculations moved from energy space to momentum space
  - Enhanced self-absorption refinement
  - Interface improvements

### Improved

- **Code Quality**: General updates and optimizations

## [v0.1.6] - 2025-05-15

### Added

- **Energy Injection**: Added missing Python bindings header for energy injection functionality

### Improved

- **Magnetar Integration**: Enhanced C-level magnetar binding for better performance

## [v0.1.5] - 2025-05-14

### Added

- **Inverse Compton Enhancements**:
  - Updated IC spectrum calculations
  - Added magnetar injection Python interface
  - Enhanced IC code with cleanup and optimization

### Changed

- **Python Support**:
  - Removed Python 3.6 support
  - Added Python 3.7 support

### Fixed

- **Unit System**: Corrected unit handling in various calculations
- **Build System**: Various build fixes and improvements

### Documentation

- Updated examples in documentation
- Enhanced README with new features
- Improved API documentation

## [v0.1.4] - 2025-05-11

### Added

- **Enhanced Build System**: Improved CMake configuration and cross-platform support

### Fixed

- **Windows Compatibility**: Resolved various Windows build issues
- **External Dependencies**: Fixed external library integration

## [v0.1.3] - 2025-05-07

### Improved

- **Build System**: Enhanced CMake build configuration
- **Documentation**: Updated README and logo integration

### Fixed

- **Cross-platform Issues**: Resolved various build issues across platforms

## [v0.1.2] - 2025-05-06

### Added

- **Enhanced Documentation**: Added comprehensive documentation system
- **Logo Integration**: Added project logo and branding

### Improved

- **Code Structure**: Enhanced code organization and commenting
- **Build System**: Improved build workflow and packaging

## [v0.1.1] - 2025-05-06

### Fixed

- **Build System**: Fixed wheel building and packaging issues
- **Documentation**: Updated README with installation instructions

### Improved

- **Layout**: Enhanced README layout and organization

## [v0.1.0] - 2025-05-05

### Added

- **Initial Release**: First public release of VegasAfterglow
- **Core Features**:
  - High-performance C++ framework with Python interface
  - Comprehensive GRB afterglow modeling
  - Forward and reverse shock dynamics
  - Synchrotron and inverse Compton radiation
  - Structured jet configurations
  - MCMC parameter fitting capabilities
- **Radiation Mechanisms**:
  - Synchrotron emission with self-absorption
  - Inverse Compton scattering with Klein-Nishina corrections
  - Synchrotron Self-Compton (SSC)
- **Physical Models**:
  - Multiple ambient medium types (ISM, wind, user-defined)
  - Various jet structures (top-hat, power-law, Gaussian, user-defined)
  - Relativistic and non-relativistic shock evolution
  - Energy and mass injection
- **Performance**: Ultra-fast light curve computation (millisecond timescales)
- **Cross-platform Support**: Linux, macOS, and Windows compatibility
- **Python Interface**: User-friendly Python bindings for easy integration

### Documentation

- Comprehensive API documentation
- Installation guides for Python and C++
- Usage examples and tutorials
- Quick start guide with Jupyter notebooks

---

## Version History Summary

| Version | Release Date | Key Features |
|---------|--------------|--------------|
| v1.0.0  | 2025-09-06   | **Major Release**: Advanced MCMC framework, adaptive mesh generation, new jet models |
| v0.2.8  | 2025-08-15   | Inverse Compton performance optimization (~10x), interface improvements |
| v0.2.7  | 2025-07-31   | Internal quantities evolution interface, performance fixes |
| v0.2.6  | 2025-07-19   | Two-component jet support |
| v0.2.5  | 2025-07-19   | Enhanced jet edge detection |
| v0.2.4  | 2025-06-28   | Magnetar spin-down documentation |
| v0.2.3  | 2025-06-27   | Unit system updates, parameter naming |
| v0.2.2  | 2025-06-23   | Default resolution improvements |
| v0.2.1  | 2025-06-22   | Reverse shock enhancements |
| v0.1.9  | 2025-06-19   | Python 3.7 support removed, reverse shock fixes |
| v0.1.8  | 2025-06-09   | Reverse shock refinements |
| v0.1.7  | 2025-05-24   | macOS fixes, deep Newtonian regime improvements |
| v0.1.6  | 2025-05-15   | Magnetar integration, energy injection bindings |
| v0.1.5  | 2025-05-14   | Inverse Compton enhancements, Python support updates |
| v0.1.4  | 2025-05-11   | Build system improvements, Windows compatibility |
| v0.1.3  | 2025-05-07   | Enhanced build system, cross-platform fixes |
| v0.1.2  | 2025-05-06   | Documentation system, logo integration |
| v0.1.1  | 2025-05-06   | Build fixes, documentation updates |
| v0.1.0  | 2025-05-05   | Initial public release |

---

For detailed information about each release, see the individual version sections above.
For the latest updates and development progress, visit our [GitHub repository](https://github.com/YihanWangAstro/VegasAfterglow).