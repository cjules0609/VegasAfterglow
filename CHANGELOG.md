# Changelog

All notable changes to VegasAfterglow will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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