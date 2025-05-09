Examples
========

.. contents:: Table of Contents
   :local:
   :depth: 2

Basic Usage
----------

Setting up a simple afterglow model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from VegasAfterglow import Model, TophatJet, ISM, Synchrotron

    # Set up the model components
    jet = TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300)
    medium = ISM(n_ism=1.0)
    radiation = Synchrotron(epsilon_e=0.1, epsilon_B=0.01, p=2.2)

    # Create the model
    model = Model(jet=jet, medium=medium, radiation=radiation)

    # Define times and frequencies for calculation
    times = np.logspace(2, 7, 100)  # 100 seconds to 10^7 seconds
    frequencies = np.array([1e9, 1e14, 1e17])  # Radio, optical, X-ray

    # Run the simulation
    results = model.calculate_light_curves(times, frequencies)

    # Plot the results
    plt.figure(figsize=(10, 6))
    for i, nu in enumerate(frequencies):
        plt.loglog(times, results[:, i], label=f'{nu:.1e} Hz')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()
    plt.title('GRB Afterglow Light Curves')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

Structured Jet Models
------------------

Gaussian Jet
^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import GaussianJet

    # Create a structured jet with Gaussian energy profile
    gaussian_jet = GaussianJet(
        theta_c=0.05,         # Core angular size (radians)
        E_iso=1e53,           # Isotropic-equivalent energy (ergs)
        Gamma0=300            # Initial Lorentz factor
    )

    # Update the model with the structured jet
    model.set_jet(gaussian_jet)
    
    # Off-axis viewing angle
    model.set_viewing_angle(0.2)  # 0.2 radians off-axis
    
    # Recalculate with the structured jet
    results_gaussian = model.calculate_light_curves(times, frequencies)

Power-Law Jet
^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import PowerLawJet

    # Create a power-law structured jet
    powerlaw_jet = PowerLawJet(
        theta_c=0.05,         # Core angular size (radians)
        E_iso=1e53,           # Isotropic-equivalent energy (ergs)
        Gamma0=300,           # Initial Lorentz factor
        k=2.0                 # Power-law index
    )

    # Update the model with the power-law jet
    model.set_jet(powerlaw_jet)
    
    # Recalculate with the power-law jet
    results_powerlaw = model.calculate_light_curves(times, frequencies)

Ambient Media Models
-----------------

Wind Medium
^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Wind

    # Create a stellar wind medium
    wind = Wind(A_star=0.1)  # A* parameter

    # Update the model with the wind medium
    model.set_medium(wind)
    
    # Recalculate with the wind medium
    results_wind = model.calculate_light_curves(times, frequencies)

User-Defined Medium
^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import UserDefinedMedium

    # Define a custom density profile function
    def custom_density(phi, theta, r):
        # Example: A medium with a density cavity
        r_cavity = 1e17  # Cavity radius in cm
        rho_0 = 1.67e-24  # Base density in g/cm³
        
        if r < r_cavity:
            return 0.1 * rho_0  # Lower density inside cavity
        else:
            return rho_0 * (r/r_cavity)**(-2)  # Wind-like outside
    
    # Create a user-defined medium
    custom_medium = UserDefinedMedium(density_func=custom_density)
    
    # Update the model
    model.set_medium(custom_medium)
    
    # Recalculate with the custom medium
    results_custom = model.calculate_light_curves(times, frequencies)

Radiation Processes
----------------

Synchrotron Self-Compton
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import SynchrotronSelfCompton

    # Create a model with synchrotron self-Compton
    ssc = SynchrotronSelfCompton(
        epsilon_e=0.1,
        epsilon_B=1e-3,  # Lower magnetization favors IC
        p=2.2,
        include_KN=True  # Include Klein-Nishina effects
    )
    
    # Update the model
    model.set_radiation(ssc)
    
    # Calculate over a broader frequency range to capture IC component
    frequencies_broad = np.logspace(9, 24, 50)  # Radio to gamma-rays
    
    # Calculate spectrum at a specific time
    t_spec = 1e4  # 10,000 seconds
    spectrum = model.calculate_spectrum(t_spec, frequencies_broad)
    
    # Plot the spectrum with components
    plt.figure(figsize=(10, 6))
    plt.loglog(frequencies_broad, spectrum, 'b-', label='Total')
    plt.loglog(frequencies_broad, model.get_synchrotron_spectrum(), 'r--', label='Synchrotron')
    plt.loglog(frequencies_broad, model.get_ic_spectrum(), 'g--', label='Inverse Compton')
    
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()
    plt.title(f'GRB Afterglow Spectrum at t = {t_spec} s')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

Advanced Features
---------------

Reverse Shock
^^^^^^^^^^^

.. code-block:: python

    # Create a model with reverse shock component
    model_with_rs = Model(
        jet=jet, 
        medium=medium, 
        radiation=radiation,
        include_reverse_shock=True
    )
    
    # Set reverse shock parameters
    model_with_rs.set_reverse_shock_parameters(
        RB=0.1,  # Magnetic field ratio between reverse and forward shock
        Re=1.0   # Electron energy ratio between reverse and forward shock
    )
    
    # Calculate light curves including reverse shock
    results_with_rs = model_with_rs.calculate_light_curves(times, frequencies)
    
    # Plot forward vs reverse shock components
    plt.figure(figsize=(10, 6))
    for i, nu in enumerate(frequencies):
        plt.loglog(times, results_with_rs[:, i], label=f'Total {nu:.1e} Hz')
        plt.loglog(times, model_with_rs.get_forward_shock_light_curve(i), '--', 
                  label=f'FS {nu:.1e} Hz')
        plt.loglog(times, model_with_rs.get_reverse_shock_light_curve(i), ':', 
                  label=f'RS {nu:.1e} Hz')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()
    plt.title('GRB Afterglow with Reverse Shock')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

MCMC Parameter Fitting
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import ObsData, Fitter, ParamDef, Scale

    # Create observation data object
    data = ObsData()

    # Add some observational data (light curves)
    t_data = np.array([1e3, 2e3, 5e3, 1e4, 2e4])  # Time in seconds
    flux_data = np.array([1e-26, 8e-27, 5e-27, 3e-27, 2e-27])  # Specific flux
    flux_err = np.array([1e-28, 8e-28, 5e-28, 3e-28, 2e-28])  # Flux error
    
    # Add a light curve at optical frequency (5e14 Hz)
    data.add_light_curve(nu=5e14, t=t_data, flux=flux_data, flux_err=flux_err)
    
    # Define parameters with priors
    params = [
        ParamDef("E_iso", 51.0, 54.0, Scale.LOG10),  # log10(E_iso/erg)
        ParamDef("theta_c", 0.01, 0.3, Scale.LINEAR),  # Core angle in radians
        ParamDef("theta_v", 0.0, 0.5, Scale.LINEAR),  # Viewing angle in radians
        ParamDef("n_ism", -3.0, 1.0, Scale.LOG10),  # log10(n/cm^-3)
        ParamDef("p", 2.1, 2.7, Scale.LINEAR),  # Electron energy index
        ParamDef("epsilon_e", -2.5, -0.5, Scale.LOG10),  # log10(epsilon_e)
        ParamDef("epsilon_B", -5.0, -0.5, Scale.LOG10),  # log10(epsilon_B)
    ]
    
    # Create the fitter with default model setup
    fitter = Fitter(data=data, params=params)
    
    # Run MCMC
    samples, log_probs = fitter.run_mcmc(
        n_walkers=32,  # Number of walkers
        n_steps=1000,  # Number of steps per walker
        n_burn=200,    # Number of burn-in steps to discard
        progress=True  # Show progress bar
    )
    
    # Plot the posterior distributions
    fitter.plot_corner()

Parameter Study
^^^^^^^^^^^^

.. code-block:: python

    # Study the effect of electron energy index p
    p_values = np.linspace(2.0, 3.0, 5)
    
    plt.figure(figsize=(10, 6))
    
    # Fix a frequency to study (optical)
    nu_index = 1  # Optical band
    
    for p in p_values:
        # Update the radiation model
        model.radiation.p = p
        
        # Calculate new light curve
        results_p = model.calculate_light_curves(times, frequencies)
        
        # Plot
        plt.loglog(times, results_p[:, nu_index], label=f'p = {p:.1f}')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()
    plt.title('Effect of Electron Energy Index (p) on Optical Light Curves')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

Sample Scripts
------------

The repository includes several example scripts in the `scripts/` directory:

1. **Basic light curve fitting**: ``scripts/fit_lc.py``
2. **MCMC parameter estimation**: ``scripts/mcmc.py``
3. **Structured jet modeling**: ``scripts/structured_jet.py``
4. **Kilonova afterglow**: ``scripts/kilonova.py``

You can run these examples directly:

.. code-block:: bash

    python scripts/structured_jet.py 