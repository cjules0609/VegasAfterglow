Examples
========

.. contents:: Table of Contents
   :local:
   :depth: 2

Basic Usage
-----------

Setting up a simple afterglow model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np
    from VegasAfterglow import ISM, TophatJet, Observer, Radiation, Model
    
    # Define the circumburst environment (constant density ISM)
    medium = ISM(n_ism=1)

    # Configure the jet structure (top-hat with opening angle, energy, and Lorentz factor)
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)

    # Set observer parameters (distance, redshift, viewing angle)
    obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0)

    # Define radiation microphysics parameters
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)

    # Combine all components into a complete afterglow model
    model = Model(jet=jet, medium=medium, observer=obs, forward_rad=rad)

    # Define time range for light curve calculation
    times = np.logspace(2, 8, 200)  

    # Define observing frequencies (radio, optical, X-ray bands in Hz)
    bands = np.array([1e9, 1e14, 1e17])  

    # Calculate the afterglow emission at each time and frequency
    results = model.specific_flux(times, bands)

    # Visualize the multi-wavelength light curves
    plt.figure(figsize=(4.8, 3.6),dpi=200)

    # Plot each frequency band 
    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
    plt.loglog(times, results['syn'][i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz')

    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()

    # Define broad frequency range (10⁵ to 10²² Hz) 
    frequencies = np.logspace(5, 22, 200)  

    # Select specific time epochs for spectral snapshots 
    epochs = np.array([1e2, 1e3, 1e4, 1e5 ,1e6, 1e7, 1e8])

    # Calculate spectra at each epoch
    results = model.spectra(frequencies, epochs)

    # Plot broadband spectra at each epoch
    plt.figure(figsize=(4.8, 3.6),dpi=200)
    colors = plt.cm.viridis(np.linspace(0,1,len(epochs)))

    for i, t in enumerate(epochs):
        exp = int(np.floor(np.log10(t)))
        base = t / 10**exp
        plt.loglog(frequencies, results['syn'][i,:], color=colors[i], label=fr'${base:.1f} \times 10^{{{exp}}}$ s')

    # Add vertical lines marking the bands from the light curve plot
    for i, band in enumerate(bands):
        exp = int(np.floor(np.log10(band)))
        base = band / 10**exp
        plt.axvline(band,ls='--',color='C'+str(i))

    plt.xlabel('frequency (Hz)')
    plt.ylabel('flux density (erg/cm²/s/Hz)')
    plt.legend(ncol=2)
    plt.title('Synchrotron Spectra')


Ambient Media Models
--------------------

Wind Medium
^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Wind

    # Create a stellar wind medium
    wind = Wind(A_star=0.1)  # A* parameter

    #..other settings
    model = Model(medium=wind, ...)

User-Defined Medium
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Medium

    # Define a custom density profile function
    def custom_density(phi, theta, r):
        #return what ever density profile you want as a function of phi, theta, and r
        

    def custom_mass(phi, theta, r):
        #return the integral of the density profile over r.
        #you may keep the consistency of the mass profile with the density profile
        #the purpose of providing the extra mass profile is to reduce the extra computations.
    
    # Create a user-defined medium
    medium = Medium(rho=custom_density, mass=custom_mass)
    
    #..other settings
    model = Model(medium=medium, ...)


Structured Jet Models
---------------------

Gaussian Jet
^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import GaussianJet

    # Create a structured jet with Gaussian energy profile
    jet = GaussianJet(
        theta_c=0.05,         # Core angular size (radians)
        E_iso=1e53,           # Isotropic-equivalent energy (ergs)
        Gamma0=300            # Initial Lorentz factor
    )

    #..other settings
    model = Model(jet=jet, ...)

Power-Law Jet
^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import PowerLawJet

    # Create a power-law structured jet
    jet = PowerLawJet(
        theta_c=0.05,         # Core angular size (radians)
        E_iso=1e53,           # Isotropic-equivalent energy (ergs)
        Gamma0=300,           # Initial Lorentz factor
        k=2.0                 # Power-law index
    )

    #..other settings
    model = Model(jet=jet, ...)

Jet with Spreading
^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import TophatJet

    # Create a power-law structured jet
    jet = TophatJet(
        theta_c=0.05,        
        E_iso=1e53,          
        Gamma0=300,         
        spreading=True       # Enable spreading
    )

    #..other settings
    model = Model(jet=jet, ...)

.. note::
    The jet spreading (Lateral Expansion) is experimental and only works for the top-hat jet, Gaussian jet, and power-law jet with a jet core.
    The spreading prescription may not work for arbitrary user-defined jet structures.

User-Defined Jet
^^^^^^^^^^^^^^^^

You may also define your own jet structure by providing the energy and lorentz factor profile.
Those two profiles are required to complete a jet structure. You may also provide the magnetization profile, enregy injection profile, and mass injection profile.
Those profiles are optional and will be set to zero function if not provided.

.. code-block:: python

    from VegasAfterglow import Ejecta

    # Define a custom energy profile function, required to complete the jet structure
    def energy(phi, theta):
        #return what ever energy PER SOLID ANGLE profile you want as a function of phi and theta

    # Define a custom lorentz factor profile function, required to complete the jet structure
    def lorentz(phi, theta):
        #return what ever lorentz factor profile you want as a function of phi and theta
    
    # Define a custom magnetization profile function, optional
    def magnetization(phi, theta):
        #return what ever magnetization profile you want as a function of phi and theta

    # Define a custom energy injection profile function, optional
    def energy_injection(phi, theta, t):
        #return what ever energy injection PER SOLID ANGLE profile you want as a function of phi, theta, and time

    # Define a custom mass injection profile function, optional
    def mass_injection(phi, theta, t):
        #returnwhat ever mass injection PER SOLID ANGLE profile you want as a function of phi, theta, and time

    # Create a user-defined jet
    jet = Ejecta(energy=energy, lorentz=lorentz, magnetization=magnetization, energy_injection=energy_injection, mass_injection=mass_injection)

    #..other settings

    #if your jet is not axisymmetric, set axisymmetric to False
    model = Model(jet=jet, ..., axisymmetric=False, resolutions=(0.3, 2., 5.))

    # the user-defined jet structure could be spiky, if the default resolution may not resolve the jet structure, if that is the case, you can try a finer resolution (phi_ppd, theta_ppd, t_ppd)
    # where phi_ppd is the number of points per degree in the phi direction, theta_ppd is the number of points per degree in the theta direction, and t_ppd is the number of points per decade in the time direction    .
    
.. note::
    Setting usere-defined structured jet in the Python level is OK for light curve and spectrum calculation. However, it is not recommended for MCMC parameter fitting.
    The reason is that setting user-defined profiles in the Python level leads to a large overhead due to the Python-C++ inter-process communication.
    Users are recommended to set up the user-defined jet structure in the C++ level for MCMC parameter fitting for better performance.
        

Radiation Processes
-------------------

Synchrotron Self-Compton
^^^^^^^^^^^^^^^^^^^^^^^^    

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
-----------------

Reverse Shock
^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^

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
--------------

The repository includes several example scripts in the ``script`` directory:

1. **MCMC parameter estimation**: ``script/mcmc.py``

You can run these examples directly:

.. code-block:: bash

    python script/mcmc.py 