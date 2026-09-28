Time-Dependent Electron Cooling
================================

``ElectronCooling`` evolves an arbitrary isotropic electron distribution in
the comoving frame.  It solves

.. math::

   \frac{\partial N}{\partial t'}+
   \frac{\partial}{\partial\gamma}(\dot\gamma N)
   = Q-\frac{N}{t_{\rm esc}}

on a logarithmic Lorentz-factor grid.  The finite-volume state stores
bin-integrated particle numbers, uses a positivity-preserving implicit upwind
step, and does not renormalize the electron energy after cooling.

The current implementation is deliberately a **one-zone solver**.  Times are
comoving seconds, not observer times.  VegasAfterglow's afterglow dynamics
stores an effective emitting shell rather than Lagrangian fluid elements and
their shock-crossing histories, so this release does not present one-zone
evolution as an exact equal-arrival-time afterglow cooling calculation.

Basic use
---------

.. code-block:: python

   import numpy as np
   from VegasAfterglow import ElectronCooling, ElectronDistribution

   injected = ElectronDistribution(
       function=lambda gamma: gamma**-2.3,
       gamma_min=1,
       gamma_max=1e8,
       normalization="number",
       samples_per_decade=64,
   )

   cooling = ElectronCooling(synchrotron=True)
   history = cooling.evolve(
       injected,
       times=np.geomspace(1, 1e6, 100),       # comoving seconds
       magnetic_field=np.geomspace(10, 0.1, 100),  # gauss
   )

``history.distribution[k]`` is :math:`dN/d\gamma` at the corresponding
comoving time.  ``history.number`` and ``history.kinetic_energy`` are numerical
conservation diagnostics.  The energy is the dimensionless moment
:math:`\int(\gamma-1)N\,d\gamma`; multiply by :math:`m_ec^2` when the input
normalization is a physical particle number.

Loss channels and injection
---------------------------

Synchrotron losses use

.. math::

   \dot\gamma_{\rm syn}=-\frac{\sigma_T B'^2}{6\pi m_ec}(\gamma^2-1).

Thomson inverse-Compton cooling can be specified through ``compton_y`` or a
``photon_energy_density`` history in erg cm\ :sup:`-3`. If both are supplied,
their loss rates are added. These are
energy-independent Thomson approximations; they are not a Klein--Nishina SSC
solution.  Adiabatic cooling accepts the comoving expansion scalar
``expansion_rate=dln(V')/dt'`` in s\ :sup:`-1`.

``injection`` may be a gamma vector, an interval-by-gamma array, or a callable
``injection(gamma, midpoint_time)`` returning
:math:`Q=dN/(d\gamma\,dt')` per second. ``escape_time`` is an optional comoving
time in seconds.  Every loss channel is independently switchable.

Boundary and radiation semantics
--------------------------------

By default particles reaching the lower grid edge accumulate in the first bin,
which conserves number. Set ``accumulate_at_gamma_min=False`` to use an open
lower boundary; ``escaped_lower`` then records the removed population. Cooling
never transports particles through the upper boundary.

``history.final_electrons()`` converts the last shape back to the existing
instantaneous numerical-radiation API and defaults to number normalization.
That is physically appropriate for a closed population because synchrotron and
adiabatic losses conserve particle number while reducing its energy.  Energy
normalization would restore the radiated energy and should only be requested
when that is intentionally the desired instantaneous constraint.  If injection
or escape changes the absolute particle fraction, the present ``Radiation``
interface still applies its configured ``xi_e``; use the history arrays for
absolute normalization until hydrodynamically coupled shell histories are
implemented.
