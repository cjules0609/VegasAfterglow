Numerical Electron Distributions
================================

VegasAfterglow can calculate optically thin or self-absorbed synchrotron emission from a bounded,
instantaneous electron column distribution. The supplied distribution is the
actual emitting population :math:`d\Sigma_e/d\gamma`; it is not an injection
function and is not automatically cooled or broken at :math:`\gamma_c`.

The built-in flat distribution can be configured as follows:

.. code-block:: python

   from VegasAfterglow import FlatElectrons, Radiation

   electrons = FlatElectrons(
       gamma_min=1e2,
       gamma_max=1e7,
       normalization="energy",
   )

   radiation = Radiation(
       eps_e=0.1,
       eps_B=1e-3,
       p=2.3,
       electrons=electrons,
   )

``PowerLawElectrons`` and ``CutoffPowerLawElectrons`` provide bounded power-law
and exponentially cut off power-law shapes. Every numerical distribution must
explicitly select ``normalization="energy"`` or ``normalization="number"``.
Energy normalization enforces ``eps_e`` and reports the implied electron number
fraction; number normalization enforces ``xi_e`` and reports the implied electron
energy fraction.

Numerical electron distributions currently have the following limitations:

* emission remains optically thin by default for backward compatibility; pass ``ssa=True`` to ``Radiation``
  to apply general numerical synchrotron self-absorption;
* ``ssc=True`` and ``kn=True`` are rejected because the existing SSC and IC
  cooling algorithms assume the standard afterglow electron model;
* the distribution need not contain :math:`\gamma_m`, and :math:`\gamma_c` does
  not modify it;
* the legacy radiative blast-wave efficiency still uses ``p``, ``eps_e``, and
  ``xi_e``. Numerical electrons change the instantaneous radiation calculation,
  not that dynamics approximation.

The per-cell population and spectrum can be inspected through ``details()``:

.. code-block:: python

   details = model.details(t_min=1e2, t_max=1e7)
   dSigma_dgamma = details.fwd.electron_column_distribution[0, 0, 5](gamma)
   I_nu = details.fwd.sync_spectrum[0, 0, 5](nu_comoving)

Numerical self-absorption is enabled explicitly:

.. code-block:: python

   radiation = Radiation(
       eps_e=0.1,
       eps_B=1e-3,
       p=2.3,
       electrons=electrons,
       ssa=True,
   )

The backend computes optical depth directly from the electron column distribution and applies

.. math::

   I_\nu = I_{\nu,\mathrm{thin}}\frac{1-e^{-\tau_\nu}}{\tau_\nu}.

The corresponding per-cell diagnostic is
``details.fwd.sync_optical_depth[i, j, k](nu_comoving)``. It is dimensionless and uses comoving frequency in Hz.
No shell width or local absorption coefficient is defined. Numerical distributions do not fabricate a
``gamma_a`` or ``nu_a``; those legacy standard-model diagnostics remain NaN.

``electron_gamma_min``, ``electron_gamma_max``,
``electron_number_fraction``, and ``electron_energy_fraction`` contain the
corresponding per-cell diagnostics.

User-Defined Shapes
-------------------

``ElectronDistribution`` samples an arbitrary Python callable once during
construction. The callable receives the complete logarithmic gamma grid as a
one-dimensional NumPy array and must return a real, finite, non-negative array
of the same shape. After construction, VegasAfterglow retains only native
sampled arrays; the callable is not used by light curves, spectra, details, sky
images, observer integration, or cell loops.

.. code-block:: python

   import numpy as np

   from VegasAfterglow import ElectronDistribution, Radiation

   def thermal_tail(gamma):
       thermal = gamma**2 * np.exp(-gamma / 3e2)
       tail = 1e-3 * gamma**-2.4
       return thermal + tail

   electrons = ElectronDistribution(
       function=thermal_tail,
       gamma_min=1,
       gamma_max=1e8,
       normalization="energy",
   )

   radiation = Radiation(
       eps_e=0.1,
       eps_B=1e-3,
       p=2.3,
       electrons=electrons,
   )

The function specifies an unnormalized shape :math:`f(\gamma)`. The physical
population remains

.. math::

   \frac{d\Sigma_e}{d\gamma} = A_{\rm cell} f(\gamma),

with ``normalization`` selecting the same energy or number constraint used by
the built-in numerical distributions. It is an instantaneous emitting
population, not an injection function, probability density, or kinetic
solution. VegasAfterglow does not evolve or automatically cool this shape.

The optional ``samples_per_decade`` argument defaults to 32, matching the
convergence-tested built-in grid. Values from 4 through 512 are accepted; more
samples resolve narrow features more accurately at higher construction cost.
``electrons.gamma`` and ``electrons.shape`` return copies of the native samples
for inspection. Scalar returns are rejected rather than broadcast.
