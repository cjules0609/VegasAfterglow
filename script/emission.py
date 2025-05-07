from VegasAfterglow import Emission, ModelParams, Setups
import matplotlib.pyplot as plt
import numpy as np

Jy2cgs = 1e-23

config = Setups()
config.medium = "ism"
config.jet = "tophat"
config.spreading = False

param = ModelParams()
param.E_iso = 1e52
param.theta_c = 0.2
param.theta_v = 0
param.p = 2.2
param.eps_e = 1e-2
param.eps_B = 1e-4
param.n_ism = 1e-2
param.xi = 1

# Initialize the Emission object
emission = Emission(param, config)

# Generate the light curve (lc)
t_list = np.geomspace(1e3, 1e8, 100)
nu = 1e16
lc = emission.generate_lc(nu=nu, t=t_list)

# Generate the spectrum (spec)
t = 100
nu_list = np.geomspace(1e8, 1e18, 100)
spec = emission.generate_spec(nu=nu_list, t=t)

fig, ax = plt.subplots(1, 2, figsize=(8, 4))
ax[0].loglog(t_list, np.array(lc) * nu * Jy2cgs)
ax[0].set_xlabel(r'$t$ [s]')
ax[0].set_ylabel(r'$F$ [cgs]')

ax[1].loglog(nu_list, np.array(spec) * nu * Jy2cgs)
ax[1].set_xlabel(r'$nu$ [Hz]')
ax[1].set_ylabel(r'$F$ [cgs]')

plt.tight_layout()
plt.savefig('emission.pdf')
plt.show()
