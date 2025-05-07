from VegasAfterglow import Emission, ModelParams, Setups
import matplotlib.pyplot as plt
import numpy as np

Jy2cgs = 1e-23

t_list = np.geomspace(1e3, 1e8, 100)
nu = 1e16
lc = Emission(nu=nu, t=t_list)

t = 100
nu_list = np.geomspace(1e8, 1e18, 100)
spec = Emission(nu=nu_list, t=t)

config = Setups()
config.medium = "ism"
config.jet = "tophat"

param = ModelParams()
param.E_iso = 1e52
param.theta_c = 0.2
param.theta_v = 0
param.p = 2.2
param.eps_e = 1e-2
param.eps_B = 1e-4
# param.A_star = 1
param.n_ism = 1e-2
param.xi = 1

lc.generate(param, config)
spec.generate(param, config)

fig, ax = plt.subplots(1, 2, figsize=(8, 4))
ax[0].loglog(lc.t, np.array(lc.f_nu) * nu * Jy2cgs)
ax[0].set_xlabel(r'$t$ [s]')
ax[0].set_ylabel(r'$F$ [cgs]')

ax[1].loglog(spec.nu, np.array(spec.f_nu) * nu * Jy2cgs)
ax[1].set_xlabel(r'$nu$ [Hz]')
ax[1].set_ylabel(r'$F$ [cgs]')

plt.tight_layout()
plt.savefig('emission.pdf')
plt.show()