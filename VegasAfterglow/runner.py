# src/vegasglow/runner.py
from typing import Optional, Sequence, Tuple

import numpy as np

from .sampler import MultiThreadEmcee
from .types import FitResult, ModelParams, ObsData, ParamDef, Scale, Setups, VegasMC


class Fitter:
    """
    High-level MCMC interface for fitting an afterglow model.
    """

    def __init__(
        self, data: ObsData, config: Setups, num_workers: Optional[int] = None
    ):
        """
        Parameters
        ----------
        data : ObsData
            Observed light curves and spectra.
        config : Setups
            Model configuration (grids, environment, etc).
        """
        self.data = data
        self.config = config
        self.num_workers = num_workers
        # placeholders to be set in fit()
        self._param_defs = None
        self._to_params = None

    def validate_parameters(self, param_defs: Sequence[ParamDef]) -> None:
        """
        Validate that parameter definitions are compatible with the configuration.

        Parameters
        ----------
        param_defs : Sequence[ParamDef]
            Parameter definitions to validate

        Raises
        ------
        ValueError
            If required parameters are missing or incompatible parameters are provided
        """
        param_names = {pd.name for pd in param_defs}
        missing_params = []
        incompatible_params = []

        # Check medium-specific parameters
        if self.config.medium == "wind":
            if "A_star" not in param_names:
                missing_params.append("A_star (required for wind medium)")
            # Wind medium supports n_ism and n0 for stratified configurations
        elif self.config.medium == "ism":
            if "n_ism" not in param_names:
                missing_params.append("n_ism (required for ISM medium)")
            incompatible = {"A_star", "n0", "k_m"} & param_names
            if incompatible:
                incompatible_params.extend(
                    [f"{p} (not used with ISM medium)" for p in incompatible]
                )

        # Check jet-specific parameters
        if self.config.jet == "tophat":
            required = {"theta_c", "E_iso", "Gamma0"}
            missing = required - param_names
            if missing:
                missing_params.extend(
                    [f"{p} (required for tophat jet)" for p in missing]
                )

            incompatible = {"k_e", "k_g", "E_iso_w", "Gamma0_w"} & param_names
            if incompatible:
                incompatible_params.extend(
                    [f"{p} (not used with tophat jet)" for p in incompatible]
                )

        elif self.config.jet == "powerlaw":
            required = {"theta_c", "E_iso", "Gamma0", "k_e", "k_g"}
            missing = required - param_names
            if missing:
                missing_params.extend(
                    [f"{p} (required for powerlaw jet)" for p in missing]
                )

            incompatible = {"E_iso_w", "Gamma0_w"} & param_names
            if incompatible:
                incompatible_params.extend(
                    [f"{p} (not used with powerlaw jet)" for p in incompatible]
                )

        elif self.config.jet == "gaussian":
            required = {"theta_c", "E_iso", "Gamma0"}
            missing = required - param_names
            if missing:
                missing_params.extend(
                    [f"{p} (required for gaussian jet)" for p in missing]
                )

            incompatible = {"k_e", "k_g", "E_iso_w", "Gamma0_w"} & param_names
            if incompatible:
                incompatible_params.extend(
                    [f"{p} (not used with gaussian jet)" for p in incompatible]
                )

        elif self.config.jet == "two_component":
            required = {"theta_c", "E_iso", "Gamma0", "theta_w", "E_iso_w", "Gamma0_w"}
            missing = required - param_names
            if missing:
                missing_params.extend(
                    [f"{p} (required for two_component jet)" for p in missing]
                )

            incompatible = {"k_e", "k_g"} & param_names
            if incompatible:
                incompatible_params.extend(
                    [f"{p} (not used with two_component jet)" for p in incompatible]
                )

        elif self.config.jet == "step_powerlaw":
            required = {
                "theta_c",
                "E_iso",
                "Gamma0",
                "E_iso_w",
                "Gamma0_w",
                "k_e",
                "k_g",
            }
            missing = required - param_names
            if missing:
                missing_params.extend(
                    [f"{p} (required for step_powerlaw jet)" for p in missing]
                )

        elif self.config.jet == "powerlaw_wing":
            required = {"theta_c", "E_iso_w", "Gamma0_w", "k_e", "k_g"}
            missing = required - param_names
            if missing:
                missing_params.extend(
                    [f"{p} (required for powerlaw_wing jet)" for p in missing]
                )

            incompatible = {"E_iso", "Gamma0"} & param_names
            if incompatible:
                incompatible_params.extend(
                    [f"{p} (not used with powerlaw_wing jet)" for p in incompatible]
                )

        # Check required forward shock parameters (always needed)
        fwd_required = {"eps_e", "eps_B", "p"}
        missing_fwd = fwd_required - param_names
        if missing_fwd:
            missing_params.extend(
                [f"{p} (required for forward shock radiation)" for p in missing_fwd]
            )

        # Check reverse shock parameters
        if self.config.rvs_shock:
            rvs_required = {"p_r", "eps_e_r", "eps_B_r", "tau"}
            missing_rvs = rvs_required - param_names
            if missing_rvs:
                missing_params.extend(
                    [
                        f"{p} (required when reverse shock is enabled)"
                        for p in missing_rvs
                    ]
                )
        else:
            rvs_params = {"p_r", "eps_e_r", "eps_B_r", "xi_e_r"} & param_names
            if rvs_params:
                incompatible_params.extend(
                    [f"{p} (reverse shock is disabled)" for p in rvs_params]
                )

        # Check magnetar parameters
        if self.config.magnetar:
            mag_required = {"L0", "t0", "q"}
            missing_mag = mag_required - param_names
            if missing_mag:
                missing_params.extend(
                    [f"{p} (required when magnetar is enabled)" for p in missing_mag]
                )
        else:
            mag_params = {"L0", "t0", "q"} & param_names
            if mag_params:
                incompatible_params.extend(
                    [f"{p} (magnetar is disabled)" for p in mag_params]
                )

        # Report errors
        if missing_params or incompatible_params:
            error_msg = "Parameter validation failed:\n"
            if missing_params:
                error_msg += (
                    "Missing required parameters:\n  - "
                    + "\n  - ".join(missing_params)
                    + "\n"
                )
            if incompatible_params:
                error_msg += (
                    "Incompatible parameters for current configuration:\n  - "
                    + "\n  - ".join(incompatible_params)
                    + "\n"
                )
            error_msg += (
                f"\nCurrent configuration: medium='{self.config.medium}', "
                f"jet='{self.config.jet}', "
            )
            error_msg += (
                f"rvs_shock={self.config.rvs_shock}, magnetar={self.config.magnetar}"
            )
            raise ValueError(error_msg)

    def fit(
        self,
        param_defs: Sequence[ParamDef],
        resolution: Tuple[float, float, float] = (0.3, 1, 10),
        total_steps: int = 10_000,
        burn_frac: float = 0.3,
        thin: int = 1,
        top_k: int = 10,
    ) -> FitResult:
        """
        Run the MCMC sampler.

        Parameters
        ----------
        param_defs :
            A sequence of (name, init, lower, upper) for each free parameter.
        resolution :
            (t_grid, theta_grid, phi_grid) for the coarse MCMC stage.
        total_steps :
            Total number of MCMC steps.
        burn_frac :
            Fraction of steps to discard as burn-in.
        thin :
            Thinning factor for the returned chain.
        top_k :
            Number of top fits to save in the result.

        Returns
        -------
        FitResult
        """

        self.validate_parameters(param_defs)

        defs = list(param_defs)
        self._param_defs = defs

        labels, inits, lowers, uppers = zip(
            *(
                (
                    pd.name,
                    (
                        0.5 * (np.log10(pd.lower) + np.log10(pd.upper))
                        if pd.scale is Scale.LOG
                        else 0.5 * (pd.lower + pd.upper)
                    ),
                    (np.log10(pd.lower) if pd.scale is Scale.LOG else pd.lower),
                    (np.log10(pd.upper) if pd.scale is Scale.LOG else pd.upper),
                )
                for pd in defs
                if pd.scale is not Scale.FIXED
            )
        )
        init = np.array(inits)
        pl = np.array(lowers)
        pu = np.array(uppers)

        # Validate that all parameter names are valid ModelParams attributes
        p_test = ModelParams()
        for pd in defs:
            try:
                getattr(p_test, pd.name)
            except AttributeError:
                raise AttributeError(f"'{pd.name}' is not a valid MCMC parameter")

        # build a fast transformation closure
        def to_params(x: np.ndarray) -> ModelParams:
            p = ModelParams()
            i = 0
            for pd in defs:
                if pd.scale is Scale.FIXED:
                    # fixed param: always pd.init
                    setattr(p, pd.name, 0.5 * (pd.lower + pd.upper))
                else:
                    v = x[i]
                    if pd.scale is Scale.LOG:
                        real = 10**v
                    else:
                        real = v
                    setattr(p, pd.name, real)
                    i += 1
            return p

        self._to_params = to_params

        mcmc = MultiThreadEmcee(
            param_config=(labels, init, pl, pu),
            to_params=to_params,
            model_cls=VegasMC,
            num_workers=self.num_workers,
        )
        result: FitResult = mcmc.run(
            data=self.data,
            base_cfg=self.config,
            resolution=resolution,
            total_steps=total_steps,
            burn_frac=burn_frac,
            thin=thin,
            top_k=top_k,
        )
        return result

    def _with_resolution(self, resolution: Tuple[float, float, float]) -> Setups:
        """
        Clone self.config (without pickle) and override t/theta/phi grids.
        """
        cfg = type(self.config)()
        # copy all public attributes
        for attr in dir(self.config):
            if attr.startswith("_"):
                continue
            if hasattr(cfg, attr):
                try:
                    setattr(cfg, attr, getattr(self.config, attr))
                except Exception:
                    pass
        # override grids
        cfg.phi_resol, cfg.theta_resol, cfg.t_resol = resolution
        return cfg

    def flux_density_grid(
        self,
        best_params: np.ndarray,
        t: np.ndarray,
        nu: np.ndarray,
        resolution: Optional[Tuple[float, float, float]] = (0.3, 1, 10),
    ) -> np.ndarray:
        """
        Compute light curves at the best-fit parameters.

        Parameters
        ----------
        best_params : 1D numpy array
            The vector returned in FitResult.best_params.
        t : 1D numpy array
            Times at which to evaluate.
        nu : 1D numpy array
            Frequencies at which to evaluate.

        Returns
        -------
        array_like
            Shape (n_bands, t.size)
        """
        if self._to_params is None:
            raise RuntimeError("Call .fit(...) before .light_curves()")

        cfg_local = self._with_resolution(resolution)
        p = self._to_params(best_params)

        model = VegasMC(self.data)
        model.set(cfg_local)
        return model.flux_density_grid(p, t, nu)

    def flux(
        self,
        best_params: np.ndarray,
        t: np.ndarray,
        nu_min: float,
        nu_max: float,
        num_points: int,
        resolution: Optional[Tuple[float, float, float]] = (0.3, 1, 10),
    ) -> np.ndarray:
        """
        Compute light curves at the best-fit parameters.

        Parameters
        ----------
        best_params : 1D numpy array
            The vector returned in FitResult.best_params.
        t : 1D numpy array
            Times at which to evaluate.
        nu_min : float
            Minimum frequency at which to evaluate.
        nu_max : float
            Maximum frequency at which to evaluate.
        num_points : int
            Number of points to sample between nu_min and nu_max.

        Returns
        -------
        array_like
            Shape (n_bands, t.size)
        """

        if self._to_params is None:
            raise RuntimeError("Call .fit(...) before .light_curves()")

        cfg_local = self._with_resolution(resolution)
        p = self._to_params(best_params)

        model = VegasMC(self.data)
        model.set(cfg_local)
        return model.flux(p, t, nu_min, nu_max, num_points)
