from pathlib import Path
import sys

import numpy as np

# Allow running demo directly from repository root without installation.
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from girf import GirfApplier, GirfProvider, bw_filter, sweeps, trapezoid


def main() -> None:
    dt_in = 10e-6
    t_in = np.arange(0, 500e-3 + dt_in, dt_in)

    in_sweep, _, _ = sweeps(t_in, 50e-3, 0, 40e3, 0, 1, 5e-3, "linear")
    in_trapz = trapezoid(t_in, 5e-3, 10, 0.1e-3, 0)

    in_data = np.zeros((len(t_in), 1, 2))
    in_data[:, 0, 0] = in_sweep[:, 0]
    in_data[:, 0, 1] = in_trapz[:, 0]

    dt_out = dt_in
    t_out = t_in.copy()
    out_basis = list(range(1, 17))
    out = np.zeros((len(t_out), len(out_basis), in_data.shape[2]))

    self_out = bw_filter(in_data, t_in, 20e3, "rc", 0.6)
    cross_out = 0.1 * bw_filter(in_data, t_in, 10e3, "rc", 1.0)
    out[:, 3, :] = self_out[:, 0, :]
    out[:, 0, :] = cross_out[:, 0, :]

    out_k = dt_out * np.cumsum(out, axis=0)
    out_k = out_k + 1e-7 * np.random.randn(*out_k.shape)
    out = np.concatenate([np.zeros((1, len(out_basis), in_data.shape[2])), np.diff(out_k, axis=0) / dt_out], axis=0)

    provider = GirfProvider(
        in_channels=["Z"],
        out_basis=out_basis,
        time_in=t_in,
        in_data=in_data,
        time_out=t_out,
        out_data=out,
    )
    provider.compute_girf()
    provider.window_freq(60e3, "rc", 1 / 4)
    provider.var_smooth_freq(30e3)
    provider.window_time(10e-3, "rc", 0)

    essential = provider
    applier = GirfApplier(
        is_freq_domain_girf=essential.is_freq_domain_girf,
        in_channels=essential.in_channels,
        out_basis=essential.out_basis,
        girf=essential.girf,
        freq=essential.freq,
        girf_time=essential.girf_time,
        time=essential.time,
    )
    g_out, _, _ = applier.predict_grad(t_in, in_data[:, :, 0], t_out, ["Z"], "conv")
    print("Prediction complete", g_out.shape)


if __name__ == "__main__":
    main()
