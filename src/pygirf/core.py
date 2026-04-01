from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np
from scipy.io import loadmat, savemat

from .utils import bw_window, time2freq, variable_smoothing


@dataclass
class GirfEssential:
    is_freq_domain_girf: bool | None = None
    in_channels: list[str] = field(default_factory=list)
    out_basis: list[int] = field(default_factory=list)
    girf: np.ndarray | None = None
    freq: np.ndarray | None = None
    girf_time: np.ndarray | None = None
    time: np.ndarray | None = None

    @property
    def self_basis(self) -> np.ndarray:
        mapping = {"Z0": 1, "B0": 1, "X": 2, "Y": 3, "Z": 4}
        return np.array([mapping[ch] for ch in self.in_channels], dtype=int)

    @property
    def df(self) -> float:
        if self.freq is None and self.time is not None:
            self.freq = time2freq(self.time)[0]
        return float(self.freq[1] - self.freq[0])

    @property
    def dt(self) -> float:
        if self.time is None and self.freq is not None:
            self.time = time2freq(self.freq)[0]
        return float(self.time[1] - self.time[0])

    def convert_domain(self, domain: Literal["freq2time", "f2t", "time2freq", "t2f"] | None = None) -> None:
        if domain is None:
            domain = "freq2time" if self.is_freq_domain_girf else "time2freq"

        if domain in ("freq2time", "f2t"):
            self.time = time2freq(np.asarray(self.freq))[0]
            self.girf_time = np.real(np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(self.girf, axes=0), axis=0), axes=0))
        elif domain in ("time2freq", "t2f"):
            self.freq = time2freq(np.asarray(self.time))[0]
            self.girf = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(self.girf_time, axes=0), axis=0), axes=0)

    def save(self, filename: str, overwrite: bool = False) -> None:
        p = Path(filename)
        if p.exists() and not overwrite:
            return
        savemat(
            filename,
            {
                "isFreqDomainGirf": self.is_freq_domain_girf,
                "inChannels": np.array(self.in_channels, dtype=object),
                "outBasis": np.asarray(self.out_basis),
                "girf": self.girf,
                "freq": self.freq,
                "girfTime": self.girf_time,
                "time": self.time,
            },
        )

    def load(self, filename: str) -> None:
        mat = loadmat(filename, squeeze_me=True, struct_as_record=False)
        self.is_freq_domain_girf = bool(mat.get("isFreqDomainGirf", self.is_freq_domain_girf))
        if "inChannels" in mat:
            self.in_channels = [str(x) for x in np.atleast_1d(mat["inChannels"]).tolist()]
        if "outBasis" in mat:
            self.out_basis = np.atleast_1d(mat["outBasis"]).astype(int).tolist()
        self.girf = mat.get("girf", self.girf)
        self.freq = mat.get("freq", self.freq)
        self.girf_time = mat.get("girfTime", self.girf_time)
        self.time = mat.get("time", self.time)


@dataclass
class GirfProvider(GirfEssential):
    time_in: np.ndarray | None = None
    in_data: np.ndarray | None = None
    time_out: np.ndarray | None = None
    out_data: np.ndarray | None = None
    resample_method: str | None = None
    filter_info: list[str] = field(default_factory=list)

    @property
    def in_freq(self) -> np.ndarray:
        return np.fft.fftshift(np.fft.fft(self.in_data, axis=0), axes=0)

    @property
    def out_freq(self) -> np.ndarray:
        return np.fft.fftshift(np.fft.fft(self.out_data, axis=0), axes=0)

    @property
    def dt_in(self) -> float:
        return float(self.time_in[1] - self.time_in[0])

    @property
    def dt_out(self) -> float:
        return float(self.time_out[1] - self.time_out[0])

    def compute_girf(self):
        if self.time_in.shape[0] != self.time_out.shape[0]:
            raise ValueError("Resample input and output waveforms onto same time grid before GIRF calculation")

        nS = len(self.time_out)
        nIn = self.in_data.shape[1]
        nOut = self.out_data.shape[1]
        nWave = self.in_data.shape[2]

        freq = time2freq(self.time_out[:, 0] if self.time_out.ndim > 1 else self.time_out)[0]
        in_freq = self.in_freq
        out_freq = self.out_freq

        use_matrix = True
        if nIn > 1:
            in_waveforms = np.zeros((nIn, nWave), dtype=int)
            for ch in range(nIn):
                for w in range(nWave):
                    if np.count_nonzero(in_freq[:, ch, w]) > 4:
                        in_waveforms[ch, w] = 1
            use_matrix = np.max(np.sum(in_waveforms, axis=0)) > 1

        girf = np.zeros((nS, nOut, nIn), dtype=complex)
        if nIn == 1:
            inv = 1.0 / np.sum(np.abs(in_freq) ** 2, axis=2)
            for i_out in range(nOut):
                girf[:, i_out, 0] = np.sum(out_freq[:, i_out, :] * np.conj(in_freq[:, 0, :]), axis=1) * inv[:, 0]
        elif not use_matrix:
            for i_in in range(nIn):
                idx = np.where(in_waveforms[i_in, :] == 1)[0]
                inv = 1.0 / np.sum(np.abs(in_freq[:, i_in, idx]) ** 2, axis=1)
                for i_out in range(nOut):
                    girf[:, i_out, i_in] = np.sum(out_freq[:, i_out, idx] * np.conj(in_freq[:, i_in, idx]), axis=1) * inv
        else:
            for i_f in range(nS):
                girf[i_f, :, :] = out_freq[i_f, :, :] @ np.linalg.pinv(in_freq[i_f, :, :])

        girf[np.isnan(girf)] = 0
        self.is_freq_domain_girf = True
        self.girf = girf
        self.freq = freq
        self.convert_domain("freq2time")
        return girf, freq

    def window_freq(self, BW: float, filter_type: str = "rc", beta: float = 1.0 / 3.0):
        if self.girf is None and self.girf_time is not None:
            self.convert_domain("time2freq")
        girf, filt = bw_window(self.girf, self.freq, BW, filter_type, beta)
        freq = self.freq.copy()
        nF = len(filt)
        zero_left = np.where(filt[: nF // 2] == 0)[0]
        if zero_left.size > 0:
            i0 = zero_left[-1]
            girf = girf[i0 : nF - i0, ...]
            freq = freq[i0 : nF - i0]
        self.girf = girf
        self.freq = freq
        self.convert_domain("freq2time")
        return girf, freq

    def window_time(self, T: float, filter_type: str = "rc", beta: float = 0.0):
        if self.girf_time is None and self.girf is not None:
            self.convert_domain("freq2time")
        girf_time, filt = bw_window(self.girf_time, self.time, T, filter_type, beta)
        t = self.time.copy()
        nS = len(filt)
        zero_left = np.where(filt[: nS // 2] == 0)[0]
        if zero_left.size > 0:
            i0 = zero_left[-1]
            girf_time = girf_time[i0 + 1 : nS - i0, ...]
            t = t[i0 + 1 : nS - i0]
        self.girf_time = girf_time
        self.time = t
        self.convert_domain("time2freq")
        return girf_time, t

    def var_smooth_freq(self, f_max: float, BW_range: tuple[float, float] | None = None, vsBW: float | None = None):
        if self.girf is None and self.girf_time is not None:
            self.convert_domain("time2freq")
        girf, freq, *_ = variable_smoothing(self.girf, self.freq, f_max, BW_range, vsBW)
        self.girf = girf
        self.freq = freq
        self.convert_domain("freq2time")
        return girf, freq

    def peak_elimination(self, f_center: float, f_width: float, f_reps: int = 1, n_fit: int = 1):
        if self.girf is None and self.girf_time is not None:
            self.convert_domain("time2freq")

        girf = self.girf.copy()
        f = self.freq
        if abs(f_center) > np.max(np.abs(f)):
            return

        nS = len(f)
        nC = int(np.floor(nS / 2.0)) + 1
        for i_r in range(1, f_reps + 1):
            inds_peak = np.where(np.abs(f - i_r * f_center) < f_width / 2.0)[0]
            inds_peak_neg = np.arange(-len(inds_peak), 0) + 1 + 2 * nC - inds_peak[0] - 1
            inds_fit = np.where(np.abs(f - i_r * f_center) < 3.0 * f_width / 2.0)[0]
            inds_fit = np.setdiff1d(inds_fit, inds_peak)
            for i_k in range(girf.shape[1]):
                pre = np.polyfit(f[inds_fit], np.real(girf[inds_fit, i_k]), n_fit)
                pim = np.polyfit(f[inds_fit], np.imag(girf[inds_fit, i_k]), n_fit)
                girf[inds_peak, i_k] = np.polyval(pre, f[inds_peak]) + 1j * np.polyval(pim, f[inds_peak])
            girf[inds_peak_neg, :] = np.flip(np.conj(girf[inds_peak, :]), axis=0)

        self.girf = girf
        self.convert_domain("freq2time")

    def vis(
        self,
        plot_type: str = "GIRF",
        domain: str = "f",
        in_channel: list[str] | None = None,
        out_basis: list[int] | None = None,
        waveforms: list[int] | None = None,
        integration: bool = False,
    ):
        import matplotlib.pyplot as plt

        if waveforms is None:
            waveforms = list(range(self.in_data.shape[2]))
        if out_basis is None:
            out_basis = self.self_basis.tolist()
        if in_channel is None:
            in_channel = self.in_channels

        figs = []
        for i_in, ch in enumerate(in_channel):
            fig = plt.figure()
            fig.suptitle(f"{plot_type} {ch}, {domain}")
            ptype = plot_type.lower()
            if ptype == "girf":
                if domain == "f":
                    if self.girf is None or self.freq is None:
                        self.convert_domain("time2freq")
                    g = self.girf[:, out_basis[i_in] - 1, i_in]
                    ax1 = fig.add_subplot(2, 1, 1)
                    ax1.plot(self.freq, np.abs(g))
                    ax2 = fig.add_subplot(2, 1, 2)
                    ax2.plot(self.freq, np.unwrap(np.angle(g)) - np.unwrap(np.angle(g))[len(g) // 2])
                    ax2.set_xlim(ax1.get_xlim())
                else:
                    if self.girf_time is None or self.time is None:
                        self.convert_domain("freq2time")
                    t = self.time
                    g = self.girf_time[:, out_basis[i_in] - 1, i_in]
                    if integration:
                        t = t + self.dt / 2.0
                        g = np.cumsum(g) * self.dt
                    ax = fig.add_subplot(1, 1, 1)
                    ax.plot(t, g)
            elif ptype == "in":
                if domain == "f":
                    f = time2freq(self.time_in)[0]
                    dat = self.in_freq[:, i_in, waveforms]
                    ax1 = fig.add_subplot(2, 1, 1)
                    ax1.plot(f, np.abs(dat))
                    ax2 = fig.add_subplot(2, 1, 2)
                    ph = np.unwrap(np.angle(dat), axis=0)
                    ph = ph - ph[len(ph) // 2 : len(ph) // 2 + 1, ...]
                    ax2.plot(f, ph)
                else:
                    t = self.time_in
                    dat = self.in_data[:, i_in, waveforms]
                    if integration:
                        t = t + self.dt_in / 2.0
                        dat = np.cumsum(dat, axis=0) * self.dt_in
                    ax = fig.add_subplot(1, 1, 1)
                    ax.plot(t, dat)
            elif ptype == "out":
                if domain == "f":
                    f = time2freq(self.time_out)[0]
                    dat = self.out_freq[:, out_basis[i_in] - 1, waveforms]
                    ax1 = fig.add_subplot(2, 1, 1)
                    ax1.plot(f, np.abs(dat))
                    ax2 = fig.add_subplot(2, 1, 2)
                    ph = np.unwrap(np.angle(dat), axis=0)
                    ph = ph - ph[len(ph) // 2 : len(ph) // 2 + 1, ...]
                    ax2.plot(f, ph)
                else:
                    t = self.time_out
                    dat = self.out_data[:, out_basis[i_in] - 1, waveforms]
                    if integration:
                        t = t + self.dt_out / 2.0
                        dat = np.cumsum(dat, axis=0) * self.dt_out
                    ax = fig.add_subplot(1, 1, 1)
                    ax.plot(t, dat)
            elif ptype == "inout":
                if domain == "f":
                    f_in = time2freq(self.time_in)[0]
                    f_out = time2freq(self.time_out)[0]
                    dat_in = self.in_freq[:, i_in, waveforms]
                    dat_out = self.out_freq[:, out_basis[i_in] - 1, waveforms]
                    ax1 = fig.add_subplot(2, 1, 1)
                    ax1.plot(f_in, np.abs(dat_in), "k")
                    ax1.plot(f_out, np.abs(dat_out))
                    ax2 = fig.add_subplot(2, 1, 2)
                    ph_in = np.unwrap(np.angle(dat_in), axis=0)
                    ph_in = ph_in - ph_in[len(ph_in) // 2 : len(ph_in) // 2 + 1, ...]
                    ph_out = np.unwrap(np.angle(dat_out), axis=0)
                    ph_out = ph_out - ph_out[len(ph_out) // 2 : len(ph_out) // 2 + 1, ...]
                    ax2.plot(f_in, ph_in, "k")
                    ax2.plot(f_out, ph_out)
                else:
                    t_in = self.time_in
                    t_out = self.time_out
                    dat_in = self.in_data[:, i_in, waveforms]
                    dat_out = self.out_data[:, out_basis[i_in] - 1, waveforms]
                    if integration:
                        t_in = t_in + self.dt_in / 2.0
                        t_out = t_out + self.dt_out / 2.0
                        dat_in = np.cumsum(dat_in, axis=0) * self.dt_in
                        dat_out = np.cumsum(dat_out, axis=0) * self.dt_out
                    ax = fig.add_subplot(1, 1, 1)
                    ax.plot(t_in, dat_in, "k")
                    ax.plot(t_out, dat_out)
            elif ptype == "sensitivity":
                f = time2freq(self.time_in)[0]
                in_f = self.in_freq[:, i_in, waveforms]
                in_rss = np.sqrt(np.sum(np.abs(in_f) ** 2, axis=1))
                ax = fig.add_subplot(1, 1, 1)
                ax.plot(f, np.abs(in_f))
                ax.plot(f, np.abs(in_rss), "k")
            else:
                raise ValueError(f"Unsupported plot_type: {plot_type}")
            figs.append(fig)
        return figs


@dataclass
class GirfApplier(GirfEssential):
    gamma1H: float = 2.675222099e8

    def predict_grad(self, t_in: np.ndarray, g_in: np.ndarray, t_out: np.ndarray | None = None, pred_channels: list[str] | None = None, conv_type: str = "conv"):
        if pred_channels is None:
            pred_channels = ["X", "Y", "Z"]
        if t_out is None:
            t_out = t_in

        n_ch = len(pred_channels)
        n_in = g_in.shape[1]
        if n_in != n_ch:
            raise ValueError("Number of inputs must match number of assigned channels")

        n_out = self.girf.shape[1]
        ns_h = len(self.freq)
        n_il = g_in.shape[2] if g_in.ndim == 3 else 1
        if g_in.ndim == 2:
            g_in = g_in[:, :, None]

        H = np.zeros((ns_h, n_in, n_out), dtype=complex)
        for i, ch in enumerate(pred_channels):
            if ch not in self.in_channels:
                raise ValueError(f"Input channel {ch} cannot be found")
            ch_idx = self.in_channels.index(ch)
            H[:, i, :] = self.girf[:, :, ch_idx]

        f_h = self.freq
        df_h = f_h[1] - f_h[0]

        TH = 1.0 / df_h
        TO = max(t_in[-1], t_out[-1]) - min(t_in[0], t_out[0])
        TZF = min(TH, TO)

        dt_in = t_in[1] - t_in[0]
        n_zf = int(np.ceil(TZF / dt_in))
        g_in = np.concatenate([np.zeros((n_zf, n_in, n_il)), g_in, np.zeros((n_zf, n_in, n_il))], axis=0)
        t_in = np.concatenate([t_in[0] + dt_in * np.arange(-n_zf, 0), t_in, t_in[-1] + dt_in * np.arange(1, n_zf + 1)])
        f_in = time2freq(t_in)[0]
        IN = np.fft.fftshift(np.fft.fft(g_in, axis=0), axes=0) * dt_in
        ns_in = len(t_in)

        HIp = np.zeros((len(f_in), n_in, n_out), dtype=complex)
        for i in range(n_in):
            for j in range(n_out):
                HIp[:, i, j] = np.interp(f_in, f_h, H[:, i, j], left=0, right=0)

        if conv_type == "conv":
            h_time = np.real(np.fft.ifft(np.fft.ifftshift(HIp, axes=0), axis=0))
            t_start = 1e-3
            t_settle = TZF
            t_end_ind = int(np.ceil(t_settle / dt_in))
            t_start_ind = int(np.ceil(t_start / dt_in))
            h_time[t_end_ind : ns_in - t_start_ind, :, :] = 0
            HIp = np.fft.fftshift(np.fft.fft(h_time, axis=0), axes=0)
            HIp[np.isnan(HIp)] = 0

        OUT = np.zeros((ns_in, n_out, n_il), dtype=complex)
        for i_il in range(n_il):
            OUT[:, :, i_il] = np.sum(IN[:, :, i_il][:, :, None] * HIp, axis=1)

        g_out = np.real(np.fft.ifft(np.fft.ifftshift(OUT, axes=0), axis=0) / dt_in)
        k_out = np.cumsum(g_out, axis=0) * dt_in * self.gamma1H

        dt_out = t_out[1] - t_out[0]
        g_out_i = np.zeros((len(t_out), n_out, n_il), dtype=float)
        k_out_i = np.zeros((len(t_out), n_out, n_il), dtype=float)
        for i in range(n_out):
            for j in range(n_il):
                g_out_i[:, i, j] = np.interp(t_out, t_in, g_out[:, i, j], left=0, right=0)
                k_out_i[:, i, j] = np.interp(t_out + dt_out / 2.0, t_in + dt_in / 2.0, k_out[:, i, j], left=0, right=0)

        t_k = t_out + dt_out / 2.0
        if n_il == 1:
            g_out_i = g_out_i[:, :, 0]
            k_out_i = k_out_i[:, :, 0]
        return g_out_i, k_out_i, t_k
