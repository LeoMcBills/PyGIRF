from __future__ import annotations

from typing import Any

import numpy as np
from scipy.special import erf
from scipy.signal.windows import blackman, blackmanharris


def time2freq(t: np.ndarray) -> tuple[np.ndarray, float, float]:
    t = np.asarray(t).reshape(-1)
    nrs = t.size
    dt = t[1] - t[0]
    f_max = 1.0 / dt
    df = f_max / nrs
    f = (np.arange(nrs) - np.floor(nrs / 2.0)) * df
    return f, df, f_max


def centered_time(dt: float, nrs: int) -> np.ndarray:
    t = np.arange(0, dt * nrs, dt)
    cs = int(np.floor(nrs / 2.0))
    return t - t[cs]


def centered_phase(array: np.ndarray) -> np.ndarray:
    a = np.asarray(array)
    phase = np.unwrap(np.angle(a), axis=0)
    middle = int(np.floor(a.shape[0] / 2.0))
    return phase - phase[middle : middle + 1, ...]


def raised_cosine(f: np.ndarray, T: float, beta: float) -> np.ndarray:
    f = np.asarray(f).reshape(-1)
    filt = np.zeros_like(f, dtype=float)
    ind = np.abs(f) <= (1.0 + beta) / (2.0 * T)
    if beta != 0:
        filt[ind] = (1.0 + np.cos(np.pi * T * (np.abs(f[ind]) - (1.0 - beta) / (2.0 * T)) / beta)) / 2.0
    filt[np.abs(f) <= (1.0 - beta) / (2.0 * T)] = 1.0
    return filt


def _bw_filter_vec(f: np.ndarray, BW: float, filter_type: str, beta: float) -> np.ndarray:
    df = f[1] - f[0]
    cs = int(np.floor(len(f) / 2.0))
    filt = np.zeros_like(f, dtype=float)
    if filter_type == "ga":
        filt = np.exp(-((f / (BW / 2.0)) ** 2) * np.log(2.0))
    elif filter_type == "rc":
        filt = raised_cosine(f, 1.0 / BW, beta)
    elif filter_type == "bl":
        ns = int(np.floor(BW / df * 1000.0 / 405.1 / 2.0))
        b = blackman(2 * ns + 1)
        if len(b) <= len(f):
            filt[cs - ns : cs + ns + 1] = b
        else:
            bi = (np.arange(len(f)) - cs + ns).astype(int)
            filt = b[bi]
    elif filter_type == "bh":
        ns = int(np.floor(BW / df * 1000.0 / 342.8 / 2.0))
        bh = blackmanharris(2 * ns + 1)
        if len(bh) <= len(f):
            filt[cs - ns : cs + ns + 1] = bh
        else:
            bi = (np.arange(len(f)) - cs + ns).astype(int)
            filt = bh[bi]
    else:
        raise ValueError(f"Unsupported filter_type: {filter_type}")
    return filt


def bw_window(array: np.ndarray, f: np.ndarray, BW: float, filter_type: str = "bh", beta: float = 1.0 / 3.0) -> tuple[np.ndarray, np.ndarray]:
    if BW is None or BW == 0:
        return np.asarray(array), np.ones(np.asarray(array).shape[0], dtype=float)

    arr = np.asarray(array)
    f = np.asarray(f).reshape(-1)
    if f.size != arr.shape[0]:
        raise ValueError("Number of samples in data and frequency vector must agree")

    filt = _bw_filter_vec(f, BW, filter_type, beta)
    reshaped = filt.reshape((-1,) + (1,) * (arr.ndim - 1))
    return arr * reshaped, filt


def bw_filter(k_data: np.ndarray, t: np.ndarray, BW: float, filter_type: str = "rc", beta: float = 1.0 / 3.0) -> np.ndarray:
    if BW is None or BW == 0:
        return np.asarray(k_data)

    k = np.asarray(k_data)
    t = np.asarray(t).reshape(-1)
    if t.size != k.shape[0]:
        raise ValueError("Number of samples in input data and time vector must agree")

    is_real = np.max(np.abs(np.imag(k))) == 0

    linear_phase = np.zeros_like(k, dtype=float)
    if not is_real:
        for c1 in range(k.shape[1]):
            for c2 in range(k.shape[2]):
                tmp = np.unwrap(np.angle(k[:, c1, c2]))
                b = np.polyfit(t, tmp, 1)
                linear_phase[:, c1, c2] = b[0] * t
        k = k * np.exp(-1j * linear_phase)

    k = np.concatenate([np.flip(k, axis=0), k], axis=0)
    k = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(k, axes=0), axis=0), axes=0)

    dt = t[1] - t[0]
    nrs = k.shape[0]
    df = 1.0 / (nrs * dt)
    f = (np.arange(nrs) - np.floor(nrs / 2.0)) * df

    k, _ = bw_window(k, f, BW, filter_type, beta)
    k = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(k, axes=0), axis=0), axes=0)
    k = k[k.shape[0] // 2 :, ...]

    if not is_real:
        k = k * np.exp(1j * linear_phase)
    if is_real:
        k = np.real(k)
    return k


def variable_smoothing(SIRF: np.ndarray, f: np.ndarray, f_max: float, BW_range: tuple[float, float] | None = None, vsBW: float | None = None):
    if BW_range is None:
        BW_kmin, BW_kmax = 0.05, 10.0
    else:
        BW_kmin, BW_kmax = BW_range
    if vsBW is None:
        vsBW = 40.0

    SIRF = np.asarray(SIRF)
    f = np.asarray(f).reshape(-1)

    df = f[1] - f[0]
    Ts = 10.0 / df
    dtk = 0.01 / BW_kmax
    tk = centered_time(dtk, int(np.round(Ts / dtk)))
    fK, _, _ = time2freq(tk)

    pb = 1.0 / BW_kmax
    beta = 1.0 / 3.0
    TK = pb / (1.0 - beta)
    BW = 1.0 / TK

    k = raised_cosine(tk, BW, beta)
    K = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(k)))
    K = K / np.max(np.abs(K))

    ind0 = np.where(np.abs(K) < 1e-6)[0]
    cs = int(np.ceil(len(ind0) / 2.0)) - 1
    nSL = 3
    K = K[ind0[cs - nSL] : ind0[cs + nSL + 1] + 1]
    fK = fK[ind0[cs - nSL] : ind0[cs + nSL + 1] + 1]

    csS = int(np.floor(SIRF.shape[0] / 2.0))
    nrf1 = int(np.fix(-fK[0] / df))
    nrf2 = int(np.fix(fK[-1] / df))
    nrBW = int(np.floor(f_max / df))
    ind_SK = np.arange(csS - nrf1, csS + nrBW + nrf2 + 1)

    f_sel = f[ind_SK]
    SIRF_sel = SIRF[ind_SK, ...]

    nS = nrBW + 1
    nOut = SIRF_sel.shape[1]
    nIn = SIRF_sel.shape[2]
    S_temp = np.zeros((nS, nOut, nIn), dtype=complex)

    ft = f_sel[: 1 + nrf1 + nrf2]
    vsr = BW_kmax / BW_kmin
    sci = np.where(f_sel < vsBW)[0][-1]
    scf = (1.0 / vsr - 1.0) * np.exp(-((f_sel[nrf1 : sci + 1] / 20.0) ** 2)) + 1.0
    nrscf = len(scf)

    Kmat = np.zeros((len(ft), nS), dtype=complex)
    for l in range(nrscf):
        Ki = np.interp(ft, fK * scf[l], K, left=0, right=0)
        Ki = Ki * np.exp(1j * 2.0 * np.pi * ft / scf[l] * pb / 2.0)
        Ki = Ki / np.abs(np.sum(Ki))
        Kmat[:, l] = Ki

    Ki = np.interp(ft, fK, K, left=0, right=0)
    Ki = Ki * np.exp(1j * 2.0 * np.pi * ft * pb / 2.0)
    Ki = Ki / np.abs(np.sum(Ki))
    if nrscf < nS:
        Kmat[:, nrscf:] = np.tile(Ki[:, None], (1, nS - nrscf))

    for j in range(nS):
        chunk = SIRF_sel[j : j + nrf1 + nrf2 + 1, ...]
        S_temp[j, ...] = np.sum(chunk * Kmat[:, j][:, None, None], axis=0)

    SIRF_sm = np.concatenate([np.conj(S_temp[-1:0:-1, ...]), S_temp], axis=0)
    # Match MATLAB indexing: -f(nrf1+1+nrBW:-1:nrf1+2), f(nrf1+1:nrf1+1+nrBW)
    fs = np.concatenate([-f_sel[nrf1 + nrBW : nrf1 : -1], f_sel[nrf1 : nrf1 + nrBW + 1]])

    return SIRF_sm, fs, f_max, (BW_kmin, BW_kmax), vsBW


def sweeps(t: np.ndarray, T_acq: Any, f1: Any, f2: Any, phi0: Any, A: Any, t_start: Any, sweep_type: Any, AM: dict | None = None, smoothing: dict | None = None, slew: Any = None):
    if AM is None:
        AM = {"type": "none", "slew": None}
    if smoothing is None:
        smoothing = {"type": "none"}

    t = np.asarray(t).reshape(-1)
    sweep = np.zeros_like(t, dtype=float)
    f_t = np.zeros_like(t, dtype=float)

    if isinstance(sweep_type, (list, tuple)):
        n_segments = len(sweep_type)
        if n_segments > 1:
            sweep1, f_t1, phi_end1 = sweeps(t, T_acq[:-1], f1, f2[:-1], phi0, 1, t_start, sweep_type[:-1])
            sweep += sweep1
            f_t += f_t1
            f1 = f2[-2]
            f2 = f2[-1]
            phi0 = phi_end1
            t_start = t_start + np.sum(T_acq[:-1])
            T_acq = T_acq[-1]
            sweep_type = sweep_type[-1]
        else:
            sweep_type = sweep_type[0]

    if not isinstance(sweep_type, str):
        n = sweep_type
    else:
        mapping = {"beta": 2, "linear": 2, "gamma": 3, "delta": 4, "betaAM": 2, "gammaAM": 3, "deltaAM": 4}
        n = mapping[sweep_type]
        if sweep_type.endswith("AM"):
            AM = {"type": "gaussian", "A_max": 10, "width": 1}

    t_shift = t - t_start
    inds = np.where((t_shift >= 0) & (t_shift <= T_acq))[0]

    beta = 2.0 * np.pi * (f2 - f1) / (n * (T_acq ** (n - 1)))
    sweep[inds] = np.sin(beta * (t_shift[inds] ** n) + 2.0 * np.pi * f1 * t_shift[inds] + phi0)
    f_t[inds] = (n * beta / (2.0 * np.pi)) * (t_shift[inds] ** (n - 1)) + f1
    phi_end = np.mod(beta * (T_acq ** n) + 2.0 * np.pi * f1 * T_acq + phi0, 2.0 * np.pi)

    AM_boost = np.ones_like(sweep)
    if AM.get("type") == "gaussian":
        AM_boost[inds] = (AM["A_max"] - 1.0) * np.exp(-((f_t[inds] / AM["width"]) ** 2)) + 1.0

    AM_smooth = np.ones_like(sweep)
    if smoothing.get("type") == "erf":
        cf = f2 - smoothing["width"] / 2.0
        AM_smooth[inds] = 0.5 - 0.5 * erf((f_t[inds] - cf) / (smoothing["width"] / 3.0))
    elif smoothing.get("type") == "time-defined":
        smooth_time = smoothing["width"]
        smooth_ct = T_acq - smooth_time / 2.0
        AM_smooth[inds] = 0.5 - 0.5 * erf((t_shift[inds] - smooth_ct) / (smooth_time / 3.0))

    A_t = A * AM_boost
    if slew:
        S_limit = np.abs(slew / (2.0 * np.pi * f_t))
        S_limit[np.isnan(S_limit)] = 0
        A_t[inds] = np.sign(A) * np.minimum(np.abs(A_t[inds]), S_limit[inds])

    sweep = A_t * AM_smooth * sweep
    return sweep.reshape(-1, 1), f_t.reshape(-1, 1), phi_end


def trapezoid(t: np.ndarray, ons: np.ndarray | float, amp: np.ndarray | float, dur: np.ndarray | float, plateau: np.ndarray | float | None = None, dur2: np.ndarray | float | None = None) -> np.ndarray:
    t = np.asarray(t).reshape(-1)
    ons = np.atleast_1d(ons)
    amp = np.atleast_1d(amp)
    dur = np.atleast_1d(dur)
    n_pulses = len(amp)

    if plateau is None:
        plateau = np.zeros(n_pulses)
    else:
        plateau = np.atleast_1d(plateau)

    if dur2 is None:
        dur2 = dur.copy()
    else:
        dur2 = np.atleast_1d(dur2)

    grads = np.zeros((len(t), n_pulses), dtype=float)
    for i in range(n_pulses):
        t0 = np.where(t > ons[i])[0][0]
        ts1 = np.where(t < ons[i] + dur[i])[0][-1]
        ts2 = np.where(t < ons[i] + dur[i] + plateau[i] + dur2[i])[0][-1]
        slope1 = amp[i] / dur[i]
        slope2 = amp[i] / dur2[i]
        grads[t0 : ts1 + 1, i] = slope1 * (t[t0 : ts1 + 1] - ons[i])
        grads[ts1 + 1 : ts2 + 1, i] = amp[i] - slope2 * (t[ts1 + 1 : ts2 + 1] - (ons[i] + plateau[i] + dur[i]))
        if plateau[i] != 0:
            tp = np.where(t < ons[i] + dur[i] + plateau[i])[0][-1]
            grads[ts1 + 1 : tp + 1, i] = amp[i]
    return grads


def compute_inputs(dt_in: float, t_out: np.ndarray, params: dict, pulse_type: list[str], ip_type: str) -> np.ndarray:
    if not pulse_type or not isinstance(pulse_type, list):
        raise ValueError("pulse_type must be a list, e.g. ['blips'] or ['sweeps']")

    inputs_all = []
    for ptype in pulse_type:
        is_blips = ptype.lower() in ("blips", "blip")
        is_sweeps = ptype.lower() in ("sweeps", "sweep")

        if is_blips:
            n_pulses = len(params["dur"])
            dur = np.asarray(params["dur"])
            t0 = np.asarray(params.get("t0", [0]))
            amp = np.asarray(params.get("amp", [1]))
            fixed_slope = params.get("fixedSlope", None)
            if fixed_slope is not None:
                amp = fixed_slope * dur
            if t0.size < n_pulses:
                t0 = np.full_like(dur, t0.item())
            if amp.size < n_pulses:
                amp = np.full_like(dur, amp.item())
            plateau = np.asarray(params.get("plateau", np.zeros_like(dur)))
            dur2 = np.asarray(params.get("dur2", dur))
        elif is_sweeps:
            Tsweep = params["Tsweep"]
            if not isinstance(Tsweep, list):
                Tsweep = [Tsweep]
                params["Tsweep"] = Tsweep
                for k in ("ampSweep", "phi0", "f1", "f2", "sweepType", "t0Sweep", "slew", "AM", "smoothing"):
                    params[k] = [params[k]]
            n_pulses = len(params["Tsweep"])
        else:
            raise ValueError(f"Unsupported pulse type: {ptype}")

        if ip_type == "ideal":
            if is_blips:
                inputs = trapezoid(t_out, t0, amp, dur, plateau, dur2)
            else:
                inputs = np.zeros((len(t_out), n_pulses))
                for i in range(n_pulses):
                    s, _, _ = sweeps(
                        t_out,
                        params["Tsweep"][i],
                        params["f1"][i],
                        params["f2"][i],
                        params["phi0"][i],
                        params["ampSweep"][i],
                        params["t0Sweep"][i],
                        params["sweepType"][i],
                        params["AM"][i],
                        params["smoothing"][i],
                        params["slew"][i],
                    )
                    inputs[:, i] = s[:, 0]
        else:
            inputs = np.zeros((len(t_out), n_pulses))
            for i in range(n_pulses):
                if is_blips:
                    t0i = t0[i]
                    pulse_dur = dur[i] + plateau[i] + dur2[i]
                    t_in = np.arange(0, pulse_dur + 2 * dt_in, dt_in) - dt_in / 2
                    pulse_in = trapezoid(t_in, 0, amp[i], dur[i], plateau[i], dur2[i])[:, 0]
                else:
                    t0i = params["t0Sweep"][i]
                    pulse_dur = np.sum(params["Tsweep"][i])
                    t_in = np.arange(0, pulse_dur + 2 * dt_in, dt_in) - dt_in / 2
                    pulse_in = sweeps(
                        t_in,
                        params["Tsweep"][i],
                        params["f1"][i],
                        params["f2"][i],
                        params["phi0"][i],
                        params["ampSweep"][i],
                        0,
                        params["sweepType"][i],
                        params["AM"][i],
                        params["smoothing"][i],
                        params["slew"][i],
                    )[0][:, 0]

                if ip_type == "sampleAndHold":
                    mask = np.where((t_out >= t0i) & (t_out <= t0i + pulse_dur))[0]
                    inds_in = np.ceil((t_out[mask] - t0i) / dt_in).astype(int)
                    inputs[mask, i] = pulse_in[inds_in]
                elif ip_type == "interpolated":
                    mask = np.where((t_out >= t_in[0] + t0i) & (t_out <= t_in[-1] + t0i))[0]
                    inputs[mask, i] = np.interp(t_out[mask], t_in + t0i, pulse_in)
                else:
                    raise ValueError("ip_type must be 'ideal', 'sampleAndHold', or 'interpolated'")

        inputs_all.append(inputs)

    return np.concatenate(inputs_all, axis=1)
