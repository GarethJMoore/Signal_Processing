#!/usr/bin/env python3
"""TAS data pipeline based on CAS_M1_GMoore.ipynb."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit


@dataclass(frozen=True)
class PipelineConfig:
    data_dir: Path
    output_dir: Path
    delimiter: str = "\t"
    background_time_threshold: float = 0.0
    chirp_grad_threshold_scale: float = 0.001
    chirp_fit_maxfev: int = 10000
    save_csv: bool = False


def load_scans(data_dir: Path, delimiter: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    files = sorted(data_dir.glob("*.txt"))
    if not files:
        raise FileNotFoundError(f"No .txt files found in {data_dir}")

    scans = [np.loadtxt(f, delimiter=delimiter) for f in files]
    data_raw = np.stack(scans, axis=0)
    t = np.array(data_raw[0, 0, 1:])
    w = np.array(data_raw[0, 1:, 0])
    return data_raw, t, w


def average_scans(data_raw: np.ndarray) -> np.ndarray:
    return np.mean(data_raw[:, 1:, 1:], axis=0)


def subtract_background(data_raw: np.ndarray, t: np.ndarray, threshold: float) -> tuple[np.ndarray, np.ndarray]:
    neg_mask = t < threshold
    if not np.any(neg_mask):
        raise ValueError(
            "No negative-time columns found for background subtraction. "
            "Adjust --background-time-threshold."
        )
    data_neg = np.mean(data_raw[:, neg_mask], axis=1)
    data_sub = data_raw - data_neg[:, None]
    return data_sub, data_neg


def chirp_func(x: np.ndarray, y0: float, amplitude: float, x0: float, tau: float) -> np.ndarray:
    return y0 + amplitude * np.exp(-(x - x0) / tau)


def estimate_chirp_bounds(times: np.ndarray, scale: float) -> tuple[int, int]:
    grad = np.gradient(np.gradient(times))
    threshold = scale * (np.median(grad) + np.std(grad))
    upper = len(grad) - 1
    lower = 0
    for i in range(int(len(grad) / 2), len(grad)):
        if np.abs(grad[i]) > threshold:
            upper = i
            break
    for i in range(int(len(grad) / 2), 0, -1):
        if np.abs(grad[i]) > threshold:
            lower = i
            break
    if upper <= lower:
        raise ValueError("Failed to estimate chirp bounds; adjust chirp parameters.")
    return lower, upper


def fit_chirp(w: np.ndarray, t: np.ndarray, data: np.ndarray, scale: float, maxfev: int) -> tuple[np.ndarray, np.ndarray]:
    times = np.array([t[np.argmax(data[i, :])] for i in range(len(w))])
    lower, upper = estimate_chirp_bounds(times, scale)
    p0 = [0.0, -1.0, float(np.max(w[lower:upper])), 1.0]
    popt, pcov = curve_fit(chirp_func, w[lower:upper], times[lower:upper], p0, maxfev=maxfev)
    return popt, pcov


def apply_chirp_correction(w: np.ndarray, t: np.ndarray, data: np.ndarray, popt: np.ndarray) -> np.ndarray:
    data_corrected = np.zeros_like(data)
    for i in range(len(w)):
        zero = chirp_func(w[i], *popt)
        td = t - zero
        data_corrected[i, :] = np.interp(t, td, data[i, :], left=np.nan, right=np.nan)
    return data_corrected


def save_outputs(
    output_dir: Path,
    t: np.ndarray,
    w: np.ndarray,
    data_avg: np.ndarray,
    data_neg: np.ndarray,
    data_sub: np.ndarray,
    data_corrected: np.ndarray,
    popt: np.ndarray,
    save_csv: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    np.save(output_dir / "timedelay.npy", t)
    np.save(output_dir / "wavelength.npy", w)
    np.save(output_dir / "data_average.npy", data_avg)
    np.save(output_dir / "data_background.npy", data_neg)
    np.save(output_dir / "data_subtracted.npy", data_sub)
    np.save(output_dir / "data_chirp_corrected.npy", data_corrected)

    metadata = {
        "chirp_fit_params": {
            "y0": popt[0],
            "amplitude": popt[1],
            "x0": popt[2],
            "tau": popt[3],
        },
        "shapes": {
            "timedelay": list(t.shape),
            "wavelength": list(w.shape),
            "data_average": list(data_avg.shape),
            "data_background": list(data_neg.shape),
            "data_subtracted": list(data_sub.shape),
            "data_chirp_corrected": list(data_corrected.shape),
        },
    }
    (output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2))

    if save_csv:
        np.savetxt(output_dir / "timedelay.csv", t, delimiter=",")
        np.savetxt(output_dir / "wavelength.csv", w, delimiter=",")
        np.savetxt(output_dir / "data_chirp_corrected.csv", data_corrected, delimiter=",")


def run_pipeline(config: PipelineConfig) -> None:
    data_raw, t, w = load_scans(config.data_dir, config.delimiter)
    data_avg = average_scans(data_raw)
    data_sub, data_neg = subtract_background(data_avg, t, config.background_time_threshold)
    popt, _ = fit_chirp(w, t, data_sub, config.chirp_grad_threshold_scale, config.chirp_fit_maxfev)
    data_corrected = apply_chirp_correction(w, t, data_sub, popt)
    save_outputs(
        config.output_dir,
        t,
        w,
        data_avg,
        data_neg,
        data_sub,
        data_corrected,
        popt,
        config.save_csv,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Transient Absorption Spectroscopy data processing pipeline."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("Data"),
        help="Directory containing raw .txt scan files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Directory to write processed outputs.",
    )
    parser.add_argument(
        "--delimiter",
        default="\t",
        help="Delimiter used in the raw scan files (default: tab).",
    )
    parser.add_argument(
        "--background-time-threshold",
        type=float,
        default=0.0,
        help="Time threshold (ps) for negative-time background subtraction.",
    )
    parser.add_argument(
        "--chirp-grad-threshold-scale",
        type=float,
        default=0.001,
        help="Scaling factor for chirp gradient threshold detection.",
    )
    parser.add_argument(
        "--chirp-fit-maxfev",
        type=int,
        default=10000,
        help="Maximum function evaluations for chirp curve fitting.",
    )
    parser.add_argument(
        "--save-csv",
        action="store_true",
        help="Also emit CSV versions of key outputs.",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    config = PipelineConfig(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        delimiter=args.delimiter,
        background_time_threshold=args.background_time_threshold,
        chirp_grad_threshold_scale=args.chirp_grad_threshold_scale,
        chirp_fit_maxfev=args.chirp_fit_maxfev,
        save_csv=args.save_csv,
    )
    run_pipeline(config)


if __name__ == "__main__":
    main()
