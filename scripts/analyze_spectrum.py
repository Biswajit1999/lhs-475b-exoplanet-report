"""Reproduce the LHS 475 b spectrum, model, and pipeline audit.

Inputs are the public products in Zenodo record 7925111. The analysis keeps
two questions separate: how a supplied, already offset model fits as published,
and how its wavelength-dependent shape fits after one local vertical offset.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2 as chi2_dist

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
MODEL_DIR = DATA_DIR / "models"
PIPELINE_DIR = DATA_DIR / "pipeline_spectra"
FIG_DIR = ROOT / "figures"

MODEL_LABELS = {
    "1000xsolar": "1000× solar H₂-rich",
    "100xsolar": "100× solar H₂-rich",
    "10xsolar": "10× solar H₂-rich",
    "1xsolar": "1× solar H₂-rich",
    "ch4": "Pure CH₄",
    "clear_titanlike": "Clear Titan-like",
    "clear_venuslike": "Clear Venus-like",
    "cloudy_venuslike": "Cloudy Venus-like",
    "co2": "Pure CO₂",
    "earthlike": "Earth-like",
    "h2o": "Pure H₂O",
    "hazy_titanlike": "Hazy Titan-like",
    "marslike": "Mars-like",
}


def load_ecsv(path: Path, ncols: int) -> np.ndarray:
    """Load numeric columns from the simple commented ECSV products."""
    rows: list[list[float]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            try:
                rows.append([float(value) for value in line.split()[:ncols]])
            except ValueError:
                continue
    array = np.asarray(rows, dtype=float)
    if array.ndim != 2 or array.shape[1] != ncols or not np.all(np.isfinite(array)):
        raise ValueError(f"Invalid {ncols}-column numeric product: {path}")
    return array


def weighted_mean(values: np.ndarray, errors: np.ndarray) -> tuple[float, float]:
    if values.shape != errors.shape or values.size == 0:
        raise ValueError("values and errors must be non-empty and have equal shape")
    if np.any(~np.isfinite(values)) or np.any(~np.isfinite(errors)) or np.any(errors <= 0):
        raise ValueError("values must be finite and uncertainties finite and positive")
    weights = errors**-2
    return float(np.sum(weights * values) / np.sum(weights)), float(np.sum(weights) ** -0.5)


def fit_constant(values: np.ndarray, errors: np.ndarray) -> dict[str, float | int]:
    level, level_error = weighted_mean(values, errors)
    statistic = float(np.sum(((values - level) / errors) ** 2))
    dof = values.size - 1
    return {
        "level": level,
        "level_error": level_error,
        "chi2": statistic,
        "dof": dof,
        "reduced_chi2": statistic / dof,
        "p_value": float(chi2_dist.sf(statistic, dof)),
    }


def compare_model(observed: np.ndarray, errors: np.ndarray, model: np.ndarray) -> dict[str, float | int]:
    if observed.shape != model.shape:
        raise ValueError("model and observed spectrum grids must have equal length")
    if np.any(errors <= 0):
        raise ValueError("uncertainties must be positive")
    fixed_chi2 = float(np.sum(((observed - model) / errors) ** 2))
    fixed_dof = observed.size
    offset, offset_error = weighted_mean(observed - model, errors)
    shape_chi2 = float(np.sum(((observed - model - offset) / errors) ** 2))
    shape_dof = observed.size - 1
    return {
        "fixed_chi2": fixed_chi2,
        "fixed_dof": fixed_dof,
        "fixed_reduced_chi2": fixed_chi2 / fixed_dof,
        "fixed_p_value": float(chi2_dist.sf(fixed_chi2, fixed_dof)),
        "offset": offset,
        "offset_error": offset_error,
        "shape_chi2": shape_chi2,
        "shape_dof": shape_dof,
        "shape_reduced_chi2": shape_chi2 / shape_dof,
        "shape_p_value": float(chi2_dist.sf(shape_chi2, shape_dof)),
    }


def rebin_hr_to_lr(hr_wave: np.ndarray, hr_depth: np.ndarray, hr_error: np.ndarray,
                   lr_wave: np.ndarray, lr_width: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Inverse-variance rebin a high-resolution spectrum into LR bin edges."""
    values = np.full(lr_wave.shape, np.nan)
    errors = np.full(lr_wave.shape, np.nan)
    for index, (centre, width) in enumerate(zip(lr_wave, lr_width, strict=True)):
        mask = (hr_wave >= centre - width / 2) & (hr_wave < centre + width / 2)
        if np.any(mask):
            values[index], errors[index] = weighted_mean(hr_depth[mask], hr_error[mask])
    return values, errors


def compare_pipeline_shapes(reference: np.ndarray, reference_error: np.ndarray,
                            candidate: np.ndarray, candidate_error: np.ndarray) -> dict[str, float | int]:
    valid = np.isfinite(candidate) & np.isfinite(candidate_error)
    combined_error = np.hypot(reference_error[valid], candidate_error[valid])
    result = compare_model(reference[valid], combined_error, candidate[valid])
    residual = reference[valid] - candidate[valid] - float(result["offset"])
    return {
        "n_bins": int(np.count_nonzero(valid)),
        "offset": result["offset"], "offset_error": result["offset_error"],
        "shape_chi2": result["shape_chi2"], "shape_dof": result["shape_dof"],
        "shape_reduced_chi2": result["shape_reduced_chi2"],
        "shape_p_value": result["shape_p_value"],
        "residual_rms": float(np.sqrt(np.mean(residual**2))),
    }


def model_name(path: Path) -> str:
    key = path.stem.removeprefix("lhs_475b_").removesuffix("_model_ecsv")
    return MODEL_LABELS[key]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    spectrum = load_ecsv(DATA_DIR / "nirspec_transmission_spectrum_LR.txt", 4)
    wave, width = spectrum[:, 0], spectrum[:, 1]
    # The FIREFLy LR product uses percent; models and HR products are fractional.
    depth, depth_error = spectrum[:, 2] / 100, spectrum[:, 3] / 100
    flat = fit_constant(depth, depth_error)

    model_rows: list[dict[str, object]] = []
    model_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for path in sorted(MODEL_DIR.glob("*_model_ecsv.txt")):
        model = load_ecsv(path, 2)
        if model.shape[0] != depth.size or np.any(np.diff(model[:, 0]) <= 0):
            raise ValueError(f"Model must contain 56 ordered samples: {path}")
        label = model_name(path)
        model_rows.append({"model": label, **compare_model(depth, depth_error, model[:, 1])})
        model_curves[label] = (model[:, 0], model[:, 1])
    if len(model_rows) != len(MODEL_LABELS):
        raise ValueError(f"Expected {len(MODEL_LABELS)} models; found {len(model_rows)}")

    pipeline_rows: list[dict[str, object]] = []
    rebinned: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for path in sorted(PIPELINE_DIR.glob("transit_spectrum_HR_*_ecsv.txt")):
        pipeline = path.stem.split("_HR_")[1].split("_ecsv")[0].capitalize()
        product = load_ecsv(path, 5)
        values, errors = rebin_hr_to_lr(product[:, 1], product[:, 3], product[:, 4], wave, width)
        valid = np.isfinite(values)
        constant = fit_constant(values[valid], errors[valid])
        sensitivity = compare_pipeline_shapes(depth, depth_error, values, errors)
        pipeline_rows.append({
            "pipeline": pipeline, "n_bins": sensitivity["n_bins"],
            "mean_depth": constant["level"], "median_formal_error": float(np.nanmedian(errors)),
            "flat_chi2": constant["chi2"], "flat_dof": constant["dof"],
            "flat_reduced_chi2": constant["reduced_chi2"],
            "lr_shape_offset": sensitivity["offset"],
            "lr_shape_offset_error": sensitivity["offset_error"],
            "lr_shape_chi2": sensitivity["shape_chi2"],
            "lr_shape_dof": sensitivity["shape_dof"],
            "lr_shape_reduced_chi2": sensitivity["shape_reduced_chi2"],
            "lr_shape_p_value": sensitivity["shape_p_value"],
            "lr_shape_residual_rms": sensitivity["residual_rms"],
        })
        rebinned[pipeline] = (values, errors)

    write_csv(FIG_DIR / "model_comparison.csv", model_rows)
    write_csv(FIG_DIR / "pipeline_comparison.csv", pipeline_rows)
    summary = {
        "source": "Zenodo 7925111", "units": "fractional transit depth unless stated",
        "n_lr_points": int(depth.size), "flat_fit": flat, "model_count": len(model_rows),
        "model_estimands": {
            "fixed": "supplied upstream-offset curve; no locally fitted parameter",
            "shape": "one locally fitted vertical offset; dof=N-1",
        },
        "pipeline_note": "same two transits reduced by different pipelines; not independent observations",
        "models": model_rows, "pipelines": pipeline_rows,
    }
    (FIG_DIR / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    with (FIG_DIR / "summary_statistics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        writer.writerow(["n_wavelength_points", depth.size, "count"])
        for key, value in flat.items():
            writer.writerow([f"flat_line_{key}", value, "fractional/chi-square as applicable"])

    plt.rcParams.update({"font.size": 9, "axes.titleweight": "bold", "figure.facecolor": "white"})
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 12), constrained_layout=True)
    ax = axes[0]
    ax.errorbar(wave, depth * 1e6, yerr=depth_error * 1e6, fmt="o", ms=3.5,
                color="#153e5c", capsize=1.5, label="FIREFLy LR data")
    ax.axhline(float(flat["level"]) * 1e6, color="#c05640", ls="--", lw=1.5,
               label=f"Flat fit: χ²ν={float(flat['reduced_chi2']):.2f}")
    choices = [("1× solar H₂-rich", "#b33a3a"), ("Pure CH₄", "#dd8a1e"),
               ("Cloudy Venus-like", "#6d597a"), ("Mars-like", "#2a9d8f")]
    for label, color in choices:
        model_wave, model_depth = model_curves[label]
        ax.plot(model_wave, model_depth * 1e6, lw=1.15, color=color, alpha=0.9, label=label)
    ax.set(title="A · Published low-resolution spectrum and representative atmosphere curves",
           ylabel="Transit depth [ppm]", xlabel="Wavelength [µm]")
    ax.legend(ncol=3, fontsize=7.5)
    ax.grid(alpha=0.18)

    ax = axes[1]
    ordered = sorted(model_rows, key=lambda row: float(row["shape_reduced_chi2"]))
    ax.barh([str(row["model"]) for row in ordered],
            [float(row["shape_reduced_chi2"]) for row in ordered], color="#457b9d")
    ax.axvline(1, color="#8c8c8c", lw=1, ls=":")
    ax.set_xscale("log")
    ax.set(title="B · Shape-only comparison after one fitted vertical offset",
           xlabel="Reduced χ² (log scale; dof = 55)")
    ax.grid(axis="x", alpha=0.18, which="both")

    ax = axes[2]
    ax.errorbar(wave, (depth - float(flat["level"])) * 1e6, yerr=depth_error * 1e6,
                fmt="o", ms=3, color="#111827", alpha=0.75, label="FIREFLy LR")
    colors = {"Eureka": "#e76f51", "Firefly": "#457b9d", "Tiberius": "#2a9d8f"}
    for row in pipeline_rows:
        pipeline = str(row["pipeline"])
        values_hr, errors_hr = rebinned[pipeline]
        valid = np.isfinite(values_hr)
        ax.errorbar(wave[valid], (values_hr[valid] - float(row["mean_depth"])) * 1e6,
                    yerr=errors_hr[valid] * 1e6, fmt=".", ms=3, lw=0.7, alpha=0.55,
                    color=colors[pipeline], label=f"{pipeline} HR→LR")
    ax.axhline(0, color="#8c8c8c", lw=1, ls=":")
    ax.set(title="C · Reduction-pipeline sensitivity (each series mean-subtracted)",
           ylabel="Relative transit depth [ppm]", xlabel="Wavelength [µm]")
    ax.legend(ncol=4, fontsize=7.5)
    ax.grid(alpha=0.18)
    fig.suptitle("LHS 475 b · JWST/NIRSpec G395H evidence audit", fontsize=15, fontweight="bold")
    fig.savefig(FIG_DIR / "lhs475b_transmission_spectrum.png", dpi=220)
    plt.close(fig)

    print(f"Flat fit: chi2/dof={flat['chi2']:.2f}/{flat['dof']}, p={flat['p_value']:.3f}")
    print(f"Compared {len(model_rows)} atmosphere scenarios and {len(pipeline_rows)} reductions")


if __name__ == "__main__":
    main()
