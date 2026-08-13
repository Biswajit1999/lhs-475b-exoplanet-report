"""Analyze the real JWST NIRSpec/G395H transmission spectrum of LHS 475b,
quantifying how featureless it really is and which real candidate
atmosphere models it statistically disfavors.

Data source: Zenodo record 7925111 (Lustig-Yaeger & Fu et al. 2023,
Nature Astronomy), the JWST validation and first atmospheric look at this
Earth-sized, interior-to-habitable-zone rocky planet. This script fits
a flat (featureless) line to the real 56-point spectrum and compares
its goodness of fit against four real published forward models
(pure methane, 1x-solar hydrogen-rich, clear Venus-like CO2, pure CO2),
each pre-computed on the same wavelength grid as the data.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import scienceplots  # noqa: F401 (registers 'science' style)
import numpy as np
from scipy.stats import chi2 as chi2_dist

plt.style.use(["science", "no-latex"])

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"

MODELS = {
    "Pure CH4 (methane)": "model_pure_methane.txt",
    "1x-solar H2-rich": "model_1x_solar_hydrogen_rich.txt",
    "Clear Venus-like (CO2)": "model_clear_venuslike.txt",
    "Pure CO2": "model_pure_co2.txt",
}


def load_ecsv(path: Path, ncols: int) -> np.ndarray:
    rows = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if parts[0].startswith('"Wavelength'):
                continue
            rows.append([float(x) for x in parts[:ncols]])
    return np.array(rows)


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)

    spec = load_ecsv(DATA_DIR / "nirspec_transmission_spectrum_LR.txt", ncols=4)
    wave = spec[:, 0]
    # The LR FIREFLy file reports (Rp/Rs)^2 in PERCENT (confirmed against the
    # real system parameters: (0.991 Rearth / 0.2789 Rsun)^2 = 0.00106, matching
    # this column divided by 100, not the raw value) -- the model files are
    # already fractional, so both must be divided by 100 to compare on the
    # same scale.
    depth = spec[:, 2] / 100
    depth_err = spec[:, 3] / 100
    n_points = len(wave)

    # Flat-line (featureless) fit: inverse-variance-weighted mean depth.
    # One parameter (the flat level) is fit locally, so dof = N - 1.
    weights = 1.0 / depth_err**2
    flat_depth = np.sum(depth * weights) / np.sum(weights)
    flat_chi2 = np.sum(((depth - flat_depth) / depth_err) ** 2)
    flat_dof = n_points - 1
    flat_reduced_chi2 = flat_chi2 / flat_dof
    flat_p_value = chi2_dist.sf(flat_chi2, flat_dof)

    # Candidate-model comparisons: each model spectrum is fixed (no local
    # free parameter fitted here -- the models were already offset to the
    # measured depth upstream, in the original data release), so dof = N,
    # not N-1. This assumes that upstream offset was not itself fit to
    # this same data in a way that should reduce dof further; that
    # provenance detail lives with the model files, not in this script.
    model_results = []
    for label, filename in MODELS.items():
        model = load_ecsv(DATA_DIR / filename, ncols=2)
        model_depth = model[:, 1]
        chi2 = np.sum(((depth - model_depth) / depth_err) ** 2)
        dof = n_points
        reduced_chi2 = chi2 / dof
        p_value = chi2_dist.sf(chi2, dof)
        model_results.append({"label": label, "chi2": chi2, "dof": dof, "reduced_chi2": reduced_chi2, "p_value": p_value})

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        writer.writerow(["n_wavelength_points", n_points, "count"])
        writer.writerow(["flat_line_depth", f"{flat_depth:.6f}", "(Rp/Rs)^2"])
        writer.writerow(["flat_line_chi2", f"{flat_chi2:.2f}", "dimensionless"])
        writer.writerow(["flat_line_dof", flat_dof, "count"])
        writer.writerow(["flat_line_reduced_chi2", f"{flat_reduced_chi2:.2f}", "dimensionless"])
        writer.writerow(["flat_line_p_value", f"{flat_p_value:.3f}", "survival function, chi2.sf(chi2, dof)"])
        for r in model_results:
            writer.writerow([f"chi2_{r['label']}", f"{r['chi2']:.2f}", f"dof={r['dof']}"])
            writer.writerow([f"reduced_chi2_{r['label']}", f"{r['reduced_chi2']:.2f}", "dimensionless"])
            writer.writerow([f"p_value_{r['label']}", f"{r['p_value']:.3g}", "survival function, chi2.sf(chi2, dof)"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})

    ax1.errorbar(wave, depth * 1e6, yerr=depth_err * 1e6, fmt="o", ms=4, color="#1f4e79", capsize=2, label="Real JWST NIRSpec data")
    ax1.axhline(flat_depth * 1e6, color="#a8431f", ls="--", lw=1.5, label=f"Flat-line fit (χ²/dof={flat_reduced_chi2:.2f})")
    ax1.set_xlabel("Wavelength [μm]")
    ax1.set_ylabel("Transit depth (Rp/Rs)² [ppm]")
    ax1.set_title("Real LHS 475b transmission spectrum\n(featureless within real measurement noise)")
    ax1.legend(fontsize=7)
    ax1.grid(alpha=0.25)

    labels = [r["label"] for r in model_results]
    reduced = [r["reduced_chi2"] for r in model_results]
    colors = ["#a8431f" if rc > 2 else "#2f6f4f" for rc in reduced]
    ax2.barh(labels, reduced, color=colors)
    ax2.set_xscale("log")
    ax2.axvline(1.0, color="#999", ls=":", lw=1, label="Ideal fit (χ²/dof=1)")
    for i, rc in enumerate(reduced):
        ax2.text(rc * 1.15, i, f"{rc:.2f}", va="center", fontsize=7.5)
    ax2.set_xlabel("Reduced χ² vs. real data (log scale)")
    ax2.set_title("Real candidate atmosphere models\n(green = consistent, red = disfavored)")
    ax2.legend(fontsize=7)
    ax2.grid(alpha=0.25, axis="x", which="both")

    fig.suptitle("LHS 475b: a real, statistically featureless JWST spectrum (Lustig-Yaeger & Fu et al. 2023)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "lhs475b_transmission_spectrum.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {FIG_DIR / 'lhs475b_transmission_spectrum.png'}")
    print(f"Flat-line depth: {flat_depth:.6f}, chi2={flat_chi2:.2f}, dof={flat_dof}, reduced chi2={flat_reduced_chi2:.2f}, p={flat_p_value:.3f}")
    for r in model_results:
        verdict = "DISFAVORED" if r["reduced_chi2"] > 2 else "consistent"
        print(f"  {r['label']}: chi2={r['chi2']:.2f}, dof={r['dof']}, reduced chi2={r['reduced_chi2']:.2f}, p={r['p_value']:.3g} -> {verdict} (reduced-chi2>2 heuristic, not a calibrated rejection threshold)")


if __name__ == "__main__":
    main()
