"""Executable checks on the ECSV loader and the chi-squared/p-value
formula, and a regression guard that the pipeline still reproduces the
documented p-values (added specifically so the reduced-chi2>2 heuristic
wasn't the only statistical statement in this repo) when run on the
real downloaded spectrum and models."""

import csv

import numpy as np
from scipy.stats import chi2 as chi2_dist
import analyze_spectrum as spec


def test_load_ecsv_returns_expected_shape():
    spec_data = spec.load_ecsv(spec.DATA_DIR / "nirspec_transmission_spectrum_LR.txt", ncols=4)
    assert spec_data.shape[0] == 56
    assert spec_data.shape[1] == 4


def test_chi2_p_value_matches_known_case():
    values = np.full(10, 5.0)
    errors = np.full(10, 1.0)
    weights = 1.0 / errors**2
    flat = np.sum(values * weights) / np.sum(weights)
    chi2 = np.sum(((values - flat) / errors) ** 2)
    assert np.isclose(chi2, 0.0, atol=1e-10)
    assert np.isclose(chi2_dist.sf(chi2, len(values) - 1), 1.0)


def test_pipeline_reproduces_documented_pvalues():
    spec.FIG_DIR.mkdir(exist_ok=True)
    spec.main()
    rows = {}
    with (spec.FIG_DIR / "summary_statistics.csv").open() as f:
        for row in csv.DictReader(f):
            rows[row["quantity"]] = row["value"]
    assert int(rows["n_wavelength_points"]) == 56
    assert abs(float(rows["flat_line_chi2"]) - 50.70) < 0.1
    assert abs(float(rows["flat_line_p_value"]) - 0.640) < 0.01
    assert float(rows["p_value_Pure CO2"]) > 0.05  # consistent
    assert float(rows["p_value_Pure CH4 (methane)"]) < 0.001  # disfavored
