"""Regression and validation tests for the published LHS 475 b products."""

import csv

import analyze_spectrum as analysis
import numpy as np
import pytest


def test_lr_product_shape_and_units():
    data = analysis.load_ecsv(analysis.DATA_DIR / "nirspec_transmission_spectrum_LR.txt", 4)
    assert data.shape == (56, 4)
    assert np.all(data[:, 3] > 0)


def test_weighted_mean_rejects_invalid_uncertainty():
    with pytest.raises(ValueError, match="positive"):
        analysis.weighted_mean(np.array([1.0, 2.0]), np.array([0.1, 0.0]))


def test_model_catalog_is_complete_and_regression_values_hold():
    paths = sorted(analysis.MODEL_DIR.glob("*_model_ecsv.txt"))
    assert len(paths) == 13
    spectrum = analysis.load_ecsv(analysis.DATA_DIR / "nirspec_transmission_spectrum_LR.txt", 4)
    depth, error = spectrum[:, 2] / 100, spectrum[:, 3] / 100
    results = {}
    for path in paths:
        model = analysis.load_ecsv(path, 2)
        results[analysis.model_name(path)] = analysis.compare_model(depth, error, model[:, 1])
    assert results["1× solar H₂-rich"]["fixed_chi2"] == pytest.approx(11541.08, abs=0.02)
    assert results["Pure CH₄"]["shape_chi2"] == pytest.approx(108.90, abs=0.02)
    assert results["Mars-like"]["shape_chi2"] == pytest.approx(55.94, abs=0.02)
    assert results["Mars-like"]["shape_dof"] == 55


def test_rebin_and_pipeline_comparison_regression():
    spectrum = analysis.load_ecsv(analysis.DATA_DIR / "nirspec_transmission_spectrum_LR.txt", 4)
    depth, error = spectrum[:, 2] / 100, spectrum[:, 3] / 100
    path = analysis.PIPELINE_DIR / "transit_spectrum_HR_tiberius_ecsv.txt"
    product = analysis.load_ecsv(path, 5)
    values, errors = analysis.rebin_hr_to_lr(
        product[:, 1], product[:, 3], product[:, 4], spectrum[:, 0], spectrum[:, 1]
    )
    assert np.count_nonzero(np.isfinite(values)) == 56
    result = analysis.compare_pipeline_shapes(depth, error, values, errors)
    assert float(result["offset"]) * 1e6 == pytest.approx(37.48, abs=0.03)
    assert float(result["shape_reduced_chi2"]) == pytest.approx(0.592, abs=0.002)


def test_main_writes_machine_readable_products():
    analysis.main()
    with (analysis.FIG_DIR / "model_comparison.csv").open(encoding="utf-8") as handle:
        models = list(csv.DictReader(handle))
    with (analysis.FIG_DIR / "pipeline_comparison.csv").open(encoding="utf-8") as handle:
        pipelines = list(csv.DictReader(handle))
    assert len(models) == 13
    assert {row["pipeline"] for row in pipelines} == {"Eureka", "Firefly", "Tiberius"}
    assert (analysis.FIG_DIR / "analysis_summary.json").is_file()
