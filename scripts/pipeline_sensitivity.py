"""Conditional pipeline-difference and spectral-shape diagnostics.

The three reductions share photons; correlation is a sensitivity parameter,
not an estimated inter-pipeline covariance.
"""
import numpy as np


def paired_reduction_difference(first, second, sigma_first, sigma_second, rho=0.0):
    """Difference and standard error under an assumed same-bin correlation."""
    a, b, sa, sb = [np.asarray(x, dtype=float) for x in
                    (first, second, sigma_first, sigma_second)]
    if a.ndim != 1 or a.size < 2 or any(x.shape != a.shape for x in (b, sa, sb)):
        raise ValueError("matched vectors with at least two bins required")
    if any(not np.all(np.isfinite(x)) for x in (a, b, sa, sb)):
        raise ValueError("all values must be finite")
    if np.any(sa <= 0) or np.any(sb <= 0) or not np.isfinite(rho) or abs(rho) > 1:
        raise ValueError("positive uncertainties and -1 <= rho <= 1 required")
    variance = sa**2 + sb**2 - 2 * rho * sa * sb
    if np.any(variance <= 0):
        raise ValueError("degenerate difference variance")
    return a - b, np.sqrt(variance)


def profiled_template_amplitude(depth, template, covariance):
    """GLS template coefficient after profiling a free constant offset.

    The coefficient is only identifiable if the template is nonconstant.
    """
    y = np.asarray(depth, dtype=float)
    t = np.asarray(template, dtype=float)
    c = np.asarray(covariance, dtype=float)
    if y.ndim != 1 or y.size < 3 or t.shape != y.shape:
        raise ValueError("at least three matching spectral bins required")
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(t)):
        raise ValueError("finite data and template required")
    if c.shape != (y.size, y.size) or not np.all(np.isfinite(c)):
        raise ValueError("invalid covariance")
    if not np.allclose(c, c.T, rtol=1e-12, atol=1e-12):
        raise ValueError("covariance must be symmetric")
    l = np.linalg.cholesky(c)
    one = np.linalg.solve(l, np.ones_like(y))
    shape = np.linalg.solve(l, t)
    data = np.linalg.solve(l, y)
    shape = shape - one * (one @ shape) / (one @ one)
    data = data - one * (one @ data) / (one @ one)
    info = float(shape @ shape)
    if info <= 1e-20 * float(t @ t + 1):
        raise ValueError("template indistinguishable from constant offset")
    return float(shape @ data / info), float(info ** -0.5)
