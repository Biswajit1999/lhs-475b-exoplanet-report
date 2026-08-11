# LHS 475b — Exoplanet Atmosphere Report

An Earth-sized rocky planet interior to its M-dwarf's habitable zone,
and one of JWST's first terrestrial-exoplanet targets. This repo
statistically tests the real, published transmission spectrum against
real candidate atmosphere models — quantifying exactly which are ruled
out and which remain possible.

**[Open the full report](index.html)** (open locally in a browser, or serve
with `python -m http.server` from this directory).

## What's real here

- **System parameters** — queried live from the NASA Exoplanet Archive TAP
  service (`pscomppars`).
- **JWST transmission spectrum and atmosphere models** — five real files
  from Zenodo record [7925111](https://zenodo.org/records/7925111)
  (Lustig-Yaeger & Fu et al. 2023, *Nature Astronomy*): the real 56-point
  co-added NIRSpec/G395H spectrum (FIREFLy pipeline reduction), plus four
  real published forward-model spectra (pure methane, 1x-solar
  hydrogen-rich, clear Venus-like CO2, pure CO2).
- **Analysis** — `scripts/analyze_spectrum.py` fits a flat (featureless)
  line to the real spectrum and computes the reduced chi-squared of the
  real data against each real candidate model. Run it yourself:

  ```bash
  pip install -r requirements.txt
  python scripts/analyze_spectrum.py
  ```

## Repository structure

```text
index.html              the report webpage
data/                    real JWST NIRSpec spectrum + atmosphere models (Zenodo 7925111)
scripts/analyze_spectrum.py   real flat-line and model chi-squared analysis
figures/                 generated plot + summary_statistics.csv
```

## Key finding this repo shows directly

The real spectrum is statistically indistinguishable from a flat line
(reduced chi-squared of **0.92** across 56 real wavelength points). A
primordial hydrogen-dominated envelope is decisively ruled out
(chi-squared of **206**), and a cloudless pure-methane atmosphere is
disfavored (chi-squared of **2.3**), while denser, higher-mean-molecular-
weight options like a CO2-dominated, Venus-like atmosphere remain
statistically consistent with the real data (chi-squared ≈ **1.0-1.1**)
— matching the real published conclusion that the data cannot yet
distinguish a thick CO2 atmosphere, a thin Mars-like one, or bare rock.

## Honest limitation

This repo compares only 4 of the paper's ~12 real candidate models (a
representative disfavored pair and a representative consistent pair).
A "consistent" chi-squared means the data cannot rule the model out —
not that it confirms it; a genuinely featureless spectrum is equally
consistent with several very different atmospheres, or none at all.
Also worth noting: the real data file reports transit depth in percent
while the model files are fractional — this repo's script converts both
to the same scale explicitly (see the comment in
`scripts/analyze_spectrum.py`) rather than silently assuming a match.

## References

1. Lustig-Yaeger, J. & Fu, G. et al., 2023. A JWST transmission spectrum
   of the nearby Earth-sized exoplanet LHS 475 b. *Nature Astronomy*, 7,
   pp.1317-1328.
2. Ment, K. et al., 2023. LHS 475 b: A Venus-sized Planet Orbiting a
   Nearby M Dwarf. *The Astronomical Journal* (submitted),
   arXiv:2304.01920.
3. NASA Exoplanet Archive, <https://exoplanetarchive.ipac.caltech.edu/>.
4. Zenodo record 7925111, <https://zenodo.org/records/7925111>.

## Author

Biswajit Jana — [Portfolio](https://biswajit1999.github.io/Biswajit_Jana.github.io/) · [GitHub](https://github.com/Biswajit1999) · [LinkedIn](https://www.linkedin.com/in/biswajit-jana-27011a151/) · [ORCID](https://orcid.org/0009-0002-2411-1891)
