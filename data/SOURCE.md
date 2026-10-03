# Data provenance

All analysis inputs are unmodified numeric products from Zenodo record
[7925111](https://doi.org/10.5281/zenodo.7925111), “Products and Models for
‘A JWST transmission spectrum of the nearby Earth-sized exoplanet LHS 475 b’”
(Lustig-Yaeger, Fu et al. 2023). The archive `LHS47b_Zenodo_Files.zip` has
Zenodo-reported MD5 `fe9d81e069fbef4787506fb78a389876`.

- `nirspec_transmission_spectrum_LR.txt` is the 56-bin co-added FIREFLy
  NIRSpec/G395H spectrum from the two GO 1981 transits. Depth and uncertainty
  are expressed in percent in this one product.
- `models/` contains all 13 released PICASO/CHIMERA atmosphere curves on the
  LR wavelength grid. Model depths are fractional and already carry the
  upstream depth offset.
- `pipeline_spectra/` contains the released high-resolution Eureka!, FIREFLy,
  and Tiberius reductions. Their depths are fractional. They analyze the same
  two visits and are not independent observations.

Files were copied from the archive and renamed only by their source directory
layout. `source_manifest.json` lists normalized SHA-256 hashes, verified by
`python scripts/validate_sources.py`. The line-ending normalization makes the
check stable across Git checkouts on Windows and Linux.
