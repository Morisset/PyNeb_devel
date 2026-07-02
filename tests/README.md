# PyNeb test suite

14 test files, 241 tests. Run from the repository root with the same
invocation as CI (`.github/workflows/tests.yml`, Python 3.8 and 3.12):

```
pytest tests --doctest-modules
```

Note: `--doctest-modules` also collects docstring examples, so test files
must not contain `>>>` in their docstrings.

## Shared fixtures — [conftest.py](conftest.py)

- Session-scoped atoms built once per run: `O3`, `N2`, `S2`, `O2`
  (`pn.Atom`) and `H1`, `He2` (`pn.RecAtom`). Tests must treat them as
  read-only.
- `O3_grid`: a small 10x8 `pn.EmisGrid` (Te 8000-12000 K, Ne 1e2-1e4),
  shared read-only by the emisGrid and plotting tests.
- `smc24_path`: path to the `smc24.dat` observation file shipped in
  `pyneb/sample_scripts` (`lines_in_rows` format).
- Forces the matplotlib Agg backend before `pyneb` is imported.

## Core physics classes

**[test_atom.py](test_atom.py) — 60 tests.** The collisional `Atom` class.
Constructor variants and attributes (`elem`/`spec`/`Z`/`NLevels`, energy and
wavelength arrays); `getA` (matrix structure, known [OIII] 5007/4959 values,
lookup by wavelength); `getEnergy` ordering; `getOmega`/`getCollRates`
shapes and positivity; `getPopulations` normalization; `getEmissivity`
(absolute values for 5007/4959/4363, line ratios, scalar/array/product
shapes, temperature dependence); `getCritDensity`; `getTemDen` roundtrips
(N2 temperature, S2 density) with `to_eval`, wave and level syntaxes;
`getTransition`; `printIonic` smoke tests.

**[test_recatom.py](test_recatom.py) — 18 tests.** The `RecAtom` class
(H I, He II): constructor attributes, `getWave` for Halpha/Hbeta/He II 4686,
emissivity values including the Balmer decrement (~2.85) and the
Hgamma/Hbeta ratio, temperature dependence, density insensitivity, and the
module-level `getRecEmissivity` helper's agreement with `RecAtom`.

**[test_continuum.py](test_continuum.py) — 11 tests.** The `Continuum`
class: continuum shape/positivity/magnitude, the Balmer-jump discontinuity
sign, normalized output shapes for temperature arrays, `BJ_HI` decreasing
with Te, a `T_BJ` roundtrip recovering 1e4 K to 0.1%, the closed-form
`T_BJ_Liu`, the two-photon continuum (zero below Ly-alpha, suppressed at
high density), Gaunt factors of order unity, free-free increasing with
helium, and NaN returns for invalid `make_cont_Ercolano` inputs.

## Observations and analysis

**[test_observation.py](test_observation.py) — 39 tests.** `EmissionLine`
(constructors by elem/spec/wave and by label, blends, invalid labels,
relative vs absolute errors, `correctIntens`, `addObs`, `printLine`) and
`Observation`: reading all four ASCII formats (`lines_in_cols`,
`lines_in_cols2`, `lines_in_rows` via the shipped `smc24.dat`,
`lines_in_rows_err_cols`), line access (`getLine`, `getIntens`, `getError`,
sorting, unique atoms), editing (`addLine`, `removeLine`, `fillObs`,
`addObs`, `addSum`), extinction (`def_EBV` + `correctData` restoring the
Balmer ratio to 2.85), Monte Carlo replication for one and two objects,
and `setAllErrors`.

**[test_diagnostics.py](test_diagnostics.py) — 16 tests.** The
`Diagnostics` class: the diagnostics catalog (>=100 entries, known labels),
`addDiag`/`addAll`/`delDiag`/`getDiagFromLabel`, `getUniqueAtoms`, and
`getCrossTemDen` (an [OIII]+[SII] roundtrip and an [NII]+[SII]
physical-range check).

**[test_icf.py](test_icf.py) — 12 tests.** The `ICF` class: catalog queries
(`getAvailableICFs`, filtered by type), metadata (`getExpression`,
`getType`, `getComment`, `getReference`, `getURL`), and `getElemAbundance`
with hand-checkable inputs (direct O sum = 3e-5, ICF factor exactly 3 for
TPP77_14, ICF >= 1 sanity), array inputs, NaN for absent ions, and
`addICF`/`delICF`.

**[test_emisgrid.py](test_emisgrid.py) — 13 tests.** The `EmisGrid` class:
grid shapes and Te/Ne ranges, `getGrid` by wavelength/levels/`to_eval`,
physical invariants (5007/4959 constant ~2.98, 4363/5007 strictly
increasing with Te), save/restore roundtrip and recompute-on-mismatch,
`getEmisGridDict` compute-then-cache behavior, and smoke tests for
`plotImage`/`plotContours`/`plotLineRatio`.

## Extinction, data management, utilities

**[test_extinction.py](test_extinction.py) — 21 tests.** The `RedCorr`
class: defaults (R_V = 3.1), E(B-V) <-> cHbeta roundtrip, `getCorr` (unity
for no correction or zero reddening, wavelength dependence, array input),
`getCorrHb`, `setCorr` from observed ratios, and the catalog of >=12 laws
(CCM89, F99, Cal00, ...).

**[test_atomic_manager.py](test_atomic_manager.py) — 7 tests.** The
`pn.atomicData` singleton: `getDataFile`, `getAllAvailableFiles` (with the
`*` current-file marker), `getAllAtoms` (collisional and recombination),
predefined dictionaries (`PYNEB_23_01`), and a `setDataFile`
switch-and-restore roundtrip verified on a freshly built atom.

**[test_chianti.py](test_chianti.py) — 1 test.** Loads O III
atom+collision data in CHIANTI format from the bundled `tests/CHIANTI/`
tree and builds a 19-level atom.

**[test_physics.py](test_physics.py) — 10 tests.** `airtovac`/`vactoair`
(IDL reference value, roundtrip, vacuum > air, no-op outside the conversion
window), `gsFromAtom` ground states, `CST` constants (Hbeta wavelength,
Rydberg, Boltzmann), and `getHbEmissivity` against both the Aller fit and
H I recombination data.

**[test_helpers.py](test_helpers.py) — 14 tests.** `getLineLabel` /
`parseLineLabel` (plain lines, blends, micron wavelengths, `H1r` labels),
`isValid`, `getAtomDict`, and `save`/`restore` (kwargs and positional
paths, error on unknown names).

**[test_utils.py](test_utils.py) — 14 tests.** `int_to_roman`/
`roman_to_int` (including a full 1-3999 roundtrip), `parseAtom`, and
`ROOT_DIR`.

**[test_plotting.py](test_plotting.py) — 5 tests.** Agg-backend smoke
tests that assert content was actually drawn: `Atom.plotEmiss`,
`Atom.plotGrotrian`, `RedCorr.plot`, `DataPlot.plotA` (loads every O III
dataset), and `Diagnostics.plot` over the `O3_grid` with a synthetic
corrected observation.

## Conventions

- Tests are grouped in classes; numeric comparisons use `pytest.approx`,
  preferring physical-sanity assertions (known line ratios, monotonic
  trends) over exact snapshots.
- Never mutate the session-scoped fixtures or global state
  (`pn.atomicData`, `pn.config`) without restoring in a `finally` block.
- All disk writes go to pytest's `tmp_path`.
- `pn.log_.error` raises `PyNebError`; tests assert the exception and that
  no state changed.
