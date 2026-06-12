# PyNeb Optimization — June 2026

Summary of the performance and robustness work done on branch `FableOptimisation`
(11 commits, `3f815ba`…`574a5e8`). Scope: hot-path performance and robustness
fixes, with numerical behavior preserved and no public-API changes.
Explicitly out of scope (deferred): dead-code removal, replacement of the
vendored `chebyshev`/`polyutils`/`polytemplate` modules, consolidation of the
three `_chianti_tools*` modules, and lazy loading at import time.

## Safety net: golden-value regression tests (`3f815ba`)

Because the existing suite only checked object instantiation, numerical
regression tests were added **before** any refactor:

- `tests/generate_golden.py` — computes 85 reference quantities and stores
  them in `tests/golden/golden_refs.npz`, together with the atomic data
  files in use and the PyNeb version.
- `tests/test_golden_regression.py` — recomputes every quantity with the
  current code and compares to the stored values
  (`rtol=1e-10`; `rtol=1e-5` for the iterative `getTemDen`; structural-zero
  patterns of emissivity/collision matrices must match exactly).

Coverage: `getPopulations` (scalar, product grid, product=False 1-D and 2-D),
`getCollRates`, `getCritDensity`, `getEmissivity` (all transitions, single
line, product=False), `getTemDen` (temperature and density branches),
RecAtom emissivities (H I, He I, He II), the full `Continuum` class
(free-bound, two-photon, free-free, Balmer jump, total), and `X`/`getCorr`
for **all** extinction laws. Atoms: O III, N II, S II and Fe III (144 levels).

If the default atomic data changes on purpose, regenerate with:

```bash
python tests/generate_golden.py
```

Every optimization below was committed only after this suite passed.

## Performance changes

| Commit | Change |
|---|---|
| `e55f407` | **Duplicate linear solve removed** in `getPopulations(product=True)`: the same system was solved twice per (tem, den) grid point — once outside the `try` block and once inside it. |
| `68e043d` | **Vectorized `getPopulations`**: the O(n²) Python loops building `sum_q_up`/`sum_q_down` (both branches) and the per-element coefficient-matrix loops of the `product=False` branch replaced by masked numpy operations. |
| `b62999d` | **Batched LAPACK solves**: all population matrices are solved in a single `np.linalg.solve` call over stacked matrices (bit-identical results); a per-matrix fallback preserves the NaN-per-singular-point semantics. |
| `d44a269` | **Vectorized `getCollRates`**: broadcasting over level pairs; the detailed-balance exponent is masked before `np.exp` to avoid overflow in the unused triangle. |
| `40f2481` | **Vectorized `getEmissivity`** (all-transitions branch): deltaE and the population/density factors computed over all level pairs at once instead of rebuilding the same (n_tem, n_den) denominator per pair. |
| `8baec35` | **Continuum caching**: the Ercolano & Storey free-bound tables were re-read from disk (`np.loadtxt`) and the `interp1d` objects rebuilt on *every* call; both are now cached at class level. The per-wavelength `Delta_E` loop is a broadcast minimum. |
| `e957257` | `poly()` in `red_corr` evaluated via `np.polynomial.polynomial.polyval` (Horner) instead of a per-coefficient Python loop. |
| `dd5bf53` | **Directory-listing cache** in `_ManageAtomicData`: `getDirForFile` ran `os.listdir` per directory on every lookup (several times per Atom construction). Listings are cached as frozensets, invalidated when the path list changes, with a refresh-and-retry on miss so files dropped into a data directory mid-session are still found. |

Measured impact (macOS, Python 3.12, numpy 1.26):

| Benchmark | Before | After | Speedup |
|---|---|---|---|
| O III emissivity grid 60×60 | 20.5 ms | 3.9 ms | **5.3×** |
| `Continuum.get_continuum`, 50-tem array (warm cache) | 583 ms | 93 ms | **6.3×** |
| Fe III (144-level) populations, 10×10 grid | 232 ms | 164 ms | 1.4× |

## Robustness fixes

- `c6b48ac` — `saverestore.py`: `f.close` was missing parentheses, so files
  were **never closed**; rewritten with context managers.
- `574a5e8` — deprecated numpy aliases that crash on numpy ≥ 1.24 / 2.0
  (`np.int()`, `np.NAN`) replaced by `int()` / `np.nan`; mutable default
  arguments replaced by tuples/`None` in `getAvailableICFs`, `printAllICFs`,
  `getElemAbundance`, `my_X2` and `DataPlot.__init__`;
  `type(x) == type('')` checks became `isinstance` in the touched functions.
- Bare `except:` clauses narrowed (to `np.linalg.LinAlgError` /
  `except Exception` / `except ValueError`) in the functions modified above.

## Intentional behavior changes

Both are strict improvements over crashes/garbage and are documented in the
commit messages:

1. A **singular population matrix** with `product=True` now yields NaN for
   that grid point (consistent with the `product=False` branch) instead of
   raising out of the solve loop.
2. **Continuum wavelengths below every ionization threshold** now return a
   null free-bound continuum instead of raising `ValueError` from an empty
   `np.min`.
