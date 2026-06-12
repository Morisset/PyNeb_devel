"""
Generate golden reference values for the numerical regression tests.

Run ONCE on a known-good revision:

    python tests/generate_golden.py

It writes tests/golden/golden_refs.npz, storing the outputs of the main
numerical entry points (populations, collision rates, emissivities,
diagnostics, continuum, extinction laws) together with the atomic data
files in use, so that tests/test_golden_regression.py can detect any
numerical change introduced by a refactor.

If the default atomic data shipped with PyNeb changes on purpose, the
goldens must be regenerated with the command above.
"""

import os
import numpy as np

GOLDEN_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'golden', 'golden_refs.npz')

# Atoms whose data files are recorded at generation time and pinned at test time
PINNED_ATOMS = ['O3', 'N2', 'S2', 'Fe3', 'H1', 'He1', 'He2']

# Keys for which the structural zero pattern (exact 0.0 entries) must be preserved
ZERO_PATTERN_PREFIXES = ('em_all', 'collrates')

# Keys compared with a loose tolerance (iterative solver, converges to 1e-3)
LOOSE_PREFIXES = ('temden',)

TEM_GRID = np.array([8000., 10000., 12000.])
DEN_GRID = np.array([100., 1000., 10000., 100000.])
TEM_PAIRS = np.linspace(8000., 20000., 5)
DEN_PAIRS = np.logspace(2., 5., 5)
WL_CONT = np.array([1200., 3000., 3500., 3640., 3650., 3700., 4861., 9000., 20000.])
WAVES_RC = np.array([1300., 1500., 2000., 2500., 3000., 4000., 4861., 5007.,
                     6563., 9000., 12000., 20000.])


def get_datafiles():
    """Return {atom: 'file1|file2|...'} for the atoms used here (None entries skipped)."""
    import pyneb as pn
    files = {}
    for atom in PINNED_ATOMS:
        entry = pn.atomicData.getDataFile(atom)
        if isinstance(entry, str):
            entry = (entry,)
        files[atom] = '|'.join(f for f in entry if f is not None)
    return files


def pin_datafiles(files):
    """Set the atomic data files recorded in the goldens (inverse of get_datafiles)."""
    import pyneb as pn
    for atom in sorted(files):
        for f in files[atom].split('|'):
            if f:
                pn.atomicData.setDataFile(f)


def compute_all():
    """Compute every golden quantity; returns {key: float64 ndarray}."""
    import pyneb as pn
    pn.log_.level = 1
    res = {}

    line_waves = {'O3': 5007, 'N2': 6584, 'S2': 6731}
    for name in ['O3', 'N2', 'S2']:
        atom = pn.Atom(name[:-1], int(name[-1]))
        wave = line_waves[name]
        res[f'{name}_pop_scalar'] = atom.getPopulations(1e4, 1e2)
        res[f'{name}_pop_prod'] = atom.getPopulations(TEM_GRID, DEN_GRID, product=True)
        res[f'{name}_pop_noprod'] = atom.getPopulations(TEM_PAIRS, DEN_PAIRS, product=False)
        res[f'{name}_pop_noprod_2d'] = atom.getPopulations(TEM_PAIRS.reshape(5, 1),
                                                           DEN_PAIRS.reshape(5, 1),
                                                           product=False)
        res[f'{name}_collrates_scalar'] = atom.getCollRates(1e4)
        res[f'{name}_collrates_arr'] = atom.getCollRates(np.array([5e3, 1e4, 2e4]))
        res[f'{name}_critden'] = atom.getCritDensity(1e4)
        res[f'{name}_critden_arr'] = atom.getCritDensity(np.array([1e4, 1.5e4]))
        res[f'{name}_em_all_scalar'] = atom.getEmissivity(1e4, 1e2)
        res[f'{name}_em_all_grid'] = atom.getEmissivity(TEM_GRID[:2], DEN_GRID[:2])
        res[f'{name}_em_line'] = atom.getEmissivity(1e4, 1e2, wave=wave)
        res[f'{name}_em_noprod'] = atom.getEmissivity(TEM_PAIRS, DEN_PAIRS,
                                                      wave=wave, product=False)

    # Large matrix (144 levels): scalar calls only to keep the file small
    Fe3 = pn.Atom('Fe', 3)
    res['Fe3_pop_scalar'] = Fe3.getPopulations(1e4, 1e2)
    res['Fe3_collrates_scalar'] = Fe3.getCollRates(1e4)
    res['Fe3_em_all_scalar'] = Fe3.getEmissivity(1e4, 1e2)

    O3 = pn.Atom('O', 3)
    S2 = pn.Atom('S', 2)
    res['O3_temden_T'] = O3.getTemDen(0.00415, den=1e2,
                                      to_eval='L(4363)/(L(5007)+L(4959))')
    res['O3_temden_T_arr'] = O3.getTemDen(np.array([0.003, 0.004, 0.005]), den=1e2,
                                          to_eval='L(4363)/(L(5007)+L(4959))')
    res['S2_temden_N'] = S2.getTemDen(1.18, tem=1e4, wave1=6716, wave2=6731)

    H1 = pn.RecAtom('H', 1)
    He1 = pn.RecAtom('He', 1)
    He2 = pn.RecAtom('He', 2)
    res['H1_rec_4_2'] = H1.getEmissivity(1e4, 1e2, label='4_2')
    res['H1_rec_11_2'] = H1.getEmissivity(1e4, 1e2, label='11_2')
    res['H1_rec_grid'] = H1.getEmissivity(TEM_GRID, DEN_GRID, label='4_2', product=True)
    res['H1_rec_noprod'] = H1.getEmissivity(TEM_PAIRS, DEN_PAIRS, label='4_2',
                                            product=False)
    res['He1_rec_5876'] = He1.getEmissivity(1e4, 1e2, label='5876.0')
    res['He2_rec_4_3'] = He2.getEmissivity(1e4, 1e2, label='4_3')

    C = pn.Continuum()
    res['cont_erc_H'] = C.make_cont_Ercolano(1e4, 'H', WL_CONT)
    res['cont_erc_He1'] = C.make_cont_Ercolano(1e4, 'He1', WL_CONT)
    res['cont_erc_He2'] = C.make_cont_Ercolano(1e4, 'He2', WL_CONT)
    res['cont_erc_H_scalar_wl'] = C.make_cont_Ercolano(1e4, 'H', 4000.)
    res['cont_total_scalar'] = C.get_continuum(1e4, 1e2, He1_H=0.08, He2_H=0.02,
                                               wl=WL_CONT)
    res['cont_total_temarr'] = C.get_continuum(np.array([8e3, 1e4, 1.5e4]), 1e2,
                                               He1_H=0.08, He2_H=0.02, wl=WL_CONT)
    res['cont_2q'] = C.two_photon(1e4, 1e2, WL_CONT)
    res['cont_FF'] = C.FreeFree(1e4, WL_CONT, He1_H=0.08, He2_H=0.02)
    res['cont_BJ'] = C.BJ_HI(1e4, 1e2, 0.08, 0.02)

    for law in pn.RedCorr().getLaws():
        if law in ('K76', 'F99-like'):
            # These laws read attributes from inside X() (K76: E_BV and cHbeta;
            # F99-like: user-supplied FitzParams), so they cannot be set up in a
            # single constructor call; initialize under CCM89, then switch law.
            RC = pn.RedCorr(E_BV=1.2, R_V=3.1, law='CCM89')
            if law == 'F99-like':
                # standard F99 values for the FM90 parameters
                c2 = -0.824 + 4.717 / 3.1
                RC.FitzParams = [4.596, 0.99, 2.030 - 3.007 * c2, c2, 3.23, 0.41]
            RC.law = law
        else:
            RC = pn.RedCorr(E_BV=1.2, R_V=3.1, law=law)
        res[f'rc_X_{law}'] = RC.X(WAVES_RC)
        res[f'rc_corr_{law}'] = RC.getCorr(WAVES_RC)
    RC = pn.RedCorr(E_BV=1.2, R_V=3.1, law='F99')
    res['rc_chbeta'] = RC.cHbetaFromEbv(1.2)
    res['rc_ebv'] = RC.EbvFromCHbeta(0.5)

    return {k: np.asarray(v, dtype=np.float64) for k, v in res.items()}


def main():
    import pyneb as pn
    data = compute_all()
    files = get_datafiles()
    for atom, ff in files.items():
        data[f'__datafiles_{atom}'] = np.array(ff)
    data['__pyneb_version'] = np.array(pn.__version__)
    os.makedirs(os.path.dirname(GOLDEN_PATH), exist_ok=True)
    np.savez_compressed(GOLDEN_PATH, **data)
    print(f'Wrote {len(data)} arrays to {GOLDEN_PATH}')


if __name__ == '__main__':
    main()
