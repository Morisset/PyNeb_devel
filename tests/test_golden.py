"""
Golden-value regression tests.

The reference values stored in golden_values.json were computed once with a
given PyNeb version and the default atomic data set, and are compared to the
current results every time the suite runs. A failure means the numerical
output of the library changed: either an unintended regression, or a
deliberate change (new atomic data, algorithm update), in which case the
reference file must be regenerated and committed:

    python tests/test_golden.py --update

The stored values are plain floats in a two-level dictionary
(section -> label -> value), so a git diff of golden_values.json shows
exactly which quantities moved.
"""
import json
import os

import numpy as np
import pytest

import pyneb as pn
from pyneb.core.pynebcore import getHbEmissivity

GOLDEN_FILE = os.path.join(os.path.dirname(__file__), 'golden_values.json')

# relative tolerance per section; anything not listed uses DEFAULT_RTOL.
# Iterative solvers (getTemDen, getCrossTemDen) get a looser tolerance.
DEFAULT_RTOL = 1e-6
TOLERANCES = {'temden': 1e-4,
              'diagnostics': 1e-3}

SECTIONS = ('atom_emissivity',
            'atom_populations',
            'atom_A',
            'atom_critical_densities',
            'temden',
            'rec_emissivity',
            'redcorr',
            'continuum',
            'icf',
            'diagnostics')


def compute_golden_values():
    """
    Compute the full dictionary of golden quantities from the current library.
    Every value must be deterministic (no random numbers) and stable across
    platforms and the Python/numpy versions used in CI.
    """
    values = {section: {} for section in SECTIONS}
    tem, den = 1e4, 1e2

    O3 = pn.Atom('O', 3)
    N2 = pn.Atom('N', 2)
    S2 = pn.Atom('S', 2)
    O2 = pn.Atom('O', 2)
    H1 = pn.RecAtom('H', 1)
    He2 = pn.RecAtom('He', 2)

    # --- collisional line emissivities (erg.s-1.cm3) at Te=1e4 K, Ne=1e2 cm-3
    emis = values['atom_emissivity']
    for atom, waves in ((O3, (5007, 4959, 4363)),
                        (N2, (6584, 5755)),
                        (S2, (6716, 6731)),
                        (O2, (3726, 3729))):
        for wave in waves:
            emis['{0}_{1}'.format(atom.atom, wave)] = float(atom.getEmissivity(tem, den, wave=wave))

    # --- O III level populations at Te=1e4 K, Ne=1e2 cm-3 (first 5 levels)
    pops = O3.getPopulations(tem, den)
    values['atom_populations'] = {'O3_level{0}'.format(i + 1): float(pops[i]) for i in range(5)}

    # --- transition probabilities
    values['atom_A'] = {'O3_A_4_2': float(O3.getA(4, 2)),
                        'O3_A_4_3': float(O3.getA(4, 3)),
                        'O3_A_5_4': float(O3.getA(5, 4)),
                        'N2_A_4_3': float(N2.getA(4, 3)),
                        'S2_A_2_1': float(S2.getA(2, 1))}

    # --- critical densities at Te=1e4 K (levels 2-5)
    crit = values['atom_critical_densities']
    for atom in (O3, S2):
        cd = atom.getCritDensity(tem)
        for lev in (2, 3, 4, 5):
            crit['{0}_level{1}'.format(atom.atom, lev)] = float(cd[lev - 1])

    # --- temperature / density diagnostics from fixed line ratios
    values['temden'] = {
        'O3_tem_from_4363_5007': float(O3.getTemDen(3.0e-3, den=1e2,
                                                    to_eval='L(4363)/L(5007)')),
        'N2_tem_from_5755_6584': float(N2.getTemDen(1.0e-2, den=1e2,
                                                    to_eval='L(5755)/L(6584)')),
        'S2_den_from_6716_6731': float(S2.getTemDen(1.0, tem=1e4,
                                                    to_eval='L(6716)/L(6731)'))}

    # --- recombination emissivities
    values['rec_emissivity'] = {
        'H1_4_2': float(H1.getEmissivity(tem, den, label='4_2')),
        'H1_3_2': float(H1.getEmissivity(tem, den, label='3_2')),
        'H1_6_2': float(H1.getEmissivity(tem, den, label='6_2')),
        'He2_4686': float(He2.getEmissivity(tem, den, label='4_3')),
        'Hbeta_fit_1e4': float(getHbEmissivity(tem=1e4)),
        'Hbeta_fit_15e3': float(getHbEmissivity(tem=1.5e4))}

    # --- extinction corrections at E(B-V)=0.5
    rc = values['redcorr']
    for law in ('CCM89', 'F99'):
        RC = pn.RedCorr(law=law, E_BV=0.5)
        for wave in (3727., 5007., 6563.):
            rc['{0}_{1:.0f}'.format(law, wave)] = float(RC.getCorr(wave))

    # --- nebular continuum
    C = pn.Continuum()
    cont = C.get_continuum(tem=tem, den=den, wl=np.array([3600., 3700., 4000.]),
                           HI_label=None)
    tp = C.two_photon(tem, den, np.array([2000., 4000.]))
    values['continuum'] = {
        'cont_3600': float(cont[0]),
        'cont_3700': float(cont[1]),
        'cont_4000': float(cont[2]),
        'BJ_HI': float(C.BJ_HI(tem=tem, den=den, He1_H=0.1, He2_H=0.01)),
        'two_photon_2000': float(tp[0]),
        'two_photon_4000': float(tp[1]),
        'gff_4000': float(C.gff(1., tem, np.array([4000., 6000.]))[0])}

    # --- ionization correction factors from a fixed abundance set
    icf = pn.ICF()
    atom_abun = {'O2': 1e-5, 'O3': 2e-5, 'N2': 1e-6, 'Ne3': 1.2e-5}
    res = icf.getElemAbundance(atom_abun,
                               icf_list=['direct_O.23', 'TPP77_14', 'TPP77_15'])
    values['icf'] = {label: float(res[label]) for label in res}

    # --- crossed diagnostics: Te and Ne from two fixed line ratios
    diags = pn.Diagnostics()
    diags.addDiag(['[OIII] 4363/5007', '[SII] 6731/6716'])
    tem_x, den_x = diags.getCrossTemDen('[OIII] 4363/5007', '[SII] 6731/6716',
                                        3.0e-3, 1.05)
    values['diagnostics'] = {'crossTemDen_O3_S2_tem': float(tem_x),
                             'crossTemDen_O3_S2_den': float(den_x)}

    return values


@pytest.fixture(scope='module')
def computed():
    return compute_golden_values()


@pytest.fixture(scope='module')
def golden():
    assert os.path.exists(GOLDEN_FILE), \
        'golden_values.json is missing; generate it with: python tests/test_golden.py --update'
    with open(GOLDEN_FILE) as f:
        return json.load(f)


@pytest.mark.parametrize('section', SECTIONS)
def test_golden_values(section, computed, golden):
    assert section in golden, \
        'section {0} missing from golden_values.json; regenerate it'.format(section)
    rtol = TOLERANCES.get(section, DEFAULT_RTOL)
    assert set(computed[section]) == set(golden[section]), \
        'quantities of section {0} changed; regenerate golden_values.json'.format(section)
    for label, value in computed[section].items():
        # abs=0 disables pytest.approx's default absolute tolerance (1e-12),
        # which would swamp the relative check for values like emissivities (~1e-21)
        assert value == pytest.approx(golden[section][label], rel=rtol, abs=0), \
            '{0}/{1}: computed {2!r} vs golden {3!r}'.format(section, label, value,
                                                             golden[section][label])


if __name__ == '__main__':
    import sys
    if '--update' in sys.argv:
        new_values = compute_golden_values()
        with open(GOLDEN_FILE, 'w') as f:
            json.dump(new_values, f, indent=2, sort_keys=True)
            f.write('\n')
        n = sum(len(v) for v in new_values.values())
        print('Wrote {0} golden values to {1}'.format(n, GOLDEN_FILE))
    else:
        print(__doc__)
        print('Run with --update to (re)generate {0}'.format(GOLDEN_FILE))
