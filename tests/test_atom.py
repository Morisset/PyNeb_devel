"""
Comprehensive tests for pyneb.Atom (collisionally-excited ions).

Covers: constructor, getA, getEnergy, getOmega, getCollRates,
        getPopulations, getEmissivity, getCritDensity, getTemDen,
        getTransition, and smoke tests for printIonic.
"""
import numpy as np
import pytest
import pyneb as pn


# ---------------------------------------------------------------------------
# A — Constructor
# ---------------------------------------------------------------------------

class TestAtomConstructor:
    def test_O3_basic_attributes(self, O3):
        assert O3.atom == 'O3'
        assert O3.elem == 'O'
        assert O3.spec == 3
        assert O3.name == 'oxygen'
        assert O3.Z == 8
        assert O3.is_valid is True

    def test_N2_basic_attributes(self, N2):
        assert N2.atom == 'N2'
        assert N2.elem == 'N'
        assert N2.spec == 2
        assert N2.name == 'nitrogen'
        assert N2.Z == 7

    def test_atom_keyword(self):
        a = pn.Atom(atom='O3')
        assert a.atom == 'O3'
        assert a.elem == 'O'
        assert a.spec == 3

    def test_spec_as_string(self):
        a = pn.Atom('N', '2')
        assert isinstance(a.spec, int)
        assert a.spec == 2
        assert a.atom == 'N2'

    def test_NLevels_is_integer(self, O3):
        assert isinstance(O3.NLevels, (int, np.integer))

    def test_NLevels_constrained(self):
        a = pn.Atom('O', 3, NLevels=3)
        assert a.NLevels == 3

    def test_energy_eV_ground_state_sentinel(self, O3):
        assert O3.energy_eV[0, 0] == np.inf

    def test_wave_Ang_shape(self, O3):
        assert O3.wave_Ang.shape == (O3.atomNLevels, O3.atomNLevels)

    def test_lineList_is_array(self, O3):
        assert isinstance(O3.lineList, np.ndarray)
        assert len(O3.lineList) > 0


# ---------------------------------------------------------------------------
# B — Einstein A coefficients
# ---------------------------------------------------------------------------

class TestGetA:
    def test_full_matrix_shape(self, O3):
        A = O3.getA()
        assert A.shape == (O3.atomNLevels, O3.atomNLevels)

    def test_upper_triangle_zero(self, O3):
        A = O3.getA()
        # No spontaneous emission upward: upper triangle must be zero
        assert np.all(np.triu(A, k=1) == 0)

    def test_lower_triangle_has_values(self, O3):
        A = O3.getA()
        assert np.any(np.tril(A, k=-1) > 0)

    def test_A_5007(self, O3):
        # [O III] 5007 Å: level 4 → 3
        assert O3.getA(4, 3) == pytest.approx(0.02046, rel=1e-3)

    def test_A_4959(self, O3):
        # [O III] 4959 Å: level 4 → 2
        assert O3.getA(4, 2) == pytest.approx(0.006791, rel=1e-3)

    def test_A_by_wave_matches_levels(self, O3):
        assert O3.getA(wave=5007) == pytest.approx(O3.getA(4, 3), rel=1e-6)

    def test_A_reversed_order_is_zero(self, O3):
        # Upward direction: no spontaneous emission
        assert O3.getA(2, 4) == 0.0


# ---------------------------------------------------------------------------
# C — Energy levels
# ---------------------------------------------------------------------------

class TestGetEnergy:
    def test_ground_level_zero_eV(self, O3):
        assert O3.getEnergy(1, unit='eV') == pytest.approx(0.0, abs=1e-10)

    def test_levels_monotonically_ordered(self, O3):
        energies = [O3.getEnergy(i, unit='eV') for i in range(1, O3.atomNLevels + 1)]
        assert all(e2 > e1 for e1, e2 in zip(energies, energies[1:]))

    def test_full_array_shape(self, O3):
        e = O3.getEnergy()
        assert e.ndim == 1
        assert len(e) == O3.atomNLevels


# ---------------------------------------------------------------------------
# D — Collision strengths
# ---------------------------------------------------------------------------

class TestGetOmegaGetCollRates:
    def test_getOmega_scalar(self, O3):
        omega = O3.getOmega(1e4, 4, 3)
        assert omega == pytest.approx(1.262, rel=0.02)
        assert omega > 0

    def test_getOmega_array_tem_shape(self, O3):
        tems = np.array([8e3, 1e4, 1.2e4])
        result = O3.getOmega(tems, 4, 3)
        assert result.shape == (3,)
        assert np.all(result > 0)

    def test_getCollRates_scalar_shape(self, O3):
        cr = O3.getCollRates(1e4)
        assert cr.shape == (O3.NLevels, O3.NLevels)

    def test_getCollRates_all_non_negative(self, O3):
        cr = O3.getCollRates(1e4)
        assert np.all(cr >= 0)

    def test_getCollRates_array_tem_shape(self, O3):
        tems = np.array([8e3, 1e4])
        cr = O3.getCollRates(tems)
        assert cr.shape == (O3.NLevels, O3.NLevels, 2)


# ---------------------------------------------------------------------------
# E — Level populations
# ---------------------------------------------------------------------------

class TestGetPopulations:
    def test_scalar_normalization(self, O3):
        pop = O3.getPopulations(1e4, 1e2)
        assert pop.sum() == pytest.approx(1.0, abs=1e-10)

    def test_scalar_shape(self, O3):
        pop = O3.getPopulations(1e4, 1e2)
        assert pop.shape == (O3.NLevels,)

    def test_scalar_non_negative(self, O3):
        pop = O3.getPopulations(1e4, 1e2)
        assert np.all(pop >= 0)

    def test_ground_level_dominant(self, O3):
        pop = O3.getPopulations(1e4, 1e2)
        assert pop[0] > pop[1]

    def test_array_tem_shape(self, O3):
        tems = np.array([8e3, 1e4, 1.5e4])
        pop = O3.getPopulations(tems, 1e2)
        assert pop.shape == (O3.NLevels, 3)

    def test_array_tem_each_column_sums_to_one(self, O3):
        tems = np.array([8e3, 1e4, 1.5e4])
        pop = O3.getPopulations(tems, 1e2)
        assert np.allclose(pop.sum(axis=0), 1.0)

    def test_array_tem_den_product_shape(self, O3):
        tems = np.array([8e3, 1e4, 1.5e4])
        dens = np.array([1e2, 1e3])
        pop = O3.getPopulations(tems, dens, product=True)
        assert pop.shape == (O3.NLevels, 3, 2)

    def test_product_false_shape(self, O3):
        tems = np.array([8e3, 1e4, 1.5e4])
        dens = np.array([1e2, 5e2, 1e3])
        pop = O3.getPopulations(tems, dens, product=False)
        assert pop.shape == (O3.NLevels, 3)

    def test_product_false_normalization(self, O3):
        tems = np.array([8e3, 1e4, 1.5e4])
        dens = np.array([1e2, 5e2, 1e3])
        pop = O3.getPopulations(tems, dens, product=False)
        assert np.allclose(pop.sum(axis=0), 1.0)


# ---------------------------------------------------------------------------
# F — Emissivity
# ---------------------------------------------------------------------------

class TestGetEmissivity:
    def test_O3_5007_value(self, O3):
        em = O3.getEmissivity(1e4, 1e2, wave=5007)
        assert em == pytest.approx(3.497e-21, rel=5e-3)

    def test_O3_4959_value(self, O3):
        em = O3.getEmissivity(1e4, 1e2, wave=4959)
        assert em == pytest.approx(1.172e-21, rel=5e-3)

    def test_O3_4363_value(self, O3):
        em = O3.getEmissivity(1e4, 1e2, wave=4363)
        assert em == pytest.approx(2.279e-23, rel=5e-3)

    def test_5007_4959_ratio(self, O3):
        ratio = (O3.getEmissivity(1e4, 1e2, wave=5007) /
                 O3.getEmissivity(1e4, 1e2, wave=4959))
        assert ratio == pytest.approx(2.984, rel=1e-3)

    def test_wave_vs_levels_agreement(self, O3):
        em_wave = O3.getEmissivity(1e4, 1e2, wave=5007)
        em_lev = O3.getEmissivity(1e4, 1e2, lev_i=4, lev_j=3)
        assert em_wave == pytest.approx(em_lev, rel=1e-6)

    def test_scalar_tem_scalar_den_is_scalar(self, O3):
        em = O3.getEmissivity(1e4, 1e2, wave=5007)
        assert np.ndim(em) == 0 or (np.ndim(em) == 1 and len(em) == 1)

    def test_array_tem_scalar_den_shape(self, O3):
        tems = np.array([8e3, 1e4, 1.5e4])
        em = O3.getEmissivity(tems, 1e2, wave=5007)
        assert em.shape == (3,)
        assert np.all(em > 0)

    def test_scalar_tem_array_den_shape(self, O3):
        dens = np.array([1e2, 1e3, 1e4])
        em = O3.getEmissivity(1e4, dens, wave=5007)
        assert em.shape == (3,)
        assert np.all(em > 0)

    def test_array_tem_array_den_product_true_shape(self, O3):
        tems = np.array([8e3, 1e4, 1.5e4])
        dens = np.array([1e2, 1e3, 1e4])
        em = O3.getEmissivity(tems, dens, wave=5007)
        assert em.shape == (3, 3)

    def test_array_tem_array_den_product_false_shape(self, O3):
        tems = np.array([8e3, 1e4, 1.5e4])
        dens = np.array([1e2, 1e3, 1e4])
        em = O3.getEmissivity(tems, dens, wave=5007, product=False)
        assert em.shape == (3,)
        assert np.all(em > 0)

    def test_4363_increases_with_temperature(self, O3):
        tems = np.array([5e3, 1e4, 2e4])
        em = O3.getEmissivity(tems, 1e2, wave=4363)
        assert np.all(np.diff(em) > 0)

    def test_reversed_levels_is_zero(self, O3):
        em = O3.getEmissivity(1e4, 1e2, lev_i=2, lev_j=4)
        assert em == 0.0

    def test_full_matrix_shape(self, O3):
        em = O3.getEmissivity(1e4, 1e2)
        assert em.shape == (O3.NLevels, O3.NLevels)


# ---------------------------------------------------------------------------
# G — Critical density
# ---------------------------------------------------------------------------

class TestGetCritDensity:
    def test_O3_level4_value(self, O3):
        ncrit = float(np.atleast_1d(O3.getCritDensity(1e4, level=4))[0])
        assert ncrit == pytest.approx(6.9e5, rel=0.05)

    def test_O2_level2_lower_than_O3(self, O3, O2):
        ncrit_O3 = float(np.atleast_1d(O3.getCritDensity(1e4, level=4))[0])
        ncrit_O2 = float(np.atleast_1d(O2.getCritDensity(1e4, level=2))[0])
        assert ncrit_O2 < ncrit_O3

    def test_O2_level2_value(self, O2):
        ncrit = float(np.atleast_1d(O2.getCritDensity(1e4, level=2))[0])
        assert ncrit == pytest.approx(1207., rel=0.05)


# ---------------------------------------------------------------------------
# H — Temperature / density diagnostics (getTemDen)
# ---------------------------------------------------------------------------

class TestGetTemDen:
    def test_N2_temperature_roundtrip(self, N2):
        ratio = (N2.getEmissivity(1e4, 1e2, wave=5755) /
                 N2.getEmissivity(1e4, 1e2, wave=6584))
        tem = N2.getTemDen(ratio, den=1e2, wave1=5755, wave2=6584)
        assert float(tem) == pytest.approx(1e4, rel=1e-3)

    def test_S2_density_roundtrip(self, S2):
        ratio = (S2.getEmissivity(1e4, 300., wave=6731) /
                 S2.getEmissivity(1e4, 300., wave=6716))
        den = S2.getTemDen(ratio, tem=1e4, wave1=6731, wave2=6716)
        assert float(den) == pytest.approx(300., rel=5e-3)

    def test_to_eval_matches_wave_syntax(self, N2):
        t1 = N2.getTemDen(0.01, den=1e3, wave1=5755, wave2=6584)
        t2 = N2.getTemDen(0.01, den=1e3, to_eval='L(5755)/L(6584)')
        assert float(t1) == pytest.approx(float(t2), rel=0.01)

    def test_level_syntax_matches_wave_syntax(self, N2):
        t1 = N2.getTemDen(0.01, den=1e3, wave1=5755, wave2=6584)
        t2 = N2.getTemDen(0.01, den=1e3, lev_i1=5, lev_j1=4, lev_i2=4, lev_j2=3)
        assert float(t1) == pytest.approx(float(t2), rel=0.01)

    def test_array_ratios_shape(self, N2):
        ratios = np.array([0.014, 0.020, 0.030])
        tems = N2.getTemDen(ratios, den=1e2, wave1=5755, wave2=6584)
        assert tems.shape == (3,)

    def test_higher_ratio_gives_higher_temp(self, N2):
        ratios = np.array([0.014, 0.020, 0.030])
        tems = N2.getTemDen(ratios, den=1e2, wave1=5755, wave2=6584)
        assert np.all(np.diff(tems) > 0)


# ---------------------------------------------------------------------------
# I — Transition lookup
# ---------------------------------------------------------------------------

class TestGetTransition:
    def test_O3_5007(self, O3):
        assert O3.getTransition(5007) == (4, 3)

    def test_O3_4959(self, O3):
        assert O3.getTransition(4959) == (4, 2)

    def test_O3_4363(self, O3):
        assert O3.getTransition(4363) == (5, 4)


# ---------------------------------------------------------------------------
# J — Smoke tests (no exception expected)
# ---------------------------------------------------------------------------

class TestPrintIonic:
    def test_printIonic_no_args(self, O3, capsys):
        O3.printIonic()
        captured = capsys.readouterr()
        assert 'O' in captured.out

    def test_printIonic_with_conditions(self, O3, capsys):
        O3.printIonic(tem=1e4, den=1e2)
        captured = capsys.readouterr()
        assert len(captured.out) > 0
