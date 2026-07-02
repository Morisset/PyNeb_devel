"""
Tests for pn.EmissionLine and pn.Observation: constructors, file reading in
the different formats, editing, extinction correction and Monte Carlo.
"""
import numpy as np
import pytest

import pyneb as pn


OBS_LINES_IN_COLS = """NAME  H1r_4861A H1r_4861Ae H1r_6563A H1r_6563Ae O3_5007A O3_5007Ae O3_4363A N2_6584A
obj1  100.0     0.05       570.0      0.05       500.0    0.05      5.0      120.0
obj2  100.0     0.05       300.0      0.05       430.0    0.05      4.0      80.0
"""

OBS_LABELS = {'H1r_4861A', 'H1r_6563A', 'O3_5007A', 'O3_4363A', 'N2_6584A'}


def write_obs_file(tmp_path, content=OBS_LINES_IN_COLS, name='obs.dat'):
    path = tmp_path / name
    path.write_text(content)
    return str(path)


@pytest.fixture
def obs2(tmp_path):
    """A fresh two-object Observation read from a lines_in_cols file."""
    return pn.Observation(write_obs_file(tmp_path))


class TestEmissionLine:

    def test_constructor_elem_spec_wave(self):
        line = pn.EmissionLine('O', 3, 5007, obsIntens=[1.4, 1.3])
        assert line.label == 'O3_5007A'
        assert line.atom == 'O3'
        assert line.wave == 5007
        assert line.is_valid is True
        assert line.to_eval == 'L(5007)'
        # not corrected yet
        assert np.all(line.corrIntens == 0.)

    def test_constructor_label_corrected(self):
        line = pn.EmissionLine(label='O3_5007A', obsIntens=[320.], corrected=True)
        assert line.corrIntens == pytest.approx(line.obsIntens)

    def test_error_not_relative(self):
        line = pn.EmissionLine(label='O3_5007A', obsIntens=[200.], obsError=[10.],
                               errIsRelative=False)
        assert line.obsError[0] == pytest.approx(0.05)

    def test_blend_label(self):
        line = pn.EmissionLine(label='O2_3727A+', obsIntens=[1.])
        assert line.blend is True
        assert line.is_valid is True
        # to_eval either refers to the blend label or comes from BLEND_LIST
        assert line.to_eval in ('S("3727A+")', pn.BLEND_LIST.get('O2_3727A+'))

    def test_invalid_line(self):
        line = pn.EmissionLine(label='O3_9999A', obsIntens=[1.])
        assert line.is_valid is False
        assert line.to_eval is None

    def test_correctIntens(self):
        line = pn.EmissionLine(label='O3_5007A', obsIntens=[100.])
        RC = pn.RedCorr(law='CCM89', E_BV=1.)
        line.correctIntens(RC)
        assert line.corrIntens[0] > line.obsIntens[0]

    def test_correctIntens_no_extinction(self):
        line = pn.EmissionLine(label='O3_5007A', obsIntens=[100.])
        RC = pn.RedCorr(law='CCM89', E_BV=0.)
        line.correctIntens(RC)
        assert line.corrIntens[0] == pytest.approx(100.)

    def test_addObs(self):
        line = pn.EmissionLine(label='O3_5007A', obsIntens=[1., 2.])
        line.addObs(3., 0.1)
        assert line.obsIntens.size == 3
        assert line.obsIntens[-1] == pytest.approx(3.)
        assert line.obsError[-1] == pytest.approx(0.1)

    def test_printLine(self, capsys):
        line = pn.EmissionLine(label='O3_5007A', obsIntens=[100.])
        line.printLine()
        captured = capsys.readouterr()
        assert 'O3_5007A' in captured.out
        assert 'evaluated as' in captured.out

    def test_repr(self):
        line = pn.EmissionLine(label='O3_5007A', obsIntens=[100.])
        assert 'O3_5007A' in repr(line)


class TestObservationBasic:

    def test_counts(self, obs2):
        assert obs2.n_lines == 5
        assert obs2.n_valid_lines == 5
        assert obs2.n_obs == 2
        assert obs2.names == ['obj1', 'obj2']

    def test_lineLabels(self, obs2):
        assert isinstance(obs2.lineLabels, np.ndarray)
        assert set(obs2.lineLabels) == OBS_LABELS

    def test_getLine_by_label(self, obs2):
        line = obs2.getLine(label='O3_5007A')
        assert line.obsIntens == pytest.approx([500., 430.])

    def test_getLine_by_elem_spec_wave(self, obs2):
        assert obs2.getLine('O', 3, 5007) is obs2.getLine(label='O3_5007A')

    def test_getLine_unknown(self, obs2):
        assert obs2.getLine(label='Ar4_4740A') is None

    def test_errors_from_file(self, obs2):
        assert obs2.getLine(label='O3_5007A').obsError == pytest.approx([0.05, 0.05])
        # N2_6584A has no error column: err_default = 0.10 applies
        assert obs2.getLine(label='N2_6584A').obsError == pytest.approx([0.10, 0.10])

    def test_getUniqueAtoms(self, obs2):
        assert set(obs2.getUniqueAtoms()) == {'H1r', 'N2', 'O3'}

    def test_getSortedLines_wave(self, obs2):
        waves = [l.wave for l in obs2.getSortedLines(crit='wave')]
        assert waves == sorted(waves)

    def test_getSortedLines_atom(self, obs2):
        atoms = [l.atom for l in obs2.getSortedLines(crit='atom')]
        assert atoms == sorted(atoms)

    def test_getIntens(self, obs2):
        intens = obs2.getIntens(returnObs=True)
        assert intens['O3_5007A'] == pytest.approx([500., 430.])
        assert obs2.getIntens(returnObs=True, obsName='obj2')['O3_5007A'] == pytest.approx(430.)

    def test_getError(self, obs2):
        err = obs2.getError(returnObs=True)
        assert err['O3_5007A'] == pytest.approx([0.05, 0.05])


class TestObservationEditing:

    def test_addLine(self, obs2):
        obs2.addLine(pn.EmissionLine(label='S2_6716A', obsIntens=[10., 20.]))
        assert obs2.n_lines == 6
        assert obs2.getLine(label='S2_6716A').obsIntens == pytest.approx([10., 20.])

    def test_addLine_invalid_object(self, obs2):
        with pytest.raises(pn.utils.logging.PyNebError):
            obs2.addLine('not a line')
        assert obs2.n_lines == 5

    def test_removeLine(self, obs2):
        obs2.removeLine('O3_4363A')
        assert obs2.n_lines == 4
        assert 'O3_4363A' not in obs2.lineLabels

    def test_fillObs(self, obs2):
        obs2.fillObs('Ne3_3869A')
        line = obs2.getLine(label='Ne3_3869A')
        assert line.obsIntens.size == 2
        assert np.all(np.isnan(line.obsIntens))

    def test_addObs(self, obs2):
        obs2.addObs('obj3', [110., 400., 450., 4.5, 100.])
        assert obs2.n_obs == 3
        assert obs2.names == ['obj1', 'obj2', 'obj3']

    def test_addObs_wrong_length(self, obs2):
        with pytest.raises(pn.utils.logging.PyNebError):
            obs2.addObs('obj3', [1., 2.])
        assert obs2.n_obs == 2

    def test_addObs_duplicate_name(self, obs2):
        with pytest.raises(pn.utils.logging.PyNebError):
            obs2.addObs('obj1', [110., 400., 450., 4.5, 100.])
        assert obs2.n_obs == 2

    def test_addSum(self, tmp_path):
        content = """NAME  O2_3726A O2_3729A
obj1  100.0    140.0
"""
        obs = pn.Observation(write_obs_file(tmp_path, content))
        obs.addSum(('O2_3726A', 'O2_3729A'), 'O2_3727A+', to_eval='S("3727+")')
        line = obs.getLine(label='O2_3727A+')
        assert line.obsIntens == pytest.approx([240.])


class TestObservationExtinction:

    def test_def_EBV_correctData(self, obs2):
        obs2.def_EBV()  # default: H1r_6563A / H1r_4861A = 2.85
        obs2.correctData()
        intens = obs2.getIntens()
        ratio = intens['H1r_6563A'] / intens['H1r_4861A']
        assert ratio == pytest.approx([2.85, 2.85], rel=1e-3)
        # both objects have observed ratios above 2.85, so E(B-V) > 0
        assert np.all(obs2.extinction.E_BV > 0)

    def test_corrected_larger_than_observed(self, obs2):
        obs2.def_EBV()
        obs2.correctData()
        # positive extinction: corrected intensity of a blue line increases
        corr = obs2.getIntens()['O3_4363A']
        obs = obs2.getIntens(returnObs=True)['O3_4363A']
        assert np.all(corr > obs)

    def test_corrected_constructor(self, tmp_path):
        obs = pn.Observation(write_obs_file(tmp_path), corrected=True)
        assert obs.extinction.law == 'No correction'
        for label in obs.lineLabels:
            assert obs.getIntens()[label] == pytest.approx(obs.getIntens(returnObs=True)[label])


class TestObservationReadFormats:

    def test_lines_in_rows(self, smc24_path):
        obs = pn.Observation(smc24_path, fileFormat='lines_in_rows', err_default=0.05)
        assert obs.n_obs == 1
        assert obs.names == ['SMC_24']
        assert obs.extinction.cHbeta[0] == pytest.approx(0.047)
        assert 'O3_5007A' in obs.lineLabels
        assert obs.getIntens(returnObs=True)['O3_5007A'][0] == pytest.approx(435.09)

    def test_lines_in_rows_correct(self, smc24_path):
        obs = pn.Observation(smc24_path, fileFormat='lines_in_rows', err_default=0.05)
        obs.correctData()
        assert obs.getIntens()['O3_5007A'][0] > obs.getIntens(returnObs=True)['O3_5007A'][0]

    def test_lines_in_cols2(self, tmp_path):
        obs = pn.Observation(write_obs_file(tmp_path), fileFormat='lines_in_cols2')
        assert obs.n_lines == 5
        assert obs.getLine(label='O3_5007A').obsIntens == pytest.approx([500., 430.])
        assert obs.getLine(label='O3_5007A').obsError == pytest.approx([0.05, 0.05])

    def test_lines_in_rows_err_cols(self, tmp_path):
        content = """LINE  IC418  errIC418
H1r_4861A 100.0 5.0
O3_5007A  200.0 10.0
N2_6584A  50.0  5.0
"""
        obs = pn.Observation(write_obs_file(tmp_path, content),
                             fileFormat='lines_in_rows_err_cols', errIsRelative=False)
        assert obs.names == ['IC418']
        assert obs.getLine(label='O3_5007A').obsIntens == pytest.approx([200.])
        assert obs.getLine(label='O3_5007A').obsError == pytest.approx([0.05])


class TestMonteCarlo:

    def test_addMonteCarloObs(self, tmp_path):
        content = """NAME  H1r_4861A H1r_4861Ae O3_5007A O3_5007Ae
obj1  100.0     0.05       500.0    0.05
"""
        obs = pn.Observation(write_obs_file(tmp_path, content))
        obs.addMonteCarloObs(N=5, random_seed=42)
        assert obs.n_obs == 6
        assert obs.MC_added is True
        assert obs.N_MC == 5
        assert obs.names[0] == 'obj1'
        assert obs.names[1] == 'obj1-MC-1'
        assert obs.n_obs_origin == 1
        # first entry of each line keeps the original value
        assert obs.getLine(label='O3_5007A').obsIntens[0] == pytest.approx(500.)
        # fake observations scatter around the original value
        assert obs.getLine(label='O3_5007A').obsIntens[1:] == pytest.approx(500., rel=0.5)

    def test_addMonteCarloObs_two_objects(self, tmp_path):
        obs = pn.Observation(write_obs_file(tmp_path))
        obs.addMonteCarloObs(N=3, random_seed=42)
        assert obs.n_obs == 8
        assert len(obs.names) == 8
        assert obs.names[0] == 'obj1'
        assert obs.names[4] == 'obj2'
        assert obs.n_obs_origin == 2


class TestSetAllErrors:

    def test_setAllErrors(self, obs2):
        obs2.setAllErrors(0.2)
        for label, err in obs2.getError(returnObs=True).items():
            assert err == pytest.approx([0.2, 0.2])
