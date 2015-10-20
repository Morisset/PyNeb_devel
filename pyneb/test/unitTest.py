import unittest
import pyneb as pn
import numpy as np

def rel_error(calculated, expected):
	return (abs((expected - calculated) * 1.0) / (expected * 1.0))

class UnitTest(unittest.TestCase):

	def getTransitions (self):
		atom = pn.Atom("O", 3)
		tran01 = atom.getTransition(4363)
		tran02 = atom.getTransition(4959)
		tran03 = atom.getTransition(5010)
		print tran01, tran02, tran03

	def testInit(self):
		atomO2 = pn.Atom("O", 2)
		assert (atomO2.elem == "O")
		assert (atomO2.spec == 2)
		assert (atomO2.atom == "O2")

	def testEmissivity(self, error=5e-2):
		atom = pn.Atom("O", 3)

		result = atom.getEmissivity(20000, 1e2)
		emi = [9.749E-22, 4.107E-28, 8.131E-22, 1.012E-24, 4.165E-21, 1.204E-20, 1.087E-22, 
			2.906E-25, 3.979E-22, 3.850E-22, 9.559E-22, 5.517E-27]
		a = []
		for i in result:
			for j in i:
				if (j > 0) :
					a.append(j)

		for aux in range(len(a) - 1):
			result = (a[aux] / emi[aux])
			assert (rel_error(a[aux], emi[aux]) < error)
	
	def testEmissivityShape(self):
		atom = pn.Atom("O", 3)
		result = atom.getEmissivity(20000, 1e2)
		assert (result.shape == (6, 6))
		result = atom.getEmissivity(20000, [1e2, 1e3])
		assert (result.shape == (6, 6, 2))
		result = atom.getEmissivity([20000, 10000], 1e2)
		assert (result.shape == (6, 6, 2))
		result = atom.getEmissivity([20000, 10000], [1e2, 1e3])		
		assert (result.shape == (6, 6, 2, 2))
		result = atom.getEmissivity([20000, 10000], [1e2, 1e3], product = True)		
		assert (result.shape == (6, 6, 2, 2))
		result = atom.getEmissivity([20000, 10000], [1e2, 1e3], product = False)		
		assert (result is None)
		
		result = atom.getEmissivity(20000, 1e2, wave = 5007)
		assert (result.shape == ())
		result = atom.getEmissivity(20000, [1e2, 1e3], wave = 5007)
		assert (result.shape == (2,))
		result = atom.getEmissivity([20000, 10000], 1e2, wave = 5007)
		assert (result.shape == (2,)) # the bug I'm looking for
		result = atom.getEmissivity([20000, 10000], [1e2, 1e3], wave = 5007)		
		assert (result.shape == (2, 2))
		result = atom.getEmissivity([20000, 10000], [1e2, 1e3], product = True, wave = 5007)		
		assert (result.shape == (2, 2))
		

	def testHbEmissivity(self, max_error=5e-4):
		H1 = pn.RecAtom('H', 1)
		tem = np.array([.5, 1, 1.5, 2., 3., 5.]) * 1e4
		den = 1e2
		Hb_expected = np.array([2.255e-25, 1.258e-25, 8.724e-26, 6.683e-26, 4.560e-26, 2.796e-26])
		Hb_computed = H1.getEmissivity(tem, den, label='4_2')
		for i in range(len(tem)):
			assert(rel_error(H1.getEmissivity(tem[i], den, label='4_2'), Hb_expected[i]) < max_error)
			assert(rel_error((Hb_computed[i]), Hb_expected[i]) < max_error)
		
	def testIonAbundance(self, max_error=5e-3):
		atom = pn.Atom("O", 3)
		assert(rel_error(atom.getIonAbundance(10000, 10000, 100, 4, 3), 3.562e-3) < max_error)
		assert(rel_error(atom.getIonAbundance(10000, 10000, 100, 4, 3), 3.562e-3) < max_error)
		assert(rel_error(atom.getIonAbundance(10000, 20000, 100, 4, 3), 5.552e-4) < max_error)
		atom = pn.Atom("S", 2)
		assert(rel_error(atom.getIonAbundance(200, 20000, 100, 3, 1), 2.202e-6) < max_error)
		atom = pn.Atom("O", 2)
		assert(rel_error(atom.getIonAbundance(200, 20000, 1000, 3, 1), 1.608e-5) < max_error)
		atom = pn.Atom("O", 3)
		assert(rel_error(atom.getIonAbundance(10000, 10000, 100, wave=5007), 3.562e-3) < max_error)
		assert(rel_error(atom.getIonAbundance(10000, 10000, 100, -5, 1, 5007), 3.562e-3) < max_error)
		assert(rel_error(atom.getIonAbundance(10000, 20000, 100, wave=5007), 5.552e-4) < max_error)

	def testGetA(self, max_error=5e-3):
		atom = pn.Atom("O", 3)
		assert(rel_error(atom.getA(1, 1) + 1, 1.0) < max_error)
		assert(rel_error(atom.getA(2, 1), 2.664e-05) < max_error)
		assert(rel_error(atom.getA(3, 1), 3.094e-11) < max_error)
		assert(rel_error(atom.getA(4, 1), 1.69e-06) < max_error)
		assert(rel_error(atom.getA(3, 2), 9.695e-05) < max_error)
		assert(rel_error(atom.getA(4, 2), 0.006995) < max_error)
		assert(rel_error(atom.getA(5, 2), 0.2268) < max_error)

	def testGetOmega(self, max_error=5e-3):
		atom = pn.Atom("O", 3)
		assert(rel_error(atom.getOmega(10000, 2, 1), 0.543264060967) < max_error)
		assert(rel_error(atom.getOmega(10000, 3, 1), 0.27035364411) < max_error)
		assert(rel_error(atom.getOmega(10000, 4, 1), 0.25398979945) < max_error)
		assert(rel_error(atom.getOmega(10000, 5, 1), 0.0325890388813) < max_error)
		assert(rel_error(atom.getOmega(10000, 3, 2), 1.28634364203) < max_error)
		assert(rel_error(atom.getOmega(10000, 4, 2), 0.761941062781) < max_error)
		assert(rel_error(atom.getOmega(10000, 5, 2), 0.0977670995235) < max_error)
		assert(rel_error(atom.getOmega(10000, wave=5007), 1.2698350039427169) < max_error)
		
	def testGetTemDen(self, max_error=5e-3):
		atom = pn.Atom('O', 3)
		assert(rel_error(atom.getTemDen(100., den=10000., wave1=5007, wave2=4363), 11444.2) < max_error)
		assert(rel_error(atom.getTemDen(1., tem=10000., wave1=88e4, wave2=52e4), 174.379) < max_error)
		atom = pn.Atom('O', 2)
		assert(rel_error(atom.getTemDen(1.5, tem=10000., wave1=3726, wave2=3729), 1095.59) < max_error)


if __name__ == '__main__':

	unittest.main()
