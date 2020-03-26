import unittest
import strapy as ts
import numpy as np


class TestIdealIsolator(unittest.TestCase):
    def test_ideal(self):
        """Test that 45° polarised light reflected off of an ideal mirror
        through an isolator does not reach the source.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.IdealIsolator, 'isolator',
                            ('n1', 'n2'))
        model.add_component(ts.components.Mirror, 'mirror', 'n3')

        model.add_component(ts.components.Stack, 's01', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 's23', ('n2', 'n3'))

        model.add_detector('out', 'n0', ('amplitude'))

        model.components['laser'].amplitude[0] = 1 / np.sqrt(2)
        model.components['laser'].amplitude[1] = 1 / np.sqrt(2)

        model.components['mirror'].rP = 1
        model.components['mirror'].rS = 1

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0],
                               1 / np.sqrt(2))
        self.assertAlmostEqual(model.detectors['out'].amplitudes[1],
                               1 / np.sqrt(2))
        self.assertAlmostEqual(model.detectors['out'].amplitudes[2], 0)
        self.assertAlmostEqual(model.detectors['out'].amplitudes[3], 0)


if __name__ == '__main__':
    unittest.main()
