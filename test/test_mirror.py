import unittest
import strapy as ts
import numpy as np


class TestMirror(unittest.TestCase):
    def test_ideal(self):
        """Test that 45° polarised light reflected off of an ideal mirror
        results in the same amplitude of 45° polarised light being returned.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Mirror, 'mirror', 'n1')

        model.add_component(ts.components.Stack, 'stack', ('n0', 'n1'))

        model.add_detector('out', 'n0', ('amplitude'))

        model.components['laser'].amplitude[0] = 1 / np.sqrt(2)
        model.components['laser'].amplitude[1] = 1 / np.sqrt(2)

        model.components['mirror'].rP = 1
        model.components['mirror'].rS = 1

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[2],
                               -1 / np.sqrt(2))
        self.assertAlmostEqual(model.detectors['out'].amplitudes[3],
                               -1 / np.sqrt(2))

    def test_90_perc(self):
        """Test that S polarised light reflected off of a mirror with 90%
        reflectivity results in 90% of the input intensity being returned.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Mirror, 'mirror', 'n1')

        model.add_component(ts.components.Stack, 'stack', ('n0', 'n1'))

        model.add_detector('out', 'n0', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.components['mirror'].rP = np.sqrt(0.9)
        model.components['mirror'].rS = np.sqrt(0.9)

        model.build()
        model.evaluate()

        # 1 + 0.9 as intensity detector monitors light propagating in both
        # directions.
        self.assertAlmostEqual(model.detectors['out'].intensity,
                               1 + 0.9)


if __name__ == '__main__':
    unittest.main()
