import unittest
import scampy as ts
import numpy as np


class TestSource(unittest.TestCase):
    def test_S(self):
        """Test passing S polarised light through a stack to a dump and
        detector.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Dump, 'dump', 'n1')

        model.add_component(ts.components.Stack, 'stack', ('n0', 'n1'))

        model.add_detector('out', 'n1', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0], 0)
        self.assertAlmostEqual(model.detectors['out'].amplitudes[1], 1)

    def test_P(self):
        """Test passing P polarised light through a stack to a dump and
        detector.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Dump, 'dump', 'n1')

        model.add_component(ts.components.Stack, 'stack', ('n0', 'n1'))

        model.add_detector('out', 'n1', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 1
        model.components['laser'].amplitude[1] = 0

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0], 1)
        self.assertAlmostEqual(model.detectors['out'].amplitudes[1], 0)

    def test_RHC(self):
        """Test passing right hand circularly  polarised light through a stack
        to a dump and detector.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Dump, 'dump', 'n1')

        model.add_component(ts.components.Stack, 'stack', ('n0', 'n1'))

        model.add_detector('out', 'n1', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 1
        model.components['laser'].amplitude[1] = -1j
        model.components['laser'].amplitude = \
            model.components['laser'].amplitude / np.sqrt(2)

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0],
                               1 / np.sqrt(2))
        self.assertAlmostEqual(model.detectors['out'].amplitudes[1],
                               -1j / np.sqrt(2))
        self.assertAlmostEqual(model.detectors['out'].intensity, 1)

    def test_LHC(self):
        """Test passing left hand circularly  polarised light through a stack
        to a dump and detector.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Dump, 'dump', 'n1')

        model.add_component(ts.components.Stack, 'stack', ('n0', 'n1'))

        model.add_detector('out', 'n1', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 1
        model.components['laser'].amplitude[1] = 1j
        model.components['laser'].amplitude = \
            model.components['laser'].amplitude / np.sqrt(2)

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0],
                               1 / np.sqrt(2))
        self.assertAlmostEqual(model.detectors['out'].amplitudes[1],
                               1j / np.sqrt(2))
        self.assertAlmostEqual(model.detectors['out'].intensity, 1)


if __name__ == '__main__':
    unittest.main()
