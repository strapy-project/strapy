import unittest
import scampy as ts
import numpy as np


class TestPolariser(unittest.TestCase):
    def test_S(self):
        """Test that 45° polarised light passing through an S polarising linear
        polariser results in pure S polarised light.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))
        model.add_component(ts.components.Dump, 'outDump', 'n3')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 'sOut', ('n2', 'n3'))

        model.add_detector('out', 'n3', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = np.sqrt(2)
        model.components['laser'].amplitude[1] = np.sqrt(2)

        model.components['inPol'].rotation = 0
        model.components['inPol'].extinction = 0
        model.components['inPol'].loss = 0
        model.components['inPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0], 0)
        self.assertAlmostEqual(
            model.detectors['out'].amplitudes[1],
            np.sqrt(2))

    def test_P(self):
        """Test that 45° polarised light passing through an P polarising linear
        polariser results in pure P polarised light.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))
        model.add_component(ts.components.Dump, 'outDump', 'n3')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 'sOut', ('n2', 'n3'))

        model.add_detector('out', 'n3', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = np.sqrt(2)
        model.components['laser'].amplitude[1] = np.sqrt(2)

        model.components['inPol'].rotation = np.pi / 2
        model.components['inPol'].extinction = 0
        model.components['inPol'].loss = 0
        model.components['inPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(
            model.detectors['out'].amplitudes[0],
            np.sqrt(2))
        self.assertAlmostEqual(model.detectors['out'].amplitudes[1], 0)

    def test_45(self):
        """Test that vertically (S) polarised light passing through an 45°
        polarising linear polariser results in pure 45° polarised light.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))
        model.add_component(ts.components.Dump, 'outDump', 'n3')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 'sOut', ('n2', 'n3'))

        model.add_detector('out', 'n3', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.components['inPol'].rotation = np.pi / 4
        model.components['inPol'].extinction = 0
        model.components['inPol'].loss = 0
        model.components['inPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0],
                               model.detectors['out'].amplitudes[1])

    def test_loss(self):
        """Test that vertically (S) polarised light passing through an S
        polarising linear polariser with 10% intensity loss results in the
        output intensity.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))
        model.add_component(ts.components.Dump, 'outDump', 'n3')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 'sOut', ('n2', 'n3'))

        model.add_detector('out', 'n3', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.components['inPol'].rotation = 0
        model.components['inPol'].extinction = 0
        model.components['inPol'].loss = 0.1
        model.components['inPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].intensity, 0.9)

    def test_extinction(self):
        """Test that 45° polarised light passing through an S polarising
        linear polariser with an extinction ratio of 0.01 results in the
        correct output intensities for each polarisation.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))
        model.add_component(ts.components.Dump, 'outDump', 'n3')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 'sOut', ('n2', 'n3'))

        model.add_detector('out', 'n3', ('amplitude',
                                         'S intensity', 'P intensity'))

        model.components['laser'].amplitude[0] = 1
        model.components['laser'].amplitude[1] = 1

        model.components['inPol'].rotation = 0
        model.components['inPol'].extinction = 0.01
        model.components['inPol'].loss = 0
        model.components['inPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].P_intensity
                               / model.detectors['out'].S_intensity, 0.01)

    def test_S_reversed(self):
        """Test that 45° polarised light passing 'backwards' through an S
        polarising linear polariser results in pure S polarised light.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Dump, 'outDump', 'n0')
        model.add_component(ts.components.Source, 'laser', 'n3')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))

        model.add_component(ts.components.Stack, 'sOut', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 'sIn', ('n2', 'n3'))

        model.add_detector('out', 'n0', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = np.sqrt(2)
        model.components['laser'].amplitude[1] = np.sqrt(2)

        model.components['inPol'].rotation = 0
        model.components['inPol'].extinction = 0
        model.components['inPol'].loss = 0
        model.components['inPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0], 0)
        self.assertAlmostEqual(
            model.detectors['out'].amplitudes[1],
            np.sqrt(2))


if __name__ == '__main__':
    unittest.main()
