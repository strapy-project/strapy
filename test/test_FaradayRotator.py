import unittest
import scampy as ts
import numpy as np


class TestFaradayRotator(unittest.TestCase):
    def test_rot(self):
        """Test that S polarised light passing through a Faraday rotator with a
        rotation of 45° results in light polarised at 45°.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))
        model.add_component(ts.components.FaradayRotator, 'fr', ('n3', 'n4'))
        model.add_component(ts.components.Polariser, 'outPol', ('n5', 'n6'))
        model.add_component(ts.components.Dump, 'outDump', 'n7')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 's1', ('n2', 'n3'))
        model.add_component(ts.components.Stack, 's2', ('n4', 'n5'))
        model.add_component(ts.components.Stack, 'sOut', ('n6', 'n7'))

        model.add_detector('out', 'n7', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.components['fr'].rotation = np.pi / 4
        model.components['fr'].update()

        model.components['inPol'].rotation = 0
        model.components['inPol'].extinction = 0
        model.components['inPol'].loss = 0
        model.components['inPol'].update()

        model.components['outPol'].rotation = np.pi / 4
        model.components['outPol'].extinction = 0
        model.components['outPol'].loss = 0
        model.components['outPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].intensity, 1)

    def test_isolator(self):
        """Test that S polarised light passing through a Faraday isolator
        consisting of an S polarising polariser, a 45° Faraday rotator and a
        45° rotated polariser is correctly isolated.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))
        model.add_component(ts.components.FaradayRotator, 'fr', ('n3', 'n4'))
        model.add_component(ts.components.Polariser, 'outPol', ('n5', 'n6'))
        model.add_component(ts.components.Mirror, 'mirror', 'n7')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 's1', ('n2', 'n3'))
        model.add_component(ts.components.Stack, 's2', ('n4', 'n5'))
        model.add_component(ts.components.Stack, 'sOut', ('n6', 'n7'))

        model.add_detector('out', 'n0', ('amplitude', 'intensity'))
        model.add_detector('mir', 'n7', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.components['fr'].rotation = np.pi / 4
        model.components['fr'].update()

        model.components['inPol'].rotation = 0
        model.components['inPol'].extinction = 0
        model.components['inPol'].loss = 0
        model.components['inPol'].update()

        model.components['outPol'].rotation = np.pi / 4
        model.components['outPol'].extinction = 0
        model.components['outPol'].loss = 0
        model.components['outPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(np.abs(model.detectors['mir'].amplitudes[0]),
                               1 / np.sqrt(2))
        self.assertAlmostEqual(np.abs(model.detectors['mir'].amplitudes[1]),
                               1 / np.sqrt(2))
        self.assertAlmostEqual(
            np.sign(
                model.detectors['mir'].amplitudes[1] *
                model.detectors['mir'].amplitudes[0]),
            1)
        self.assertAlmostEqual(model.detectors['out'].amplitudes[2], 0)
        self.assertAlmostEqual(model.detectors['out'].amplitudes[3], 0)


if __name__ == '__main__':
    unittest.main()
