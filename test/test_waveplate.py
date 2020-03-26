import unittest
import strapy as ts
import numpy as np


class TestWaveplate(unittest.TestCase):
    def test_HWP_SP(self):
        """Test that S polarised light passing through a half wave plate at 45°
        degrees to the vertical results in P polarised light.

        Polarisers are used to filter inputs and outputs; the polariser unit
        tests should be run first to ensure any issues can reasonably be
        isolated to the waveplate.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))
        model.add_component(ts.components.Waveplate, 'hwp', ('n3', 'n4'))
        model.add_component(ts.components.Polariser, 'outPol', ('n5', 'n6'))
        model.add_component(ts.components.Dump, 'outDump', 'n7')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 's1', ('n2', 'n3'))
        model.add_component(ts.components.Stack, 's2', ('n4', 'n5'))
        model.add_component(ts.components.Stack, 'sOut', ('n6', 'n7'))

        model.add_detector('out', 'n7', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.components['hwp'].rotation = np.pi / 4
        model.components['hwp'].retardance = np.pi
        model.components['hwp'].update()

        model.components['inPol'].rotation = 0
        model.components['inPol'].extinction = 0
        model.components['inPol'].loss = 0
        model.components['inPol'].update()

        model.components['outPol'].rotation = np.pi / 2
        model.components['outPol'].extinction = 0
        model.components['outPol'].loss = 0
        model.components['outPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].intensity, 1)

    def test_HWP_PS(self):
        """Test that P polarised light passing through a half wave plate at 45°
        degrees to the vertical results in S polarised light.

        Polarisers are used to filter inputs and outputs; the polariser unit
        tests should be run first to ensure any issues can reasonably be
        isolated to the waveplate.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))
        model.add_component(ts.components.Waveplate, 'hwp', ('n3', 'n4'))
        model.add_component(ts.components.Polariser, 'outPol', ('n5', 'n6'))
        model.add_component(ts.components.Dump, 'outDump', 'n7')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 's1', ('n2', 'n3'))
        model.add_component(ts.components.Stack, 's2', ('n4', 'n5'))
        model.add_component(ts.components.Stack, 'sOut', ('n6', 'n7'))

        model.add_detector('out', 'n7', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 1
        model.components['laser'].amplitude[1] = 0

        model.components['hwp'].rotation = np.pi / 4
        model.components['hwp'].retardance = np.pi
        model.components['hwp'].update()

        model.components['inPol'].rotation = np.pi / 2
        model.components['inPol'].extinction = 0
        model.components['inPol'].loss = 0
        model.components['inPol'].update()

        model.components['outPol'].rotation = 0
        model.components['outPol'].extinction = 0
        model.components['outPol'].loss = 0
        model.components['outPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].intensity, 1)

    def test_HWP_4545(self):
        """Test that 45° polarised light passing through a half wave plate at
        0° degrees to the vertical results in -45° polarised light.

        Polarisers are used to filter inputs and outputs; the polariser unit
        tests should be run first to ensure any issues can reasonably be
        isolated to the waveplate.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Polariser, 'inPol', ('n1', 'n2'))
        model.add_component(ts.components.Waveplate, 'hwp', ('n3', 'n4'))
        model.add_component(ts.components.Polariser, 'outPol', ('n5', 'n6'))
        model.add_component(ts.components.Dump, 'outDump', 'n7')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 's1', ('n2', 'n3'))
        model.add_component(ts.components.Stack, 's2', ('n4', 'n5'))
        model.add_component(ts.components.Stack, 'sOut', ('n6', 'n7'))

        model.add_detector('out', 'n7', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 1 / np.sqrt(2)
        model.components['laser'].amplitude[1] = 1 / np.sqrt(2)

        model.components['hwp'].rotation = 0
        model.components['hwp'].retardance = np.pi
        model.components['hwp'].update()

        model.components['inPol'].rotation = np.pi / 4
        model.components['inPol'].extinction = 0
        model.components['inPol'].loss = 0
        model.components['inPol'].update()

        model.components['outPol'].rotation = -np.pi / 4
        model.components['outPol'].extinction = 0
        model.components['outPol'].loss = 0
        model.components['outPol'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].intensity, 1)

    def test_QWP(self):
        """Test that S polarised light passing through a quarter wave plate at
        45° degrees to the vertical, reflecting off a mirror then passing back
        through the same waveplate results in P polarised light.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Waveplate, 'qwp', ('n1', 'n2'))
        model.add_component(ts.components.Mirror, 'outMirror', 'n3')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 's1', ('n2', 'n3'))

        model.add_detector('out', 'n0', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.components['outMirror'].rP = 1
        model.components['outMirror'].rS = 1

        model.components['qwp'].rotation = np.pi / 4
        model.components['qwp'].retardance = np.pi / 2
        model.components['qwp'].update()

        model.build()
        model.evaluate()

        back_intensity_P = np.abs(model.detectors['out'].amplitudes[2])**2
        back_intensity_S = np.abs(model.detectors['out'].amplitudes[3])**2

        self.assertAlmostEqual(back_intensity_P, 1)
        self.assertAlmostEqual(back_intensity_S, 0)

    def test_QWP_SS(self):
        """ Test that S polarised light passing through a quarter wave plate at
        0° degrees to the vertical, reflecting off a mirror then passing back
        through the same waveplate results in S polarised light.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Waveplate, 'qwp', ('n1', 'n2'))
        model.add_component(ts.components.Mirror, 'outMirror', 'n3')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 's1', ('n2', 'n3'))

        model.add_detector('out', 'n0', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.components['qwp'].rotation = 0
        model.components['qwp'].retardance = np.pi / 2
        model.components['qwp'].update()

        model.build()
        model.evaluate()

        back_intensity_P = np.abs(model.detectors['out'].amplitudes[2])**2
        back_intensity_S = np.abs(model.detectors['out'].amplitudes[3])**2

        self.assertAlmostEqual(back_intensity_P, 0)
        self.assertAlmostEqual(back_intensity_S, 1)


if __name__ == '__main__':
    unittest.main()
