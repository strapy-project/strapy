import unittest
import strapy as ts
import numpy as np
import pyctmm


class TestStack(unittest.TestCase):
    def test_set_length(self):
        """Test that the set_length method produces the correct phase changes.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Dump, 'dump', 'n1')

        model.add_component(ts.components.Stack, 'stack', ('n0', 'n1'))

        model.add_detector('out', 'n1', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 1
        model.components['laser'].amplitude[1] = 1

        model.components['stack'].set_length(0)

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0], 1)
        self.assertAlmostEqual(model.detectors['out'].amplitudes[1], 1)

        model.components['stack'].set_length(0.25)

        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0],
                               np.exp(1j * 2 * np.pi * 0.25))
        self.assertAlmostEqual(model.detectors['out'].amplitudes[1],
                               np.exp(1j * 2 * np.pi * 0.25))

        model.components['stack'].set_length(0.634)

        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0],
                               np.exp(1j * 2 * np.pi * 0.634))
        self.assertAlmostEqual(model.detectors['out'].amplitudes[1],
                               np.exp(1j * 2 * np.pi * 0.634))

    def test_set_pyctmm_free_space(self):
        """Test that the set_ctmm method produces the correct results.

        This does not aim to fully test the pyctmm library, but tests a few
        simple cases:
            phase change through free space
            reflection from air-glass interface
            destructive interference at air-glass interface (half-wave layer)
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Dump, 'dump', 'n1')

        model.add_component(ts.components.Stack, 'stack', ('n0', 'n1'))

        model.add_detector('out', 'n1', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        cstack = pyctmm.create_stack(1, model.wavelength, 0)

        pyctmm.set_ind(cstack, 0, 1, 0)
        pyctmm.set_d(cstack, 0, model.wavelength / 4)

        model.components['stack'].set_pyctmm(cstack)

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['out'].amplitudes[0], 0)
        self.assertAlmostEqual(model.detectors['out'].amplitudes[1],
                               np.exp(-1j * np.pi / 2))
        self.assertAlmostEqual(model.detectors['out'].amplitudes[2], 0)
        self.assertAlmostEqual(model.detectors['out'].amplitudes[3], 0)

        del cstack
        del model

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Dump, 'dump', 'n1')

        model.add_component(ts.components.Stack, 'stack', ('n0', 'n1'))

        model.add_detector('in', 'n0', ('amplitude', 'intensity'))
        model.add_detector('out', 'n1', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        cstack = pyctmm.create_stack(2, model.wavelength, 0)

        pyctmm.set_ind(cstack, 0, 1, 0)
        pyctmm.set_ind(cstack, 1, 1.5, 0)
        pyctmm.set_d(cstack, 0, 0)
        pyctmm.set_d(cstack, 1, 0)

        model.components['stack'].set_pyctmm(cstack)

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['in'].intensity - 1,
                               ((1 - 1.5) / (1 + 1.5))**2)

        del cstack
        del model

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.Dump, 'dump', 'n1')

        model.add_component(ts.components.Stack, 'stack', ('n0', 'n1'))

        model.add_detector('in', 'n0', ('amplitude', 'intensity'))
        model.add_detector('out', 'n1', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        cstack = pyctmm.create_stack(3, model.wavelength, 0)

        pyctmm.set_ind(cstack, 0, 1, 0)
        pyctmm.set_ind(cstack, 1, 1.515, 0)
        pyctmm.set_ind(cstack, 2, 1, 0)
        pyctmm.set_d(cstack, 0, 0)
        pyctmm.set_d(cstack, 1, (model.wavelength / 1.515) / 2)
        pyctmm.set_d(cstack, 2, 0)

        model.components['stack'].set_pyctmm(cstack)

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['in'].intensity - 1,
                               0)
        self.assertAlmostEqual(model.detectors['out'].intensity,
                               1)


if __name__ == '__main__':
    unittest.main()
