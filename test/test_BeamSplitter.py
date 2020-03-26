import unittest
import strapy as ts
import numpy as np


class TestBeamSplitter(unittest.TestCase):
    def test_S(self):
        """Test that S polarised light passing through a beamsplitter produces
        two equal intensity S polarised beams.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.BeamSplitter, 'bs',
                            ('n1', 'n2', 'n3', 'n4'))
        model.add_component(ts.components.Dump, 'd1', 'n5')
        model.add_component(ts.components.Dump, 'd2', 'n6')
        model.add_component(ts.components.Dump, 'd3', 'n7')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 's1', ('n2', 'n5'))
        model.add_component(ts.components.Stack, 's2', ('n3', 'n6'))
        model.add_component(ts.components.Stack, 's3', ('n4', 'n7'))

        model.add_detector('pd1', 'n5', ('amplitude', 'intensity'))
        model.add_detector('pd2', 'n6', ('amplitude', 'intensity'))
        model.add_detector('pd3', 'n7', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['pd1'].amplitudes[0], 0)
        self.assertAlmostEqual(
            model.detectors['pd1'].amplitudes[1],
            1 / np.sqrt(2))
        self.assertAlmostEqual(model.detectors['pd2'].amplitudes[0], 0)
        self.assertAlmostEqual(
            model.detectors['pd2'].amplitudes[1],
            1 / np.sqrt(2))
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[0], 0)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[1], 0)

        self.assertAlmostEqual(model.detectors['pd1'].intensity, 0.5)
        self.assertAlmostEqual(model.detectors['pd2'].intensity, 0.5)
        self.assertAlmostEqual(model.detectors['pd3'].intensity, 0)


if __name__ == '__main__':
    unittest.main()
