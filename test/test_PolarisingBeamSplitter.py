import unittest
import scampy as ts
import numpy as np


class TestPolarisingBeamSplitter(unittest.TestCase):
    def test_45(self):
        """Test that 45° polarised light passing through an unrotated
        polarising beamsplitter produces S and P polarised beams.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.PolarisingBeamSplitter, 'pbs',
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

        model.components['laser'].amplitude[0] = np.sqrt(2)
        model.components['laser'].amplitude[1] = np.sqrt(2)

        model.components['pbs'].rExtinction = 0
        model.components['pbs'].tExtinction = 0
        model.components['pbs'].rLoss = 0
        model.components['pbs'].tLoss = 0
        model.components['pbs'].theta0 = 0
        model.components['pbs'].theta1 = 0
        model.components['pbs'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['pd1'].amplitudes[0], 0)
        self.assertAlmostEqual(
            model.detectors['pd1'].amplitudes[1],
            np.sqrt(2))
        self.assertAlmostEqual(
            model.detectors['pd2'].amplitudes[0],
            np.sqrt(2))
        self.assertAlmostEqual(model.detectors['pd2'].amplitudes[1], 0)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[0], 0)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[1], 0)

    def test_S(self):
        """Test that S polarised light passing through an polarising beam
        splitter is fully reflected.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.PolarisingBeamSplitter, 'pbs',
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

        model.components['pbs'].rExtinction = 0
        model.components['pbs'].tExtinction = 0
        model.components['pbs'].rLoss = 0
        model.components['pbs'].tLoss = 0
        model.components['pbs'].theta0 = 0
        model.components['pbs'].theta1 = 0
        model.components['pbs'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['pd1'].amplitudes[0], 0)
        self.assertAlmostEqual(model.detectors['pd1'].amplitudes[1], 1)
        self.assertAlmostEqual(model.detectors['pd2'].amplitudes[0], 0)
        self.assertAlmostEqual(model.detectors['pd2'].amplitudes[1], 0)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[0], 0)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[1], 0)

    def test_loss(self):
        """Test that 45° polarised light passing through a polarising
        beamsplitter with 10% intensity loss for the reflected beam and 20%
        loss for the transmitted beam results in the correct output intensity.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.PolarisingBeamSplitter, 'pbs',
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

        model.components['laser'].amplitude[0] = 1
        model.components['laser'].amplitude[1] = 1

        model.components['pbs'].rExtinction = 0
        model.components['pbs'].tExtinction = 0
        model.components['pbs'].sLoss = 0.9
        model.components['pbs'].pLoss = 0.8
        model.components['pbs'].theta0 = 0
        model.components['pbs'].theta1 = 0
        model.components['pbs'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['pd1'].intensity, 0.9)
        self.assertAlmostEqual(model.detectors['pd2'].intensity, 0.8)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[0], 0)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[1], 0)

    def test_extinction(self):
        """ Test that 45° polarised light passing through a polarising
        beamsplitter with a transmission extinction of 0.01 and a reflection
        extinction of 0.02 results in the correct intensity outputs.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.PolarisingBeamSplitter, 'pbs',
                            ('n1', 'n2', 'n3', 'n4'))
        model.add_component(ts.components.Dump, 'd1', 'n5')
        model.add_component(ts.components.Dump, 'd2', 'n6')
        model.add_component(ts.components.Dump, 'd3', 'n7')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
        model.add_component(ts.components.Stack, 's1', ('n2', 'n5'))
        model.add_component(ts.components.Stack, 's2', ('n3', 'n6'))
        model.add_component(ts.components.Stack, 's3', ('n4', 'n7'))

        model.add_detector('pd1', 'n5', ('amplitude', 'intensity',
                                         'S intensity', 'P intensity'))
        model.add_detector('pd2', 'n6', ('amplitude', 'intensity',
                                         'S intensity', 'P intensity'))
        model.add_detector('pd3', 'n7', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 1
        model.components['laser'].amplitude[1] = 1

        model.components['pbs'].rExtinction = 0.02
        model.components['pbs'].tExtinction = 0.01
        model.components['pbs'].rLoss = 0
        model.components['pbs'].tLoss = 0
        model.components['pbs'].theta0 = 0
        model.components['pbs'].theta1 = 0
        model.components['pbs'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['pd1'].intensity
                               + model.detectors['pd2'].intensity, 2)
        self.assertAlmostEqual(model.detectors['pd1'].P_intensity
                               / model.detectors['pd1'].S_intensity,
                               model.components['pbs'].rExtinction)
        self.assertAlmostEqual(model.detectors['pd2'].S_intensity
                               / model.detectors['pd2'].P_intensity,
                               model.components['pbs'].tExtinction)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[0], 0)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[1], 0)

    def test_theta0(self):
        """ Test that S polarised light passing through a polarising beam
        splitter rotated by 45° about the axis from the zeroth to the second
        node results in an even distribution of light in both outputs.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.PolarisingBeamSplitter, 'pbs',
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

        model.components['pbs'].rExtinction = 0
        model.components['pbs'].tExtinction = 0
        model.components['pbs'].rLoss = 0
        model.components['pbs'].tLoss = 0
        model.components['pbs'].theta0 = np.pi / 4
        model.components['pbs'].theta1 = 0
        model.components['pbs'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['pd1'].intensity, 0.5)
        self.assertAlmostEqual(model.detectors['pd2'].intensity, 0.5)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[0], 0)
        self.assertAlmostEqual(model.detectors['pd3'].amplitudes[1], 0)

    def test_theta1(self):
        """ Test that S polarised light passing through a polarising beam
        splitter (from the first node) rotated by 45° about the axis from the
        first to the third node results in an even distribution of light in
        both outputs.
        """

        model = ts.Model()
        model.wavelength = 633e-9

        model.add_component(ts.components.Source, 'laser', 'n0')
        model.add_component(ts.components.PolarisingBeamSplitter, 'pbs',
                            ('n1', 'n2', 'n3', 'n4'))
        model.add_component(ts.components.Dump, 'd1', 'n5')
        model.add_component(ts.components.Dump, 'd2', 'n6')
        model.add_component(ts.components.Dump, 'd3', 'n7')

        model.add_component(ts.components.Stack, 'sIn', ('n0', 'n2'))
        model.add_component(ts.components.Stack, 's1', ('n1', 'n5'))
        model.add_component(ts.components.Stack, 's2', ('n3', 'n6'))
        model.add_component(ts.components.Stack, 's3', ('n4', 'n7'))

        model.add_detector('pd1', 'n5', ('amplitude', 'intensity'))
        model.add_detector('pd2', 'n6', ('amplitude', 'intensity'))
        model.add_detector('pd3', 'n7', ('amplitude', 'intensity'))

        model.components['laser'].amplitude[0] = 0
        model.components['laser'].amplitude[1] = 1

        model.components['pbs'].rExtinction = 0
        model.components['pbs'].tExtinction = 0
        model.components['pbs'].rLoss = 0
        model.components['pbs'].tLoss = 0
        model.components['pbs'].theta0 = 0
        model.components['pbs'].theta1 = np.pi / 4
        model.components['pbs'].update()

        model.build()
        model.evaluate()

        self.assertAlmostEqual(model.detectors['pd1'].intensity, 0.5)
        self.assertAlmostEqual(model.detectors['pd2'].intensity, 0)
        self.assertAlmostEqual(model.detectors['pd3'].intensity, 0.5)


if __name__ == '__main__':
    unittest.main()
