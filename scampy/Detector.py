import numpy as np


class Detector:
    """A detector for logging electric field amplitude and related properties.

    An instance of the `Detector` class is created for each detector added to a
    model. The detector will log at least the electric field amplitude at the
    specified node, and optionally calculate other properties. Detectors can
    currently monitor:

            * `amplitude` - the amplitudes of the forward and backward
                    propagating S and P polarised electric fields.
            * `intensity` - the total intensity of light passing through the
                    detector, including both polarisations propagating in both
                    directions.
            * `S intensity` - the intensity of S polarised light passing
                    through the detector in both directions.
            * `P intensity` - the intensity of P polarised light passing
                    through the detector in both directions.

    Note that calculated intensity is proportional to both the amplitude and
    the dielectric properties of the medium; the intensity for a detector is
    only directly comparable to other intensities calculated in the same
    medium.

    Attributes
    ----------
    name : str
            Unique name of the detector.
    node : str
            Node to be monitored by detector.
    properties : tuple of str
            Optical properties to be logged - see description for current
            options.
    node_index : int
            Index of the node the detector is to monitor in the solution
            vector.
    """

    def __init__(self, name, node, properties):
        self.name = name
        self.node = node
        self.properties = properties
        self.node_index = None

        if 'amplitude' in properties:
            self.AMP = True
        else:
            self.AMP = False

        if 'intensity' in properties:
            self.INT = True
        else:
            self.INT = False

        if 'S intensity' in properties:
            self.S_INT = True
        else:
            self.S_INT = False

        if 'P intensity' in properties:
            self.P_INT = True
        else:
            self.P_INT = False

    def update(self, solution_vector):
        """Update the detected values from solution vector.

        Parameters
        ----------
        solution_vector : ndarray
                solution to network matrix equation."""

        self.amplitudes = \
            solution_vector[self.node_index:self.node_index + 4].flatten()

        if self.INT:
            self.intensity = np.sum(np.abs(self.amplitudes)**2)
        if self.S_INT:
            self.S_intensity = np.abs(self.amplitudes[1])**2 \
                + np.abs(self.amplitudes[3])**2
        if self.P_INT:
            self.P_intensity = np.abs(self.amplitudes[0])**2 \
                + np.abs(self.amplitudes[2])**2
