"""The components module holds all optical components currently defined in
scampy. Each scampy optical component is defined in a separate class, and must
inherit at least from the `_Component` class, and optionally from the
`_ScatterComponent` or `_TransferComponent` class.

When defining new components, if not using the `_ScatterComponent` or
`_TransferComponent` classes, care must be taken to ensure the correct
stack attachment logic is employed in the `_TransferComponent.initEquation`
method. Adding new components that do not inherit from the `_ScatterComponent`
or `_TransferComponent` is therefore discouraged, to avoid cluttering the stack
attachment logic.

Coupling of counter propagating electric field components is resolved as
follows:

    * Sources always emit **into** stacks; sources set the 'a' components of
        their single node.
    * Beamsplitters always take inputs **from** stacks; the input to any port
        of the beamsplitter is the 'a' component at the given node.
    * General multi-port devices will always take inputs **from** stacks; the
        input to any port is the 'a' component at the given node.
    * Dumps always block the component going **into** the stack; the 'b'
        component is set to 0.

All direction logic therefore takes place in the stack equation initialisation.
This ensures there is no change to the physics depending on how components are
ordered. This model allows light to propigate in both directions
simultaneously, so there should be no concept of 'forwards'.

Eventually this should be automated by checking all nodes connect to at least
one stack and inserting a identity matrix connection stack if this is not the
case. In the long term placing individual components in their own file is
probably a better approach as well.
"""


import sympy as sp
import numpy as np
import pyctmm


def rotationMatrix44(theta):
    """Returns rotation matrix for vector ordered (a0P, a0S, a1P, a1S).

    The calculated rotation matrix rotates both the forward and backward
    propagating polarisation vectors at a node.

    Parameters
    ----------
    theta : double
            Rotation angle, measured clockwise when looking from the first to
            the second node.

    Returns
    -------
    rMat : numpy.ndarray
            Rotation matrix.
    """
    rMat = np.zeros((4, 4), dtype=np.float64)

    rMat[0][0] = np.cos(theta)
    rMat[0][1] = -np.sin(theta)
    rMat[1][0] = np.sin(theta)
    rMat[1][1] = np.cos(theta)

    rMat[2][2] = np.cos(theta)
    rMat[2][3] = -np.sin(theta)
    rMat[3][2] = np.sin(theta)
    rMat[3][3] = np.cos(theta)

    return rMat


def rotationMatrix88(theta0, theta1):
    """Calculates rotation matrix for vector ordered (a0P, a0S, a1P, a1S, a2P,
    a2S, a3P, a3S).

    Intended for coordinate rotations at four port devices (for example
    beam splitters). Note that this does not rotate the optical component, it
    just allows the S and P polarised states of the component to be misaligned
    relative to the S and P polarised states of the rest of the model.

    Parameters
    ----------
    theta0 : double
            Rotation angle about zeroth to second node axis, measured clockwise
            when looking from the zeroth node to the second.
    theta1 : double
            Rotation angle about first to third node axis, measured clockwise
            looking from the first node to the third.

    Returns
    -------
    rMat : numpy.ndarray
            Rotation matrix.
    """
    rMat = np.zeros((8, 8), dtype=np.float64)

    rMat[0][0] = np.cos(theta0)
    rMat[0][1] = -np.sin(theta0)
    rMat[1][0] = np.sin(theta0)
    rMat[1][1] = np.cos(theta0)

    rMat[2][2] = np.cos(theta1)
    rMat[2][3] = -np.sin(theta1)
    rMat[3][2] = np.sin(theta1)
    rMat[3][3] = np.cos(theta1)

    rMat[4][4] = np.cos(theta0)
    rMat[4][5] = -np.sin(theta0)
    rMat[5][4] = np.sin(theta0)
    rMat[5][5] = np.cos(theta0)

    rMat[6][6] = np.cos(theta1)
    rMat[6][7] = -np.sin(theta1)
    rMat[7][6] = np.sin(theta1)
    rMat[7][7] = np.cos(theta1)

    return rMat


class _Component():
    """General component class for inheritance of common properties.

    Not for external use.
    """

    def __init__(self, name, nodes, model):
        """ constructor """
        self.name = name
        self.nodes = nodes
        self.equation = sp.Eq(0, 0, evaluate=False)
        self.symbols = ()
        self.model = model
        self.symbolIdxs = []


class _ScatterComponent():
    """General class for scattering matrix based components.

    Inherits from `_Component`.
    """

    def __init__(self, name, nodes, model, node_number):
        _Component.__init__(self, name, nodes, model)
        self.node_number = node_number


class Source(_Component):
    """Light source.

    Nodes: 1

    Attributes
    ----------
    name : str
            Unique name of the component.
    nodes : str
            Node to which the component is attached.
    model : scampy.Model()
            Model in which component is to be included.
    amplitude : list of complex
            Jones amplitude vector of the emitted light. Defaults to S
            polarised light. Should be set by user for other polarisations.
    node_number : int
            Number of nodes component attaches to. Should not be changed.
    """

    def __init__(self, name, nodes, model):
        _Component.__init__(self, name, nodes, model)
        self.amplitude = [0, 1]  # defaults to S polarised light
        self.node_number = 1

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(2, 1)
        rhs = sp.zeros(2, 1)

        lhs[0] = nodes[self.nodes[0]].symbols[0]
        lhs[1] = nodes[self.nodes[0]].symbols[1]

        rhs[0] = sp.symbols(self.name + '_AP', real=False)
        rhs[1] = sp.symbols(self.name + '_AS', real=False)

        self.symbols = (rhs[0], rhs[1])

        self.equation = sp.Eq(lhs, rhs)

    def setVals(self):
        """Returns numerical values needed to solve the network matrix.

        Should not need to be called by the user.
        """
        return self.amplitude


class BeamSplitter(_ScatterComponent):
    """Symmetrical beam splitter component with predefined scattering matrix.

    Transmission and reflection properties are identical for all incident ports
    of the beam splitter. Can be used for non-polarising or polarising beam
    splitters by setting the amplitude coefficients for each polarisation.

    Nodes: 4

    Attributes
    ----------
    name : str
            Unique name of the component.
    nodes : str
            Nodes to which the component is attached.
    model : scampy.Model()
            Model in which component is to be included.
    rP : complex
            Ampltiude reflectivity coefficient for P polarised light.
    rS : complex
            Ampltiude reflectivity coefficient for S polarised light.
    tP : complex
            Ampltiude transmission coefficient for P polarised light.
    tS : complex
            Ampltiude transmission coefficient for S polarised light.
    """

    def __init__(self, name, nodes, model):
        _ScatterComponent.__init__(self, name, nodes, model, 4)

        self.rP = np.sqrt(0.5)
        self.rS = np.sqrt(0.5)
        self.tP = np.sqrt(0.5)
        self.tS = np.sqrt(0.5)

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(8, 1)
        rhs = sp.zeros(8, 1)
        scatteringMatrix = sp.zeros(8, 8)

        rhs[0] = nodes[self.nodes[0]].symbols[2]
        rhs[1] = nodes[self.nodes[0]].symbols[3]
        rhs[2] = nodes[self.nodes[1]].symbols[2]
        rhs[3] = nodes[self.nodes[1]].symbols[3]
        rhs[4] = nodes[self.nodes[2]].symbols[2]
        rhs[5] = nodes[self.nodes[2]].symbols[3]
        rhs[6] = nodes[self.nodes[3]].symbols[2]
        rhs[7] = nodes[self.nodes[3]].symbols[3]

        lhs[0] = nodes[self.nodes[0]].symbols[0]
        lhs[1] = nodes[self.nodes[0]].symbols[1]
        lhs[2] = nodes[self.nodes[1]].symbols[0]
        lhs[3] = nodes[self.nodes[1]].symbols[1]
        lhs[4] = nodes[self.nodes[2]].symbols[0]
        lhs[5] = nodes[self.nodes[2]].symbols[1]
        lhs[6] = nodes[self.nodes[3]].symbols[0]
        lhs[7] = nodes[self.nodes[3]].symbols[1]

        scatteringMatrix[0, 2] = sp.symbols(self.name + '_rP')
        scatteringMatrix[0, 4] = sp.symbols(self.name + '_tP')

        scatteringMatrix[1, 3] = sp.symbols(self.name + '_rS')
        scatteringMatrix[1, 5] = sp.symbols(self.name + '_tS')

        scatteringMatrix[2, 0] = sp.symbols(self.name + '_rP')
        scatteringMatrix[2, 6] = sp.symbols(self.name + '_tP')

        scatteringMatrix[3, 1] = sp.symbols(self.name + '_rS')
        scatteringMatrix[3, 7] = sp.symbols(self.name + '_tS')

        scatteringMatrix[4, 0] = sp.symbols(self.name + '_tP')
        scatteringMatrix[4, 6] = sp.symbols(self.name + '_rP')

        scatteringMatrix[5, 1] = sp.symbols(self.name + '_tS')
        scatteringMatrix[5, 7] = sp.symbols(self.name + '_rS')

        scatteringMatrix[6, 2] = sp.symbols(self.name + '_tP')
        scatteringMatrix[6, 4] = sp.symbols(self.name + '_rP')

        scatteringMatrix[7, 3] = sp.symbols(self.name + '_tS')
        scatteringMatrix[7, 5] = sp.symbols(self.name + '_rS')

        self.symbols = (sp.symbols(self.name + '_rP'),
                        sp.symbols(self.name + '_rS'),
                        sp.symbols(self.name + '_tP'),
                        sp.symbols(self.name + '_tS'))

        self.equation = sp.Eq(scatteringMatrix * lhs, rhs)

    def setVals(self):
        """Returns numerical values needed to solve the network matrix.

        Should not need to be called by the user.
        """
        return [self.rP, self.rS, self.tP, self.tS]


class PolarisingBeamSplitter(_ScatterComponent):
    """Symmetrical polarising beam splitter component.

    A polarising beam splitter that reflects S polarised light and transmits P
    polarised light. Transmission and reflection properties are identical for
    all incident ports of the beam splitter.

    Nodes: 4

    Attributes
    ----------
    name : str
            Unique name of the component.
    nodes : str
            Nodes to which the component is attached.
    model : scampy.Model()
            Model in which component is to be included.
    rExtinction : double
            Polarisation extinction coefficient for the reflected beam. Defined
            as the intensity ratio of reflected P polarised light to reflected
            S polarised light.
    tExtinction : double
            Polarisation extinction coefficient for the transmitted beam.
            Defined as the intensity ratio of transmitted S polarised light to
            transmitted P polarised light.
    sLoss : double
            Defined as the sum of the intensities of the relfection and
            transmitted S polarised components.
    pLoss : double
            Defined as the sum of the intensities of the relfection and
            transmitted P polarised components.
    theta0 : double
            Rotation angle about the node 0 to node 2 axis, measured clockwise
            from the S polarised axis, looking from node 0 to 2.
    theta1 : double
            Rotation angle about the node 1 to node 3 axis, measured clockwise
            from the S polarised axis, looking from node 1 to 3."""

    def __init__(self, name, nodes, model):
        _ScatterComponent.__init__(self, name, nodes, model, 4)

        self.numeric_matrix = np.zeros((8, 8), dtype=np.complex)

        self.rExtinction = 0
        self.tExtinction = 0
        self.sLoss = 1
        self.pLoss = 1

        self.theta0 = 0
        self.theta1 = 0

        self.rP = self.rExtinction \
            * (self.sLoss - self.pLoss * self.tExtinction) \
            / (1 - self.tExtinction * self.rExtinction)
        self.rS = (self.sLoss - self.pLoss * self.tExtinction) \
            / (1 - self.tExtinction * self.rExtinction)
        self.tS = (self.pLoss * self.tExtinction
                   - self.sLoss * self.tExtinction * self.rExtinction) \
            / (1 - self.tExtinction * self.rExtinction)
        self.tP = (self.pLoss - self.sLoss * self.rExtinction) \
            / (1 - self.tExtinction * self.rExtinction)

        self.rP = np.sqrt(self.rP)
        self.rS = np.sqrt(self.rS)
        self.tP = np.sqrt(self.tP)
        self.tS = np.sqrt(self.tS)

        self.numeric_matrix[0][2] = self.rP
        self.numeric_matrix[0][4] = self.tP
        self.numeric_matrix[1][3] = self.rS
        self.numeric_matrix[1][5] = self.tS
        self.numeric_matrix[2][0] = self.rP
        self.numeric_matrix[2][6] = self.tP
        self.numeric_matrix[3][1] = self.rS
        self.numeric_matrix[3][7] = self.tS
        self.numeric_matrix[4][0] = self.tP
        self.numeric_matrix[4][6] = self.rP
        self.numeric_matrix[5][1] = self.tS
        self.numeric_matrix[5][7] = self.rS
        self.numeric_matrix[6][2] = self.tP
        self.numeric_matrix[6][4] = self.rP
        self.numeric_matrix[7][3] = self.tS
        self.numeric_matrix[7][5] = self.rS

        self.numeric_matrix = rotationMatrix88(-self.theta0, -self.theta1) \
            @ self.numeric_matrix @ rotationMatrix88(self.theta0, self.theta1)

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(8, 1)
        rhs = sp.zeros(8, 1)
        scatteringMatrix = sp.zeros(8, 8)

        rhs[0] = nodes[self.nodes[0]].symbols[2]
        rhs[1] = nodes[self.nodes[0]].symbols[3]
        rhs[2] = nodes[self.nodes[1]].symbols[2]
        rhs[3] = nodes[self.nodes[1]].symbols[3]
        rhs[4] = nodes[self.nodes[2]].symbols[2]
        rhs[5] = nodes[self.nodes[2]].symbols[3]
        rhs[6] = nodes[self.nodes[3]].symbols[2]
        rhs[7] = nodes[self.nodes[3]].symbols[3]

        lhs[0] = nodes[self.nodes[0]].symbols[0]
        lhs[1] = nodes[self.nodes[0]].symbols[1]
        lhs[2] = nodes[self.nodes[1]].symbols[0]
        lhs[3] = nodes[self.nodes[1]].symbols[1]
        lhs[4] = nodes[self.nodes[2]].symbols[0]
        lhs[5] = nodes[self.nodes[2]].symbols[1]
        lhs[6] = nodes[self.nodes[3]].symbols[0]
        lhs[7] = nodes[self.nodes[3]].symbols[1]

        scatteringMatrix[0, 0] = sp.symbols(self.name + '00')
        scatteringMatrix[0, 1] = sp.symbols(self.name + '01')
        scatteringMatrix[0, 2] = sp.symbols(self.name + '02')
        scatteringMatrix[0, 3] = sp.symbols(self.name + '03')
        scatteringMatrix[0, 4] = sp.symbols(self.name + '04')
        scatteringMatrix[0, 5] = sp.symbols(self.name + '05')
        scatteringMatrix[0, 6] = sp.symbols(self.name + '06')
        scatteringMatrix[0, 7] = sp.symbols(self.name + '07')

        scatteringMatrix[1, 0] = sp.symbols(self.name + '10')
        scatteringMatrix[1, 1] = sp.symbols(self.name + '11')
        scatteringMatrix[1, 2] = sp.symbols(self.name + '12')
        scatteringMatrix[1, 3] = sp.symbols(self.name + '13')
        scatteringMatrix[1, 4] = sp.symbols(self.name + '14')
        scatteringMatrix[1, 5] = sp.symbols(self.name + '15')
        scatteringMatrix[1, 6] = sp.symbols(self.name + '16')
        scatteringMatrix[1, 7] = sp.symbols(self.name + '17')

        scatteringMatrix[2, 0] = sp.symbols(self.name + '20')
        scatteringMatrix[2, 1] = sp.symbols(self.name + '21')
        scatteringMatrix[2, 2] = sp.symbols(self.name + '22')
        scatteringMatrix[2, 3] = sp.symbols(self.name + '23')
        scatteringMatrix[2, 4] = sp.symbols(self.name + '24')
        scatteringMatrix[2, 5] = sp.symbols(self.name + '25')
        scatteringMatrix[2, 6] = sp.symbols(self.name + '26')
        scatteringMatrix[2, 7] = sp.symbols(self.name + '27')

        scatteringMatrix[3, 0] = sp.symbols(self.name + '30')
        scatteringMatrix[3, 1] = sp.symbols(self.name + '31')
        scatteringMatrix[3, 2] = sp.symbols(self.name + '32')
        scatteringMatrix[3, 3] = sp.symbols(self.name + '33')
        scatteringMatrix[3, 4] = sp.symbols(self.name + '34')
        scatteringMatrix[3, 5] = sp.symbols(self.name + '35')
        scatteringMatrix[3, 6] = sp.symbols(self.name + '36')
        scatteringMatrix[3, 7] = sp.symbols(self.name + '37')

        scatteringMatrix[4, 0] = sp.symbols(self.name + '40')
        scatteringMatrix[4, 1] = sp.symbols(self.name + '41')
        scatteringMatrix[4, 2] = sp.symbols(self.name + '42')
        scatteringMatrix[4, 3] = sp.symbols(self.name + '43')
        scatteringMatrix[4, 4] = sp.symbols(self.name + '44')
        scatteringMatrix[4, 5] = sp.symbols(self.name + '45')
        scatteringMatrix[4, 6] = sp.symbols(self.name + '46')
        scatteringMatrix[4, 7] = sp.symbols(self.name + '47')

        scatteringMatrix[5, 0] = sp.symbols(self.name + '50')
        scatteringMatrix[5, 1] = sp.symbols(self.name + '51')
        scatteringMatrix[5, 2] = sp.symbols(self.name + '52')
        scatteringMatrix[5, 3] = sp.symbols(self.name + '53')
        scatteringMatrix[5, 4] = sp.symbols(self.name + '54')
        scatteringMatrix[5, 5] = sp.symbols(self.name + '55')
        scatteringMatrix[5, 6] = sp.symbols(self.name + '56')
        scatteringMatrix[5, 7] = sp.symbols(self.name + '57')

        scatteringMatrix[6, 0] = sp.symbols(self.name + '60')
        scatteringMatrix[6, 1] = sp.symbols(self.name + '61')
        scatteringMatrix[6, 2] = sp.symbols(self.name + '62')
        scatteringMatrix[6, 3] = sp.symbols(self.name + '63')
        scatteringMatrix[6, 4] = sp.symbols(self.name + '64')
        scatteringMatrix[6, 5] = sp.symbols(self.name + '65')
        scatteringMatrix[6, 6] = sp.symbols(self.name + '66')
        scatteringMatrix[6, 7] = sp.symbols(self.name + '67')

        scatteringMatrix[7, 0] = sp.symbols(self.name + '70')
        scatteringMatrix[7, 1] = sp.symbols(self.name + '71')
        scatteringMatrix[7, 2] = sp.symbols(self.name + '72')
        scatteringMatrix[7, 3] = sp.symbols(self.name + '73')
        scatteringMatrix[7, 4] = sp.symbols(self.name + '74')
        scatteringMatrix[7, 5] = sp.symbols(self.name + '75')
        scatteringMatrix[7, 6] = sp.symbols(self.name + '76')
        scatteringMatrix[7, 7] = sp.symbols(self.name + '77')

        self.symbols = tuple(sorted(list(scatteringMatrix.free_symbols),
                                    key=str))

        self.equation = sp.Eq(scatteringMatrix * lhs, rhs)

    def setVals(self):
        """Returns numerical values needed to solve the network matrix.

        Should not need to be called by the user.
        """
        return self.numeric_matrix.flatten()

    def update(self):
        """Updates numeric values of matrix from user set optical parameters.

        Must be called manually when values have been changed.
        """
        self.numeric_matrix = np.zeros((8, 8), dtype=np.complex)

        self.rP = self.rExtinction \
            * (self.sLoss - self.pLoss * self.tExtinction) \
            / (1 - self.tExtinction * self.rExtinction)
        self.rS = (self.sLoss - self.pLoss * self.tExtinction) \
            / (1 - self.tExtinction * self.rExtinction)
        self.tS = (self.pLoss * self.tExtinction
                   - self.sLoss * self.tExtinction * self.rExtinction) \
            / (1 - self.tExtinction * self.rExtinction)
        self.tP = (self.pLoss - self.sLoss * self.rExtinction) \
            / (1 - self.tExtinction * self.rExtinction)

        self.rP = np.sqrt(self.rP)
        self.rS = np.sqrt(self.rS)
        self.tP = np.sqrt(self.tP)
        self.tS = np.sqrt(self.tS)

        self.numeric_matrix[0][2] = self.rP
        self.numeric_matrix[0][4] = self.tP
        self.numeric_matrix[1][3] = self.rS
        self.numeric_matrix[1][5] = self.tS
        self.numeric_matrix[2][0] = self.rP
        self.numeric_matrix[2][6] = self.tP
        self.numeric_matrix[3][1] = self.rS
        self.numeric_matrix[3][7] = self.tS
        self.numeric_matrix[4][0] = self.tP
        self.numeric_matrix[4][6] = self.rP
        self.numeric_matrix[5][1] = self.tS
        self.numeric_matrix[5][7] = self.rS
        self.numeric_matrix[6][2] = self.tP
        self.numeric_matrix[6][4] = self.rP
        self.numeric_matrix[7][3] = self.tS
        self.numeric_matrix[7][5] = self.rS

        self.numeric_matrix = rotationMatrix88(-self.theta0, -self.theta1) \
            @ self.numeric_matrix @ rotationMatrix88(self.theta0, self.theta1)

        self.model.updated.append(self.name)


class _TransferComponent(_Component):
    """General class for transfer matrix based 2-port components.

    Can be used to build other transfer matrix based components. All
    'direction' matching logic takes place in this class, ensuring that 'a'
    and 'b' components of the electric field are matched appropriately between
    components.

    Nodes: 2
    """

    def __init__(self, name, nodes, model):
        _Component.__init__(self, name, nodes, model)
        self.node_number = 2
        self.stack_matrix = np.identity(4, dtype=np.complex)
        self.left_swapped = False
        self.right_swapped = False

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(4, 1)
        rhs = sp.zeros(4, 1)
        stackMatrix = sp.zeros(4, 4)

        # check left side of stack alignment
        if nodes[self.nodes[0]].components[0] == self:
            attachedLeft = nodes[self.nodes[0]].components[1]
        else:
            attachedLeft = nodes[self.nodes[0]].components[0]

        # check right side of stack alignment
        if nodes[self.nodes[1]].components[0] == self:
            attachedRight = nodes[self.nodes[1]].components[1]
        else:
            attachedRight = nodes[self.nodes[1]].components[0]

        if isinstance(attachedLeft, Dump) or \
                issubclass(attachedLeft.__class__, _ScatterComponent):
            lhs[0] = nodes[self.nodes[0]].symbols[2]
            lhs[1] = nodes[self.nodes[0]].symbols[0]
            lhs[2] = nodes[self.nodes[0]].symbols[3]
            lhs[3] = nodes[self.nodes[0]].symbols[1]
        else:
            lhs[0] = nodes[self.nodes[0]].symbols[0]
            lhs[1] = nodes[self.nodes[0]].symbols[2]
            lhs[2] = nodes[self.nodes[0]].symbols[1]
            lhs[3] = nodes[self.nodes[0]].symbols[3]

        if isinstance(attachedRight, Dump) or \
                isinstance(attachedRight, Stack) or \
                issubclass(attachedRight.__class__, _ScatterComponent):
            rhs[0] = nodes[self.nodes[1]].symbols[0]
            rhs[1] = nodes[self.nodes[1]].symbols[2]
            rhs[2] = nodes[self.nodes[1]].symbols[1]
            rhs[3] = nodes[self.nodes[1]].symbols[3]
        else:
            rhs[0] = nodes[self.nodes[1]].symbols[2]
            rhs[1] = nodes[self.nodes[1]].symbols[0]
            rhs[2] = nodes[self.nodes[1]].symbols[3]
            rhs[3] = nodes[self.nodes[1]].symbols[1]

        if isinstance(attachedLeft, Stack):
            if attachedLeft.nodes[0] == self.nodes[0]:
                if not(attachedLeft.left_swapped):
                    lhs[0] = nodes[self.nodes[0]].symbols[2]
                    lhs[1] = nodes[self.nodes[0]].symbols[0]
                    lhs[2] = nodes[self.nodes[0]].symbols[3]
                    lhs[3] = nodes[self.nodes[0]].symbols[1]
                    self.left_swapped = True

        if isinstance(attachedRight, Stack):
            if attachedRight.nodes[1] == self.nodes[1]:
                if not(attachedRight.right_swapped):
                    rhs[0] = nodes[self.nodes[1]].symbols[2]
                    rhs[1] = nodes[self.nodes[1]].symbols[0]
                    rhs[2] = nodes[self.nodes[1]].symbols[3]
                    rhs[3] = nodes[self.nodes[1]].symbols[1]
                    self.right_swapped = True

        stackMatrix[0, 0] = sp.symbols(self.name + '00')
        stackMatrix[0, 1] = sp.symbols(self.name + '01')
        stackMatrix[0, 2] = sp.symbols(self.name + '02')
        stackMatrix[0, 3] = sp.symbols(self.name + '03')

        stackMatrix[1, 0] = sp.symbols(self.name + '10')
        stackMatrix[1, 1] = sp.symbols(self.name + '11')
        stackMatrix[1, 2] = sp.symbols(self.name + '12')
        stackMatrix[1, 3] = sp.symbols(self.name + '13')

        stackMatrix[2, 0] = sp.symbols(self.name + '20')
        stackMatrix[2, 1] = sp.symbols(self.name + '21')
        stackMatrix[2, 2] = sp.symbols(self.name + '22')
        stackMatrix[2, 3] = sp.symbols(self.name + '23')

        stackMatrix[3, 0] = sp.symbols(self.name + '30')
        stackMatrix[3, 1] = sp.symbols(self.name + '31')
        stackMatrix[3, 2] = sp.symbols(self.name + '32')
        stackMatrix[3, 3] = sp.symbols(self.name + '33')

        self.symbols = tuple(sorted(list(stackMatrix.free_symbols), key=str))

        self.equation = sp.Eq(stackMatrix * lhs, rhs)

    def setVals(self):
        """Returns numerical values needed to solve the network matrix.

        Should not need to be called by the user.
        """

        return self.stack_matrix.flatten()


class Stack(_TransferComponent):
    """
            Stack of one or more layers for linking components.
            Nodes: 2
    """

    def __init__(self, name, nodes, model):
        _TransferComponent.__init__(self, name, nodes, model)

    def set_length(self, length):
        """Sets stack transfer matrix to a single layer of thickness length,
        in units of wavelength.

        Parameters
        ----------
        length : double
                Optical thickness of stack in units of wavelength.
        """
        self.stack_matrix = np.identity(4, dtype=np.complex)

        self.stack_matrix[0][0] = np.exp(1j * length * 2 * np.pi)
        self.stack_matrix[1][1] = np.exp(-1j * length * 2 * np.pi)
        self.stack_matrix[2][2] = np.exp(1j * length * 2 * np.pi)
        self.stack_matrix[3][3] = np.exp(-1j * length * 2 * np.pi)

        self.model.updated.append(self.name)

    def set_pyctmm(self, cstack):
        """Sets stack transfer matrix to that of a pyctmm stack.

        pyctmm allows multilayer stacks to be defined, this can be used for
        optical paths consisting of several parallel interfaces between
        materials, or to model thin film stacks, for example anti-reflection
        coatings. Metallic materials can also be modelled, allowing realistic
        properties of metallic mirrors to be modeled.

        The pyctmm stack will be evaluated when passed to this function, so all
        desired properties of the stack should have already been set. If the
        stack properties are changed, this function should be called again.

        Parameters
        ----------
        cstack : pyctmm stack
                A pyctmm with pre-set layer thicknesses and refractive indexes.
        """

        pyctmm.evaluate(cstack)
        self.stack_matrix = pyctmm.get_matrix(cstack)

        self.model.updated.append(self.name)


class FaradayRotator(_ScatterComponent):
    """Faraday rotator with variable rotation angle.

    Nodes: 2

    Attributes
    ----------
    rotation : double
            Angle by which the incided polarisation state is rotated, measured
            clockwise looking from the zeroth node to the first.
    """

    def __init__(self, name, nodes, model):
        _ScatterComponent.__init__(self, name, nodes, model, 2)
        self.rotation = 0

        self.numeric_matrix = np.zeros((4, 4), dtype=complex)

        self.numeric_matrix[0][2] = np.cos(self.rotation)
        self.numeric_matrix[0][3] = -np.sin(self.rotation)
        self.numeric_matrix[1][2] = np.sin(self.rotation)
        self.numeric_matrix[1][3] = np.cos(self.rotation)
        self.numeric_matrix[2][0] = np.cos(self.rotation)
        self.numeric_matrix[2][1] = -np.sin(self.rotation)
        self.numeric_matrix[3][0] = np.sin(self.rotation)
        self.numeric_matrix[3][1] = np.cos(self.rotation)

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(4, 1)
        rhs = sp.zeros(4, 1)
        scatteringMatrix = sp.zeros(4, 4)

        rhs[0] = nodes[self.nodes[0]].symbols[2]
        rhs[1] = nodes[self.nodes[0]].symbols[3]
        rhs[2] = nodes[self.nodes[1]].symbols[2]
        rhs[3] = nodes[self.nodes[1]].symbols[3]

        lhs[0] = nodes[self.nodes[0]].symbols[0]
        lhs[1] = nodes[self.nodes[0]].symbols[1]
        lhs[2] = nodes[self.nodes[1]].symbols[0]
        lhs[3] = nodes[self.nodes[1]].symbols[1]

        scatteringMatrix[0, 0] = sp.symbols(self.name + '00')
        scatteringMatrix[0, 1] = sp.symbols(self.name + '01')
        scatteringMatrix[0, 2] = sp.symbols(self.name + '02')
        scatteringMatrix[0, 3] = sp.symbols(self.name + '03')

        scatteringMatrix[1, 0] = sp.symbols(self.name + '10')
        scatteringMatrix[1, 1] = sp.symbols(self.name + '11')
        scatteringMatrix[1, 2] = sp.symbols(self.name + '12')
        scatteringMatrix[1, 3] = sp.symbols(self.name + '13')

        scatteringMatrix[2, 0] = sp.symbols(self.name + '20')
        scatteringMatrix[2, 1] = sp.symbols(self.name + '21')
        scatteringMatrix[2, 2] = sp.symbols(self.name + '22')
        scatteringMatrix[2, 3] = sp.symbols(self.name + '23')

        scatteringMatrix[3, 0] = sp.symbols(self.name + '30')
        scatteringMatrix[3, 1] = sp.symbols(self.name + '31')
        scatteringMatrix[3, 2] = sp.symbols(self.name + '32')
        scatteringMatrix[3, 3] = sp.symbols(self.name + '33')

        self.symbols = tuple(sorted(list(scatteringMatrix.free_symbols),
                                    key=str))

        self.equation = sp.Eq(scatteringMatrix * lhs, rhs)

    def setVals(self):
        """Returns numerical values needed to solve the network matrix.

        Should not need to be called by the user.
        """
        return self.numeric_matrix.flatten()

    def update(self):
        """Updates numeric values of matrix from user set optical parameters.

        Must be called manually when values have been changed.
        """
        self.numeric_matrix = np.zeros((4, 4), dtype=np.complex)

        self.numeric_matrix[0][2] = np.cos(self.rotation)
        self.numeric_matrix[0][3] = np.sin(self.rotation)
        self.numeric_matrix[1][2] = -np.sin(self.rotation)
        self.numeric_matrix[1][3] = np.cos(self.rotation)
        self.numeric_matrix[2][0] = np.cos(self.rotation)
        self.numeric_matrix[2][1] = np.sin(self.rotation)
        self.numeric_matrix[3][0] = -np.sin(self.rotation)
        self.numeric_matrix[3][1] = np.cos(self.rotation)

        self.model.updated.append(self.name)


class IdealIsolator(_ScatterComponent):
    """Ideal optical isolator, prevents light propagating from the first node
to the zeroth.

Nodes: 2

Attributes
----------
isolationCoefficient : double
    Fraction of intensity of light propagating from node 1 to 0 that is blocked
    by the isolator.
"""

    def __init__(self, name, nodes, model):
        _ScatterComponent.__init__(self, name, nodes, model, 2)
        self.isolationCoefficient = 1

        self.numeric_matrix = np.zeros((4, 4), dtype=complex)

        self.numeric_matrix[0][2] = np.sqrt(1 - self.isolationCoefficient)
        self.numeric_matrix[1][3] = np.sqrt(1 - self.isolationCoefficient)
        self.numeric_matrix[2][0] = 1
        self.numeric_matrix[3][1] = 1

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(4, 1)
        rhs = sp.zeros(4, 1)
        scatteringMatrix = sp.zeros(4, 4)

        rhs[0] = nodes[self.nodes[0]].symbols[2]
        rhs[1] = nodes[self.nodes[0]].symbols[3]
        rhs[2] = nodes[self.nodes[1]].symbols[2]
        rhs[3] = nodes[self.nodes[1]].symbols[3]

        lhs[0] = nodes[self.nodes[0]].symbols[0]
        lhs[1] = nodes[self.nodes[0]].symbols[1]
        lhs[2] = nodes[self.nodes[1]].symbols[0]
        lhs[3] = nodes[self.nodes[1]].symbols[1]

        scatteringMatrix[0, 0] = sp.symbols(self.name + '00')
        scatteringMatrix[0, 1] = sp.symbols(self.name + '01')
        scatteringMatrix[0, 2] = sp.symbols(self.name + '02')
        scatteringMatrix[0, 3] = sp.symbols(self.name + '03')

        scatteringMatrix[1, 0] = sp.symbols(self.name + '10')
        scatteringMatrix[1, 1] = sp.symbols(self.name + '11')
        scatteringMatrix[1, 2] = sp.symbols(self.name + '12')
        scatteringMatrix[1, 3] = sp.symbols(self.name + '13')

        scatteringMatrix[2, 0] = sp.symbols(self.name + '20')
        scatteringMatrix[2, 1] = sp.symbols(self.name + '21')
        scatteringMatrix[2, 2] = sp.symbols(self.name + '22')
        scatteringMatrix[2, 3] = sp.symbols(self.name + '23')

        scatteringMatrix[3, 0] = sp.symbols(self.name + '30')
        scatteringMatrix[3, 1] = sp.symbols(self.name + '31')
        scatteringMatrix[3, 2] = sp.symbols(self.name + '32')
        scatteringMatrix[3, 3] = sp.symbols(self.name + '33')

        self.symbols = tuple(sorted(list(scatteringMatrix.free_symbols),
                                    key=str))

        self.equation = sp.Eq(scatteringMatrix * lhs, rhs)

    def setVals(self):
        """Returns numerical values needed to solve the network matrix.

        Should not need to be called by the user.
        """
        return self.numeric_matrix.flatten()

    def update(self):
        """Updates numeric values of matrix from user set optical parameters.

        Must be called manually when values have been changed.
        """
        self.numeric_matrix = np.zeros((4, 4), dtype=np.complex)

        self.numeric_matrix[0][2] = np.sqrt(1 - self.isolationCoefficient)
        self.numeric_matrix[1][3] = np.sqrt(1 - self.isolationCoefficient)
        self.numeric_matrix[2][0] = 1
        self.numeric_matrix[3][1] = 1

        self.model.updated.append(self.name)


class Reflector(_ScatterComponent):
    """Reflector with specified reflectivity coefficients, and separated inputs
    and outputs. Rotation allows for change of polarisation axes.

    Allows for polarisation mixing on reflection.

    Nodes: 2

    Attributes
    ----------
    rpp : complex
            Reflectivity coefficient for P to P polarised light.
    rss : complex
            Reflectivity coefficient for S to S polarised light.
    rps : complex
            Reflectivity coefficient for P to S polarised light.
    rsp : complex
            Reflectivity coefficient for S to P polarised light.
    rotation : double
            Rotation of polarisation coordinate system relative the rest of the
            model, measured clockwise from the S polarised axis looking from
            the zeroth node to the first.
    """

    def __init__(self, name, nodes, model):
        _ScatterComponent.__init__(self, name, nodes, model, 2)

        self.rotation = 0

        self.rpp = 1
        self.rss = 1
        self.rsp = 0
        self.rps = 0

        self.numeric_matrix = np.zeros((4, 4), dtype=complex)

        self.numeric_matrix[0][2] = self.rpp
        self.numeric_matrix[0][3] = self.rsp
        self.numeric_matrix[1][2] = self.rps
        self.numeric_matrix[1][3] = self.rss
        self.numeric_matrix[2:, :2] = self.numeric_matrix[:2, 2:]

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(4, 1)
        rhs = sp.zeros(4, 1)
        scatteringMatrix = sp.zeros(4, 4)

        rhs[0] = nodes[self.nodes[0]].symbols[2]
        rhs[1] = nodes[self.nodes[0]].symbols[3]
        rhs[2] = nodes[self.nodes[1]].symbols[2]
        rhs[3] = nodes[self.nodes[1]].symbols[3]

        lhs[0] = nodes[self.nodes[0]].symbols[0]
        lhs[1] = nodes[self.nodes[0]].symbols[1]
        lhs[2] = nodes[self.nodes[1]].symbols[0]
        lhs[3] = nodes[self.nodes[1]].symbols[1]

        scatteringMatrix[0, 0] = sp.symbols(self.name + '00')
        scatteringMatrix[0, 1] = sp.symbols(self.name + '01')
        scatteringMatrix[0, 2] = sp.symbols(self.name + '02')
        scatteringMatrix[0, 3] = sp.symbols(self.name + '03')

        scatteringMatrix[1, 0] = sp.symbols(self.name + '10')
        scatteringMatrix[1, 1] = sp.symbols(self.name + '11')
        scatteringMatrix[1, 2] = sp.symbols(self.name + '12')
        scatteringMatrix[1, 3] = sp.symbols(self.name + '13')

        scatteringMatrix[2, 0] = sp.symbols(self.name + '20')
        scatteringMatrix[2, 1] = sp.symbols(self.name + '21')
        scatteringMatrix[2, 2] = sp.symbols(self.name + '22')
        scatteringMatrix[2, 3] = sp.symbols(self.name + '23')

        scatteringMatrix[3, 0] = sp.symbols(self.name + '30')
        scatteringMatrix[3, 1] = sp.symbols(self.name + '31')
        scatteringMatrix[3, 2] = sp.symbols(self.name + '32')
        scatteringMatrix[3, 3] = sp.symbols(self.name + '33')

        self.symbols = tuple(sorted(list(scatteringMatrix.free_symbols),
                                    key=str))

        self.equation = sp.Eq(scatteringMatrix * lhs, rhs)

    def setVals(self):
        """Returns numerical values needed to solve the network matrix.

        Should not need to be called by the user.
        """
        return self.numeric_matrix.flatten()

    def update(self):
        """Updates numeric values of matrix from user set optical parameters.

        Must be called manually when values have been changed.
        """
        self.numeric_matrix = np.zeros((4, 4), dtype=np.complex)

        self.numeric_matrix[0][2] = self.rpp
        self.numeric_matrix[0][3] = self.rsp
        self.numeric_matrix[1][2] = self.rps
        self.numeric_matrix[1][3] = self.rss
        self.numeric_matrix[2:, :2] = self.numeric_matrix[:2, 2:]

        self.numeric_matrix = rotationMatrix44(-self.rotation) \
            @ self.numeric_matrix \
            @ rotationMatrix44(self.rotation)

        self.model.updated.append(self.name)


class Waveplate(_ScatterComponent):
    """Waveplate with variable retardance and rotation angle.

    Nodes: 2

    Attributes
    ----------
    rotation : double
            Rotation of waveplate, measured clockwise from the S polarised
            axis, looking from the zeroth node to the first.
    retardance : double
            Phase retardance of the waveplate - the phase difference introduced
            between polarisation components passing through the fast and slow
            axes of the waveplate.
    """

    def __init__(self, name, nodes, model):
        _ScatterComponent.__init__(self, name, nodes, model, 2)
        self.rotation = 0
        self.retardance = 0

        self.numeric_matrix = np.zeros((4, 4), dtype=complex)

        self.numeric_matrix[0][2] = np.exp(-1j * self.retardance / 2)
        self.numeric_matrix[1][3] = np.exp(1j * self.retardance / 2)
        self.numeric_matrix[2][0] = np.exp(-1j * self.retardance / 2)
        self.numeric_matrix[3][1] = np.exp(1j * self.retardance / 2)

        self.numeric_matrix = rotationMatrix44(-self.rotation) \
            @ self.numeric_matrix \
            @ rotationMatrix44(self.rotation)

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(4, 1)
        rhs = sp.zeros(4, 1)
        scatteringMatrix = sp.zeros(4, 4)

        rhs[0] = nodes[self.nodes[0]].symbols[2]
        rhs[1] = nodes[self.nodes[0]].symbols[3]
        rhs[2] = nodes[self.nodes[1]].symbols[2]
        rhs[3] = nodes[self.nodes[1]].symbols[3]

        lhs[0] = nodes[self.nodes[0]].symbols[0]
        lhs[1] = nodes[self.nodes[0]].symbols[1]
        lhs[2] = nodes[self.nodes[1]].symbols[0]
        lhs[3] = nodes[self.nodes[1]].symbols[1]

        scatteringMatrix[0, 0] = sp.symbols(self.name + '00')
        scatteringMatrix[0, 1] = sp.symbols(self.name + '01')
        scatteringMatrix[0, 2] = sp.symbols(self.name + '02')
        scatteringMatrix[0, 3] = sp.symbols(self.name + '03')

        scatteringMatrix[1, 0] = sp.symbols(self.name + '10')
        scatteringMatrix[1, 1] = sp.symbols(self.name + '11')
        scatteringMatrix[1, 2] = sp.symbols(self.name + '12')
        scatteringMatrix[1, 3] = sp.symbols(self.name + '13')

        scatteringMatrix[2, 0] = sp.symbols(self.name + '20')
        scatteringMatrix[2, 1] = sp.symbols(self.name + '21')
        scatteringMatrix[2, 2] = sp.symbols(self.name + '22')
        scatteringMatrix[2, 3] = sp.symbols(self.name + '23')

        scatteringMatrix[3, 0] = sp.symbols(self.name + '30')
        scatteringMatrix[3, 1] = sp.symbols(self.name + '31')
        scatteringMatrix[3, 2] = sp.symbols(self.name + '32')
        scatteringMatrix[3, 3] = sp.symbols(self.name + '33')

        self.symbols = tuple(sorted(list(scatteringMatrix.free_symbols),
                                    key=str))

        self.equation = sp.Eq(scatteringMatrix * lhs, rhs)

    def setVals(self):
        """Returns numerical values needed to solve the network matrix.

        Should not need to be called by the user.
        """
        return self.numeric_matrix.flatten()

    def update(self):
        """Updates numeric values of matrix from user set optical parameters.

        Must be called manually when values have been changed.
        """
        self.numeric_matrix = np.zeros((4, 4), dtype=np.complex)

        self.numeric_matrix[0][2] = np.exp(-1j * self.retardance / 2)
        self.numeric_matrix[1][3] = np.exp(1j * self.retardance / 2)
        self.numeric_matrix[2][0] = np.exp(-1j * self.retardance / 2)
        self.numeric_matrix[3][1] = np.exp(1j * self.retardance / 2)

        self.numeric_matrix = rotationMatrix44(-self.rotation) \
            @ self.numeric_matrix \
            @ rotationMatrix44(self.rotation)

        self.model.updated.append(self.name)


class Polariser(_ScatterComponent):
    """Linear polariser with variable rotation angle, extinction ratio and
    loss for transmission of in-plane polarised light.

    Nodes: 2

    Attributes
    ----------
    rotation : double
            Rotation of the transmission axis of the polariser, measured
            clockwise from the S polarised axis, looking from the zeroth node
            to the first.
    extinction : double
            Extinction coefficient, defined as the intensity ratio of
            transmitted light polarised perpendicular to the transmission axis
            of the polariser to the intensity of transmitted light polarised
            parallel to the transmission axis.
    loss : double
            Transmission loss, defined as the intensity ratio of transmitted
            light polarised parallel to the transmission axis to the intensity
            of the component of the incident light polarised along that axis.
    """

    def __init__(self, name, nodes, model):
        _ScatterComponent.__init__(self, name, nodes, model, 2)
        self.rotation = 0
        self.extinction = 0
        self.loss = 0

        self.numeric_matrix = np.zeros((4, 4))

        self.numeric_matrix[1][3] = np.sqrt(1 - self.loss)
        self.numeric_matrix[3][1] = self.numeric_matrix[1][3]
        self.numeric_matrix[0][2] = np.sqrt(self.extinction) \
            * self.numeric_matrix[1][3]
        self.numeric_matrix[2][0] = self.numeric_matrix[0][2]

        self.numeric_matrix = rotationMatrix44(-self.rotation) \
            @ self.numeric_matrix \
            @ rotationMatrix44(self.rotation)

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(4, 1)
        rhs = sp.zeros(4, 1)
        scatteringMatrix = sp.zeros(4, 4)

        rhs[0] = nodes[self.nodes[0]].symbols[2]
        rhs[1] = nodes[self.nodes[0]].symbols[3]
        rhs[2] = nodes[self.nodes[1]].symbols[2]
        rhs[3] = nodes[self.nodes[1]].symbols[3]

        lhs[0] = nodes[self.nodes[0]].symbols[0]
        lhs[1] = nodes[self.nodes[0]].symbols[1]
        lhs[2] = nodes[self.nodes[1]].symbols[0]
        lhs[3] = nodes[self.nodes[1]].symbols[1]

        scatteringMatrix[0, 0] = sp.symbols(self.name + '00')
        scatteringMatrix[0, 1] = sp.symbols(self.name + '01')
        scatteringMatrix[0, 2] = sp.symbols(self.name + '02')
        scatteringMatrix[0, 3] = sp.symbols(self.name + '03')

        scatteringMatrix[1, 0] = sp.symbols(self.name + '10')
        scatteringMatrix[1, 1] = sp.symbols(self.name + '11')
        scatteringMatrix[1, 2] = sp.symbols(self.name + '12')
        scatteringMatrix[1, 3] = sp.symbols(self.name + '13')

        scatteringMatrix[2, 0] = sp.symbols(self.name + '20')
        scatteringMatrix[2, 1] = sp.symbols(self.name + '21')
        scatteringMatrix[2, 2] = sp.symbols(self.name + '22')
        scatteringMatrix[2, 3] = sp.symbols(self.name + '23')

        scatteringMatrix[3, 0] = sp.symbols(self.name + '30')
        scatteringMatrix[3, 1] = sp.symbols(self.name + '31')
        scatteringMatrix[3, 2] = sp.symbols(self.name + '32')
        scatteringMatrix[3, 3] = sp.symbols(self.name + '33')

        self.symbols = tuple(sorted(list(scatteringMatrix.free_symbols),
                                    key=str))

        self.equation = sp.Eq(scatteringMatrix * lhs, rhs)

    def setVals(self):
        """Returns numerical values needed to solve the network matrix.

        Should not need to be called by the user.
        """
        return self.numeric_matrix.flatten()

    def update(self):
        """Updates numeric values of matrix from user set optical parameters.

        Must be called manually when values have been changed.
        """
        self.numeric_matrix = np.zeros((4, 4))

        self.numeric_matrix[1][3] = np.sqrt(1 - self.loss)
        self.numeric_matrix[3][1] = self.numeric_matrix[1][3]
        self.numeric_matrix[0][2] = np.sqrt(self.extinction) \
            * self.numeric_matrix[1][3]
        self.numeric_matrix[2][0] = self.numeric_matrix[0][2]

        self.numeric_matrix = rotationMatrix44(-self.rotation) \
            @ self.numeric_matrix \
            @ rotationMatrix44(self.rotation)

        self.model.updated.append(self.name)


class Mirror(_Component):
    """Mirror for reflecting stack output back into input.

    Connects to the end of a stack and couples together the inputs and outputs
    of the stack using the specified amplitude reflection coefficients.

    Nodes: 1

    Attributes
    ----------
    rP : complex
            P polarised amplitude reflectivity coefficient.
    rS : complex
            S polarised amplitude reflectivity coefficient.
    """

    def __init__(self, name, nodes, model):
        _Component.__init__(self, name, nodes, model)
        self.node_number = 1
        self.rP = 1
        self.rS = 1

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(2, 1)
        rhs = sp.zeros(2, 1)

        lhs[0] = -sp.symbols(self.name + '_rP') \
            * nodes[self.nodes[0]].symbols[2]
        lhs[1] = -sp.symbols(self.name + '_rS') \
            * nodes[self.nodes[0]].symbols[3]

        rhs[0] = nodes[self.nodes[0]].symbols[0]
        rhs[1] = nodes[self.nodes[0]].symbols[1]

        self.symbols = (sp.symbols(self.name + '_rP'),
                        sp.symbols(self.name + '_rS'))

        self.equation = sp.Eq(lhs, rhs)

    def setVals(self):
        """Returns numerical values needed to solve the network matrix.

        Should not need to be called by the user.
        """
        return [self.rP, self.rS]


class Dump(_Component):
    """Beam dump to terminate stack.

    Sets the input into the end of an otherwise uncoupled stack to zero. Dumps
    must be attached to all free stack ends for the network matrix equation to
    be fully defined.

    Nodes: 1
    """

    def __init__(self, name, nodes, model):
        _Component.__init__(self, name, nodes, model)
        self.node_number = 1

    def initEquation(self, nodes):
        """Initialises sympy equation for component.

        Should not need to be called by the user.

        Parameters
        ----------
        nodes : list of scampy.Node
                Nodes to which the component is attached.
        """
        lhs = sp.zeros(2, 1)
        rhs = sp.zeros(2, 1)

        lhs[0] = nodes[self.nodes[0]].symbols[2]
        lhs[1] = nodes[self.nodes[0]].symbols[3]

        rhs[0] = 0
        rhs[1] = 0

        self.equation = sp.Eq(lhs, rhs)
