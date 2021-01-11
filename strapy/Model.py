from .Node import Node
from .Detector import Detector
from . import components
import sympy as sp
import numpy as np
import timeit
from scipy.sparse.linalg import lsqr


class Model:
    """Defines the optical network to be modelled.

    The Model class holds the information that defines the structure of the
    optical network to be modelled with strapy, along with functions for
    building and evaluating the model. In general only the `wavelength`
    attribute should be accessed directly by the user - all other attributes
    should be set through the member functions.

    Attributes
    ----------
    wavelength : float
            The vacuum (n = 1) wavelength of illuminating light.
    components : dict
            Dictionary of optical components that have been added to the model.
    detectors : dict
            Dictionary of detectors that have been added to the model.
    nodes : dict
            Dictionary of nodes that have been added to the model.
    symbols : list
            List of sympy symbols used by the model.
    equations : list
            List of sympy equations used by the model.
    rhsVariables : list
            List of sympy symbols that are included in the right hand side
            vector of the matrix equation.
    matrixVariables : list
            List of sympy symbols that are included in the network matrix.
    updated : list
            List of optical components that have changed since the model was
            last evaluated.
    useLambdify : bool
            True if sympy's lambdify functionality should be used to set right
            hand side vector and matrix before evaluation. If false only
            components that are included in the `updated` list are updated in
            the matrix equation before solving. Lambdify is currently a
            fallback if there is more than one sympy symbol per matrix element.
    """

    def __init__(self):
        self.wavelength = 633e-9
        self.components = {}
        self.detectors = {}
        self.nodes = {}
        self.symbols = []
        self.equations = []
        self.rhsVariables = []
        self.matrixVariables = []
        self.updated = []
        self.useLambdify = False

    def add_component(self, component, name, nodes):
        """Adds component to model and updates the node list.

        Parameters
        ----------
        component : strapy.components.component
                Type of component to be added.
        name : str
                Unique name of component.
        nodes : str or tuple of str
                Node(s) to which the component is attached.
        """

        # allows for passing string instead of tuple of strings
        if not(isinstance(nodes, tuple)):
            nodes = (nodes,)

        self.components[name] = component(name, nodes, self)
        self.updated.append(name)

        if len(nodes) != self.components[name].node_number:
            raise Exception(
                'Specified nodes ({}) not equal to required node',
                ' number ({}).'.format(
                    len(nodes),
                    self.components[name].node_number))

        # keep nodes list up to date when adding a component
        for node in nodes:
            if not(node in self.nodes):
                self.nodes[node] = Node(node)
            if len(self.nodes[node].components) > 1:
                raise Exception(
                    'Nodes may not be linked to more that two',
                    ' components. Node {} is attached to {}.'.format(
                        node, len(self.nodes[node].components) + 1))
            else:
                self.nodes[node].components.append(self.components[name])

    def add_detector(self, name, node, properties=('amplitude',)):
        """Adds detector to model and updates list of nodes.

        For details on currently implemented properties see
        :py:class:`strapy.Detector()`.

        Parameters
        ----------
        name : str
                Unique name of the detector.
        nodes : str
                Node monitored by detector.
        properties : tuple of str
                Optical properties to be logged - see
                :py:class:`strapy.Detector()` for current options.
        """

        node = (node,)
        self.detectors[name] = Detector(name, node, properties)

        if not(node[0] in self.nodes):
            self.nodes[node[0]] = Node(node[0])

    def build(self, verbose=False):
        """Builds network matrix from defined components.

        The model must be built before `evaluate()` is called. Additionally, if
        components or detectors are added to the model between runs, the model
        must be rebuilt.

        Parameters
        ----------
        verbose : bool
                If true details of the constructed optical network will be
                printed during the build process.
        """

        sourceFlag = False

        # check for at least one source in the network
        for component in self.components.values():
            if isinstance(component, components.Source):
                sourceFlag = True
                break

        if not(sourceFlag):
            raise Exception('No sources present.')

        # Threadable.
        for node in self.nodes.values():
            if len(node.components) != 2:
                raise Exception(
                    'Nodes must be linked to exactly two',
                    ' components, node {} is attached to {}.'.format(
                        node.name, len(node.components)))

            self.symbols.extend(node.symbols)

            if verbose:
                print('Node {}'.format(node.name))
                print('\tComponent 1: {}'.format(node.components[0].name))
                print('\tComponent 2: {}'.format(node.components[1].name))

        # Should be threadable.
        for component in self.components.values():
            component.initEquation(self.nodes)
            component.symbolIdxs = [[] for _ in range(len(component.symbols))]

            for line in range(len(component.equation.lhs)):
                self.equations.append(sp.Eq(component.equation.lhs[line],
                                            component.equation.rhs[line]))
            if isinstance(component, components.Source):
                self.rhsVariables.extend(component.symbols)
            else:
                self.matrixVariables.extend(component.symbols)

        networkMatrix, rhsVector = sp.linear_eq_to_matrix(self.equations,
                                                          self.symbols)
        self.matrixShape = networkMatrix.shape

        self.matrix = np.zeros(self.matrixShape, dtype=np.complex)
        self.rhs = np.zeros(rhsVector.shape, dtype=np.complex)

        # These loops may be thread safe, other than the lambdify command, and
        # the loops could therefore be parfored. Unknown if access to a sympy
        # matrix is actually thread safe, needs testing.
        for i in range(networkMatrix.shape[0]):
            for j in range(networkMatrix.shape[1]):
                # isinstance() is significantly faster than .equals(0)
                if not(isinstance(networkMatrix[i, j], sp.core.numbers.Zero)):
                    if networkMatrix[i, j].is_constant():
                        self.matrix[i, j] = sp.N(networkMatrix[i, j])
                    else:
                        for component in self.components.values():
                            for k, symbol in enumerate(component.symbols):
                                if len(networkMatrix[i, j].free_symbols) > 1:
                                    print(
                                        "Multiple symbols per matrix element",
                                        ", falling back to lambdify matrix",
                                        " setting.")
                                    self.useLambdify = True
                                    self.setMatrix = sp.lambdify(
                                        self.matrixVariables,
                                        networkMatrix,
                                        modules=["numpy"])
                                if networkMatrix[i, j] == symbol:
                                    component.symbolIdxs[k].append((i, j, 1))
                                if networkMatrix[i, j] == -symbol:
                                    component.symbolIdxs[k].append((i, j, -1))

        # Again, may be threadable.
        for i in range(len(rhsVector)):
            # isinstance() is significantly faster than .equals(0)
            if not(isinstance(rhsVector[i], sp.core.numbers.Zero)):
                if rhsVector[i].is_constant():
                    self.rhs[i, j] = sp.N(rhsVector[i])
                else:
                    for component in self.components.values():
                        for k, symbol in enumerate(component.symbols):
                            if len(rhsVector[i].free_symbols) > 1:
                                print(
                                    "Multiple symbols per rhs vector element",
                                    ", falling back to lambdify rhs vector",
                                    " setting.")
                                self.useLambdify = True
                                self.setRhs = sp.lambdify(
                                    self.rhsVariables,
                                    rhsVector,
                                    modules=["numpy"])
                            if rhsVector[i] == symbol:
                                component.symbolIdxs[k].append((i, 0, 1))
                            if rhsVector[i] == -symbol:
                                component.symbolIdxs[k].append((i, 0, -1))

        # identify the location in the solution vector which the detector
        # should detect.
        for detector in self.detectors.values():
            detector.node_index = self.symbols.index(
                self.nodes[detector.node[0]].symbols[0])

        self.matrixPassVector = np.empty((len(self.matrixVariables),),
                                         dtype=np.complex)
        self.rhsPassVector = np.empty((len(self.rhsVariables),),
                                      dtype=np.complex)

        self.sparcity = len(self.matrixPassVector) \
            / (self.matrixShape[0] * self.matrixShape[1])

        # Should be threadable.
        for component in self.components.values():
            # print(component.name, component.symbols, component.symbolIdxs)
            if isinstance(component, components.Source):
                component.set_slice = slice(self.rhsVariables.index(
                    component.symbols[0]),
                    self.rhsVariables.index(component.symbols[0])
                    + len(component.symbols))
            elif not(isinstance(component, components.Dump)):
                component.set_slice = slice(
                    self.matrixVariables.index(component.symbols[0]),
                    self.matrixVariables.index(component.symbols[0])
                    + len(component.symbols))

        if verbose:
            print('\nNetwork built, {} matrix.\n'.format(networkMatrix.shape))

    def _update(self):
        """Updates numerical values in model matrix.

        Model must have been built first.
        Should not be called externally.
        """

        # Threadable
        for key in self.updated:
            if isinstance(self.components[key], components.Source):
                self.rhsPassVector[self.components[key].set_slice] = \
                    self.components[key].setVals()
            elif not(isinstance(self.components[key], components.Dump)):
                self.matrixPassVector[self.components[key].set_slice] = \
                    self.components[key].setVals()
        self.updated.clear()

    def evaluate(self, timing=False):
        """Solve the network matrix and log optical properties at detectors.

        `build()` must have been called before the model is evaluated.

        Parameters
        ----------
        timing : bool
                If true, the times taken to set values in the right hand side
                vector and network matrix (`set_time`), solve the network
                equation (`solve_time`) and pull out the detected values
                (`detector_time`) are returned."""

        if timing:
            set_time = timeit.default_timer()

            if not(self.useLambdify):
                for key in self.updated:
                    if isinstance(self.components[key], components.Source):
                        vals = self.components[key].setVals()
                        for i, index in enumerate(
                                self.components[key].symbolIdxs):
                            self.rhs[index[0]] = vals[i] * index[2]
                    elif not(isinstance(self.components[key],
                                        components.Dump)):
                        vals = self.components[key].setVals()
                        for i, index in enumerate(
                                self.components[key].symbolIdxs):
                            self.matrix[index[0], index[1]
                                        ] = vals[i] * index[2]
                self.updated.clear()
            else:
                self._update()
                self.matrix = self.setMatrix(*self.matrixPassVector)
                self.rhs = self.setRhs(*self.rhsPassVector)

            set_time = timeit.default_timer() - set_time

            solve_time = timeit.default_timer()
            self.solution_vector = np.linalg.solve(self.matrix, self.rhs)
            solve_time = timeit.default_timer() - solve_time

            detector_time = timeit.default_timer()
            for detector in self.detectors.values():
                detector.update(self.solution_vector)
            detector_time = timeit.default_timer() - detector_time

            return (set_time, solve_time, detector_time)
        else:
            if not(self.useLambdify):
                # Should be threadable.
                for key in self.updated:
                    if isinstance(self.components[key], components.Source):
                        vals = self.components[key].setVals()
                        for i, symbol in enumerate(
                                self.components[key].symbolIdxs):
                            for index in symbol:
                                self.rhs[index[0]] = vals[i] * index[2]
                    elif not(isinstance(self.components[key],
                                        components.Dump)):
                        vals = self.components[key].setVals()
                        for i, symbol in enumerate(
                                self.components[key].symbolIdxs):
                            for index in symbol:
                                self.matrix[index[0], index[1]] = \
                                    vals[i] * index[2]
                self.updated.clear()
            else:
                self._update()
                self.matrix = self.setMatrix(*self.matrixPassVector)
                self.rhs = self.setRhs(*self.rhsPassVector)

            self.solution_vector = np.linalg.solve(self.matrix, self.rhs)

            for detector in self.detectors.values():
                detector.update(self.solution_vector)
