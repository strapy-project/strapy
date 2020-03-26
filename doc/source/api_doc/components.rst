Components
=================================

The components module holds all optical components currently defined in
strapy. Each strapy optical component is defined in a separate class, and must
inherit at least from the `_Component` class, and optionally from the
`_ScatterComponent` or `_TransferComponent` class.

When defining new components, if not using the `_ScatterComponent` or
`_TransferComponent` classes, care must be taken to ensure the correct
stack attachment logic is employed in the `_TransferComponent.initEquation`
method. Adding new components that do not inherit from the `_ScatterComponent`
or `_TransferComponent` is therefore discouraged, to avoid cluttering the 
attachment logic.

Coupling of counter propagating electric field components is resolved as
follows:

	* Sources always emit **into** stacks; sources set the 'a' components of their single node.
	* Beamsplitters always take inputs **from** stacks; the input to any port of the beamsplitter is the 'a' component at the given node.
	* General multi-port devices will always take inputs **from** stacks; the input to any port is the 'a' component at the given node.
	* Dumps always block the component going **into** the stack; the 'b' component is set to 0.

All direction logic therefore takes place in the stack equation initialisation.
This ensures there is no change to the physics depending on how components are
ordered. This model allows light to propagate in both directions
simultaneously, so there should be no concept of 'forwards'.

Eventually this should be automated by checking all nodes connect to at least
one stack and inserting a identity matrix connection stack if this is not the
case. In the long term placing individual components in their own file is
probably a better approach as well.

.. autosummary::
    strapy.components.BeamSplitter
    strapy.components.Dump
    strapy.components.FaradayRotator
    strapy.components.IdealIsolator
    strapy.components.Mirror
    strapy.components.Polariser
    strapy.components.PolarisingBeamSplitter
    strapy.components.Reflector
    strapy.components.Source
    strapy.components.Stack
    strapy.components.Waveplate

Beam splitter
------
.. autoclass:: strapy.components.BeamSplitter
    :members:

Dump
------
.. autoclass:: strapy.components.Dump
    :members:

Faraday rotator
------
.. autoclass:: strapy.components.FaradayRotator
    :members:

Ideal isolator
------
.. autoclass:: strapy.components.IdealIsolator
    :members:

Mirror
------
.. autoclass:: strapy.components.Mirror
    :members:

Polariser
------
.. autoclass:: strapy.components.Polariser
    :members:

Polarising beam splitter
------
.. autoclass:: strapy.components.PolarisingBeamSplitter
    :members:

Reflector
------
.. autoclass:: strapy.components.Reflector
    :members:

Source
------
.. autoclass:: strapy.components.Source
    :members:

Stack
------
.. autoclass:: strapy.components.Stack
    :members:

Waveplate
------
.. autoclass:: strapy.components.Waveplate
    :members:

Utility functions
-----------------
.. automethod:: strapy.components.rotationMatrix44
.. automethod:: strapy.components.rotationMatrix88

