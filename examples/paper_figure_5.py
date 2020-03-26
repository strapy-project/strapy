###
#   A simple polarising michelson inteferometer model.
#
#   This script produces Figure 5 in "Polarisation-sensitive transfer matrix
#   modelling for displacement measuring interferometry", A. Bridges, A. Yacoot,
#   T. Kissinger and R. P. Tatam.
#
#   Last updated by Angus Bridges, 24/03/2020.
###

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import strapy as ts
import pyctmm

model = ts.Model()

model.add_component(ts.components.Source, 'laser', 'n0')
model.add_component(ts.components.Stack, 'sIn', ('n0', 'n1'))
model.add_component(ts.components.PolarisingBeamSplitter, 'pbs', \
    ('n1', 'n2', 'n3', 'n4'))
model.add_component(ts.components.Stack, 'sRefA', ('n2', 'n5'))
model.add_component(ts.components.Waveplate, 'qwpRef', ('n5', 'n6'))
model.add_component(ts.components.Stack, 'sRefB', ('n6', 'n7'))
model.add_component(ts.components.Mirror, 'mRef', 'n7')
model.add_component(ts.components.Stack, 'sMesA', ('n3', 'n8'))
model.add_component(ts.components.Waveplate, 'qwpMes', ('n8', 'n9'))
model.add_component(ts.components.Stack, 'sMesB', ('n9', 'n10'))
model.add_component(ts.components.Mirror, 'mMes', 'n10')
model.add_component(ts.components.Stack, 'sOutA', ('n4', 'n11'))
model.add_component(ts.components.BeamSplitter, 'npbs', \
    ('n11', 'n12', 'n13', 'nNPBSdumpA'))
model.add_component(ts.components.Stack, 'sNPBSdump', \
    ('nNPBSdumpA', 'nNPBSdumpB'))
model.add_component(ts.components.Dump, 'dNPBS', 'nNPBSdumpB')
model.add_component(ts.components.Stack, 'sCosA', ('n12', 'n14'))
model.add_component(ts.components.Waveplate, 'qwpCos', ('n14', 'n15'))
model.add_component(ts.components.Stack, 'sCosB', ('n15', 'n16'))
model.add_component(ts.components.Polariser, 'polCos', ('n16', 'n17'))
model.add_component(ts.components.Stack, 'sCosC', ('n17', 'n18'))
model.add_component(ts.components.Dump, 'dCos', 'n18')
model.add_component(ts.components.Stack, 'sSinA', ('n13', 'n19'))
model.add_component(ts.components.Polariser, 'polSin', ('n19', 'n20'))
model.add_component(ts.components.Stack, 'sSinB', ('n20', 'n21'))
model.add_component(ts.components.Dump, 'dSin', 'n21')
model.add_detector('pd2', 'n18', properties=('amplitude', 'intensity'))
model.add_detector('pd1', 'n21', properties=('amplitude', 'intensity'))

model.components['laser'].amplitude = [1/np.sqrt(2), 1/np.sqrt(2)]

model.build()

model.components['qwpRef'].retardance = 2*np.pi/4
model.components['qwpRef'].rotation = np.pi/4
model.components['qwpRef'].update()

model.components['qwpMes'].retardance = 2*np.pi/4
model.components['qwpMes'].rotation = np.pi/4
model.components['qwpMes'].update()

model.components['qwpCos'].retardance = 2*np.pi/4
model.components['qwpCos'].rotation = 20*np.pi/180
model.components['qwpCos'].update()

model.components['polCos'].rotation = np.pi/4
model.components['polCos'].update()

model.components['polSin'].rotation = np.pi/4
model.components['polSin'].update()

stack = pyctmm.create_stack(2, model.wavelength, 0)
pyctmm.set_ind(stack, 0, 1, 0)
pyctmm.set_ind(stack, 1, 1, 0)
pyctmm.set_d(stack, 0, 0)
pyctmm.set_d(stack, 1, 0)

model.components['sMesA'].set_pyctmm(stack)

nPoints = 100
xs = np.linspace(0, 1, nPoints)
ints1 = np.empty(xs.shape, dtype=float)
ints2 = np.empty(xs.shape, dtype=float)

for i, x in enumerate(xs):
    model.components['sMesB'].set_length(x)
    model.evaluate()
    ints1[i] = model.detectors['pd1'].intensity
    ints2[i] = model.detectors['pd2'].intensity

fig = plt.figure(figsize=(6, 2))
gs = fig.add_gridspec(1, 3)

ax0 = fig.add_subplot(gs[0,:2])
ax0.plot(xs, ints1, label='PD1', color='k')
ax0.plot(xs, ints2, label='PD2', color='k', ls='--')
lgd = ax0.legend()
ax0.set_xlabel('Displacement (wavelengths)')
ax0.set_ylabel('Intensity')
ax0.set_yticks([0, 0.25, 0.5])

ax1 = fig.add_subplot(gs[0, 2])
ax1.plot(ints1, ints2, color='k')
ax1.set_aspect('equal')
ax1.set_xlabel('PD1 intensity')
ax1.set_ylabel('PD2 intensity')
ax1.set_xticks([0, 0.25, 0.5])
ax1.set_yticks([0, 0.25, 0.5])
plt.tight_layout()
plt.show()