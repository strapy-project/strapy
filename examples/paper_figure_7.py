###
#   This script produces Figure 7 in "Polarisation-sensitive transfer matrix
#   modelling for displacement measuring interferometry", A. Bridges, A. Yacoot,
#   T. Kissinger and R. P. Tatam.
#
#   Last updated by Angus Bridges, 24/03/2020.
###

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import scampy as ts
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
model.components['qwpCos'].rotation = 0#20*np.pi/180
model.components['qwpCos'].update()

model.components['polCos'].rotation = np.pi/4
model.components['polCos'].update()

model.components['polSin'].rotation = np.pi/4
model.components['polSin'].update()

stack = pyctmm.create_stack(2, model.wavelength, 0)
pyctmm.set_ind(stack, 0, 3, 0)
pyctmm.set_ind(stack, 1, 1., 0)
pyctmm.set_d(stack, 0, 0)
pyctmm.set_d(stack, 1, 0)

model.components['sMesA'].set_pyctmm(stack)

nPoints = 100
xs = np.linspace(0, 1, nPoints)
ints1 = np.empty(xs.shape, dtype=float)
ints2 = np.empty(xs.shape, dtype=float)
ints1analytic = np.empty(xs.shape, dtype=float)
ints2analytic = np.empty(xs.shape, dtype=float)

iAir = 1
iGlass = 3

rag = (iAir - iGlass)/(iGlass + iAir)
rga = (iGlass - iAir)/(iGlass + iAir)
tga = (2*iGlass)/(iGlass + iAir)
tag = (2*iAir)/(iGlass + iAir)

k = 2*np.pi/633e-9
xOff = 0

for i, x in enumerate(xs):
    fOutS = (1/np.sqrt(2))*tga*tag*(np.exp(1j*2*k*(x + xOff)*633e-9) - \
        (rag*rag*np.exp(1j*6*k*(x + xOff)*633e-9)) \
            /(rag*rag*np.exp(1j*4*k*(x + xOff)*633e-9) - 1))
    fOutP = (1/np.sqrt(2))*np.exp(1j*2*k*0*633e-9)

    fOutS1st = (1/np.sqrt(2))*tga*tag*(np.exp(1j*2*k*(x + xOff)*633e-9) + \
        rag*rag*np.exp(1j*k*(4*0 + 6)*(x + xOff)*633e-9))

    model.components['sMesB'].set_length(-x)
    model.evaluate()
    ints2[i] = model.detectors['pd1'].intensity
    ints1[i] = model.detectors['pd2'].intensity

    ints1analytic[i] = np.abs(0.5*fOutS1st + 0.5*fOutP*np.exp(1j*np.pi/2))**2
    ints2analytic[i] = np.abs(0.5*fOutS1st + 0.5*fOutP)**2

fig = plt.figure(figsize=(6, 2/0.75))
gs = fig.add_gridspec(1, 3)

red = (1, 0.5, 0.5)

ax0 = fig.add_subplot(gs[0,:2])
ax0.plot(xs, ints1, label='PD1 (model)', color='k')
ax0.plot(xs, ints2, label='PD2 (model)', color='k', ls='--')
ax0.plot(xs, ints1analytic, label='PD1 (1^st order)', color=red, ls=':')
ax0.plot(xs, ints2analytic, label='PD2 (1^st order)', color=red, ls='-.')
ax0.set_xlabel('Displacement (wavelengths)')
ax0.set_ylabel('Intensity')
ax0.set_yticks([0, 0.25, 0.5])

ax1 = fig.add_subplot(gs[0, 2])
ax1.plot(ints1, ints2, color='k', label='model')
ax1.plot(ints1analytic[:int(len(ints1analytic)/2)], \
    ints2analytic[:int(len(ints1analytic)/2)], color=red, ls='--', \
        label='1^st order')
ax1.set_aspect('equal')
ax1.set_xlabel('PD1 intensity')
ax1.set_ylabel('PD2 intensity')
ax1.set_xticks([0, 0.25, 0.5])
ax1.set_yticks([0, 0.25, 0.5])

plt.tight_layout(rect=(0,0,1,0.75))
lgd = ax0.legend(ncol=2, bbox_to_anchor=(1.01,1.6), loc='upper right')
ax1.legend(ncol=1, bbox_to_anchor=(1.01,1.6), loc='upper right')
plt.show()