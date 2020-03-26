import numpy as np
import strapy
import pyctmm

try:
    from matplotlib import pyplot as plt
except:
    print("\nmatplotlib is not installed, install now (y/n)?")
    installFlag = input()

    if installFlag == 'y':
        import subprocess
        import sys
        subprocess.check_call([sys.executable, "-m", "pip", "install",
            "matplotlib"])
    try:
        from matplotlib import pyplot as plt
    except:
        print("matplotlib is still not installed, exiting.\n")
        exit()

model = strapy.Model()

model.add_component(strapy.components.Source, 'laser', 'n0')
model.add_component(strapy.components.BeamSplitter, 'BS', ('n1', 'n2', 'n3', 'n4'))
model.add_component(strapy.components.Mirror, 'Mmes', 'n5')
model.add_component(strapy.components.Mirror, 'Mref', 'n6')
model.add_component(strapy.components.Polariser, 'pol', ('n7', 'n8'))
model.add_component(strapy.components.Stack, 's01', ('n0', 'n1'))
model.add_component(strapy.components.Stack, 's35', ('n3', 'n5'))
model.add_component(strapy.components.Stack, 's26', ('n2', 'n6'))
model.add_component(strapy.components.Stack, 's47', ('n4', 'n7'))
model.add_component(strapy.components.Stack, 's89', ('n8', 'n9'))
model.add_component(strapy.components.Dump, 'd9', 'n9')

model.add_detector('PD1', 'n9', ('amplitude', 'intensity'))

model.components['laser'].amplitude[0] = 0
model.components['laser'].amplitude[1] = 1

model.components['BS'].rP = np.sqrt(0.5)
model.components['BS'].rS = np.sqrt(0.5)
model.components['BS'].tP = np.sqrt(0.5)
model.components['BS'].tS = np.sqrt(0.5)

model.components['Mmes'].rP = 1
model.components['Mmes'].rS = 1

model.components['Mref'].rP = 1
model.components['Mref'].rS = 1

model.components['pol'].rotation = 0
model.components['pol'].extinction = 0
model.components['pol'].loss = 0

model.components['pol'].update()

stack = pyctmm.create_stack(3, 633e-9, 0)

pyctmm.set_ind(stack, 0, 1, 0)
pyctmm.set_ind(stack, 1, 0.2, -3)
pyctmm.set_ind(stack, 2, 1, 0)

pyctmm.set_d(stack, 0, 0)
pyctmm.set_d(stack, 1, 10e-9)
pyctmm.set_d(stack, 2, 0)

model.components['s47'].set_pyctmm(stack)

model.build()

xs = np.linspace(0, 1, 100)
ys = np.linspace(0, 1, 100)
ints = np.empty((len(xs), len(ys)), dtype=float)

for i in range(len(xs)):
    for j in range(len(ys)):
        model.components['s35'].set_length(xs[i])
        model.components['s26'].set_length(ys[j])
        model.evaluate()

        ints[i, j] = model.detectors['PD1'].intensity

ints = (ints - np.min(ints))/(np.max(ints) - np.min(ints))

plt.imshow(ints, extent=[0,1,0,1])
cbar = plt.colorbar()
cbar.set_label("Normalised intensity")
plt.xlabel('$M_{mes}$ displacement (wavelengths)')
plt.ylabel('$M_{ref}$ displacement (wavelengths)')
plt.savefig("michelson_cavity_output.svg")
plt.show()

