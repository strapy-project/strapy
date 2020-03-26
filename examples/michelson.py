import numpy as np
import strapy

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

model.build()

xs = np.linspace(0, 1, 100)
ints = np.empty(xs.shape, dtype=float)

for i, x in enumerate(xs):
    model.components['s35'].set_length(x)
    model.evaluate()

    ints[i] = model.detectors['PD1'].intensity

plt.plot(xs, ints)
plt.xlabel('Displacement (wavelengths)')
plt.ylabel('Intensity at detector')
plt.savefig("michelson_output.svg")
plt.show()

