import strapy as ts
import numpy as np
from tabulate import tabulate

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

model.components['pbs'].rExtinction = 0.02
model.components['pbs'].tExtinction = 0.01
model.components['pbs'].rLoss = 0
model.components['pbs'].tLoss = 0
model.components['pbs'].theta0 = 0
model.components['pbs'].theta1 = 0
model.components['pbs'].update()

model.build()

nRuns = 10000
set_times = np.empty((nRuns,), dtype=float)
solve_times = np.empty((nRuns,), dtype=float)
detector_times = np.empty((nRuns,), dtype=float)

for i in range(nRuns):
    model.components['sIn'].set_length(i / nRuns)
    set_times[i], solve_times[i], detector_times[i] = \
        model.evaluate(timing=True)

print()
print(
    tabulate(
        [['set time', np.average(set_times),
          np.std(set_times) / np.sqrt(nRuns)],
         ['solve time', np.average(solve_times),
          np.std(solve_times) / np.sqrt(nRuns)],
         ['detector time', np.average(detector_times),
          np.std(detector_times) / np.sqrt(nRuns)]],
        headers=['', 'average', 'standard error']))
print()
