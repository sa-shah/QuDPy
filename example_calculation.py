from qutip import *
import numpy as np
from qudpy.Classes import *
import qudpy.plot_functions as pf
import ufss  # diagram generation

# Setting up the required double sided diagrams for tests
# DiagramGenerator class, or DG for short
DG = ufss.DiagramGenerator
# initialize the module
R3rd = DG()  # DG takes a single key-word argument, which has the default value detection_type = 'polarization'
# DiagramAutomation needs to know the phase-matching/-cycling condition
R3rd.set_phase_discrimination([(0, 1), (1, 0), (1, 0)])  # setting condition for rephasing diagrams R1,2,3
# Set the pulse durations for pulses 0,1,2 and the local oscillator
d0 = 2
d1 = 4
d2 = 4
dlo = 6
# Set the pulse intervals
t0 = np.array([-d0 / 2, d0 / 2])
t1 = np.array([-d1 / 2, d1 / 2])
t2 = np.array([-d2 / 2, d2 / 2])
tlo = np.array([-dlo / 2, dlo / 2])
all_pulse_intervals = [t0, t1, t2, tlo]

# These pulse intervals are given to DG by setting DG's efield_times attribute
R3rd.efield_times = all_pulse_intervals
time_ordered_diagrams_rephasing = R3rd.get_diagrams([0, 100, 200, 200])
[R3, R1, R2] = time_ordered_diagrams_rephasing
rephasing = [R1, R2, R3]
print('the rephasing diagrams are R1, R2 and R3 ')
# R3rd.display_diagrams(rephasing)

sys1 = System()
sys1.c_ops = []  # [sys1.a*0.5]

# Uncomment the line below if you need to check the a single trial of the system evolution
# states = sys1.diagram_donkey([0, 5, 10, 20], rephasing, r=10)

time_delays = [10, 5, 10]
scan_id = [0, 2]
response_list = []

for k in range(3):
    states, t1, t2, dipole = sys1.coherence2d(time_delays, rephasing[k], scan_id, parallel=True)
    response_list.append(np.real(dipole))

# Uncomment the lines below to plot the time domain dipole signal
# pf.multiplot(response_list, [0, 10, 0, 15], ['t1', 't2'], ['R1', 'R2', 'R3 imag', 'R3 real'], 'log')

spectra_list, extent, f1, f2 = sys1.spectra(response_list)
pf.multiplot(np.real(spectra_list), extent, ['w1', 'w2'], ['R1', 'R2', 'R3'], 'log')
pf.plot(np.real(spectra_list[0]), extent, ['w1', 'w2'], 'R1', 'log')

# Test2: As a second test consider the interaction given in lines below
# time_delays = [0, 0, 10, 0, 0, 5, 0, 0, 15]
# Rtest = (('Bu', 0), ('Ku', 0), ('Kd', 0), ('Ku', 1), ('Bd', 1), ('Bu', 1), ('Bd', 2), ('Kd', 2), ('Ku', 2))
# scan_id = [2, 8]
