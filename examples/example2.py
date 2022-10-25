from qutip import *
import numpy as np
from qudpy.Classes import *
import qudpy.plot_functions as pf
import ufss  # diagram generation

w = 1.0
w0 = 1.0
g = 1  # initial coupling strength for cavity mode-spin interaction for demonstration purposes only
gc = np.sqrt(w * w0)/2  # critical coupling strength

kappa = 0.05
gamma = 0.15
M = 18    # number of fock basis for cavity mode. Use larger value for stronger couplings
N = 6    # number of spins
j = N/2  # Total J for spins.
n = N+1  # dimensionality of total spin operator.

print("critical coupling strength = ", gc)

a = tensor(destroy(M), qeye(n))
Jp = tensor(qeye(M), jmat(j, '+'))
Jm = tensor(qeye(M), jmat(j, '-'))
Jz = tensor(qeye(M), jmat(j, 'z'))

H0 = w * a.dag() * a + w0 * Jz  # default basic Hamiltonian
H1 = (a + a.dag()) * (Jp + Jm)/np.sqrt(N)  # intra-cavity interaction term
H = H0 + g * H1

print("dimensionality of Hilbert-space: ", H.shape)

# average number thermal photons in the bath coupling to the resonator
n_th = 0.25
c_ops = [np.sqrt(kappa * (n_th + 1)) * a, np.sqrt(kappa * n_th) * a.dag(),np.sqrt(gamma)*Jz]

#c_ops = [sqrt(kappa) * a, sqrt(gamma) * Jm]

g_vec = np.linspace(0.01, 1, 20)

# Ground state for the Hamiltonian: H = H0 + g * H1
rho_ss_list = [steadystate(H0 + g * H1, c_ops) for g in g_vec]

# calculate the expectation value of the number of photons in the cavity
n_ss_vec = expect(a.dag() * a, rho_ss_list)
n2_ss_vec = expect(a.dag() * a*a.dag() * a, rho_ss_list)
Jz_ss_vec = expect(Jz, rho_ss_list)

a_ss_vec = expect(a, rho_ss_list)



fig, axes = plt.subplots(1, 3, sharex=True, figsize=(12,4))

axes[0].plot(g_vec, n_ss_vec, 'r', linewidth=2, label="cavity occupation")
axes[0].set_ylim(0, max(n_ss_vec))
axes[0].set_ylabel("<n>", fontsize=16)
axes[0].set_xlabel("g", fontsize=16)

axes[1].plot(g_vec, Jz_ss_vec, 'r', linewidth=2, label="<Jz>")
axes[1].set_ylim(-j, j)
axes[1].set_ylabel(r"$\langle J_z\rangle$", fontsize=16)
axes[1].set_xlabel("g", fontsize=16)

axes[2].plot(g_vec, abs(a_ss_vec), 'r', linewidth=2, label="<a>")

fig.tight_layout()

rho_ss_sublist = rho_ss_list[::5]

xvec = np.linspace(-3, 3, 200)

fig_grid = (3, len(rho_ss_sublist))
fig = plt.figure(figsize=(3 * len(rho_ss_sublist), 9))

for idx, rho_ss in enumerate(rho_ss_sublist):
    # trace out the cavity density matrix
    rho_ss_cavity = ptrace(rho_ss, 0)

    # calculate its wigner function
    W = wigner(rho_ss_cavity, xvec, xvec)

    # plot its wigner function
    ax = plt.subplot2grid(fig_grid, (0, idx))
    ax.contourf(xvec, xvec, W, 100)

    # plot its fock-state distribution
    ax = plt.subplot2grid(fig_grid, (1, idx))
    ax.bar(np.arange(0, M), np.real(rho_ss_cavity.diag()), color="blue", alpha=0.6)
    ax.set_ylim(0, 1)
    ax.set_xlim(0, M)


# plot the cavity occupation probability in the ground state
ax = plt.subplot2grid(fig_grid, (2, 0), colspan=fig_grid[1])
ax.plot(g_vec, n_ss_vec, 'r', linewidth=2, label="cavity occupation")
ax.set_xlim(0, max(g_vec))
ax.set_ylim(0, max(n_ss_vec) * 1.2)
ax.set_ylabel("Cavity gnd occ. prob.", fontsize=16)
ax.set_xlabel("interaction strength", fontsize=16)

for g in g_vec[::4]:
    ax.plot([g, g], [0, max(n_ss_vec) * 1.2], 'b:', linewidth=2.5)


#Setting up the required double sided diagrams for tests
# DiagramGenerator class, or DG for short
DG = ufss.DiagramGenerator
# initialize the module
R3rd = DG()  # DG takes a single key-word argument, which has the default value detection_type = 'polarization'
# DiagramAutomation needs to know the phase-matching/-cycling condition
R3rd.set_phase_discrimination([(0, 1), (1, 0), (1, 0)])  # setting phase-matching condition for rephasing diagrams R1,2,3
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
# Creating diagrams for pulse arrival times 0, 100, 200 and detection time 300.
[R3, R1, R2] = R3rd.get_diagrams([0, 100, 200, 300])
rephasing = [R1, R2, R3]
print('the rephasing diagrams are R1, R2 and R3 ', rephasing)

# setting conditions for and generating non-rephasing diagrams R4, R5 and R6
R3rd.set_phase_discrimination([(1, 0), (0, 1), (1, 0)])
[R6, R4, R5] = R3rd.get_diagrams([0, 100, 200, 200])
nonrephasing = [R4, R5, R6]
print('the non-rephasing diagrams are R4, R5 and R6', nonrephasing)


w  = 1.0
w0 = 1.0


gc = np.sqrt(w * w0)/2  # critical coupling strength. For this case it is 1/2
g = 0.5  # coupling strength. For this case it is 50% of the critical coupling
kappa = 0.05
gamma = 0.15
#M = 10
#N = 2
j = N/2    # unfort. this is a real.
n = N+1    # needs to be integer for qeye()

print("critical coupling strength = ",gc)

a = tensor(destroy(M), qeye(n))
Jp = tensor(qeye(M), jmat(j, '+'))
Jm = tensor(qeye(M), jmat(j, '-'))
Jz = tensor(qeye(M), jmat(j, 'z'))

H0 = w * a.dag() * a + w0 * Jz
H1 = (a + a.dag()) * (Jp + Jm)/np.sqrt(N)
H = H0 + g * H1

print("dimensionality of Hilbert-space: ",H.shape)

en,T = H.eigenstates()

Hd=H.transform(T)

# average number thermal photons in the bath coupling to the resonator
n_th = 0.25
c1 = np.sqrt(kappa * (n_th + 1)) * a
c2 = np.sqrt(kappa * n_th) * a.dag()
c3 = np.sqrt(gamma)*Jz



c1d =c1.transform(T)
c2d =c2.transform(T)
c3d =c3.transform(T)


# A and A.dag() are now the excitation/de-excitation operators due to the coupling to the external photon-field
# and we need x-form this to the eigenbasis
A = a.transform(T)

mud = A + A.dag()



print("finding steady-states of the DM")
#define density matrix as steady-state in the eigenbasis
rhoSS=steadystate(Hd, [c1d,c2d,c3d])


# properties of the SS system (need to use transformed operators )
n_ss = expect(A.dag() * A, rhoSS)
n2_ss = expect(A.dag()*A*A.dag()*A,rhoSS)

print("Steady State Properties")
print("<n> = ",n_ss)
print("<n^2> = ",n2_ss)
print("<(n-<n>)^2> = ",n2_ss-n_ss*n_ss)


sys2 = System(H=Hd, a=A, u=mud, c_ops=[c1d, c2d, c3d], rho=rhoSS, diagonalize=False)
sys2.hbar=1
print("system has been intialized")

states = sys2.diagram_donkey([0, 75, 80, 155], [R1], r=10)

print("finished setup stage")


# generating 2Dcoherence response for rephasing diagrams
time_delays = [100, 10, 100]
scan_id = [0, 2]
response_list = []
states_list = []
diagrams = rephasing+nonrephasing
for k in range(6):
    states, t1, t2, dipole = sys2.coherence2d(time_delays, diagrams[k], scan_id, r=1/2, parallel=True)
    print('diagram ', k, ' done')
    response_list.append(1j*dipole)
    states_list.append(states)
spectra_list, extent, f1, f2 = sys2.spectra(np.imag(response_list))




qsave(response_list, 'dipole_6spin_6cavity_res_half_range_100_10_100_g_050')
