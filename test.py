import cell_pattern
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection

N0 = 1
N = 32
R = 2.0
nout = 200
timesteps = 1000
kdash = 5.0/1.0
k = 1.5
T = 0.2
A = 1.0
p = 0.001
potential_well_scaling = 10


def plot(filename, sim):
    plt.clf()
    xy = sim.get_position()
    u = sim.get_orientation()
    sigma_s = sim.get_sigma_s()
    k = sim.get_k()
    n = sim.size()

    ww = np.ones(n)*sigma_s*k
    hh = np.ones(n)*sigma_s
    aa = np.arctan2(u[:, 1], u[:, 0])*360/(2*3.14)

    ax = plt.gca()
    ec = EllipseCollection(ww, hh, aa, units='x', offsets=xy, transOffset=ax.transData)
    ax.add_collection(ec)
    ax.autoscale_view()
    ax.quiver(xy[:, 0], xy[:, 1], u[:, 0], u[:, 1])
    ax.set_xlabel('X')
    ax.set_ylabel('y')
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    plt.savefig(filename)


sim = cell_pattern.Simulation()
sim.set_num_internal(N0)
sim.set_max_internal(N)
sim.set_proliferation_rate(p)
sim.set_potential_well_scaling(potential_well_scaling)
sim.set_temperature(T)
sim.set_aging(A)
sim.set_circle_radius(R)
sim.set_kdash(kdash)
sim.set_k(k)
sim.initialise()
plot("initialise.pdf", sim)

fig, ax = plt.subplots(figsize=(10, 10))
for i in range(nout):
    sim.integrate(timesteps)
    plot("snapshot%05d.png" % i, sim)
