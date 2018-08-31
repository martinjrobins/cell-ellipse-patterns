import cell_pattern
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection

N = 30
R = 2.6
nout = 100
timesteps = 100
kdash = 1.0/5.0
k = 2.0
T = 1.0


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
    ax.set_xlabel('X')
    ax.set_ylabel('y')
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    plt.savefig(filename)


sim = cell_pattern.Simulation()
sim.set_num_internal(N)
sim.set_temperature(T)
sim.set_circle_radius(R)
sim.set_kdash(kdash)
sim.set_k(k)
sim.initialise()
plot("initialise.pdf", sim)

fig, ax = plt.subplots(figsize=(10, 10))
for i in range(nout):
    sim.integrate(timesteps)
    plot("snapshot%05d.png" % i, sim)
