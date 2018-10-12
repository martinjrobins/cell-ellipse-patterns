import cell_pattern
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection
import multiprocessing


class Param():
    def __init__(self, filename):
        self.N0 = 15
        self.N = 15
        self.R = 4.0
        self.R_init = 2.0
        self.nout = 50
        self.timesteps = 1000
        self.kdash = 5.0/1.0
        self.k = 1.5
        self.T = 20.0
        self.A = 0.001
        self.p = 0.001
        self.potential_well_scaling = 50
        self.orientation_well_scaling = 0
        self.filename = filename


params = [Param('weak'), Param('medium'), Param('strong')]
params[0].potential_well_scaling = 10
params[1].potential_well_scaling = 500
params[2].potential_well_scaling = 2000
#params[0].R = 3.0
#params[0].R_init = 2.0
#params[1].R = 1.6
#params[1].R_init = 1.6
#params[2].R = 1.9
#params[2].R_init = 1.9


def plot(filename, sim, i):
    fig = plt.figure(i, clear=True, figsize=(10, 10))
    xy = sim.get_position()
    u = sim.get_orientation()
    sigma_s = sim.get_sigma_s()
    k = sim.get_k()
    n = sim.size()

    ww = np.ones(n)*sigma_s*k
    hh = np.ones(n)*sigma_s
    aa = np.arctan2(u[:, 1], u[:, 0])*360/(2*3.14)

    ax = fig.subplots()
    ec = EllipseCollection(ww, hh, aa, units='x', offsets=xy, transOffset=ax.transData)
    ax.add_collection(ec)
    ax.autoscale_view()
    ax.quiver(xy[:, 0], xy[:, 1], u[:, 0], u[:, 1])
    ax.set_xlabel('X')
    ax.set_ylabel('y')
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    fig.savefig(filename)


def run(i):
    print('running', i)
    sim = cell_pattern.Simulation()
    sim.set_num_internal(params[i].N0)
    sim.set_max_internal(params[i].N)
    sim.set_proliferation_rate(params[i].p)
    sim.set_potential_well_scaling(params[i].potential_well_scaling)
    sim.set_orientation_well_scaling(params[i].orientation_well_scaling)
    sim.set_temperature(params[i].T)
    sim.set_aging(params[i].A)
    sim.set_circle_radius(params[i].R)
    sim.set_radius_init(params[i].R_init)
    sim.set_kdash(params[i].kdash)
    sim.set_k(params[i].k)
    sim.initialise()
    print('starting plot', i)
    plot(params[i].filename+'_initialise.pdf', sim, i)
    print('ending plot', i)

    #fig, ax = plt.subplots(figsize=(10, 10))
    for j in range(params[i].nout):
        print('\titeration', j, 'of run', i)
        sim.integrate(params[i].timesteps)
        plot(params[i].filename+'_snapshot%05d.png' % j, sim, i)


p = multiprocessing.Pool(len(params))
p.map(run, range(len(params)))
