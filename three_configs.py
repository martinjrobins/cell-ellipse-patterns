import cell_pattern
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection
import multiprocessing
import os


class Param():
    def __init__(self, filename):
        self.N0 = 15
        self.N = 15
        self.R = 4.0
        self.R_init = 2.0
        self.nout = 10
        self.timesteps = 2000
        self.kdash = 5.0/1.0
        self.k = 1.2
        self.T = 20.0
        self.A = 0.001
        self.p = 0.001
        self.seed = 1
        self.potential_well_scaling = 50
        self.orientation_well_scaling = 0
        self.filename = filename


def generate_parameter_sweeps():
    Tf = np.exp(-0.001*(20*1000))*20.0
    print('final Tf = ', Tf)

    params = []

    dirName = 'vary_pot_well'
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    samples = 10

    for pot_well_scale in 2**np.arange(3, 16):
        print('pot_well = ', pot_well_scale)
        for seed in range(samples):
            params.append(Param(dirName + '/vary_pot_well_%d_' % pot_well_scale))
            params[-1].potential_well_scaling = pot_well_scale
            params[-1].seed = seed

    for pot_well_scale in [500.0]:
        dirName = 'vary_temp_pot_%d' % pot_well_scale
        if not os.path.exists(dirName):
            os.mkdir(dirName)
        for max_temp in 2**np.arange(-20.0, 20.0, 3.0):
            print('max_temp = ', max_temp, 'timesteps = ', params[-1].timesteps, ' pot_well_scale = ', pot_well_scale)
            for seed in range(samples):
                params.append(Param(dirName + '/vary_temp_pot_%d_%f_' % (pot_well_scale, max_temp)))
                params[-1].T = max_temp
                params[-1].seed = seed
                params[-1].potential_well_scaling = pot_well_scale
                params[-1].timesteps = int(-(1.0 / (params[-1].A * params[-1].nout)) * np.log(Tf/max_temp))

        dirName = 'vary_polarity_pot_%d' % pot_well_scale
        if not os.path.exists(dirName):
            os.mkdir(dirName)
        for kdash in np.arange(1.0, 25.0, 3.0)/5.0:
            print('kdash= ', kdash, ' pot_well_scale = ', pot_well_scale)
            for seed in range(samples):
                params.append(Param(dirName + '/vary_polarity_pot_%d_%f_' % (pot_well_scale, kdash)))
                params[-1].kdash = kdash
                params[-1].potential_well_scaling = pot_well_scale
                params[-1].seed = seed
    return params


def generate_three_configs():
    params = [Param('weak'), Param('medium'), Param('strong')]
    params[0].potential_well_scaling = 10
    params[1].potential_well_scaling = 500
    params[2].potential_well_scaling = 2000
    return params

# params[0].R = 3.0
# params[0].R_init = 2.0
# params[1].R = 1.6
# params[1].R_init = 1.6
# params[2].R = 1.9
# params[2].R_init = 1.9


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
    fig.savefig(filename+'.png')
    plt.close(fig)

    data = np.concatenate([xy, u], axis=1)
    np.savetxt(filename+'.csv', data, header='x y ux uy')


def plot_multi(filename, sims):
    fig, axs = plt.subplots(1, len(sims), clear=True, figsize=(12, 4))
    for i, sim in enumerate(sims):
        xy = sim.get_position()
        u = sim.get_orientation()
        sigma_s = sim.get_sigma_s()
        k = sim.get_k()
        n = sim.size()

        ww = np.ones(n)*sigma_s*k
        hh = np.ones(n)*sigma_s
        aa = np.arctan2(u[:, 1], u[:, 0])*360/(2*3.14)

        ax = axs[i]
        ec = EllipseCollection(ww, hh, aa, units='x', offsets=xy, transOffset=ax.transData)
        ax.add_collection(ec)
        ax.autoscale_view()
        ax.quiver(xy[:, 0], xy[:, 1], u[:, 0], u[:, 1])
        ax.set_xlabel('X')
        ax.set_ylabel('y')
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)

    axs[0].set_title('weak support\npotential well scaling = %d' % params[0].potential_well_scaling)
    axs[1].set_title('medium support\npotential well scaling = %d' % params[1].potential_well_scaling)
    axs[2].set_title('strong support\npotential well scaling = %d' % params[2].potential_well_scaling)
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
    sim.seed(params[i].seed)
    sim.initialise()
    # print('starting plot', i)
    # plot(params[i].filename+'_initialise.pdf', sim, i)
    # print('ending plot', i)

    # fig, ax = plt.subplots(figsize=(10, 10))
    for j in range(params[i].nout):
        # print('\titeration', j, 'of run', i)
        sim.integrate(params[i].timesteps)
        # plot(params[i].filename+'_snapshot%05d.png' % j, sim, i)

    plot(params[i].filename+'_seed_%d_final' % params[i].seed, sim, i)


def run_all(params):
    p = multiprocessing.Pool(8)
    p.map(run, range(len(params)))


def run_and_plot_multi(params):
    print('running')
    sims = [cell_pattern.Simulation() for i in range(len(params))]
    for i, sim in enumerate(sims):
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
        sim.seed(params[i].seed)
        sim.initialise()
        print('starting plot', i)
        plot(params[i].filename+'_initialise', sim, i)
        print('ending plot', i)

    # fig, ax = plt.subplots(figsize=(10, 10))
    for j in range(params[i].nout):
        print('\titeration', j)
        [sim.integrate(params[i].timesteps) for i, sim in enumerate(sims)]
        plot_multi('all_snapshot%05d.png' % j, sims)


params = generate_parameter_sweeps()
run_all(params)
