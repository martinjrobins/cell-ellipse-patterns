#include "Aboria.h"
#include "GayBernePotential.h"

using namespace Aboria;

ABORIA_VARIABLE(orientation, vdouble2, "orientation");
ABORIA_VARIABLE(fixed, uint8_t, "fixed");

using Particles_t = Particles<std::tuple<orientation, fixed>, 2>;
using position = Particles_t::position;

/*
class LatticeSimulation {
  static constexpr double pi = 3.14;
  double m_dx{0.5};
  double m_rotate_scale{2.0 * pi / 25.0};
  double m_temperature{1.0};
  double m_circle_radius{5.0};
  std::default_random_engine m_gen;

public:
  Particles_t particles;

  void initialise() {
    int num_square = std::pow(2 * m_circle_radius, 2);
    particles.resize(num_square);
    GayBernePotential potential(m_sigma_s, m_k, m_kdash, m_mu, m_nu,
                                m_epsilon_0);

    void set_rotate_scale(const double arg) { m_rotate_scale = arg; }
    void set_temperature(const double arg) { m_temperature = arg; }
    void seed(const int arg) { m_gen.seed(arg); }

    void integrate(const int timesteps) {
      GayBernePotential potential(m_sigma_s, m_k, m_kdash, m_mu, m_nu,
                                  m_epsilon_0);
      std::uniform_real_distribution<double> uniform(0, 1);
      const double cut_off2 = std::pow(potential.cut_off(), 2);
      int accepts = 0;
      const double circle_radius2 = std::pow(m_circle_radius, 2);
      for (int ts = 0; ts < timesteps; ++ts) {
        for (int ii = m_num_external; ii < particles.size(); ++ii) {
          // pick a random particle
          const int index =
              uniform(m_gen) * (particles.size() - m_num_external) +
              m_num_external;

          auto i = particles[index];
          const auto ri = get<position>(i);
          const auto ui = get<orientation>(i);

          // get a candidate position and orientation
          const auto rand_inc =
              m_translate_scale *
              vdouble2(uniform(m_gen) - 0.5, uniform(m_gen) - 0.5);
          const auto rand_inc2 = m_rotate_scale * (uniform(m_gen) - 0.5);
          auto candidate_pos = ri + rand_inc;
          const double theta = std::atan2(ui[1], ui[0]);
          const double candidate_theta = theta + rand_inc2;

          const auto candidate_u =
              vdouble2(std::cos(candidate_theta), std::sin(candidate_theta));

          auto canditate_radius2 = candidate_pos.squaredNorm();

          // if candidate outside domain reject
          if (candidate_pos.squaredNorm() > circle_radius2) {
            continue;
          }

          // calculate the difference in potential
          double Udiff = 0;
          for (int jj = 0; jj < particles.size(); ++jj) {
            const auto rj = get<position>(particles)[jj];
            const auto uj = get<orientation>(particles)[jj];
            if ((ri - rj).squaredNorm() > cut_off2) {
              continue;
            }
            if (get<id>(particles)[jj] == get<id>(i)) {
              continue;
            }
            Udiff -= potential(ri, ui, rj, uj);
            Udiff += potential(candidate_pos, candidate_u, rj, uj);
          }

          // either accept or reject step
          const double acceptance_ratio = std::exp(-Udiff / m_temperature);
          if ((Udiff <= 0) || (uniform(m_gen) < acceptance_ratio)) {
            // std::cout << "accepted" << std::endl;
            get<position>(i) = candidate_pos;
            get<orientation>(i) = candidate_u;
            accepts++;
          } else {
            // std::cout << "rejected move from r =" << ri << " to " <<
            // candidate_pos
            //          << " and u = " << ui << " to " << candidate_u <<
            //          std::endl;
          }
        }
      }

      std::cout << "finished monte carlo steps ratio of accepts to total is "
                << static_cast<double>(accepts) / (timesteps * particles.size())
                << std::endl;
    }
  };
*/

class Simulation {
  static constexpr double pi = 3.14;
  double m_epsilon_0{1.0};
  double m_sigma_s{0.5};
  double m_k{3.0};
  double m_kdash{1.0 / 5.0};
  double m_mu{1.0};
  double m_nu{3.0};
  double epsilon_0{1.0};
  double m_translate_scale{0.5 / 20.0};
  double m_rotate_scale{2.0 * pi / 25.0};
  double m_temperature{1.0};
  double m_circle_radius{5.0};
  double m_proliferation_rate{0.001};
  int m_num_external;
  int m_num_internal;
  std::default_random_engine m_gen;

public:
  Particles_t particles;

  void initialise() {
    m_num_external = std::ceil(2 * pi * m_circle_radius / (m_k * m_sigma_s));
    particles.resize(m_num_internal + m_num_external);
    GayBernePotential potential(m_sigma_s, m_k, m_kdash, m_mu, m_nu,
                                m_epsilon_0);

    for (int i = 0; i < m_num_external; ++i) {
      const double theta = i * 2.0 * pi / m_num_external;
      get<position>(particles)[i] =
          m_circle_radius * vdouble2(std::cos(theta), std::sin(theta));
      get<orientation>(particles)[i] =
          vdouble2(-std::sin(theta), std::cos(theta));
      get<fixed>(particles)[i] = true;
    }
    std::uniform_real_distribution<double> uniform(-1, 1);
    for (int i = m_num_external; i < particles.size(); ++i) {
      bool overlap;
      do {
        do {
          get<position>(particles)[i] =
              m_circle_radius * vdouble2(uniform(m_gen), uniform(m_gen));
        } while (get<position>(particles)[i].norm() > m_circle_radius);

        const double theta = 2 * pi * uniform(m_gen);
        get<orientation>(particles)[i] =
            vdouble2(std::cos(theta), std::sin(theta));

        overlap = false;

        for (int j = 0; j < i; ++j) {
          auto dx = get<position>(particles)[j] - get<position>(particles)[i];
          double r = dx.norm();
          auto rhat = dx / r;
          auto u1 = get<orientation>(particles)[i];
          auto u2 = get<orientation>(particles)[j];
          double udotu = u1.dot(u2);
          double u1dotr = u1.dot(rhat);
          double u2dotr = u2.dot(rhat);
          if (1.5 * r < potential.sigma(udotu, u1dotr, u2dotr)) {
            overlap = true;
          }
        }
      } while (overlap);

      get<fixed>(particles)[i] = false;
    }
  }

  void set_proliferation_rate(const double arg) { m_proliferation_rate = arg; }
  void set_circle_radius(const double arg) { m_circle_radius = arg; }
  void set_sigma_s(const double arg) { m_sigma_s = arg; }
  double get_sigma_s() { return m_sigma_s; }
  void set_k(const double arg) { m_k = arg; }
  double get_k() { return m_k; }
  void set_kdash(const double arg) { m_kdash = arg; }
  void set_mu(const double arg) { m_mu = arg; }
  void set_nu(const double arg) { m_nu = arg; }
  void set_num_internal(const int arg) { m_num_internal = arg; }
  void set_translate_scale(const double arg) { m_translate_scale = arg; }
  void set_rotate_scale(const double arg) { m_rotate_scale = arg; }
  void set_temperature(const double arg) { m_temperature = arg; }
  void seed(const int arg) { m_gen.seed(arg); }

  void integrate(const int timesteps) {
    GayBernePotential potential(m_sigma_s, m_k, m_kdash, m_mu, m_nu,
                                m_epsilon_0);
    std::uniform_real_distribution<double> uniform(0, 1);
    const double cut_off2 = std::pow(potential.cut_off(), 2);
    int accepts = 0;
    const double circle_radius2 = std::pow(m_circle_radius, 2);
    for (int ts = 0; ts < timesteps; ++ts) {
      for (int ii = m_num_external; ii < particles.size(); ++ii) {
        // pick a random particle
        const int index = uniform(m_gen) * (particles.size() - m_num_external) +
                          m_num_external;

        auto i = particles[index];
        const auto ri = get<position>(i);
        const auto ui = get<orientation>(i);

        // check if this cell will divide, otherwise do a move
        std::poisson_distribution<> poisson(m_proliferation_rate /
                                            particles.size());
        if (poisson(m_gen) > 0) {
          particles.push_back(ri + ui * m_sigma_s * m_k / 2);
          get<orientation>(particles[particles.size() - 1]) = ui;
        } else {
          // get a candidate position and orientation
          const auto rand_inc =
              m_translate_scale *
              vdouble2(uniform(m_gen) - 0.5, uniform(m_gen) - 0.5);
          const auto rand_inc2 = m_rotate_scale * (uniform(m_gen) - 0.5);
          auto candidate_pos = ri + rand_inc;
          const double theta = std::atan2(ui[1], ui[0]);
          const double candidate_theta = theta + rand_inc2;

          const auto candidate_u =
              vdouble2(std::cos(candidate_theta), std::sin(candidate_theta));

          auto canditate_radius2 = candidate_pos.squaredNorm();

          // if candidate outside domain reject
          if (candidate_pos.squaredNorm() > circle_radius2) {
            continue;
          }

          // calculate the difference in potential
          double Udiff = 0;
          for (int jj = 0; jj < particles.size(); ++jj) {
            const auto rj = get<position>(particles)[jj];
            const auto uj = get<orientation>(particles)[jj];
            if ((ri - rj).squaredNorm() > cut_off2) {
              continue;
            }
            if (get<id>(particles)[jj] == get<id>(i)) {
              continue;
            }
            Udiff -= potential(ri, ui, rj, uj);
            Udiff += potential(candidate_pos, candidate_u, rj, uj);
          }

          // either accept or reject step
          const double acceptance_ratio = std::exp(-Udiff / m_temperature);
          if ((Udiff <= 0) || (uniform(m_gen) < acceptance_ratio)) {
            // std::cout << "accepted" << std::endl;
            get<position>(i) = candidate_pos;
            get<orientation>(i) = candidate_u;
            accepts++;
          } else {
            // std::cout << "rejected move from r =" << ri << " to " <<
            // candidate_pos
            //          << " and u = " << ui << " to " << candidate_u <<
            //          std::endl;
          }
        }
      }
    }

    std::cout << "finished monte carlo steps ratio of accepts to total is "
              << static_cast<double>(accepts) / (timesteps * particles.size())
              << std::endl;
  }
};

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(cell_pattern, m) {
  py::class_<Simulation>(m, "Simulation")
      .def(py::init<>())
      .def("get_position",
           [](Simulation &sim) {
             return py::array_t<double>(
                 std::vector<ptrdiff_t>{
                     static_cast<ptrdiff_t>(sim.particles.size()), 2},
                 reinterpret_cast<double *>(
                     get<position>(sim.particles).data()));
           })
      .def("get_orientation",
           [](Simulation &sim) {
             return py::array_t<double>(
                 std::vector<ptrdiff_t>{
                     static_cast<ptrdiff_t>(sim.particles.size()), 2},
                 reinterpret_cast<double *>(
                     get<orientation>(sim.particles).data()));
           })
      .def("get_fixed",
           [](Simulation &sim) {
             return py::array_t<uint8_t>(
                 std::vector<ptrdiff_t>{
                     static_cast<ptrdiff_t>(sim.particles.size())},
                 reinterpret_cast<uint8_t *>(get<fixed>(sim.particles).data()));
           })
      .def("size", [](Simulation &sim) { return sim.particles.size(); })
      .def("set_sigma_s", &Simulation::set_sigma_s)
      .def("set_circle_radius", &Simulation::set_circle_radius)
      .def("get_sigma_s", &Simulation::get_sigma_s)
      .def("set_k", &Simulation::set_k)
      .def("get_k", &Simulation::get_k)
      .def("set_kdash", &Simulation::set_kdash)
      .def("set_mu", &Simulation::set_mu)
      .def("set_nu", &Simulation::set_nu)
      .def("set_num_internal", &Simulation::set_num_internal)
      .def("set_translate_scale", &Simulation::set_translate_scale)
      .def("set_rotate_scale", &Simulation::set_rotate_scale)
      .def("set_temperature", &Simulation::set_temperature)
      .def("set_proliferation_rate", &Simulation::set_proliferation_rate)
      .def("seed", &Simulation::seed)
      .def("initialise", &Simulation::initialise)
      .def("integrate", &Simulation::integrate);
}
