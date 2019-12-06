#include "PEP.h"
#include "latticeExporter.h"
int main(int argc, char* argv[]) {
  {
    Ranq2 myran(2);
    int Lx = 400;
    int Ly = Lx;
    int Lz = Lx;
    double phi = 0.01;
    double alpha = 0.01;
    int n_step = 1000;
    std::vector<Par_3> par_arr;
    HyperCubicLattice lattice(Lx, Ly, Lz);
    lattice.create_particles_randomly(par_arr, myran, phi);
    std::vector<size_t> seq;
    seq.reserve(par_arr.size());
    for (int i = 0; i < par_arr.size(); i++) {
      seq.push_back(i);
    }

    char logfile[100];
    char xyfile[100];
    snprintf(logfile, 100, "%d_%d.log", Lx, BLOCK_SIZE_X);

    snprintf(xyfile, 100, "%d_%d.extxyz", Lx, BLOCK_SIZE_X);
    LogExporter log(logfile, 0, n_step, 1000, par_arr.size());
    XYExporter xy(xyfile, 0, n_step, 1000, Lx, Ly);

    xy.dump(0, par_arr);
    unsigned long long beta = static_cast<unsigned long long>(alpha * UINT64_MAX);
    for (int i = 1; i <= n_step; i++) {
      shuffle(seq, myran);
      update_in_rand_seq(par_arr, lattice, myran, seq, beta);
      log.record(i);
      if (i % 100 == 0) {
        std::cout << i << std::endl;
      }
      xy.dump(i, par_arr);
    }
  }

  {
    Ranq2 myran(2);
    int Lx = 400;
    int Ly = Lx;
    int Lz = Lx;
    double phi = 0.01;
    double alpha = 0.01;
    int n_step = 1000;
    std::vector<Particle_3> par_arr;
    HyperCubicLattice lattice(Lx, Ly, Lz);
    lattice.create_particles_randomly(par_arr, myran, phi);
    std::vector<size_t> seq;
    seq.reserve(par_arr.size());
    for (int i = 0; i < par_arr.size(); i++) {
      seq.push_back(i);
    }

    char logfile[100];
    char xyfile[100];
    snprintf(logfile, 100, "%d_%d_2.log", Lx, BLOCK_SIZE_X);

    snprintf(xyfile, 100, "%d_%d_2.extxyz", Lx, BLOCK_SIZE_X);
    LogExporter log(logfile, 0, n_step, 1000, par_arr.size());
    XYExporter xy(xyfile, 0, n_step, 1000, Lx, Ly);

    xy.dump(0, par_arr);
    unsigned long long beta = static_cast<unsigned long long>(alpha * UINT64_MAX);
    for (int i = 1; i <= n_step; i++) {
      shuffle(seq, myran);
      update_in_rand_seq(par_arr, lattice, myran, seq, beta);
      log.record(i);
      if (i % 100 == 0) {
        std::cout << i << std::endl;
      }
      xy.dump(i, par_arr);
    }
  }
}