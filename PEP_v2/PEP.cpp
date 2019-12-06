#include <iostream>
#include "PEP.h"
#include "latticeExporter.h"

int Par_2::size_ori;
int Par_3::size_ori;
const char Particle_3::ori_arr[6][3] = { {1, 0, 0}, {-1, 0, 0},
                                         {0, 1, 0}, {0, -1, 0},
                                         {0, 0, 1}, {0, 0, -1} };
HyperCubicLattice::HyperCubicLattice(int Lx0, int Ly0, int Lz0) 
  : Lx_(Lx0), Ly_(Ly0), Lz_(Lz0), LxLy_(static_cast<size_t>(Lx_)* Ly_),
    tot_lattices_(LxLy_* Lz_) {
  n_ = new unsigned char[tot_lattices_]();

  size_t nx_b = Lx_ / BLOCK_SIZE_X;
  size_t ny_b = Ly_ / BLOCK_SIZE_Y;
  size_t nz_b = Lz_ / BLOCK_SIZE_Z;
  if (nx_b * BLOCK_SIZE_X < Lx_) {
    nx_b += 1;
  }
  if (ny_b * BLOCK_SIZE_Y < Ly_) {
    ny_b += 1;
  }
  if (nz_b * BLOCK_SIZE_Z < Lz_) {
    nz_b += 1;
  }
  Nx_b_ = nx_b;
  NxNy_b_ = nx_b * ny_b;
}

void test_rand(){
  int Lx = 1000;
  int Ly = 1000;
  double phi = 0.01;
  double alpha = 0.01;
  int n_step = 200000;

  {
    Ran myran(2);

    std::vector<Par_2> par_arr;
    SquareLattice lattice(Lx, Ly);
    lattice.create_particles_randomly(par_arr, myran, phi);
    std::vector<size_t> seq;
    seq.reserve(par_arr.size());
    for (int i = 0; i < par_arr.size(); i++) {
      seq.push_back(i);
    }

    char logfile[100] = "Ran.log";
    //snprintf(xyfile, 100, "PEP_Lx%d_Ly%d_p%g_a%g_NM%d.extxyz", Lx, Ly, phi, alpha, N_MAX);
    LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
    //XYExporter xy(xyfile, 0, n_step, 1000, Lx, Ly);

    //xy.dump(0, par_arr);
    for (int i = 1; i <= n_step; i++) {
      shuffle(seq, myran);
      update_in_rand_seq(par_arr, lattice, myran, seq, alpha);
      log.record(i);
      //xy.dump(i, par_arr);
    }
  }

  {
    Ranq1 myran(2);

    std::vector<Par_2> par_arr;
    SquareLattice lattice(Lx, Ly);
    lattice.create_particles_randomly(par_arr, myran, phi);
    std::vector<size_t> seq;
    seq.reserve(par_arr.size());
    for (int i = 0; i < par_arr.size(); i++) {
      seq.push_back(i);
    }

    char logfile[100] = "Ranq1.log";
    //snprintf(xyfile, 100, "PEP_Lx%d_Ly%d_p%g_a%g_NM%d.extxyz", Lx, Ly, phi, alpha, N_MAX);
    LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
    //XYExporter xy(xyfile, 0, n_step, 1000, Lx, Ly);

    //xy.dump(0, par_arr);
    for (int i = 1; i <= n_step; i++) {
      shuffle(seq, myran);
      update_in_rand_seq(par_arr, lattice, myran, seq, alpha);
      log.record(i);
      //xy.dump(i, par_arr);
    }
  }

  {
    Ranq2 myran(2);

    std::vector<Par_2> par_arr;
    SquareLattice lattice(Lx, Ly);
    lattice.create_particles_randomly(par_arr, myran, phi);
    std::vector<size_t> seq;
    seq.reserve(par_arr.size());
    for (int i = 0; i < par_arr.size(); i++) {
      seq.push_back(i);
    }

    char logfile[100] = "Ranq2.log";
    //snprintf(xyfile, 100, "PEP_Lx%d_Ly%d_p%g_a%g_NM%d.extxyz", Lx, Ly, phi, alpha, N_MAX);
    LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
    //XYExporter xy(xyfile, 0, n_step, 1000, Lx, Ly);

    //xy.dump(0, par_arr);
    for (int i = 1; i <= n_step; i++) {
      shuffle(seq, myran);
      update_in_rand_seq(par_arr, lattice, myran, seq, alpha);
      log.record(i);
      //xy.dump(i, par_arr);
    }
  }

  {
    Ranfib myran(2);

    std::vector<Par_2> par_arr;
    SquareLattice lattice(Lx, Ly);
    lattice.create_particles_randomly(par_arr, myran, phi);
    std::vector<size_t> seq;
    seq.reserve(par_arr.size());
    for (int i = 0; i < par_arr.size(); i++) {
      seq.push_back(i);
    }

    char logfile[100] = "Ranfib.log";
    //snprintf(xyfile, 100, "PEP_Lx%d_Ly%d_p%g_a%g_NM%d.extxyz", Lx, Ly, phi, alpha, N_MAX);
    LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
    //XYExporter xy(xyfile, 0, n_step, 1000, Lx, Ly);

    //xy.dump(0, par_arr);
    for (int i = 1; i <= n_step; i++) {
      shuffle(seq, myran);
      update_in_rand_seq(par_arr, lattice, myran, seq, alpha);
      log.record(i);
      //xy.dump(i, par_arr);
    }
  }


  {
    Ranq1 myran(2);

    std::vector<Par_2> par_arr;
    SquareLattice lattice(Lx, Ly);
    lattice.create_particles_randomly(par_arr, myran, phi);
    std::vector<size_t> seq;
    seq.reserve(par_arr.size());
    for (int i = 0; i < par_arr.size(); i++) {
      seq.push_back(i);
    }

    char logfile[100] = "Ranq1_int64.log";
    //snprintf(xyfile, 100, "PEP_Lx%d_Ly%d_p%g_a%g_NM%d.extxyz", Lx, Ly, phi, alpha, N_MAX);
    LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
    //XYExporter xy(xyfile, 0, n_step, 1000, Lx, Ly);

    //xy.dump(0, par_arr);
    unsigned long long beta = static_cast<unsigned long long>(alpha * UINT64_MAX);
    for (int i = 1; i <= n_step; i++) {
      shuffle(seq, myran);
      update_in_rand_seq(par_arr, lattice, myran, seq, beta);
      log.record(i);
      //xy.dump(i, par_arr);
    }
  }

  {
    Ranq2 myran(2);

    std::vector<Par_2> par_arr;
    SquareLattice lattice(Lx, Ly);
    lattice.create_particles_randomly(par_arr, myran, phi);
    std::vector<size_t> seq;
    seq.reserve(par_arr.size());
    for (int i = 0; i < par_arr.size(); i++) {
      seq.push_back(i);
    }

    char logfile[100] = "Ranq2_int64.log";
    //snprintf(xyfile, 100, "PEP_Lx%d_Ly%d_p%g_a%g_NM%d.extxyz", Lx, Ly, phi, alpha, N_MAX);
    LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
    //XYExporter xy(xyfile, 0, n_step, 1000, Lx, Ly);

    //xy.dump(0, par_arr);
    unsigned long long beta = static_cast<unsigned long long>(alpha * UINT64_MAX);
    for (int i = 1; i <= n_step; i++) {
      shuffle(seq, myran);
      update_in_rand_seq(par_arr, lattice, myran, seq, beta);
      log.record(i);
      //xy.dump(i, par_arr);
    }
  }
}

void test_3D(int L) {
  int Lx = L;
  int Ly = Lx;
  int Lz = Lx;
  double phi = 0.01;
  double alpha = 0.01;
  int n_step = 1000;
  {
    Ranq2 myran(2);

    std::vector<Par_3> par_arr;
    CubicLattice lattice(Lx, Ly, Lz);
    lattice.create_particles_randomly(par_arr, myran, phi);
    std::vector<size_t> seq;
    seq.reserve(par_arr.size());
    for (int i = 0; i < par_arr.size(); i++) {
      seq.push_back(i);
    }

    char logfile[100];
    char xyfile[100];
    snprintf(logfile, 100, "%d_0.log", Lx);

    snprintf(xyfile, 100, "%d_0.extxyz", Lx);
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
}
