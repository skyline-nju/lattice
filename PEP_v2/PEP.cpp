#include <iostream>
#include "PEP.h"
#include "latticeExporter.h"

int main(int agrc, char* argv[]) {
  int Lx = 500;
  int Ly = 500;
  double phi = 0.1;
  double alpha = 0.01;
  int n_step = 50000;
  Ran myran(2);
  
  std::vector<Par_2> par_arr;
  SquareLattice_2 lattice(Lx, Ly);
  lattice.create_particles_randomly(par_arr, myran, phi);
  std::vector<size_t> seq;
  seq.reserve(par_arr.size());
  for (int i = 0; i < par_arr.size(); i++) {
    seq.push_back(i);
  }

  char logfile[100];
  char xyfile[100];
  snprintf(logfile, 100, "PEP_Lx%d_Ly%d_p%g_a%g_NM%d.log", Lx, Ly, phi, alpha, N_MAX);
  snprintf(xyfile, 100, "PEP_Lx%d_Ly%d_p%g_a%g_NM%d.extxyz", Lx, Ly, phi, alpha, N_MAX);
  LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
  XYExporter xy(xyfile, 0, n_step, 1000, Lx, Ly);

  xy.dump(0, par_arr);
  for (int i = 1; i <= n_step; i++) {
    shuffle(seq, myran);
    update_in_rand_seq(par_arr, lattice, myran, seq, alpha);
    log.record(i);
    xy.dump(i, par_arr);
  }

}


