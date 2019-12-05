#include <iostream>
#include "PEP.h"
#include "latticeExporter.h"

template <>
const short SquarePar_2<short>::ori[4][2] = { {1, 0}, {0, 1}, {-1, 0}, {0, -1} };

int main(int agrc, char* argv[]) {
  int Lx = 1000;
  int Ly = 4000;
  double phi = 0.01;
  double alpha = 0.01;
  int n_step = 10000;

  std::vector<SquarePar_2<short>> par_arr;
  //SimpleSquareLattice_2 lattice(Lx, Ly);
  PackedSquareLattice_2 lattice(Lx, Ly, 4, 3);
  Ranq1 myran(2);
  create_rand(phi, par_arr, lattice, myran);
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
    run_in_rand_seq(par_arr, lattice, myran, seq, alpha);
    log.record(i);
    xy.dump(i, par_arr);
  }
}