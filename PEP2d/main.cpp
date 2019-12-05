#include <iostream>
#include "PEP2D.h"
#include "latticeExporter_2.h"

int main(int agrc, char* argv[]) {
  int Lx = 1000;
  int Ly = 1000;
  double phi = 0.01;
  double alpha = 0.01;
  int n_step = 50000;

  std::vector<PEP_2::LatticeCoord> par_arr;
  PEP_2::SquareLattice<unsigned char> lattice(Lx, Ly);
  Ranq1 myran(2);
  lattice.create_rand(phi, myran, par_arr);
  std::vector<int> seq;
  seq.reserve(par_arr.size());
  for (int i = 0; i < par_arr.size(); i++) {
    seq.push_back(i);
  }

  char logfile[100];
  char xyfile[100];
  snprintf(logfile, 100, "PEP2_Lx%d_Ly%d_p%g_a%g_NM%d.log", Lx, Ly, phi, alpha, N_MAX);
  snprintf(xyfile, 100, "PEP2_Lx%d_Ly%d_p%g_a%g_NM%d.extxyz", Lx, Ly, phi, alpha, N_MAX);
  LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
  XYExporter xy(xyfile, 0, n_step, 1000, Lx, Ly);

  xy.dump(0, par_arr);
  for (int i = 1; i <= n_step; i++) {
    lattice.update_rand_seq(par_arr, myran, seq, alpha);
    log.record(i);
    xy.dump(i, par_arr);
  }
}