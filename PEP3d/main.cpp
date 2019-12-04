#include <iostream>
#include "PEP3D.h"
#include "latticeExporter_3.h"

int main(int agrc, char* argv[]) {
#ifdef _MSC_VER
  int Lx = 100;
  int Ly = 100;
  int Lz = 100;
  double phi = 0.01;
  double alpha = 0.01;
  int n_step = 10000;
#else
  int L = atoi(argv[1]);
  int Lx = L;
  int Ly = L;
  int Lz = L;
  double phi = atof(argv[2]);
  double alpha = atof(argv[3]);
  int n_step = atoi(argv[4]);
#endif

  std::vector<PEP_3::LatticeCoord> par_arr;
  PEP_3::SquareLattice<unsigned char> lattice(Lx, Ly, Lz);
  Ranq1 myran(2);
  lattice.create_rand(phi, myran, par_arr);
  std::vector<int> seq;
  seq.reserve(par_arr.size());
  for (int i = 0; i < par_arr.size(); i++) {
    seq.push_back(i);
  }

  char logfile[100];
  char xyzfile[100];
  snprintf(logfile, 100, "PEP3_Lx%d_Ly%d_p%g_a%g_NM%d.log", Lx, Ly, phi, alpha, N_MAX);
  snprintf(xyzfile, 100, "PEP3_Lx%d_Ly%d_p%g_a%g_NM%d.extxyz", Lx, Ly, phi, alpha, N_MAX);
  LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
  XYZExporter xyz(xyzfile, 0, n_step, 10000, Lx, Ly, Lz);

  xyz.dump(0, par_arr);
  for (int i = 1; i <= n_step; i++) {
    lattice.update_rand_seq(par_arr, myran, seq, alpha);
    log.record(i);
    xyz.dump(i, par_arr);
    if (i % 100 == 0) {
      std::cout << "t = " << i << std::endl;
    }
  }
}