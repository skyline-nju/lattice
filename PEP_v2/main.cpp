#include "PEP.h"
#include "latticeExporter.h"

void run_3D(int argc, char** argv) {
#ifdef _MSC_VER
  int Lx = 200;
  int Ly = 200;
  int Lz = 10;
  double phi = 0.08;
  double alpha = 0.001;
  int n_step = 1000000;
#else
  int Lx = atoi(argv[1]);
  int Ly = atoi(argv[2]);
  int Lz = atoi(argv[3]);
  double phi = atof(argv[4]);
  double alpha = atof(argv[5]);
  int n_step = atoi(argv[6]);
#endif
  {
    Ranq2 myran(4);

    std::vector<Par_3> par_arr;
    CubicLattice lattice(Lx, Ly, Lz);
    lattice.create_particles_randomly(par_arr, myran, phi);
    std::vector<size_t> seq;
    seq.reserve(par_arr.size());
    for (int i = 0; i < par_arr.size(); i++) {
      seq.push_back(i);
    }

    char logfile[100];
    char xyzfile[100];
    //char profile[100];
    snprintf(logfile, 100, "NM%d_Lx%d_Ly%d_Lz%d_a%g_p%g.log", N_MAX, Lx, Ly, Lz, alpha, phi);
    snprintf(xyzfile, 100, "NM%d_Lx%d_Ly%d_Lz%d_a%g_p%g.extxyz", N_MAX, Lx, Ly, Lz, alpha, phi);
    //snprintf(profile, 100, "NM%d_Lx%d_Ly%d_a%g_p%g.bin", N_MAX, Lx, Ly, alpha, phi);
    LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
    XYZExporter xyz(xyzfile, 0, n_step, 10000, Lx, Ly, Lz);
    //WettingProfileExporter_2 pro(profile, 0, n_step, 1000, Lx, Ly);
    xyz.dump(0, par_arr);
    unsigned long long beta = static_cast<unsigned long long>(alpha * UINT64_MAX);
    for (int i = 1; i <= n_step; i++) {
      shuffle(seq, myran);
      update_in_rand_seq(par_arr, lattice, myran, seq, beta);
      log.record(i);
      xyz.dump(i, par_arr);
    }
  }
}

void run_3D_biased(int argc, char** argv) {
#ifdef _MSC_VER
  int Lx = 200;
  int Ly = 200;
  int Lz = 10;
  double phi = 0.05;
  double alpha = 0.01;
  int n_step = 1000000;
  double biased_ratio = 0;
#else
  int Lx = atoi(argv[1]);
  int Ly = atoi(argv[2]);
  int Lz = atoi(argv[3]);
  double phi = atof(argv[4]);
  double alpha = atof(argv[5]);
  int n_step = atoi(argv[6]);
  double biased_ratio = atof(argv[7]);
#endif

  int snap_sep = 10000;
  double tumble_rates[5];
  double tot_rates = 4. + 2 * biased_ratio;
  tumble_rates[0] = 1. / tot_rates;
  tumble_rates[1] = tumble_rates[0] * 2;
  tumble_rates[2] = tumble_rates[0] * 3;
  tumble_rates[3] = tumble_rates[0] * 4;
  tumble_rates[4] = tumble_rates[3] + biased_ratio / tot_rates;
  Ranq2 myran(5);

  std::vector<Par_3> par_arr;
  CubicLattice lattice(Lx, Ly, Lz);
  lattice.create_particles_randomly(par_arr, myran, phi);
  std::vector<size_t> seq;
  seq.reserve(par_arr.size());
  for (int i = 0; i < par_arr.size(); i++) {
    seq.push_back(i);
  }

  char logfile[100];
  char xyzfile[100];
  //char profile[100];
  snprintf(logfile, 100, "NM%d_Lx%d_Ly%d_Lz%d_p%g_a%g_b%g.log", N_MAX, Lx, Ly, Lz, phi, alpha, biased_ratio);
  snprintf(xyzfile, 100, "NM%d_Lx%d_Ly%d_Lz%d_p%g_a%g_b%g.extxyz", N_MAX, Lx, Ly, Lz, phi, alpha, biased_ratio);
  //snprintf(profile, 100, "NM%d_Lx%d_Ly%d_a%g_p%g.bin", N_MAX, Lx, Ly, alpha, phi);
  LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
  XYZExporter xyz(xyzfile, 0, n_step, snap_sep, Lx, Ly, Lz);
  //WettingProfileExporter_2 pro(profile, 0, n_step, 1000, Lx, Ly);
  xyz.dump(0, par_arr);
  unsigned long long beta = static_cast<unsigned long long>(alpha * UINT64_MAX);
  for (int i = 1; i <= n_step; i++) {
    shuffle(seq, myran);
    update_in_rand_seq(par_arr, lattice, myran, seq, beta, tumble_rates);
    log.record(i);
    xyz.dump(i, par_arr);
  }
}

void run_2D(int argc, char** argv) {
#ifdef _MSC_VER
  int Lx = 200;
  int Ly = 200;
  double phi = 0.5;
  double alpha = 0.01;
  int n_step = 1000000;
#else
  int Lx = atoi(argv[1]);
  int Ly = atoi(argv[2]);
  double phi = atof(argv[3]);
  double alpha = atof(argv[4]);
  int n_step = atoi(argv[5]);
#endif

  Ranq2 myran(5);
  std::vector<Par_2> par_arr;
  SquareLattice lattice(Lx, Ly);
  lattice.create_particles_randomly(par_arr, myran, phi);
  std::vector<size_t> seq;
  seq.reserve(par_arr.size());
  for (int i = 0; i < par_arr.size(); i++) {
    seq.push_back(i);
  }

  char logfile[100];
  char xyfile[100];
  //char profile[100];
  snprintf(logfile, 100, "NM%d_Lx%d_Ly%d_p%g_a%g.log", N_MAX, Lx, Ly, phi, alpha);
  snprintf(xyfile, 100, "NM%d_Lx%d_Ly%d_p%g_a%g.extxyz", N_MAX, Lx, Ly, phi, alpha);
  //snprintf(profile, 100, "NM%d_Lx%d_Ly%d_a%g_p%g.bin", N_MAX, Lx, Ly, alpha, phi);
  LogExporter log(logfile, 0, n_step, 10000, par_arr.size());
  XYExporter xy(xyfile, 0, n_step, 10000, Lx, Ly);
  //WettingProfileExporter_2 pro(profile, 0, n_step, 1000, Lx, Ly);
  xy.dump(0, par_arr);
  unsigned long long beta = static_cast<unsigned long long>(alpha * UINT64_MAX);
  for (int i = 1; i <= n_step; i++) {
    shuffle(seq, myran);
    update_in_rand_seq(par_arr, lattice, myran, seq, beta);
    log.record(i);
    xy.dump(i, par_arr);
    //if (i % 1000 == 0) {
    //  pro.dump(i, lattice);
    //}
  }
}

int main(int argc, char* argv[]) {
  run_3D_biased(argc, argv);
  //run_2D(argc, argv);
}