#include <iostream>
#include <fstream>
#include <iomanip>
#include "AIM.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 128;
  int Ly = 64;
  double beta = 2.5;
  double eps = 0.9;
  double rho0 = 5;
  double D = 1;
  double alpha = 0;
  int n_step = 100000;
  int dn_out = 1000;
  unsigned long long seed = 4001;
  std::string ini_condi = "resume";

  run(Lx, Ly, rho0, beta, eps, D, n_step, dn_out, seed, alpha, ini_condi);
}