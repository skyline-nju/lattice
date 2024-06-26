#include <iostream>
#include <fstream>
#include <iomanip>
#include "QS_q6.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 128;
  int Ly = 128;
  double rho_thresh = 10;
  double phi = rho_thresh;

  double Dt = 0.03;
  double Dr = 0.1;
  double v0 = 1;

  double eta = -2;

  int n_step = 10000000;
  int dn_out = 10000;
  unsigned long long seed = 1000;


  run_q6(Lx, Ly, phi, rho_thresh, Dt, Dr, v0,
    eta, n_step, dn_out, seed);
}