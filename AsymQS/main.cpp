#include <iostream>
#include <fstream>
#include <iomanip>
#include "AsymQS.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 128;
  int Ly = 16;
  double rho_thresh = 80;
  double phi = rho_thresh;

  double Dt = 0.3;
  double Dr = 0.01;
  double v0 = 1;

  double eta = 2;
  double alpha = 1;

  int n_step = 10000000;
  int dn_out = 1000;
  int seed = 1001;


  run(Lx, Ly, phi, rho_thresh, Dt, Dr, v0,
    eta, alpha, n_step, dn_out, seed);
}