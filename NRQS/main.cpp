#include <iostream>
#include <fstream>
#include <iomanip>
#include "NRQS_q6.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx =256;
  int Ly =256;
  double rho_thresh = 10;
  double phiA = rho_thresh;
  double phiB = rho_thresh;

  double Dt = 0.1;
  double Dr = 0.1;
  double v0 = 1;

  double etaAA = -1;
  double etaBB = -2;
  double etaAB = 1;
  double etaBA = -etaAB;

  int n_step = 1000000;
  int dn_out = 5000;
  unsigned long long seed = 1000;


  run_q6(Lx, Ly, phiA, phiB, rho_thresh, Dt, Dr, v0,
    etaAA, etaAB, etaBA, etaBB, n_step, dn_out, seed);
}