#include <iostream>
#include <fstream>
#include <iomanip>
#include "NRQS.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 64;
  int Ly = 8;
  double rho_thresh = 20;
  double phiA = rho_thresh;
  double phiB = rho_thresh;

  double Dt = 0.2;
  double Dr = 0.05;
  double v0 = 1;

  double etaAA = 0;
  double etaBB = 0;
  double etaAB = 1;
  double etaBA = -etaAB;

  int n_step = 10000000;
  int dn_out = 5000;
  unsigned long long seed = 1001;


  run(Lx, Ly, phiA, phiB, rho_thresh, Dt, Dr, v0,
    etaAA, etaAB, etaBA, etaBB, n_step, dn_out, seed);
}