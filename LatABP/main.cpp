#include <iostream>
#include <fstream>
#include <iomanip>
#include "ABP.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 200;
  int Ly = 200;
  double phi = 0.5;

  double rate_P = 1;
  double rate_D = 2. / 6;
  double rate_R = 0.032 / 2;


  int n_step = 1000000;
  int dn_out = 10000;
  int seed = 1001;


  run(Lx, Ly, phi, rate_P, rate_D, rate_R,
    n_step, dn_out, seed);
}