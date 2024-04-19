#include "NRQS_q4.h"
#include <cmath>
#include <fstream>
#include <iomanip>

lattice_4::lattice_4(const Vec_2<int>& l, double Dt, double Dr,
  double phiA, double phiB, double rho_thresh, double v0,
  double etaAA, double etaAB, double etaBA, double etaBB)
  :l_(l), n_sites_(l.x* l.y), Dt_(Dt), Dr_(Dr), v0_(v0),
   phiA_(phiA), phiB_(phiB), rho_thresh_(rho_thresh) {
  eta_[0][0] = etaAA / rho_thresh_;
  eta_[0][1] = etaAB / rho_thresh_;
  eta_[1][0] = etaBA / rho_thresh_;
  eta_[1][1] = etaBB / rho_thresh_;
  self_couplings_on_ = !(etaAA == 0 && etaBB == 0);

  ori_[0][0] = Vec_2<int>(1, 0);
  ori_[0][1] = Vec_2<int>(0, 1);
  ori_[0][2] = Vec_2<int>(-1, 0);
  ori_[0][3] = Vec_2<int>(0, -1);

  ori_[1][3] = Vec_2<int>(1, 0);
  ori_[1][0] = Vec_2<int>(0, 1);
  ori_[1][1] = Vec_2<int>(-1, 0);
  ori_[1][2] = Vec_2<int>(0, -1);

  ori_[2][2] = Vec_2<int>(1, 0);
  ori_[2][3] = Vec_2<int>(0, 1);
  ori_[2][0] = Vec_2<int>(-1, 0);
  ori_[2][1] = Vec_2<int>(0, -1);

  ori_[3][1] = Vec_2<int>(1, 0);
  ori_[3][2] = Vec_2<int>(0, 1);
  ori_[3][3] = Vec_2<int>(-1, 0);
  ori_[3][0] = Vec_2<int>(0, -1);

  double max_vel;
  if (self_couplings_on_) {
    max_vel = 4 * v0;
  } else {
    max_vel = 2 * v0;
  }

  prob_arr_[0] = max_vel + Dt_;
  prob_arr_[1] = prob_arr_[0] + Dt_;
  prob_arr_[2] = prob_arr_[1] + Dt_;
  prob_arr_[3] = prob_arr_[2] + Dt_;
  prob_arr_[4] = prob_arr_[3] + Dr_;
  prob_arr_[5] = prob_arr_[4] + Dr_;

  for (int i = 0; i < 6; i++) {
    std::cout << prob_arr_[i] << std::endl;
  }

  delta_t_ = 1. / prob_arr_[5];

  sigma_ = new unsigned short[2 * n_sites_ * 4] {};

  n_par_A_ = int(l_.x * l_.y * phiA_);
  n_par_ = n_par_A_ + int(l_.x * l_.y * phiB_);
  p_arr_.reserve(n_par_);
  std::cout << "phiA = " << phiA_ << "\n";
  std::cout << "phiB = " << phiB_ << "\n";
  std::cout << "L = " << l_ << "\n";
  std::cout << "total particles = " << n_par_ << std::endl;
}

lattice_4::~lattice_4() {
  delete[] sigma_;
}

void lattice_4::get_rho(const Vec_2<short>& pos, double& rho_A, double& rho_B) const {
  int j = (int(pos.x) + pos.y * l_.x) * 8;
  rho_A = sigma_[j] + sigma_[j + 1] + sigma_[j + 2] + sigma_[j + 3];
  rho_B = sigma_[j + 4] + sigma_[j + 5] + sigma_[j + 6] + sigma_[j + 7];
}

double lattice_4::get_v(const Par_4& p) const {
  double rho_A, rho_B;
  get_rho(p.pos, rho_A, rho_B);
  double drho[2] = { rho_A - rho_thresh_, rho_B - rho_thresh_ };
  int s1 = p.species;
  if (self_couplings_on_) {
    return v0_ * (1 + tanh(eta_[s1][0] * drho[0])) * (1 + tanh(eta_[s1][1] * drho[1]));
  } else {
    int s2 = 1 - p.species;
    return v0_ * (1 + tanh(eta_[s1][s2] * drho[s2]));
  }
}

void lattice_4::hop_rot(Par_4& p, double rand_val) const {
  if (rand_val < prob_arr_[0]) {
    double v = get_v(p);
    if (rand_val < Dt_ + v) {
      hop(p, ori_[p.spin][0]);
    }
  } else if (rand_val < prob_arr_[1]) {
    hop(p, ori_[p.spin][1]);
  } else if (rand_val < prob_arr_[2]) {
    hop(p, ori_[p.spin][2]);
  } else if (rand_val < prob_arr_[3]) {
    hop(p, ori_[p.spin][3]);
  } else if (rand_val < prob_arr_[4]) {
    rot(p, 1);
  } else {
    rot(p, -1);
  }
}

void lattice_4::hop(Par_4& p, const Vec_2<int>& ori) const {
  del_particle(p);
  if (ori.x == 0) {
    p.pos.y += ori.y;
    tangle_1(p.pos.y, l_.y);
  } else {
    p.pos.x += ori.x;
    tangle_1(p.pos.x, l_.x);
  }
  add_particle(p);
}

void lattice_4::rot(Par_4& p, int ds) const {
  int j = (int(p.pos.x) + p.pos.y * l_.x) * 8 + p.species * 4;
  int i_old = j + p.spin;
  p.spin += ds;
  if (p.spin < 0) {
    p.spin += 4;
  } else if (p.spin >= 4) {
    p.spin -= 4;
  }
  int i_new = j + p.spin;
  sigma_[i_old] -= 1;
  sigma_[i_new] += 1;
}


double lattice_4::cal_rho0() const {
  double rho_sum = 0;
  size_t end = n_sites_ * 4 * 2;
  for (int i = 0; i < end; i++) {
    rho_sum += sigma_[i];
  }
  return rho_sum / n_sites_;
}

void lattice_4::output_snap(std::ofstream& fout) {
  fout.write((const char*)&sigma_[0], sizeof(unsigned short) * 8 * n_sites_);
}

void run_q4(int Lx, int Ly,
  double phiA, double phiB, double rho_thresh,
  double Dt, double Dr, double v0,
  double etaAA, double etaAB, double etaBA, double etaBB,
  int n_step, int dn_out, int seed) {
  Vec_2<int> l(Lx, Ly);
  Ranq2 myran(seed);
  lattice_4 domain(l, Dt, Dr, phiA, phiB, rho_thresh, v0, etaAA, etaAB, etaBA, etaBB);
  domain.ini_rand(myran);
  //domain.ini_rand(myran, 1);
  char folder[255] = "data";
  char basename[255];
  snprintf(basename, 255, "L%d_%d_Dr%g_Dt%g_e%g_%g_J%g_%g_v%g_r%g_s%d",
    Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho_thresh, seed);

  //char order_para_file[255];
  char snap_file[255];

  int t_start = 0;
  //snprintf(order_para_file, 256, "%s/%s_t%d.dat", folder, basename, t_start);
  //std::ofstream fout(order_para_file);

  snprintf(snap_file, 255, "%s/%s_dt%d_t%d.bin", folder, basename, dn_out, t_start);
  std::ofstream fsnap(snap_file, std::ios::binary);

  std::cout << sizeof(Par_4) << std::endl;

  for (int i = 1; i <= n_step; i++) {
    domain.one_step(myran);
    if (i % 5000 == 0) {
      std::cout << i << "\t" << domain.cal_rho0() << "\n";
      if (i % dn_out == 0) {
        domain.output_snap(fsnap);
      }
    }
  }
}



