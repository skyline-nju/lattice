#include "QS_q6.h"
#include <cmath>
#include <fstream>
#include <iomanip>

lattice_6::lattice_6(const Vec_2<int>& l, double Dt, double Dr,
                     double phi, double rho_thresh, double v0, double eta)
  :l_(l), n_sites_(l.x* l.y), Dt_(Dt), Dr_(Dr), v0_(v0),
  phi_(phi), rho_thresh_(rho_thresh), eta_(eta /rho_thresh) {
  ori_[0][0] = Vec_2<int>(1, 0);
  ori_[0][1] = Vec_2<int>(1, 1);
  ori_[0][2] = Vec_2<int>(0, 1);
  ori_[0][3] = Vec_2<int>(-1, 0);
  ori_[0][4] = Vec_2<int>(-1, -1);
  ori_[0][5] = Vec_2<int>(0, -1);

  for (int i = 1; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      int k = j + i;
      if (k >= 6) {
        k -= 6;
      }
      ori_[i][j] = ori_[0][k];
    }
  }

  prob_arr_[0] = 2 * v0_ + Dt_;
  //prob_arr_[0] = 4 + Dt_;
  prob_arr_[1] = prob_arr_[0] + Dt_;
  prob_arr_[2] = prob_arr_[1] + Dt_;
  prob_arr_[3] = prob_arr_[2] + Dt_;
  prob_arr_[4] = prob_arr_[3] + Dt_;
  prob_arr_[5] = prob_arr_[4] + Dt_;
  prob_arr_[6] = prob_arr_[5] + Dr_;
  prob_arr_[7] = prob_arr_[6] + Dr_;

  for (int i = 0; i < 8; i++) {
    std::cout << prob_arr_[i] << std::endl;
  }

  delta_t_ = 1. / prob_arr_[7];

  sigma_ = new unsigned char[size_t(n_sites_) * 6] {};
  rho_ = new unsigned short[n_sites_] {};

  n_par_ = int(l_.x * l_.y * phi_);
  p_arr_.reserve(n_par_);
  std::cout << "phi = " << phi_ << "\n";
  std::cout << "L = " << l_ << "\n";
  std::cout << "total particles = " << n_par_ << std::endl;
}

lattice_6::~lattice_6() {
  delete[] sigma_;
  delete[] rho_;
}

void lattice_6::add_particle(const Par_6& p) const {
  int site_idx = get_site_idx(p);
  int sigma_idx = site_idx * 6 + p.spin;
  rho_[site_idx] += 1;
  sigma_[sigma_idx] += 1;
}

void lattice_6::del_particle(const Par_6& p) const {
  int site_idx = get_site_idx(p);
  int sigma_idx = site_idx * 6 + p.spin;
  rho_[site_idx] -= 1;
  sigma_[sigma_idx] -= 1;
}

double lattice_6::get_v_back(const Par_6& p) const {
  Vec_2<int> my_ori = ori_[p.spin][0];
  int x_back = p.pos.x - my_ori.x;
  int y_back = p.pos.y - my_ori.y;
  tangle_1(x_back, l_.x);
  tangle_1(y_back, l_.y);
  int idx_back0 = x_back + y_back * l_.x;

  my_ori = ori_[p.spin][1];
  x_back = p.pos.x - my_ori.x;
  y_back = p.pos.y - my_ori.y;
  tangle_1(x_back, l_.x);
  tangle_1(y_back, l_.y);
  int idx_back1 = x_back + y_back * l_.x;

  my_ori = ori_[p.spin][5];
  x_back = p.pos.x - my_ori.x;
  y_back = p.pos.y - my_ori.y;
  tangle_1(x_back, l_.x);
  tangle_1(y_back, l_.y);
  int idx_back2 = x_back + y_back * l_.x;

  double rho_back = (rho_[idx_back0] + rho_[idx_back1] + rho_[idx_back2]) / 3;

  return v0_ * (1 + tanh(eta_ * (rho_back - rho_thresh_)));
}

double lattice_6::get_v_front(const Par_6& p) const {
  Vec_2<int> my_ori = ori_[p.spin][0];
  int x_front = p.pos.x + my_ori.x;
  int y_front = p.pos.y + my_ori.y;
  tangle_1(x_front, l_.x);
  tangle_1(y_front, l_.y);
  int idx_front = x_front + y_front * l_.x;
  return v0_ * (1 + tanh(eta_ * (rho_[idx_front] - rho_thresh_)));
}


void lattice_6::hop_rot(Par_6& p, double rand_val) const {
  if (rand_val < prob_arr_[0]) {
    double v = get_v(p);
    // double v = get_v_back(p);
    //double v = get_v_front(p);
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
    hop(p, ori_[p.spin][4]);
  } else if (rand_val < prob_arr_[5]) {
    hop(p, ori_[p.spin][5]);
  } else if (rand_val < prob_arr_[6]) {
    rot(p, 1);
  } else {
    rot(p, -1);
  }
}

void lattice_6::hop(Par_6& p, const Vec_2<int>& ori) const {
  del_particle(p);
  p.pos.x += ori.x;
  p.pos.y += ori.y;
  tangle_1(p.pos.x, l_.x);
  tangle_1(p.pos.y, l_.y);
  add_particle(p);
}

void lattice_6::rot(Par_6& p, int ds) const {
  int j = get_site_idx(p) * 6;
  int i_old = j + p.spin;
  p.spin += ds;
  if (p.spin < 0) {
    p.spin += 6;
  } else if (p.spin >= 6) {
    p.spin -= 6;
  }
  int i_new = j + p.spin;
  sigma_[i_old] -= 1;
  sigma_[i_new] += 1;
}


double lattice_6::cal_rho0() const {
  double rho_sum = 0;
  size_t end = size_t(n_sites_) * 6;
  for (int i = 0; i < end; i++) {
    rho_sum += sigma_[i];
  }
  return rho_sum / n_sites_;
}

void lattice_6::output_snap(std::ofstream& fout) {
  fout.write((const char*)&sigma_[0], sizeof(unsigned char) * 6 * n_sites_);
}

void run_q6(int Lx, int Ly, double phi, double rho_thresh,
            double Dt, double Dr, double v0, double eta,
            int n_step, int dn_out, int seed) {
  Vec_2<int> l(Lx, Ly);
  Ranq2 myran(seed);
  lattice_6 domain(l, Dt, Dr, phi, rho_thresh, v0, eta);
  domain.ini_rand(myran);
  //domain.ini_rand(myran, 1);
  char folder[255] = "data";
  char basename[255];
  snprintf(basename, 255, "L%d_%d_Dr%.3f_Dt%.3f_e%g_r%g_p%g_v%g_s%d",
    Lx, Ly, Dr, Dt, eta, rho_thresh, phi, v0, seed);

  //char order_para_file[255];
  char snap_file[255];

  int t_start = 0;
  //snprintf(order_para_file, 256, "%s/%s_t%d.dat", folder, basename, t_start);
  //std::ofstream fout(order_para_file);

  snprintf(snap_file, 255, "%s/%s_dt%d_t%d.bin", folder, basename, dn_out, t_start);
  std::ofstream fsnap(snap_file, std::ios::binary);

  std::cout << sizeof(Par_6) << std::endl;

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



