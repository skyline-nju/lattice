#include "ABP.h"
#include "lattice_io.h"
#include <cmath>
#include <fstream>
#include <iomanip>

hexa_lattice_2::hexa_lattice_2(const Vec_2<int>& l, double phi,
                               double rate_P, double rate_D, double rate_R)
  :l_(l), n_sites_(l.x * l.y), phi_(phi), rate_P_(rate_P), rate_D_(rate_D), rate_R_(rate_R) {
  
  double p0[8]{};
  p0[0] = rate_D_;
  p0[1] = p0[0] + rate_D_;
  p0[2] = p0[1] + rate_D_;
  p0[3] = p0[2] + rate_D_;
  p0[4] = p0[3] + rate_D_;
  p0[5] = p0[4] + rate_D_;
  p0[6] = p0[5] + rate_R_;
  p0[7] = p0[6] + rate_R_;

  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 8; i++) {
      prob_arr_[j][i] = p0[i];
      if (i >= j) {
        prob_arr_[j][i] += rate_P_;
      }
    }
  }
  
  n_par_ = int(round(phi_ * n_sites_));
  delta_t_ = 1 / (rate_P_ + 6 * rate_D_ + 2 * rate_R_);
  site_state_.reserve(n_sites_);
  for (int i = 0; i < n_sites_; i++) {
    site_state_.push_back(0);
  }
  p_arr_.reserve(n_par_);
  std::cout << "rho_0 = " << phi_ << "\n";
  std::cout << "L = " << l_ << "\n";
  std::cout << "total particles = " << n_par_ << std::endl;
}

int hexa_lattice_2::get_site_idx(const Par_2& p, int x_offset, int y_offset) const {
  int y_new = p.pos.y + y_offset;
  if (y_offset != 0) {
    tangle_1(y_new, l_.y);
  }
  int idx = y_new * l_.x;

  if (x_offset == 0) {
    idx += p.pos.x;
  } else {
    int x_new = p.pos.x + x_offset;
    tangle_1(x_new, l_.x);
    idx += x_new;
  }
  return idx;
}


void hexa_lattice_2::add_particle(const Par_2& p) {
  int j = get_site_idx(p);
  site_state_[j] = 1;
}

void hexa_lattice_2::del_particle(const Par_2& p){
  int j = get_site_idx(p);
  site_state_[j] = 0;
}

void hexa_lattice_2::hop(Par_2& p, const Vec_2<int>& dR) {
  int x_new = p.pos.x + dR.x;
  tangle_1(x_new, l_.x);
  int y_new = p.pos.y + dR.y;
  tangle_1(y_new, l_.y);

  int idx_new = x_new * l_.x + y_new;
  if (site_state_[idx_new] == 0) {
    int idx_old = p.pos.x * l_.x + p.pos.y;
    site_state_[idx_new] = 1;
    site_state_[idx_old] = 0;
    p.pos.x = x_new;
    p.pos.y = y_new;
  }
}

void hexa_lattice_2::hop_rot(Par_2& p, double rand_val) {
  int s = p.spin;
  if (rand_val < prob_arr_[s][0]) {
    hop(p, Vec_2<int>(1, 0));
  } else if (rand_val < prob_arr_[s][1]) {
    hop(p, Vec_2<int>(1, 1));
  } else if (rand_val < prob_arr_[s][2]) {
    hop(p, Vec_2<int>(0, 1));
  } else if (rand_val < prob_arr_[s][3]) {
    hop(p, Vec_2<int>(-1, 0));
  } else if (rand_val < prob_arr_[s][4]) {
    hop(p, Vec_2<int>(-1, -1));
  } else if (rand_val < prob_arr_[s][5]) {
    hop(p, Vec_2<int>(0, -1));
  } else if (rand_val < prob_arr_[s][6]) {
    p.spin += 1;
    if (p.spin >= 6) {
      p.spin -= 6;
    }
  } else {
    p.spin -= 1;
    if (p.spin < 0) {
      p.spin += 6;
    }
  }
}


void run(int Lx, int Ly, double phi,
         double rate_P, double rate_D, double rate_R,
         int n_step, int dn_out, int seed){
  Vec_2<int> l(Lx, Ly);
  Ranq2 myran(seed);
  hexa_lattice_2 domain(l, phi, rate_P, rate_D, rate_R);
  domain.ini_rand(myran);
  //domain.ini_rand(myran, 1);
  char folder[255] = "data";
  char basename[255];
  snprintf(basename, 255, "L%d_%d_r%.3f_v%g_Dt%.3f_Dr%.3f_s%d",
    Lx, Ly, phi, rate_P, rate_D, rate_R, seed);

  char log_file[255];
  char snap_file[255];

  int t_start = 0;
  snprintf(log_file, 255, "%s/log_%s.dat", folder, basename);
  LogExporter log(log_file, 0, n_step, 10000, domain.get_n_par());

  snprintf(snap_file, 255, "%s/%s.extxyz", folder, basename);
  XYExporter xy_exporter(snap_file, 0, n_step, dn_out, Lx, Ly);

  xy_exporter.dump(0, domain.get_p_arr());


  for (int i = 1; i <= n_step; i++) {
    domain.one_step(myran);
    log.record(i);
    xy_exporter.dump(i, domain.get_p_arr());
  }
}



