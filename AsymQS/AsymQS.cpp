#include "AsymQS.h"
#include <cmath>
#include <fstream>
#include <iomanip>

lattice_2::lattice_2(const Vec_2<int>& l, double Dt, double Dr,
                     double phi, double rho_thresh, double v0,
                     double eta, double alpha)
  :l_(l), n_sites_(l.x * l.y), Dt_(Dt), Dr_(Dr), v0_(v0),
   phi_(phi), rho_thresh_(rho_thresh), eta_(eta/rho_thresh), alpha_(alpha) {
  prob_arr_[0] = Dr_;
  if (l.y > 1) {
    prob_arr_[1] = prob_arr_[0] + Dt_;
    prob_arr_[2] = prob_arr_[1] + Dt_;
    prob_arr_[3] = prob_arr_[2] + Dt_;
  } else {
    prob_arr_[1] = prob_arr_[0];
    prob_arr_[2] = prob_arr_[1];
    prob_arr_[3] = prob_arr_[2] + Dt_;
  }
  delta_t_ = 1 / (prob_arr_[3] + Dt_ + 2 * v0);
  
  rho_ = new unsigned short[n_sites_] {};
  m_ = new short[n_sites_] {};
  n_par_ = int(l_.x * l_.y * phi_);
  p_arr_.reserve(n_par_);
  std::cout << "rho_0 = " << phi_ << "\n";
  std::cout << "L = " << l_ << "\n";
  std::cout << "total particles = " << n_par_ << std::endl;
}

lattice_2::~lattice_2() {
  delete[] rho_;
  delete[] m_;
}


void lattice_2::add_particle(const Par_2& p) const {
  int j = get_field_idx(p);
  rho_[j] += 1;
  m_[j] += p.spin;
}

void lattice_2::del_particle(const Par_2& p) const {
  int j = get_field_idx(p);
  rho_[j] -= 1;
  m_[j] -= p.spin;
}

int lattice_2::get_field_idx(const Par_2& p, int offset) const {
  int idx = p.pos.y * l_.x;
  if (offset == 0) {
    idx += p.pos.x;
  } else {
    int x_new = p.pos.x + offset;
    tangle_1(x_new, l_.x);
    idx += x_new;
  }
  return idx;
}

double lattice_2::get_v(const Par_2& p) const {
  int j0 = get_field_idx(p);
  int j1 = get_field_idx(p, p.spin);
  double rho_m = (1 - alpha_) * rho_[j0] + alpha_ * rho_[j1];
  return v0_ * (1 + tanh(eta_ * (rho_m - rho_thresh_)));
}

void lattice_2::flip_hop(Par_2& p, double rand_val) const {
  if (rand_val < prob_arr_[0]) {
    int j0 = get_field_idx(p);
    p.spin = -p.spin;
    m_[j0] += 2 * p.spin;
  } else {
    double v = get_v(p);
    if (rand_val < prob_arr_[3] + Dt_ + v) {
      hop(p, rand_val);
    }
  }
}

void lattice_2::hop(Par_2 & p, double rand_val) const {
  del_particle(p);
  if (rand_val < prob_arr_[1]) {
    p.pos.y += 1;
    if (p.pos.y >= l_.y)
      p.pos.y = 0;
  } else if (rand_val < prob_arr_[2]) {
    p.pos.y -= 1;
    if (p.pos.y < 0)
      p.pos.y += l_.y;
  } else if (rand_val < prob_arr_[3]) {
    p.pos.x -= p.spin;
    tangle_1(p.pos.x, l_.x);
  } else {
    p.pos.x += p.spin;
    tangle_1(p.pos.x, l_.x);
  }
  add_particle(p);
}


double lattice_2::cal_m_mean() const {
  double m_sum = 0;
  for (int i = 0; i < n_sites_; i++) {
    m_sum += m_[i];
  }
  return m_sum / n_par_;
}

double lattice_2::cal_rho0() const {
  double rho_sum = 0;
  for (int i = 0; i < n_sites_; i++) {
    rho_sum += rho_[i];
  }
  return rho_sum / n_sites_;
}

void lattice_2::output_snap(std::ofstream& fout) {
  fout.write((const char*)&rho_[0], sizeof(unsigned short) * n_sites_);
  fout.write((const char*)&m_[0], sizeof(short) * n_sites_);
}

void run(int Lx, int Ly,
         double phi, double rho_thresh,
         double Dt, double Dr, double v0,
         double eta, double alpha,
         int n_step, int dn_out, int seed) {
  Vec_2<int> l(Lx, Ly);
  Ranq2 myran(seed);
  lattice_2 domain(l, Dt, Dr, phi, rho_thresh, v0, eta, alpha);
  domain.ini_rand(myran, 1);
  //domain.ini_rand(myran, 1);
  char folder[255] = "data";
  char basename[255];
  snprintf(basename, 255, "L%d_%d_Dr%g_Dt%g_e%g_a%g_v%g_r%g_s%d",
    Lx, Ly, Dr, Dt, eta, alpha, v0, rho_thresh, seed);

  //char order_para_file[255];
  char snap_file[255];

  int t_start = 0;
  //snprintf(order_para_file, 256, "%s/%s_t%d.dat", folder, basename, t_start);
  //std::ofstream fout(order_para_file);

  snprintf(snap_file, 255, "%s/%s_dt%d_t%d.bin", folder, basename, dn_out, t_start);
  std::ofstream fsnap(snap_file, std::ios::binary);


  for (int i = 1; i <= n_step; i++) {
    domain.one_step(myran);
    if (i % 1000 == 0) {
      std::cout << i << "\t" << domain.cal_rho0() << "\t" << domain.cal_m_mean() << "\n";
      if (i % dn_out == 0) {
        domain.output_snap(fsnap);
      }
    }
  }
}



