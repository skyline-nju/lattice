#include "NRQS.h"
#include <cmath>
#include <fstream>
#include <iomanip>

lattice_2::lattice_2(const Vec_2<int>& l, double Dt, double Dr,
                     double phiA, double phiB, double rho_thresh, double v0,
                     double etaAA, double etaAB, double etaBA, double etaBB)
  :l_(l), n_sites_(l.x * l.y), Dt_(Dt), Dr_(Dr), v0_(v0), phiA_(phiA), phiB_(phiB), rho_thresh_(rho_thresh) {
  eta_[0][0] = etaAA / rho_thresh_;
  eta_[0][1] = etaAB / rho_thresh_;
  eta_[1][0] = etaBA / rho_thresh_;
  eta_[1][1] = etaBB / rho_thresh_;
  self_inhibition_on_ = !(etaAA == 0 && etaBB == 0);

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

  
  if (self_inhibition_on_) {
    delta_t_ = 1 / (prob_arr_[3] + Dt_ + 4 * v0);
  } else {
    delta_t_ = 1 / (prob_arr_[3] + Dt_ + 2 * v0);
  }
  
  rho_ = new unsigned short[n_sites_ * 2] {};
  m_ = new short[n_sites_ * 2] {};
  n_par_A_ = int(l_.x * l_.y * phiA_);
  n_par_ = n_par_A_ + int(l_.x * l_.y * phiB_);
  p_arr_.reserve(n_par_);
  std::cout << "phiA = " << phiA_ << "\n";
  std::cout << "phiB = " << phiB_ << "\n";
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

double lattice_2::get_v(const Par_2& p) const {
  int j = (int(p.pos.x) + p.pos.y * l_.x) * 2;
  double drho[2] = { rho_[j] - rho_thresh_, rho_[j + 1] - rho_thresh_ };
  int s1 = p.species;
  if (self_inhibition_on_) {
    return v0_ * (1 + tanh(eta_[s1][0] * drho[0])) * (1 + tanh(eta_[s1][1] * drho[1]));
  } else {
    int s2 = 1 - p.species;
    return v0_ * (1 + tanh(eta_[s1][s2] * drho[s2]));
  }
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
  for (int i = 0; i < n_sites_ * 2; i++) {
    m_sum += m_[i];
  }
  return m_sum / n_par_;
}

double lattice_2::cal_rho0() const {
  double rho_sum = 0;
  for (int i = 0; i < n_sites_ * 2; i++) {
    rho_sum += rho_[i];
  }
  return rho_sum / n_sites_;
}

void lattice_2::output_snap(std::ofstream& fout) {
  fout.write((const char*)&rho_[0], sizeof(unsigned short) * n_sites_ * 2);
  fout.write((const char*)&m_[0], sizeof(short) * n_sites_ * 2);
}

void run(int Lx, int Ly,
         double phiA, double phiB, double rho_thresh,
         double Dt, double Dr, double v0,
         double etaAA, double etaAB, double etaBA, double etaBB,
         int n_step, int dn_out, int seed) {
  Vec_2<int> l(Lx, Ly);
  Ranq2 myran(seed);
  lattice_2 domain(l, Dt, Dr, phiA, phiB, rho_thresh, v0, etaAA, etaAB, etaBA, etaBB);
  domain.ini_rand(myran);
  //domain.ini_rand(myran, 1);
  char folder[255] = "data";
  char basename[255];
  snprintf(basename, 255, "L%d_%d_Dr%g_Dt%g_e%g_%g_J%g_%g_v%g_r%g_s%llu",
    Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho_thresh, seed);

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



