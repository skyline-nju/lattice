#pragma once
#include "vect.h"
#include "comn.h"

class Par_6 {
public:
  Par_6(int spin0, int x, int y)
    : spin(spin0), pos(x, y) {}

  template <typename TRan>
  Par_6(TRan& myran, const Vec_2<int>& l);

  template <typename TRan>
  Par_6(TRan& myran, const Vec_2<int>& l, int spin0);

  Vec_2<short> pos;
  char spin;
};

template<typename TRan>
Par_6::Par_6(TRan& myran, const Vec_2<int>& l) {
  spin = int(myran.doub() * 6);
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

template<typename TRan>
Par_6::Par_6(TRan& myran, const Vec_2<int>& l, int spin0)
  : spin(spin0) {
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

class lattice_6 {
public:
  lattice_6(const Vec_2<int>& l, double Dt, double Dr,
    double phi, double rho_thresh, double v0, double eta);

  ~lattice_6();

  template <typename TRan>
  void ini_rand(TRan& myran);

  template <typename TRan>
  void ini_rand(TRan& myran, int spin0);

  int get_site_idx(const Par_6& p) const {
    return int(p.pos.x) + p.pos.y * l_.x;
  }

  void add_particle(const Par_6& p) const;

  void del_particle(const Par_6& p) const;

  double get_rho(const Par_6& p) const {
    return rho_[get_site_idx(p)];
  }

  double get_v(const Par_6& p) const {
    return v0_ * (1 + tanh(eta_ * (get_rho(p) - rho_thresh_)));
  }

  void hop_rot(Par_6& p, double rand_val) const;

  void hop(Par_6& p, const Vec_2<int>& ori) const;

  void rot(Par_6& p, int ds) const;

  double cal_rho0() const;

  template <typename TRan>
  void one_step(TRan& myran);

  void output_snap(std::ofstream& fout);

private:
  Vec_2<int> l_;
  int n_sites_;
  unsigned char* sigma_;
  unsigned short* rho_;

  double Dt_;
  double Dr_;
  double v0_;
  double eta_;

  double delta_t_;
  Vec_2<int> ori_[6][6]{};
  double prob_arr_[8]{};

  double phi_;
  double rho_thresh_;
  int n_par_;
  std::vector<Par_6> p_arr_;
};

template <typename TRan>
void lattice_6::ini_rand(TRan& myran) {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran, l_);
    add_particle(p_arr_.back());
  }
}

template <typename TRan>
void lattice_6::ini_rand(TRan& myran, int spin0) {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran, l_, spin0);
    add_particle(p_arr_.back());
  }
}

template <typename TRan>
void lattice_6::one_step(TRan& myran) {
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    const double rand_val = myran.doub() / delta_t_;
    hop_rot(p_arr_[i], rand_val);
  }
}

void run_q6(int Lx, int Ly,
            double phi, double rho_thresh,
            double Dt, double Dr,
            double v0, double eta,
            int n_step, int dn_out, int seed);


