#pragma once
#include "vect.h"
#include "comn.h"

class Par_2 {
public:
  Par_2(int spin0, int x, int y, int species0): spin(spin0), pos(x, y), species(species0) {}

  template <typename TRan>
  Par_2(TRan &myran, const Vec_2<int> &l, int species0);

  template <typename TRan>
  Par_2(TRan& myran, const Vec_2<int>& l, int spin0, int species0);

  Vec_2<short> pos;
  char spin;
  char species;
};

template<typename TRan>
Par_2::Par_2(TRan & myran, const Vec_2<int> &l, int species0): species(species0) {
  if (myran.doub() < 0.5) {
    spin = 1;
  } else {
    spin = -1;
  }
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

template<typename TRan>
Par_2::Par_2(TRan& myran, const Vec_2<int>& l, int spin0, int species0): spin(spin0), species(species0) {
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

class lattice_2 {
public:
  lattice_2(const Vec_2<int> &l, double Dt, double Dr,
            double phi_A, double phi_B, double rho_thresh, double v0,
            double etaAA, double etaAB, double etaBA, double etaBB);

  ~lattice_2();

  template <typename TRan>
  void ini_rand(TRan &myran);

  template <typename TRan>
  void ini_rand(TRan &myran, int spin0);

  int get_field_idx(const Par_2& p) const {
    return (int(p.pos.x) + p.pos.y * l_.x) * 2 + p.species;
  }

  void add_particle(const Par_2 &p) const;

  void del_particle(const Par_2 &p) const;

  double get_v(const Par_2 &p) const;

  void flip_hop(Par_2& p, double rand_val) const;

  void hop(Par_2 &p, double rand_val) const;

  double cal_m_mean() const;

  double cal_rho0() const;

  template <typename TRan>
  void one_step(TRan& myran);

  void output_snap(std::ofstream &fout);

private:
  Vec_2<int> l_;
  int n_sites_;
  unsigned short* rho_;
  short *m_;

  double Dt_;
  double Dr_;
  double v0_;
  double eta_[2][2]{};

  bool self_couplings_on_;
  double delta_t_;
  double prob_arr_[4]{};

  double phiA_;
  double phiB_;
  double rho_thresh_;
  int n_par_;
  int n_par_A_;
  std::vector<Par_2> p_arr_;
};

template <typename TRan>
void lattice_2::ini_rand(TRan& myran) {
  for (int i = 0; i < n_par_; i++) {
    if (i < n_par_A_) {
      p_arr_.emplace_back(myran, l_, 0);
    } else {
      p_arr_.emplace_back(myran, l_, 1);
    }
    add_particle(p_arr_.back());
  }
}

template <typename TRan>
void lattice_2::ini_rand(TRan& myran, int spin0) {
  for (int i = 0; i < n_par_; i++) {
    if (i < n_par_A_) {
      p_arr_.emplace_back(myran, l_, spin0, 0);
    } else {
      p_arr_.emplace_back(myran, l_, spin0, 1);
    }
    add_particle(p_arr_.back());
  }
}

template <typename TRan>
void lattice_2::one_step(TRan& myran) {
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    const double rand_val = myran.doub() / delta_t_;
    flip_hop(p_arr_[i], rand_val);
  }
}

void run(int Lx, int Ly,
         double phiA, double phiB, double rho_thresh,
         double Dt, double Dr, double v0,
         double etaAA, double etaAB, double etaBA, double etaBB,
         int n_step, int dn_out, int seed);

