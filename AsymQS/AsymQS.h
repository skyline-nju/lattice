#pragma once
#include "vect.h"
#include "comn.h"

class Par_2 {
public:
  Par_2(int spin0, int x, int y, int species0): spin(spin0), pos(x, y) {}

  template <typename TRan>
  Par_2(TRan &myran, const Vec_2<int> &l);

  template <typename TRan>
  Par_2(TRan& myran, const Vec_2<int>& l, int spin0);

  Vec_2<short> pos;
  char spin;
};

template<typename TRan>
Par_2::Par_2(TRan & myran, const Vec_2<int> &l) {
  if (myran.doub() < 0.5) {
    spin = 1;
  } else {
    spin = -1;
  }
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

template<typename TRan>
Par_2::Par_2(TRan& myran, const Vec_2<int>& l, int spin0): spin(spin0) {
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

class lattice_2 {
public:
  lattice_2(const Vec_2<int> &l, double Dt, double Dr,
            double phi, double rho_thresh, double v0,
            double eta, double alpha);

  ~lattice_2();

  template <typename TRan>
  void ini_rand(TRan &myran);

  template <typename TRan>
  void ini_rand(TRan &myran, int spin0);

  int get_field_idx(const Par_2& p, int offset = 0) const;

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
  double eta_;
  double alpha_;

  double delta_t_;
  double prob_arr_[4]{};

  double phi_;
  double rho_thresh_;
  int n_par_;
  std::vector<Par_2> p_arr_;
};

template <typename TRan>
void lattice_2::ini_rand(TRan& myran) {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran, l_);
    add_particle(p_arr_.back());
  }
}

template <typename TRan>
void lattice_2::ini_rand(TRan& myran, int spin0) {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran, l_, spin0);
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
         double phi, double rho_thresh,
         double Dt, double Dr, double v0,
         double eta, double alpha,
         int n_step, int dn_out, int seed);

