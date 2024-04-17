#pragma once
#include "vect.h"
#include "comn.h"

class Par_4 {
public:
  Par_4(int spin0, int x, int y, int species0)
    : spin(spin0), pos(x, y), species(species0) {}

  template <typename TRan>
  Par_4(TRan& myran, const Vec_2<int>& l, int species0);

  template <typename TRan>
  Par_4(TRan& myran, const Vec_2<int>& l, int spin0, int species0);

  Vec_2<short> pos;
  char spin;
  char species;
};

template<typename TRan>
Par_4::Par_4(TRan& myran, const Vec_2<int>& l, int species0)
  : species(species0) {
  spin = int(myran.doub() * 4);
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

template<typename TRan>
Par_4::Par_4(TRan& myran, const Vec_2<int>& l, int spin0, int species0)
  : spin(spin0), species(species0) {
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

class lattice_4 {
public:
  lattice_4(const Vec_2<int>& l, double Dt, double Dr,
    double phi_A, double phi_B, double rho_thresh, double v0,
    double etaAA, double etaAB, double etaBA, double etaBB);

  ~lattice_4();

  template <typename TRan>
  void ini_rand(TRan& myran);

  template <typename TRan>
  void ini_rand(TRan& myran, int spin0);

  int get_sigma_idx(const Par_4& p) const {
    return (int(p.pos.x) + p.pos.y * l_.x) * 8 + p.species * 4 + p.spin;
  }

  void add_particle(const Par_4& p) const { sigma_[get_sigma_idx(p)] += 1; }

  void del_particle(const Par_4& p) const { sigma_[get_sigma_idx(p)] -= 1; }

  void get_rho(const Vec_2<short>& pos, double& rho_A, double& rho_B) const;

  double get_v(const Par_4& p) const;

  void hop_rot(Par_4& p, double rand_val) const;

  void hop(Par_4& p, const Vec_2<int>& ori) const;

  void rot(Par_4& p, int ds) const;

  double cal_rho0() const;

  template <typename TRan>
  void one_step(TRan& myran);

  void output_snap(std::ofstream& fout);

private:
  Vec_2<int> l_;
  int n_sites_;
  unsigned short* sigma_;
  //unsigned short* rho_;

  double Dt_;
  double Dr_;
  double v0_;
  double eta_[2][2]{};

  bool self_inhibition_on_;
  double delta_t_;
  Vec_2<int> ori_[4][4]{};
  double prob_arr_[6]{};

  double phiA_;
  double phiB_;
  double rho_thresh_;
  int n_par_;
  int n_par_A_;
  std::vector<Par_4> p_arr_;
};

template <typename TRan>
void lattice_4::ini_rand(TRan& myran) {
  for (int i = 0; i < n_par_; i++) {
    if (i < n_par_A_) {
      p_arr_.emplace_back(myran, l_, 0);
    } else {
      p_arr_.emplace_back(myran, l_, 1);
    }
    add_particle(p_arr_.back());
  }

  for (int i = 0; i < n_par_; i++) {
    int j = get_sigma_idx(p_arr_[i]);
    if (sigma_[j] == 0) {
      std::cout << "j = " << j << std::endl;
      exit(1);
    }
  }
}

template <typename TRan>
void lattice_4::ini_rand(TRan& myran, int spin0) {
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
void lattice_4::one_step(TRan& myran) {
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    const double rand_val = myran.doub() / delta_t_;
    hop_rot(p_arr_[i], rand_val);
  }
}

void run_q4(int Lx, int Ly,
  double phiA, double phiB, double rho_thresh,
  double Dt, double Dr, double v0,
  double etaAA, double etaAB, double etaBA, double etaBB,
  int n_step, int dn_out, int seed);


