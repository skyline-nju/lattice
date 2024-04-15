#pragma once
#include "vect.h"
#include "comn.h"

class Par_2 {
public:
  Par_2(int spin0, int x, int y): spin(spin0), pos(x, y) {}

  template <typename TRan>
  Par_2(TRan &myran, const Vec_2<int> &l);

  template <typename TRan>
  Par_2(TRan& myran, const Vec_2<int>& l, int spin0);

  template <typename TRan>
  Par_2(TRan& myran, const Vec_2<int>& l, const Vec_2<int>& R);

  Vec_2<short> pos;
  char spin;
};

template<typename TRan>
Par_2::Par_2(TRan & myran, const Vec_2<int> &l) {
  spin = int(myran.doub() * 6);
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

template<typename TRan>
Par_2::Par_2(TRan& myran, const Vec_2<int>& l, int spin0): spin(spin0) {
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

template<typename TRan>
Par_2::Par_2(TRan& myran, const Vec_2<int>& l, const Vec_2<int>& R) {
  pos.x = R.x;
  pos.y = R.y;
  spin = int(myran.doub() * 6);
}

class hexa_lattice_2 {
public:
  hexa_lattice_2(const Vec_2<int> &l, double phi,
                 double rate_P, double rate_D, double rate_R);

  template <typename TRan>
  void ini_rand(TRan &myran);

  int get_site_idx(const Par_2& p, int x_offset = 0, int y_offset = 0) const;

  void add_particle(const Par_2 &p);

  void del_particle(const Par_2 &p);

  void hop(Par_2& p, const Vec_2<int>& dR);

  void hop_rot(Par_2& p, double rand_val);

  const std::vector<Par_2>& get_p_arr() const { return p_arr_;};

  int get_n_par() const { return n_par_; }

  template <typename TRan>
  void one_step(TRan& myran);

private:
  Vec_2<int> l_;
  int n_sites_;
 
  double rate_P_;
  double rate_D_;
  double rate_R_;

  double delta_t_;

  double phi_;
  int n_par_;

  double prob_arr_[6][8]{};
  Vec_2<int> ori_arr_[6][6]{};

  std::vector<unsigned char> site_state_;
  std::vector<Par_2> p_arr_;
};

template <typename TRan>
void hexa_lattice_2::ini_rand(TRan& myran) {
  int* site_idx = new int[n_sites_] {};
  for (int i = 0; i < n_sites_; i++) {
    site_idx[i] = i;
  }
  shuffle(site_idx, n_sites_, myran);
  for (int i = 0; i < n_par_; i++) {
    Vec_2<int> R;
    R.x = site_idx[i] % l_.x;
    R.y = site_idx[i] / l_.x;
    p_arr_.emplace_back(myran, l_, R);
    if (site_state_[site_idx[i]] != 0) {
      std::cout << "Error, site (" << R.x << ", " << R.y << ") has been occupied!" << std::endl;
      exit(1);
    } else {
      site_state_[site_idx[i]] = 1;
    }
  }
  delete[] site_idx;
  std::cout << "initialize " << n_par_ << " particles randomly" << std::endl;
}

template <typename TRan>
void hexa_lattice_2::one_step(TRan& myran) {
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    const double rand_val = myran.doub() / delta_t_;
    hop_rot(p_arr_[i], rand_val);
  }
}

void run(int Lx, int Ly, double phi,
         double rate_P, double rate_D, double rate_R,
         int n_step, int dn_out, int seed);

