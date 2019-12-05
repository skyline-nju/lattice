#pragma once
#include <iostream>
#include <vector>
#include "rand.h"
#define N_MAX 1

namespace PEP_2 {

struct LatticeCoord {
  int x;
  int y;
  int i;
};

const int ori[4][2] = { {1, 0}, {0, 1}, {-1, 0}, {0, -1} };


template <typename T>
class SquareLattice {
public:
  SquareLattice(int Lx, int Ly);

  int get_idx0(int x, int y) const { return N_MAX * (x + y * Lx_); }
  int get_idx(int x, int y, int i) const { return get_idx0(x, y) + i; }
  int get_idx(const LatticeCoord& coord) const { return N_MAX * (coord.x + coord.y * Lx_) + coord.i; }

  int get_n(int idx0) const;
  int get_n(int x, int y) const { return get_n(get_idx0(x, y)); }

  template <typename TRan>
  void create_rand(double phi, TRan& myran, std::vector<LatticeCoord>& par);

  bool add_new(int x, int y, T s);
  bool add_new(LatticeCoord& coord, T s);

  template <typename TRan>
  void tumble(int idx0, TRan& myran);
  template <typename TRan>
  void tumble(const LatticeCoord& coord, TRan& myran);

  void jump(LatticeCoord& coord);

  template <typename TRan>
  void update_rand_seq(std::vector<LatticeCoord>& par, TRan& myran, std::vector<int>& seq, double alpha);

protected:
  int Lx_;
  int Ly_;
  int n_tot_;
  T* lattice_;
};
template<typename T>
SquareLattice<T>::SquareLattice(int Lx, int Ly) :Lx_(Lx), Ly_(Ly) {
  n_tot_ = Lx_ * Ly_ * N_MAX;
  lattice_ = new T[n_tot_];
  for (int i = 0; i < n_tot_; i++) {
    lattice_[i] = 4;
  }
}

template<typename T>
int SquareLattice<T>::get_n(int idx0) const {
#if N_MAX == 1
  return lattice_[idx0] < 4;
#elif N_MAX == 2
  return static_cast<int>(lattice_[idx0] < 4)
    + static_cast<int>(lattice_[idx0 + 1] < 4);
#elif N_MAX == 3
  return static_cast<int>(lattice_[idx0] < 4)
    + static_cast<int>(lattice_[idx0 + 1] < 4)
    + static_cast<int>(lattice_[idx0 + 2] < 4);
#elif N_MAX == 4
  return static_cast<int>(lattice_[idx0] < 4)
    + static_cast<int>(lattice_[idx0 + 1] < 4)
    + static_cast<int>(lattice_[idx0 + 2] < 4)
    + static_cast<int>(lattice_[idx0 + 3] < 4);
#endif
}


template<typename T>
bool SquareLattice<T>::add_new(int x, int y, T s) {
  int idx0 = get_idx0(x, y);
  int n = get_n(idx0);
  if (n < N_MAX) {
    lattice_[idx0 + n] = s;
    return true;
  } else {
    return false;
  }
}

template<typename T>
bool SquareLattice<T>::add_new(LatticeCoord& coord, T s) {
  int idx0 = get_idx0(coord.x, coord.y);
  int n = get_n(idx0);
  if (n < N_MAX) {
    coord.i = n;
    lattice_[idx0 + n] = s;
    return true;
  } else {
    return false;
  }
}

template<typename T>
void SquareLattice<T>::jump(LatticeCoord& coord) {
  int idx = get_idx(coord);
  int s = lattice_[idx];
  int x_new = coord.x + ori[s][0];
  int y_new = coord.y + ori[s][1];

  if (x_new < 0) {
    x_new = Lx_ - 1;
  } else if (x_new >= Lx_) {
    x_new = 0;
  }

  if (y_new < 0) {
    y_new = Ly_ - 1;
  } else if (y_new >= Ly_) {
    y_new = 0;
  }

  int idx0_new = get_idx0(x_new, y_new);
  int i_new = -1;
#if N_MAX == 1
  if (lattice_[idx0_new] == 4) {
    i_new = 0;
  }
#elif N_MAX == 2
  if (lattice_[idx0_new] == 4) {
    i_new = 0;
  } else if (lattice_[idx0_new + 1] == 4) {
    i_new = 1;
  }
#elif N_MAX == 3
  if (lattice_[idx0_new] == 4) {
    i_new = 0;
  } else if (lattice_[idx0_new + 1] == 4) {
    i_new = 1;
  } else if (lattice_[idx0_new + 2] == 4) {
    i_new = 2;
  }
#elif N_MAX == 4
  if (lattice_[idx0_new] == 4) {
    i_new = 0;
  } else if (lattice_[idx0_new + 1] == 4) {
    i_new = 1;
  } else if (lattice_[idx0_new + 2] == 4) {
    i_new = 2;
  } else if (lattice_[idx0_new + 3] == 4) {
    i_new = 3;
  }
#endif

  if (i_new >= 0) {
    coord.x = x_new;
    coord.y = y_new;
    coord.i = i_new;
    lattice_[idx] = 4;
    lattice_[idx0_new + i_new] = s;
  }
}

template<typename T>
template<typename TRan>
void SquareLattice<T>::create_rand(double phi, TRan& myran, std::vector<LatticeCoord>& par) {
  int n_par = phi * Lx_ * Ly_;
  par.reserve(n_par);
  for (int j = 0; j < n_par; j++) {
    int state = static_cast<int>(myran.doub() * 4);
    while (true) {
      LatticeCoord coord;
      coord.x = static_cast<int>(myran.doub() * Lx_);
      coord.y = static_cast<int>(myran.doub() * Ly_);
      if (add_new(coord, state)) {
        par.push_back(coord);
        break;
      }
    }
  }
}

template<typename T>
template<typename TRan>
void SquareLattice<T>::tumble(int idx0, TRan& myran) {
  lattice_[idx0] = static_cast<T>(myran.doub() * 4);
}

template<typename T>
template<typename TRan>
void SquareLattice<T>::tumble(const LatticeCoord& coord, TRan& myran) {
  lattice_[get_idx(coord)] = static_cast<T>(myran.doub() * 4);
}

template<typename T>
template<typename TRan>
void SquareLattice<T>::update_rand_seq(std::vector<LatticeCoord>& par, TRan& myran, std::vector<int>& seq, double alpha) {
  shuffle(seq, myran);
  auto end = seq.cend();
  for (auto it = seq.cbegin(); it != end; ++it) {
    if (myran.doub() < alpha) {
      tumble(par[*it], myran);
    }
    jump(par[*it]);
  }
}
}
