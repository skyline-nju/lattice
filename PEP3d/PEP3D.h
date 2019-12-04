#pragma once
#include <iostream>
#include <vector>
#include "rand.h"
#define N_MAX 2

namespace PEP_3 {

struct LatticeCoord {
  int x;
  int y;
  int z;
  int i;
};

const int ori[6][3] = {
  {1, 0, 0},
  {-1, 0, 0},
  {0, 1, 0},
  {0, -1, 0},
  {0, 0, 1},
  {0, 0, -1}
};

template <typename T>
class SquareLattice {
public:
  SquareLattice(int Lx, int Ly, int Lz);

  size_t get_idx0(int x, int y, int z) const { return N_MAX * (x + y * Lx_ + z * LxLy_); }
  size_t get_idx(int x, int y, int z, int i) const { return get_idx0(x, y, z) + i; }
  size_t get_idx0(const LatticeCoord& coord) const { return N_MAX * (coord.x + coord.y * Lx_ + coord.z * LxLy_); }
  size_t get_idx(const LatticeCoord& coord) const { return N_MAX * (coord.x + coord.y * Lx_ + coord.z* LxLy_) + coord.i; }

  int get_n(size_t idx0) const;
  int get_n(int x, int y, int z) const { return get_n(get_idx0(x, y, z)); }

  template <typename TRan>
  void create_rand(double phi, TRan& myran, std::vector<LatticeCoord>& par);

  bool add_new(int x, int y, int z, T s);
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
  int Lz_;
  size_t LxLy_;
  size_t n_tot_;
  T* lattice_;
};
template<typename T>
SquareLattice<T>::SquareLattice(int Lx, int Ly, int Lz) :Lx_(Lx), Ly_(Ly), Lz_(Lz), LxLy_(Lx * Ly) {
  n_tot_ = Lx_ * LxLy_ * N_MAX;
  lattice_ = new T[n_tot_];
  for (size_t i = 0; i < n_tot_; i++) {
    lattice_[i] = 6;
  }
}

template<typename T>
int SquareLattice<T>::get_n(size_t idx0) const {
#if N_MAX == 1
  return lattice_[idx0] < 6;
#elif N_MAX == 2
  return static_cast<int>(lattice_[idx0] < 6)
    + static_cast<int>(lattice_[idx0 + 1] < 6);
#elif N_MAX == 3
  return static_cast<int>(lattice_[idx0] < 6)
    + static_cast<int>(lattice_[idx0 + 1] < 6)
    + static_cast<int>(lattice_[idx0 + 2] < 6);
#elif N_MAX == 4
  return static_cast<int>(lattice_[idx0] < 6)
    + static_cast<int>(lattice_[idx0 + 1] < 6)
    + static_cast<int>(lattice_[idx0 + 2] < 6)
    + static_cast<int>(lattice_[idx0 + 3] < 6);
#endif
}


template<typename T>
bool SquareLattice<T>::add_new(int x, int y, int z, T s) {
  size_t idx0 = get_idx0(x, y, z);
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
  size_t idx0 = get_idx0(coord);
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
  size_t idx = get_idx(coord);
  T s = lattice_[idx];
  int x_new = coord.x + ori[s][0];
  int y_new = coord.y + ori[s][1];
  int z_new = coord.z + ori[s][2];

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
  if (z_new < 0) {
    z_new = Lz_ - 1;
  } else if (z_new >= Lz_) {
    z_new = 0;
  }
  size_t idx0_new = get_idx0(x_new, y_new, z_new);
  int i_new = -1;
#if N_MAX == 1
  if (lattice_[idx0_new] == 6) {
    i_new = 0;
  }
#elif N_MAX == 2
  if (lattice_[idx0_new] == 6) {
    i_new = 0;
  } else if (lattice_[idx0_new + 1] == 6) {
    i_new = 1;
  }
#elif N_MAX == 3
  if (lattice_[idx0_new] == 6) {
    i_new = 0;
  } else if (lattice_[idx0_new + 1] == 6) {
    i_new = 1;
  } else if (lattice_[idx0_new + 2] == 6) {
    i_new = 2;
  }
#elif N_MAX == 4
  if (lattice_[idx0_new] == 6) {
    i_new = 0;
  } else if (lattice_[idx0_new + 1] == 6) {
    i_new = 1;
  } else if (lattice_[idx0_new + 2] == 6) {
    i_new = 2;
  } else if (lattice_[idx0_new + 3] == 6) {
    i_new = 3;
  }
#endif

  if (i_new >= 0) {
    coord.x = x_new;
    coord.y = y_new;
    coord.z = z_new;
    coord.i = i_new;
    lattice_[idx] = 6;
    lattice_[idx0_new + i_new] = s;
  }
}

template<typename T>
template<typename TRan>
void SquareLattice<T>::create_rand(double phi, TRan& myran, std::vector<LatticeCoord>& par) {
  int n_par = static_cast<int>(phi * n_tot_);
  par.reserve(n_par);
  for (int j = 0; j < n_par; j++) {
    int state = static_cast<int>(myran.doub() * 6);
    while (true) {
      LatticeCoord coord;
      coord.x = static_cast<int>(myran.doub() * Lx_);
      coord.y = static_cast<int>(myran.doub() * Ly_);
      coord.z = static_cast<int>(myran.doub() * Lz_);
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
  lattice_[idx0] = static_cast<T>(myran.doub() * 6);
}

template<typename T>
template<typename TRan>
void SquareLattice<T>::tumble(const LatticeCoord& coord, TRan& myran) {
  lattice_[get_idx(coord)] = static_cast<T>(myran.doub() * 6);
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
