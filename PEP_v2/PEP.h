#pragma once
#include "config.h"
#include "rand.h"

struct Par_2 {
public:
  Par_2() = default;
  Par_2(short a, short b, unsigned char s) : x(a), y(b), ori(s) {}
  short x;
  short y;
  unsigned char ori;
};

class SquareLattice_2 {
public:
  SquareLattice_2(int Lx0, int Ly0)
    : Lx_(Lx0), Ly_(Ly0), tot_lattices_(static_cast<size_t>(Lx_) * Ly_) { 
    n_ = new unsigned char[tot_lattices_](); }

  ~SquareLattice_2() { delete[] n_; }

  template <typename TRan>
  void create_particles_randomly(std::vector<Par_2>& p_arr, TRan& myran, double phi);

  template <typename TRan>
  void tumble(Par_2& p, TRan& myran, double alpha) const;  
  template <typename TPar>
  void move(TPar& p);

  template <typename T>
  void boundary_condi(T &x, T &y) const;

private:
  int Lx_;
  int Ly_;
  const size_t tot_lattices_;
  unsigned char* n_;
  const int n_oris_ = 4;
  const short ori_arr_[4][2] = { {1, 0}, {0, 1}, {-1, 0}, {0, -1} };
};

template <typename TRan>
void SquareLattice_2::create_particles_randomly(std::vector<Par_2>& p_arr, TRan& myran, double phi) {
  size_t n_par = phi * tot_lattices_;
  p_arr.reserve(n_par);
  for (int j = 0; j < n_par; j++) {
    unsigned char ori = static_cast<unsigned char>(myran.doub() * n_oris_);
    short x;
    short y;
    while (true) {
      x = static_cast<short>(myran.doub() * Lx_);
      y = static_cast<short>(myran.doub() * Ly_);
      size_t idx = x + y * Lx_;
      if (n_[idx] < N_MAX) {
        n_[idx]++;
        break;
      }
    }
    p_arr.emplace_back(x, y, ori);
  }
}

template <typename T>
void SquareLattice_2::boundary_condi(T &x, T&y) const {
  if (x < 0) {
    x = Lx_ - 1;
  } else if (x >= Lx_) {
    x = 0;
  }

  if (y < 0) {
    y = Ly_ - 1;
  } else if (y >= Ly_) {
    y = 0;
  }
}

template<typename TRan>
void SquareLattice_2::tumble(Par_2& p, TRan& myran, double alpha) const {
  if (myran.doub() < alpha) {
    p.ori = static_cast<unsigned char>(myran.doub() * n_oris_);
  }
}

template <typename TPar>
void SquareLattice_2::move(TPar& p) {
  short x = p.x + ori_arr_[p.ori][0];
  short y = p.y + ori_arr_[p.ori][1];
  boundary_condi(x, y);
  size_t idx1 = x + y * static_cast<size_t>(Lx_);
  if (n_[idx1] < N_MAX) {
    size_t idx0 = p.x + p.y * static_cast<size_t>(Lx_);
    n_[idx0]--;
    n_[idx1]++;
    p.x = x;
    p.y = y;
  }
}

template <typename TPar, typename TLattice, typename TRan>
void update_in_rand_seq(std::vector<TPar>& p_arr, TLattice& lat, TRan& myran,
  std::vector<size_t>& idx_arr, double alpha) {
  auto end = idx_arr.cend();
  for (auto it = idx_arr.cbegin(); it != end; ++it) {
    auto& p = p_arr[*it];
    lat.tumble(p, myran, alpha);
    lat.move(p);
  }
}