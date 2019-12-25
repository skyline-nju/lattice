#pragma once
#include "config.h"
#include "rand.h"
#include <iostream>

struct Par_2 {
public:
  Par_2() = default;
  Par_2(short a, short b, unsigned char s) : x(a), y(b), ori(s) {}
  short x;
  short y;
  unsigned char ori;

  static int size_ori;
};

struct Par_3 {
public:
  Par_3() = default;
  Par_3(short a, short b, short c, unsigned char s) : x(a), y(b), z(c), ori(s) {}
  short x;
  short y;
  short z;
  unsigned char ori;

  static int size_ori;
};

class SquareLattice {
public:
  SquareLattice(int Lx0, int Ly0)
    : Lx_(Lx0), Ly_(Ly0), tot_lattices_(static_cast<size_t>(Lx_) * Ly_) { 
    n_ = new unsigned char[tot_lattices_](); }

  ~SquareLattice() { delete[] n_; }

  template <typename TRan>
  void create_particles_randomly(std::vector<Par_2>& p_arr, TRan& myran, double phi);

  template <typename TPar>
  void move(TPar& p);

  template <typename TPar, typename TRan>
  void move_rand_reflect(TPar& p, TRan& myran) {}

  template <typename T>
  void boundary_condi(T &x, T &y) const;

  int count_empty_sites(int row) const;

  template <typename IntType>
  void get_wetting_profile(IntType* num, int row) const;

  bool near_bottom_wall(const Par_2& p) const { return p.y == 0; }
private:
  int Lx_;
  int Ly_;
  const size_t tot_lattices_;
  unsigned char* n_;
  const int n_oris_ = 4;
  const short ori_arr_[4][2] = { {1, 0}, {-1, 0}, {0, 1}, {0, -1} };
};


class CubicLattice {
public:
  CubicLattice(int Lx0, int Ly0, int Lz0)
    : Lx_(Lx0), Ly_(Ly0), Lz_(Lz0), LxLy_(static_cast<size_t>(Lx_)* Ly_),
    tot_lattices_(LxLy_* Lz_) {
    n_ = new unsigned char[tot_lattices_]();
  }

  ~CubicLattice() { delete[] n_; }

  int get_n_ori() const { return n_oris_; }

  size_t get_idx(short x, short y, short z) const { return x + y * Lx_ + z * LxLy_; }
  size_t get_idx(const Par_3& p) const { return p.x + p.y * Lx_ + p.z * LxLy_; }

  template <typename TRan>
  void create_particles_randomly(std::vector<Par_3>& p_arr, TRan& myran, double phi);

  template <typename TPar>
  void move(TPar& p);

  template <typename TPar, typename TRan>
  void move_rand_reflect(TPar& p, TRan& myran);

  template <typename T>
  void boundary_condi(T& x, T& y, T& z) const;

  bool near_bottom_wall(const Par_3& p) const { return p.z == 0; }

private:
  int Lx_;
  int Ly_;
  int Lz_;
  size_t LxLy_;
  const size_t tot_lattices_;
  unsigned char* n_;
  const int n_oris_ = 6;
  const short ori_arr_[6][3] = { {1, 0, 0}, {-1, 0, 0},
                                 {0, 1, 0}, {0, -1, 0},
                                 {0, 0, 1}, {0, 0, -1} };
};

template <typename TRan>
void SquareLattice::create_particles_randomly(std::vector<Par_2>& p_arr, TRan& myran, double phi) {
  size_t n_par = phi * tot_lattices_;
  p_arr.reserve(n_par);
  for (int j = 0; j < n_par; j++) {
    unsigned char ori = static_cast<unsigned char>(myran.doub() * n_oris_);
    short x;
    short y;
    while (true) {
      x = static_cast<short>(myran.doub() * Lx_);
      y = static_cast<short>(myran.doub() * Ly_);
      size_t idx = x + y * static_cast<size_t>(Lx_);
      if (n_[idx] < N_MAX) {
        n_[idx]++;
        break;
      }
    }
    p_arr.emplace_back(x, y, ori);
  }
  Par_2::size_ori = n_oris_;
}

template <typename T>
void SquareLattice::boundary_condi(T &x, T&y) const {
  if (x < 0) {
    x = Lx_ - 1;
  } else if (x >= Lx_) {
    x = 0;
  }

  if (y < 0) {
#ifdef WALL_TOP_BOTTOM
    y = 0;
#else
    y = Ly_ - 1;
#endif
  } else if (y >= Ly_) {
#ifdef WALL_TOP_BOTTOM
    y = Ly_ - 1;
#else
    y = 0;
#endif
  }
}

template <typename TPar>
void SquareLattice::move(TPar& p) {
  short x = p.x + ori_arr_[p.ori][0];
  short y = p.y + ori_arr_[p.ori][1];
  boundary_condi(x, y);
  size_t idx1 = x + y * static_cast<size_t>(Lx_);
  if (n_[idx1] < N_MAX) {
    n_[idx1]++;
    size_t idx0 = p.x + p.y * static_cast<size_t>(Lx_);
    n_[idx0]--;
    p.x = x;
    p.y = y;
  }
}

template <typename IntType>
void SquareLattice::get_wetting_profile(IntType* num, int row) const{
  for (int i = 0; i < Lx_; i++) {
    num[i] = 0;
  }
  int drow;
  if (row == 0) {
    drow = 1;
  } else if (row == Ly_ - 1) {
    drow = -1;
  } else {
    std::cout << "Error! The value of row must be 0 or Ly-1" << std::endl;
    exit(1);
  }

  for (int i = 0; i < Lx_; i++) {
    size_t row_cur = row;
    while (true) {
      size_t idx = i + row_cur * Lx_;
      if (n_[idx] > 0) {
        num[i] += n_[idx];
        row_cur += drow;
      } else {
        break;
      }
    }
  }
}

template <typename TRan>
void CubicLattice::create_particles_randomly(std::vector<Par_3>& p_arr, TRan& myran, double phi) {
  size_t n_par = phi * tot_lattices_;
  p_arr.reserve(n_par);
  for (int j = 0; j < n_par; j++) {
    unsigned char ori = static_cast<unsigned char>(myran.doub() * n_oris_);
    short x;
    short y;
    short z;
    while (true) {
      x = static_cast<short>(myran.doub() * Lx_);
      y = static_cast<short>(myran.doub() * Ly_);
      z = static_cast<short>(myran.doub() * Lz_);
      size_t idx = x + y * Lx_ + z * LxLy_;
      if (n_[idx] < N_MAX) {
        n_[idx]++;
        break;
      }
    }
    p_arr.emplace_back(x, y, z, ori);
  }
  Par_3::size_ori = n_oris_;
}

template <typename T>
void CubicLattice::boundary_condi(T& x, T& y, T& z) const {
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

  if (z < 0) {
#ifdef WALL_TOP_BOTTOM
    z = 0;
#else
    z = Lz_ - 1;
#endif
  } else if (z >= Lz_) {
#ifdef WALL_TOP_BOTTOM
    z = Lz_ - 1;
#else
    z = 0;
#endif
  }
}

template <typename TPar>
void CubicLattice::move(TPar& p) {
  short x = p.x + ori_arr_[p.ori][0];
  short y = p.y + ori_arr_[p.ori][1];
  short z = p.z + ori_arr_[p.ori][2];
  boundary_condi(x, y, z);
  size_t idx1 = get_idx(x, y, z);
  if (n_[idx1] < N_MAX) {
    n_[idx1]++;
    size_t idx0 = get_idx(p);
    n_[idx0]--;
    p.x = x;
    p.y = y;
    p.z = z;
  }

}

template<typename TPar, typename TRan>
void CubicLattice::move_rand_reflect(TPar& p, TRan& myran) {

  short z = p.z + ori_arr_[p.ori][2];
  bool flag_reflect = false;
  if (z < 0) {
    z = 0;
  } else if (z >= Lz_) {
    z = Lz_ - 1;
    flag_reflect = true;
  }

  if (flag_reflect) {
    size_t idx1;
    short x, y;
    do {
      x = static_cast<short>(Lx_ * myran.doub());
      y = static_cast<short>(Ly_ * myran.doub());
      idx1 = get_idx(x, y, z);
    } while (n_[idx1] >= N_MAX);
    n_[idx1]++;
    size_t idx0 = get_idx(p);
    n_[idx0]--;
    p.x = x;
    p.y = y;
    p.z = z;
    p.ori = 5;
  } else {
    short x = p.x + ori_arr_[p.ori][0];
    short y = p.y + ori_arr_[p.ori][1];
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
    size_t idx1 = get_idx(x, y, z);
    if (n_[idx1] < N_MAX) {
      n_[idx1]++;
      size_t idx0 = get_idx(p);
      n_[idx0]--;
      p.x = x;
      p.y = y;
      p.z = z;
    }
  }


}

template <typename TPar, typename TRan>
void tumble(TPar& p, TRan& myran) {
  p.ori = static_cast<unsigned char>(myran.doub() * TPar::size_ori);
}

template <typename TRan>
void tumble(Par_2& p, TRan& myran, const double* rates) {
  double rval = myran.doub();
  if (rval < rates[0]) {
    p.ori = 0;
  } else if (rval < rates[1]) {
    p.ori = 1;
  } else if (rval < rates[2]) {
    p.ori = 2;
  } else {
    p.ori = 3;
  }
}

template <typename TRan>
void tumble(Par_3& p, TRan& myran, const double* rates) {
  double rval = myran.doub();
  if (rval < rates[0]) {
    p.ori = 0;
  } else if (rval < rates[1]) {
    p.ori = 1;
  } else if (rval < rates[2]) {
    p.ori = 2;
  } else if (rval < rates[3]){
    p.ori = 3;
  } else if (rval < rates[4]) {
    p.ori = 4;
  } else {
    p.ori = 5;
  }
}

template <typename TPar, typename TLattice, typename TRan>
void update_in_rand_seq(std::vector<TPar>& p_arr, TLattice& lat, TRan& myran,
  std::vector<size_t>& idx_arr, double alpha) {
  auto end = idx_arr.cend();
  for (auto it = idx_arr.cbegin(); it != end; ++it) {
    auto& p = p_arr[*it];
    if (myran.doub() < alpha) {
      tumble(p, myran);
    }
    lat.move(p);
  }
}

template <typename TPar, typename TLat, typename TRan>
void update_in_rand_seq(std::vector<TPar>& p_arr, TLat& lat, TRan& myran,
                        std::vector<size_t>& idx_arr,
                        unsigned long long alpha) {
  auto end = idx_arr.cend();
  for (auto it = idx_arr.cbegin(); it != end; ++it) {
    auto& p = p_arr[*it];
    if (myran.int64() < alpha) {
      tumble(p, myran);
    }
#ifndef WALL_BOTTOM
    lat.move(p);
#else
    lat.move_rand_reflect(p, myran);
#endif
  }
}

template <typename TPar, typename TLat, typename TRan>
void update_in_rand_seq(std::vector<TPar>& p_arr, TLat& lat, TRan& myran,
                        std::vector<size_t>& idx_arr,
                        unsigned long long alpha, const double* rates) {
  auto end = idx_arr.cend();
  for (auto it = idx_arr.cbegin(); it != end; ++it) {
    auto& p = p_arr[*it];
    if (myran.int64() < alpha) {
      if (lat.near_bottom_wall(p)) {
        tumble(p, myran, rates);
      } else {
        tumble(p, myran);
      }
    }
#ifndef WALL_BOTTOM
    lat.move(p);
#else
    lat.move_rand_reflect(p, myran);
#endif
  }
}

void test_rand();

void test_3D(int Lx);
