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

class SquareLattice_2 {
public:
  SquareLattice_2(int Lx0, int Ly0)
    : Lx_(Lx0), Ly_(Ly0), tot_lattices_(static_cast<size_t>(Lx_) * Ly_) { 
    n_ = new unsigned char[tot_lattices_](); }

  ~SquareLattice_2() { delete[] n_; }

  template <typename TRan>
  void create_particles_randomly(std::vector<Par_2>& p_arr, TRan& myran, double phi);

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


class SquareLattice_3 {
public:
  SquareLattice_3(int Lx0, int Ly0, int Lz0)
    : Lx_(Lx0), Ly_(Ly0), Lz_(Lz0), LxLy_(static_cast<size_t>(Lx_)* Ly_),
    tot_lattices_(LxLy_* Lz_) {
    n_ = new unsigned char[tot_lattices_]();
  }

  ~SquareLattice_3() { delete[] n_; }

  int get_n_ori() const { return n_oris_; }

  size_t get_idx(short x, short y, short z) const { return x + y * Lx_ + z * LxLy_; }
  size_t get_idx(const Par_3& p) const { return p.x + p.y * Lx_ + p.z * LxLy_; }

  template <typename TRan>
  void create_particles_randomly(std::vector<Par_3>& p_arr, TRan& myran, double phi);

  template <typename TPar>
  void move(TPar& p);

  template <typename T>
  void boundary_condi(T& x, T& y, T& z) const;

  const short ori_arr_[6][3] = { {1, 0, 0}, {-1, 0, 0},
                                 {0, 1, 0}, {0, -1, 0},
                                 {0, 0, 1}, {0, 0, -1} };
  unsigned char* n_;
private:
  int Lx_;
  int Ly_;
  int Lz_;
  size_t LxLy_;
  const size_t tot_lattices_;
  const int n_oris_ = 6;
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

template <typename TRan>
void SquareLattice_3::create_particles_randomly(std::vector<Par_3>& p_arr, TRan& myran, double phi) {
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
void SquareLattice_3::boundary_condi(T& x, T& y, T& z) const {
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
    z = Lz_ - 1;
  } else if (z >= Lz_) {
    z = 0;
  }
}

template <typename TPar>
void SquareLattice_3::move(TPar& p) {
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

template <typename TPar, typename TRan>
void tumble(TPar& p, TRan& myran) {
  p.ori = static_cast<unsigned char>(myran.doub() * TPar::size_ori);
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

template <typename TPar, typename TLattice, typename TRan>
void update_in_rand_seq(std::vector<TPar>& p_arr, TLattice& lat, TRan& myran,
  std::vector<size_t>& idx_arr, unsigned long long alpha) {
  auto end = idx_arr.cend();
  for (auto it = idx_arr.cbegin(); it != end; ++it) {
    auto& p = p_arr[*it];
    if (myran.int64() < alpha) {
      tumble(p, myran);
    }
    lat.move(p);
  }
}

template <typename TPar, typename TRan>
void update_in_rand_seq(std::vector<TPar>& p_arr, SquareLattice_3& lat, TRan& myran,
  std::vector<size_t>& idx_arr, unsigned long long alpha) {
  auto end = idx_arr.cend();
  for (auto it = idx_arr.cbegin(); it != end; ++it) {
    auto& p = p_arr[*it];
    if (myran.int64() < alpha) {
      tumble(p, myran);
    }
    //lat.move(p);
    short x = p.x + lat.ori_arr_[p.ori][0];
    short y = p.y + lat.ori_arr_[p.ori][1];
    short z = p.z + lat.ori_arr_[p.ori][2];
    lat.boundary_condi(x, y, z);
    size_t idx1 = lat.get_idx(x, y, z);
    if (lat.n_[idx1] < N_MAX) {
      lat.n_[idx1]++;
      size_t idx0 = lat.get_idx(p);
      lat.n_[idx0]--;
      p.x = x;
      p.y = y;
      p.z = z;
    }
  }
}

template <typename TPar, typename TRan>
void update_in_rand_seq1(std::vector<TPar>& p_arr, SquareLattice_3& lat, TRan& myran,
  std::vector<size_t>& idx_arr, unsigned long long alpha) {
  auto end = idx_arr.cend();
  for (auto it = idx_arr.cbegin(); it != end; ++it) {
    auto& p = p_arr[*it];
    if (myran.int64() < alpha) {
      tumble(p, myran);
    }
    //lat.move(p);
    short x = p.x + lat.ori_arr_[p.ori][0];
    short y = p.y + lat.ori_arr_[p.ori][1];
    short z = p.z + lat.ori_arr_[p.ori][2];
    lat.boundary_condi(x, y, z);
    size_t idx1 = lat.get_idx(x, y, z);
    if (lat.n_[idx1] < N_MAX) {
      size_t idx0 = lat.get_idx(p);
      lat.n_[idx1]++;
      lat.n_[idx0]--;
      p.x = x;
      p.y = y;
      p.z = z;
    }
  }
}

template <typename TPar, typename TRan>
void update_in_rand_seq2(std::vector<TPar>& p_arr, SquareLattice_3& lat, TRan& myran,
  std::vector<size_t>& idx_arr, unsigned long long alpha) {
  auto end = idx_arr.cend();
  for (auto it = idx_arr.cbegin(); it != end; ++it) {
    auto& p = p_arr[*it];
    if (myran.int64() < alpha) {
      tumble(p, myran);
    }
    //lat.move(p);
    short x = p.x + lat.ori_arr_[p.ori][0];
    short y = p.y + lat.ori_arr_[p.ori][1];
    short z = p.z + lat.ori_arr_[p.ori][2];
    lat.boundary_condi(x, y, z);
    size_t idx1 = lat.get_idx(x, y, z);
    if (lat.n_[idx1] < N_MAX) {
      size_t idx0 = lat.get_idx(p);
      lat.n_[idx0]--;
      lat.n_[idx1]++;
      p.x = x;
      p.y = y;
      p.z = z;
    }
  }
}

template <typename TPar, typename TLat, typename TRan>
void update_in_rand_seq3(std::vector<TPar>& p_arr, TLat& lat, TRan& myran,
  std::vector<size_t>& idx_arr, unsigned long long alpha) {
  auto end = idx_arr.cend();
  for (auto it = idx_arr.cbegin(); it != end; ++it) {
    auto& p = p_arr[*it];
    if (myran.int64() < alpha) {
      tumble(p, myran);
    }
    lat.move(p);
  }
}