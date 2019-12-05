#pragma once
#include "config.h"
#include "rand.h"

template <typename T>
struct Pos_2 {
  Pos_2() = default;
  Pos_2(T a, T b) : x(a), y(b) {}

  T x;
  T y;
};

template <typename T>
struct SquarePar_2 {
  SquarePar_2() = default;
  template <typename T2>
  SquarePar_2(const Pos_2<T2>& pos, unsigned char s0) : x(pos.x), y(pos.y), s(s0) {}
  Pos_2<T> jump() const { return Pos_2<T>(x + ori[s][0], y + ori[s][1]); }
  void reset_pos(const Pos_2<T>& pos) { x = pos.x; y = pos.y; }

  T x;
  T y;
  unsigned char s;

  const static int size_s = 4;
  const static short ori[4][2];
};

class SimpleSquareLattice_2 {
public:
  SimpleSquareLattice_2(int Lx, int Ly) : Lx_(Lx), Ly_(Ly) { n_par_ = new unsigned char[Lx * Ly](); }
  ~SimpleSquareLattice_2() { delete[] n_par_; }

  size_t get_vol() const { return Lx_ * Ly_; }

  template <typename T>
  unsigned char get_np(T x, T y) const { return n_par_[x + y * Lx_]; }

  template <typename TPar, typename TPos>
  void try_move(TPar& p, const TPos& pos_new);

  template <typename TRan, typename TPos>
  void add_one_randomly(TRan& myran, TPos& pos);

  int Lx_;
  int Ly_;

protected:
  unsigned char* n_par_;
};

class PackedSquareLattice_2 {
public:
  PackedSquareLattice_2(int Lx, int Ly, int nj, int base)
    : Lx_(Lx), Ly_(Ly), ncols_(Lx_), nrows_(Ly_ / nj), nj_(nj), base_(base) {
    n_par_ = new unsigned char[ncols_ * nrows_]();
  }

  ~PackedSquareLattice_2() { delete[] n_par_; }

  size_t get_vol() const { return Lx_ * Ly_; }

  template <typename TPar, typename TPos>
  void try_move(TPar& p, const TPos& pos_new);

  template <typename TRan, typename TPos>
  void add_one_randomly(TRan& myran, TPos& pos);

  int Lx_;
  int Ly_;

protected:
  int ncols_;
  int nrows_;
  int nj_;
  int base_;
  unsigned char* n_par_;

  const int base_arr_[4] = { 1, 3, 9, 27 };
};


template<typename TPar, typename TPos>
void SimpleSquareLattice_2::try_move(TPar& p, const TPos& pos_new) {
  size_t idx_new = pos_new.x + pos_new.y * Lx_;
  if (n_par_[idx_new] < N_MAX) {
    size_t idx_old = p.x + p.y * Lx_;
    n_par_[idx_new]++;
    n_par_[idx_old]--;
    p.reset_pos(pos_new);
  }
}

template<typename TRan, typename TPos>
void SimpleSquareLattice_2::add_one_randomly(TRan& myran, TPos& pos) {
  while (true) {
    pos.x = myran.doub() * Lx_;
    pos.y = myran.doub() * Ly_;
    size_t idx = pos.x + pos.y * Lx_;
    if (n_par_[idx] < N_MAX) {
      n_par_[idx]++;
      break;
    }
  }
}

template<typename TPar, typename TPos>
void PackedSquareLattice_2::try_move(TPar& p, const TPos& pos_new) {
  auto x_lat_new = pos_new.x;
  auto y_lat_new = pos_new.y / nj_;
  auto j_new = pos_new.y - y_lat_new * nj_;
  size_t idx_new = x_lat_new + y_lat_new * ncols_;
  unsigned char np_new = (n_par_[idx_new] / base_arr_[j_new]) % base_;
  if (np_new < N_MAX) {
    auto x_lat = p.x;
    auto y_lat = p.y / nj_;
    auto j = p.y - y_lat * nj_;
    size_t idx = x_lat + y_lat * ncols_;
    n_par_[idx_new] += base_arr_[j_new];
    n_par_[idx] -= base_arr_[j];
    p.reset_pos(pos_new);
  }
}

template<typename TRan, typename TPos>
void PackedSquareLattice_2::add_one_randomly(TRan& myran, TPos& pos) {
  while (true) {
    pos.x = myran.doub() * Lx_;
    pos.y = myran.doub() * Ly_;
    auto x_lat = pos.x;
    auto y_lat = pos.y / nj_;
    auto j = pos.y - y_lat * nj_;
    size_t idx = x_lat + y_lat * ncols_;
    unsigned char np = (n_par_[idx] / base_arr_[j]) % base_;
    if (np < N_MAX) {
      n_par_[idx] += base_arr_[j];
      break;
    }
  }
}

template<typename TPos, typename TLat>
void boundary_condi(TPos& p, const TLat& lat) {
  if (p.x < 0) {
    p.x = lat.Lx_ - 1;
  } else if (p.x >= lat.Lx_) {
    p.x = 0;
  }

  if (p.y < 0) {
    p.y = lat.Ly_ - 1;
  } else if (p.y >= lat.Ly_) {
    p.y = 0;
  }
}

template <typename TPar, typename TLattice, typename TRan>
void create_rand(double phi, std::vector<TPar>& p_arr, TLattice& lat, TRan& myran) {
  size_t n_par = phi * lat.get_vol();
  p_arr.reserve(n_par);
  for (int j = 0; j < n_par; j++) {
    unsigned char s_new = static_cast<unsigned char>(myran.doub() * TPar::size_s);
    Pos_2<short> pos_new;
    lat.add_one_randomly(myran, pos_new);
    p_arr.emplace_back(pos_new, s_new);
  }
}

template <typename TPar, typename TLattice, typename TRan>
void run_in_rand_seq(std::vector<TPar>& p_arr, TLattice& lat, TRan& myran,
  std::vector<size_t>& idx_arr, double alpha) {
  shuffle(idx_arr, myran);
  auto end = idx_arr.cend();
  for (auto it = idx_arr.cbegin(); it != end; ++it) {
    auto& p = p_arr[*it];
    if (myran.doub() < alpha) {
      p.s = static_cast<unsigned char>(myran.doub() * TPar::size_s);
    }
    auto pos_new = p.jump();
    boundary_condi(pos_new, lat);
    lat.try_move(p, pos_new);
  }
}
