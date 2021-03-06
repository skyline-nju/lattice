#pragma once
#include "PEP.h"
#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "config.h"
/**
 * @brief Basic class for exporting data.
 *
 * Define the timming to dump data.
 */
class ExporterBase {
public:
  ExporterBase() : n_step_(0) {}

  ExporterBase(int start, int n_step, int sep) : start_(start), n_step_(n_step) {
    set_lin_frame(start, n_step, sep);
  }

  void set_lin_frame(int start, int n_step, int sep);

  bool need_export(const int i_step);

protected:
  int n_step_;    // total steps to run
  int start_ = 0; // The first step 
private:
  std::vector<int> frames_arr_; // frames that need to export
  std::vector<int>::iterator frame_iter_;
};

/**
 * @brief Exporter to output log
 *
 * Output the parameters after the initialization.
 * Output the beginning and endding time of the simulation.
 * Record time every certain time steps.
 */
class LogExporter : public ExporterBase {
public:
  LogExporter(const std::string& outfile, int start, int n_step, int sep, int np);

  ~LogExporter();

  void record(int i_step);

  std::ofstream fout;
private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  int n_par_;
  int step_count_ = 0;
};


class XYExporter : public ExporterBase {
public:
  XYExporter(const std::string outfile, int start, int n_step, int sep, int Lx, int Ly)
    : ExporterBase(start, n_step, sep), fout_(outfile), Lx_(Lx), Ly_(Ly), sep_(sep) {}

  template <typename TPar>
  void dump(int i, const std::vector<TPar>& par_arr);
private:
  std::ofstream fout_;
  int Lx_;
  int Ly_;
  int sep_;
};

template <typename TPar>
void XYExporter::dump(int i_step, const std::vector<TPar>& par_arr) {
  //if (need_export(i_step)) {
  if (i_step % sep_ == 0) {
    int n_par = par_arr.size();
    fout_ << n_par << "\n";
    // comment line
    fout_ << "Lattice=\"" << Lx_ << " 0 0 0 " << Ly_ << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2 "
      << "Time=" << i_step;
    for (int j = 0; j < n_par; j++) {
      fout_ << "\n" << "N\t"
        << par_arr[j].x << "\t" << par_arr[j].y;
    }
    fout_ << std::endl;
  }
}

class XYZExporter : public ExporterBase {
public:
  XYZExporter(const std::string outfile, int start, int n_step, int sep, int Lx, int Ly, int Lz)
    : ExporterBase(start, n_step, sep), fout_(outfile), Lx_(Lx), Ly_(Ly), Lz_(Lz), sep_(sep) {}

  template <typename TPar>
  void dump(int i, const std::vector<TPar>& p_arr);
private:
  std::ofstream fout_;
  int Lx_;
  int Ly_;
  int Lz_;
  int sep_;
};

template <typename TPar>
void XYZExporter::dump(int i_step, const std::vector<TPar>& par_arr) {
  //if (need_export(i_step)) {
  if (i_step % sep_ == 0) {
    int n_par = par_arr.size();
    fout_ << n_par << "\n";
    // comment line
    fout_ << "Lattice=\"" << Lx_ << " 0 0 0 " << Ly_ << " 0 0 0 " << Lz_
      << "\" Properties=species:S:1:pos:R:3 "
      << "Time=" << i_step;
    for (int j = 0; j < n_par; j++) {
      fout_ << "\n" << "N\t"
        << par_arr[j].x << "\t" << par_arr[j].y << "\t" << par_arr[j].z;
    }
    fout_ << std::endl;
  }
}

class WettingProfileExporter_2 : public ExporterBase {
public:
  WettingProfileExporter_2(const std::string outfile, int start, int n_step, int sep, int Lx, int Ly)
    : ExporterBase(start, n_step, sep), fout_(outfile, std::ios::binary), Lx_(Lx), Ly_(Ly), sep_(sep) {
    buf_ = new unsigned char[Lx_];
  }

  ~WettingProfileExporter_2() {
    delete[] buf_;
    fout_.close();
  }
  template <typename TLat>
  void dump(int i_step, const TLat& lattice);
private:
  std::ofstream fout_;
  int Lx_;
  int Ly_;
  int sep_;
  unsigned char* buf_;
};

template <typename TLat>
void WettingProfileExporter_2::dump(int i_step, const TLat& lattice) {
  int t = i_step;
  fout_.write((char*)&t, sizeof(int));
  lattice.get_wetting_profile(buf_, 0);
  fout_.write((char*)buf_, sizeof(unsigned char) * Lx_);
  lattice.get_wetting_profile(buf_, Ly_ - 1);
  fout_.write((char*)buf_, sizeof(unsigned char) * Lx_);
}