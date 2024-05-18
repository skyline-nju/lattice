#ifndef COMN_H
#define COMN_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <ctime>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

const double PI = 3.14159265358979;

// check whether file or folder exists
bool exist(const char *filename);

#ifdef _MSC_VER
const std::string delimiter("\\");
#else
const std::string delimiter("/");
#endif

// create folder
void mkdir(const char *folder);

// split string by a delimiter
std::vector<std::string> split(
  const std::string &str, const std::string &dlm);

template <class T>
void str_to_num(const std::string str, T &num) {
  std::stringstream ss;
  ss << str;
  ss >> num;
}

template <class T>
void num_to_str(const T &num, std::string str) {
  std::stringstream ss;
  ss << num;
  ss >> str;
}

template<typename T1, typename T2>
void tangle_1(T1& x, const T2 L) {
  if (x < T2(0)) {
    x += L;
  } else if (x >= L) {
    x -= L;
  }
}

#endif

