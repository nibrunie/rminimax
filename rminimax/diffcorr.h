#ifndef RMINIMAX_DIFFCORR
#define RMINIMAX_DIFFCORR

#include <cstdio>
#include <cstdlib>
#include <functional>
#include <gmp.h>
#include <iostream>
#include <mpfr.h>
#include <mpreal.h>
#include <utility>
#include <vector>

bool adiffcorr(std::vector<mpfr::mpreal> &num, std::vector<mpfr::mpreal> &den,
               mpfr::mpreal &errDC,
               std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nb,
               std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &db,
               std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dbounds,
               std::function<mpfr::mpreal(mpfr::mpreal)> &nfix,
               std::function<mpfr::mpreal(mpfr::mpreal)> &dfix,
               std::vector<mpfr::mpreal> const &x,
               std::vector<mpfr::mpreal> &xsmall,
               std::function<mpfr::mpreal(mpfr::mpreal)> &f,
               std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log = false);

bool diffcorr(std::vector<mpfr::mpreal> &num, std::vector<mpfr::mpreal> &den,
              mpfr::mpreal &errDC,
              std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nb,
              std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &db,
              std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dbounds,
              std::function<mpfr::mpreal(mpfr::mpreal)> &nfix,
              std::function<mpfr::mpreal(mpfr::mpreal)> &dfix,
              std::vector<mpfr::mpreal> const &x,
              std::function<mpfr::mpreal(mpfr::mpreal)> &f,
              std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log = false);

bool diffcorr(std::vector<mpfr::mpreal> &num, std::vector<mpfr::mpreal> &den,
              mpfr::mpreal &errDC,
              std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nb,
              std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &db,
              std::vector<mpfr::mpreal> const &x,
              std::function<mpfr::mpreal(mpfr::mpreal)> &f,
              std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log = false);

bool diffcorr(std::vector<mpfr::mpreal> &num, std::vector<mpfr::mpreal> &den,
              mpfr::mpreal &errDC, std::pair<size_t, size_t> &type,
              std::vector<mpfr::mpreal> const &x,
              std::function<mpfr::mpreal(mpfr::mpreal)> &f,
              std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log = false);

#endif