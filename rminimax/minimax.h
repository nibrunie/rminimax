#ifndef RMINIMAX_MINIMAX
#define RMINIMAX_MINIMAX

#include <cstdio>
#include <cstdlib>
#include <functional>
#include <gmp.h>
#include <iostream>
#include <mpreal.h>
#include <utility>
#include <vector>

enum AlgorithmType { REMEZ_FIRST, REMEZ_SECOND };

void chgvar(std::vector<mpfr::mpreal> &out, std::vector<mpfr::mpreal> const &in,
            std::pair<mpfr::mpreal, mpfr::mpreal> const &type);

void infnorm(std::pair<mpfr::mpreal, mpfr::mpreal> &norm,
             std::function<mpfr::mpreal(mpfr::mpreal)> &f,
             std::pair<mpfr::mpreal, mpfr::mpreal> const &dom);

std::pair<bool, std::vector<mpfr::mpreal>>
minimax(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
        std::vector<mpfr::mpreal> &den,
        std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
        std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nbasis,
        std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dbasis,
        std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dbounds,
        std::function<mpfr::mpreal(mpfr::mpreal)> &nfix,
        std::function<mpfr::mpreal(mpfr::mpreal)> &dfix,
        std::vector<mpfr::mpreal> &x,
        std::function<mpfr::mpreal(mpfr::mpreal)> &f,
        std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log = false,
        mp_prec_t prec = 165ul);

std::pair<bool, std::vector<mpfr::mpreal>>
minimax(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
        std::vector<mpfr::mpreal> &den,
        std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
        std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nbasis,
        std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dbasis,
        std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dbounds,
        std::vector<mpfr::mpreal> &x,
        std::function<mpfr::mpreal(mpfr::mpreal)> &f,
        std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log = false,
        mp_prec_t prec = 165ul);

std::pair<bool, std::vector<mpfr::mpreal>>
minimax(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
        std::vector<mpfr::mpreal> &den,
        std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
        std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nbasis,
        std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dbasis,
        std::function<mpfr::mpreal(mpfr::mpreal)> &f,
        std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log = false,
        mp_prec_t prec = 165ul);

std::pair<bool, std::vector<mpfr::mpreal>>
minimax(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
        std::vector<mpfr::mpreal> &den,
        std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
        std::function<mpfr::mpreal(mpfr::mpreal)> &f,
        std::function<mpfr::mpreal(mpfr::mpreal)> &w,
        std::pair<size_t, size_t> &type, AlgorithmType atype, bool log = false,
        mp_prec_t prec = 165ul);

std::pair<bool, std::vector<mpfr::mpreal>>
minimax(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
        std::vector<mpfr::mpreal> &den,
        std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
        std::function<mpfr::mpreal(mpfr::mpreal)> &f,
        std::function<mpfr::mpreal(mpfr::mpreal)> &w,
        std::pair<size_t, size_t> &type, bool log = false,
        mp_prec_t prec = 165ul);

bool remez_2(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
             std::vector<mpfr::mpreal> &den,
             std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
             std::vector<mpfr::mpreal> &x,
             std::function<mpfr::mpreal(mpfr::mpreal)> &f,
             std::function<mpfr::mpreal(mpfr::mpreal)> &w,
             std::pair<size_t, size_t> &type, size_t Nmax,
             mpfr::mpreal const &convThreshold, bool log, mp_prec_t prec);

#endif