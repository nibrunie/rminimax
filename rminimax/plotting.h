#ifndef RMINIMAX_PLOTTING
#define RMINIMAX_PLOTTING

#include <functional>
#include <gmp.h>
#include <mpreal.h>
#include <string>
#include <vector>

void plotFunc(std::string &filename,
              std::function<mpfr::mpreal(mpfr::mpreal)> &f,
              mpfr::mpreal const &a, mpfr::mpreal const &b,
              mp_prec_t prec = 165ul);

void plotFuncDouble(std::string &filename,
                    std::function<mpfr::mpreal(mpfr::mpreal)> &f,
                    mpfr::mpreal const &a, mpfr::mpreal const &b,
                    mp_prec_t prec = 165ul);

void plotFunc(std::string &filename,
              std::function<mpfr::mpreal(mpfr::mpreal)> &f,
              std::vector<mpfr::mpreal> const &x, mp_prec_t prec = 165ul);

void plotFuncs(std::string &filename,
               std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &fs,
               mpfr::mpreal const &a, mpfr::mpreal const &b,
               mp_prec_t prec = 165ul);

void plotVals(std::string &filename,
              std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &p,
              mp_prec_t prec = 165ul);

void plotFuncEtVals(std::string &filename,
                    std::function<mpfr::mpreal(mpfr::mpreal)> &f,
                    std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &p,
                    mpfr::mpreal const &a, mpfr::mpreal const &b,
                    mp_prec_t prec = 165ul);

#endif
