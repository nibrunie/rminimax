#ifndef RMINIMAX_BARY
#define RMINIMAX_BARY

#include <mpreal.h>
#include <utility>
#include <vector>

void baryweights(std::vector<mpfr::mpreal> &w, std::vector<mpfr::mpreal> &x,
                 std::vector<mpfr::mpreal> &fx,
                 std::pair<std::size_t, std::size_t> &deg,
                 mp_prec_t prec = 165u);

void baryeval(mpfr::mpreal &result, const mpfr::mpreal &in,
              std::vector<mpfr::mpreal> &x, std::vector<mpfr::mpreal> &fx,
              std::vector<mpfr::mpreal> &w, mp_prec_t prec = 165ul);

#endif