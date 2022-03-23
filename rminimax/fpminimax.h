#ifndef RMINIMAX_FPMINIMAX
#define RMINIMAX_FPMINIMAX

#include <functional>
#include <gmp.h>
#include <gmpxx.h>
#include <mpfr.h>
#include <mpreal.h>
#include <vector>

std::pair<mpz_class, mp_exp_t> mpfr_decomp(mpfr::mpreal const &val);

void fpminimax(
    std::vector<mpfr::mpreal> &fpnum, std::vector<mpfr::mpreal> &fpden,
    std::function<mpfr::mpreal(mpfr::mpreal)> const &r,
    std::function<mpfr::mpreal(mpfr::mpreal)> const &w,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> const &nbasis,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> const &dbasis,
    std::vector<mpfr::mpreal> const &num, std::vector<mpfr::mpreal> const &den,
    std::vector<mp_prec_t> const &nump, std::vector<mp_prec_t> const &denp,
    std::pair<mpfr::mpreal, mpfr::mpreal> const &dom, std::size_t idx,
    mp_prec_t prec);

void fpminimax_kernel(std::vector<mpfr::mpreal> &lll_coeffs,
                      std::vector<std::vector<mpfr::mpreal>> &svp_coeffs,
                      std::vector<mpfr::mpreal> const &x,
                      std::function<mpfr::mpreal(mpfr::mpreal)> &f,
                      std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &b,
                      mp_prec_t prec);

void factorize(
    mpfr::mpreal &scale,
    std::vector<std::pair<std::vector<mpfr::mpreal>, std::size_t>> &factors,
    std::vector<mpfr::mpreal> const &coeffs, mp_prec_t prec);

#endif