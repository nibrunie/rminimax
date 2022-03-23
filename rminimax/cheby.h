#ifndef RMINIMAX_CHEBY
#define RMINIMAX_CHEBY

#include <mpreal.h>
#include <utility>
#include <vector>

enum ChebyshevKind { FIRST, SECOND };

void roots(std::vector<mpfr::mpreal> &r, std::vector<mpfr::mpreal> &a,
           std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
           ChebyshevKind kind = FIRST, bool balance = true);

void chebpts(std::vector<mpfr::mpreal> &r, std::size_t n);

void equipts(std::vector<mpfr::mpreal> &r, std::size_t n);

void logpts(std::vector<mpfr::mpreal> &r, std::size_t n);

void chebcoeffs(std::vector<mpfr::mpreal> &c, std::vector<mpfr::mpreal> &fv);

void diffcoeffs(std::vector<mpfr::mpreal> &dc, std::vector<mpfr::mpreal> &c,
                ChebyshevKind kind = SECOND);

void cos(std::vector<mpfr::mpreal> &out, std::vector<mpfr::mpreal> &in);

#endif