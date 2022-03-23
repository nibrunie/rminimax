#include "cheby.h"
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

using mpfr::mpreal;
using std::size_t;

typedef Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> MatrixXm;
typedef Eigen::Matrix<std::complex<mpfr::mpreal>, Eigen::Dynamic, 1> VectorXcm;

void balance(MatrixXm &A) {
  size_t n = A.rows();

  mpfr::mpreal rnorm;
  mpfr::mpreal cnorm;
  bool converged = false;
  mpfr::mpreal one = mpfr::mpreal(1.0);

  mpfr::mpreal g, f, s;
  while (!converged) {
    converged = true;
    for (size_t i{0u}; i < n; ++i) {
      rnorm = cnorm = 0.0;
      for (size_t j{0u}; j < n; ++j) {
        if (i == j)
          continue;
        cnorm += mpfr::fabs(A(j, i));
        rnorm += mpfr::fabs(A(i, j));
      }
      if ((cnorm == 0) || (rnorm == 0))
        continue;

      g = rnorm >> 1u;
      f = 1.0;
      s = cnorm + rnorm;

      while (mpfr::isfinite(cnorm) && cnorm < g) {
        f <<= 1u;
        cnorm <<= 2u;
      }

      g = rnorm << 1u;

      while (mpfr::isfinite(cnorm) && cnorm > g) {
        f >>= 1u;
        cnorm >>= 2u;
      }

      if ((rnorm + cnorm) < s * f * 0.95) {
        converged = false;
        g = one / f;
        A.row(i) *= g;
        A.col(i) *= f;
      }
    }
  }
}

MatrixXm colleague(std::vector<mpfr::mpreal> const &a, ChebyshevKind kind,
                   bool bal) {
  std::vector<mpfr::mpreal> c = a;
  size_t n = a.size() - 1;
  MatrixXm C(n, n);

  for (size_t i{0u}; i < n; ++i)
    for (size_t j{0u}; j < n; ++j)
      C(i, j) = 0;

  mpfr::mpreal denom = -1;
  denom /= c[n];
  denom >>= 1;
  for (size_t i{0u}; i < a.size() - 1; ++i)
    c[i] *= denom;
  c[n - 2] += 0.5;

  for (size_t i{0u}; i < n - 1; ++i)
    C(i, i + 1) = C(i + 1, i) = 0.5;
  switch (kind) {
  case FIRST:
    C(n - 2, n - 1) = 1;
    break;
  default:
    C(n - 2, n - 1) = 0.5;
    break;
  }

  for (size_t i{0u}; i < n; ++i)
    C(i, 0) = c[n - i - 1];

  if (bal)
    balance(C);

  return C;
}

void roots(std::vector<mpfr::mpreal> &r, std::vector<mpfr::mpreal> &a,
           std::pair<mpfr::mpreal, mpfr::mpreal> const &dom, ChebyshevKind kind,
           bool balance) {
  MatrixXm C = colleague(a, kind, balance);
  Eigen::EigenSolver<MatrixXm> es(C);
  VectorXcm eigs = es.eigenvalues();

  mpfr::mpreal threshold = 10;
  threshold = mpfr::pow(threshold, -20, MPFR_RNDN);
  for (size_t i{0u}; i < eigs.size(); ++i)
    if (mpfr::abs(eigs(i).imag()) < threshold)
      if (dom.first <= eigs(i).real() && dom.second >= eigs(i).real())
        r.push_back(eigs(i).real());

  std::sort(begin(r), end(r));
}

void equipts(std::vector<mpfr::mpreal> &r, std::size_t n) {
  mpfr::mpreal pi = mpfr::const_pi();
  r.resize(n);
  for (size_t i{0u}; i < n; ++i) {
    r[n - 1 - i] = pi * i;
    r[n - 1 - i] /= (n - 1);
  }
}

void chebpts(std::vector<mpfr::mpreal> &r, std::size_t n) {
  equipts(r, n);
  for (size_t i{0u}; i < n; ++i)
    r[i] = mpfr::cos(r[i]);
}

void logpts(std::vector<mpfr::mpreal> &r, std::size_t n) {
  r.clear();
  if (n % 2 != 0)
    r.push_back(0);
  for (std::size_t i{0u}; i < n / 2; ++i) {
    r.push_back(pow(10.0, -4.0 * i / n));
    r.push_back(-pow(10.0, -4.0 * i / n));
  }
  std::sort(begin(r), end(r));
}

void clenshaw(mpfr::mpreal &r, std::vector<mpfr::mpreal> &a, mpfr::mpreal &x,
              ChebyshevKind kind) {
  mpfr::mpreal bn1, bn2, bn;

  int n = (int)a.size() - 1;
  bn2 = 0;
  bn1 = a[n];
  for (int k{n - 1}; k >= 1; --k) {
    bn = x * 2;
    bn = bn * bn1 - bn2 + a[k];
    bn2 = bn1;
    bn1 = bn;
  }

  switch (kind) {
  case FIRST:
    r = x * bn1 - bn2 + a[0];
    break;
  default:
    r = (x << 1) * bn1 - bn2 + a[0];
    break;
  }
}

void chebcoeffs(std::vector<mpfr::mpreal> &c, std::vector<mpfr::mpreal> &fv) {
  // TODO: look at Eigen's FFT to see if that
  // can be used instead here
  size_t n = fv.size();
  std::vector<mpfr::mpreal> v(n);
  c.resize(n);
  equipts(v, n);

  mpfr::mpreal buffer;

  fv[0] >>= 1;
  fv[n - 1] >>= 1;

  for (int i{0}; i < n; ++i) {
    buffer = mpfr::cos(v[i]);
    clenshaw(c[i], fv, buffer, FIRST);
    if (i == 0 || i == (n - 1))
      c[i] /= (n - 1);
    else {
      c[i] <<= 1;
      c[i] /= (n - 1);
    }
  }

  fv[0] <<= 1;
  fv[n - 1] <<= 1;
}

void diffcoeffs(std::vector<mpfr::mpreal> &dc, std::vector<mpfr::mpreal> &c,
                ChebyshevKind kind) {
  dc.resize(c.size() - 1);
  switch (kind) {
  case FIRST: {
    int n = (int)c.size() - 2;
    dc[n] = c[n + 1] * (2 * (n + 1));
    dc[n - 1] = c[n] * (2 * n);
    for (int i{n - 2}; i >= 0; --i) {
      dc[i] = 2 * (i + 1);
      dc[i] = fma(dc[i], c[i + 1], dc[i + 2]);
    }
    dc[0] >>= 1;
  }; break;
  default: {
    int n = c.size() - 1;
    for (int i{n}; i > 0; --i)
      dc[i - 1] = c[i] * i;
  }; break;
  }
}