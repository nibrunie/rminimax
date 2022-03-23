#include "diffcorr.h"
#include "plotting.h"
#include <set>

extern "C" {
#include <qsopt_ex/QSopt_ex.h>
}

template <typename T> void remove_duplicate(std::vector<T> &v) {
  std::set<T> s(begin(v), end(v));
  v.assign(begin(s), end(s));
}

void mpreal_to_mpq(mpq_t &ratVal, mpfr::mpreal &mprealVal) {
  mpf_t mpfVal;
  mpf_init(mpfVal);
  mpfr_get_f(mpfVal, mprealVal.mpfr_ptr(), GMP_RNDN);
  mpq_set_f(ratVal, mpfVal);
  mpf_clear(mpfVal);
}

bool dcinit(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
            std::vector<mpfr::mpreal> &den,
            std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nB,
            std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dB,
            std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dBounds,
            std::function<mpfr::mpreal(mpfr::mpreal)> &nFix,
            std::function<mpfr::mpreal(mpfr::mpreal)> &dFix,
            std::vector<mpfr::mpreal> const &x,
            std::function<mpfr::mpreal(mpfr::mpreal)> &f,
            std::function<mpfr::mpreal(mpfr::mpreal)> &w) {
  int rval{0}, status{0}, ncols{(int)nB.size() + (int)dB.size() + 1},
      nrows{2 * (int)x.size()}, nrowsT{3 * (int)x.size()};
  bool success{true};

  int *cmatcnt = new int[ncols];
  int *cmatbeg = new int[ncols];
  int *cmatind = new int[ncols * nrowsT];
  char *sense = new char[nrowsT];

  mpq_t *cmatval = new mpq_t[ncols * nrowsT];
  mpq_t *obj = new mpq_t[ncols];
  mpq_t *rhs = new mpq_t[nrowsT];
  mpq_t *lower = new mpq_t[ncols];
  mpq_t *upper = new mpq_t[ncols];

  mpq_QSprob p;
  // construct the objective function
  for (int i{0}; i < ncols; ++i) {
    mpq_init(obj[i]);
    if (i < ncols - 1)
      mpq_set_ui(obj[i], 0u, 1u);
    else
      mpq_set_ui(obj[i], 1u, 1u);
  }

  mpfr::mpreal coeff;
  // construct the constraint matrix
  for (int i{(int)nB.size()}; i < ncols - 1; ++i) {
    for (int j{0}; j < nrows; j += 2) {
      mpq_init(cmatval[i * nrowsT + j]);
      mpq_init(cmatval[i * nrowsT + j + 1]);

      coeff = f(x[j / 2]) * w(x[j / 2]) * dB[i - nB.size()](x[j / 2]);
      mpreal_to_mpq(cmatval[i * nrowsT + j], coeff);
      coeff = -f(x[j / 2]) * w(x[j / 2]) * dB[i - nB.size()](x[j / 2]);
      mpreal_to_mpq(cmatval[i * nrowsT + j + 1], coeff);
    }
    for (int j{nrows}; j < nrowsT; ++j) {
      mpq_init(cmatval[i * nrowsT + j]);
      coeff = -dB[i - nB.size()](x[j - nrows]);
      mpreal_to_mpq(cmatval[i * nrowsT + j], coeff);
    }
  }

  for (int i{0}; i < nB.size(); ++i) {
    for (int j{0}; j < nrows; j += 2) {
      mpq_init(cmatval[i * nrowsT + j]);
      mpq_init(cmatval[i * nrowsT + j + 1]);

      coeff = -nB[i](x[j / 2]) * w(x[j / 2]);
      mpreal_to_mpq(cmatval[i * nrowsT + j], coeff);
      coeff = -coeff;
      mpreal_to_mpq(cmatval[i * nrowsT + j + 1], coeff);
    }
    for (int j{nrows}; j < nrowsT; ++j) {
      mpq_init(cmatval[i * nrowsT + j]);
      mpq_set_d(cmatval[i * nrowsT + j], 0.0);
    }
  }

  for (int i{0}; i < nrows; i += 2) {
    mpq_init(cmatval[(ncols - 1) * nrowsT + i]);
    mpq_init(cmatval[(ncols - 1) * nrowsT + i + 1]);

    coeff = -1;
    mpreal_to_mpq(cmatval[(ncols - 1) * nrowsT + i], coeff);
    mpq_set(cmatval[(ncols - 1) * nrowsT + i + 1],
            cmatval[(ncols - 1) * nrowsT + i]);

    // initialise the right hand side of the constraint matrix
    mpq_init(rhs[i]);
    mpq_init(rhs[i + 1]);
    coeff = 0;
    coeff = (-f(x[i / 2]) * dFix(x[i / 2]) + nFix(x[i / 2])) * w(x[i / 2]);
    mpreal_to_mpq(rhs[i], coeff);
    coeff = -coeff;
    mpreal_to_mpq(rhs[i + 1], coeff);
  }

  for (int i{nrows}; i < nrowsT; ++i) {
    mpq_init(cmatval[(ncols - 1) * nrowsT + i]);

    mpreal_to_mpq(cmatval[(ncols - 1) * nrowsT + i], coeff);
    mpq_set_d(cmatval[(ncols - 1) * nrowsT + i], 0.0);

    // initialise the right hand side of the constraint matrix
    mpq_init(rhs[i]);
    coeff = -1e-20;
    mpreal_to_mpq(rhs[i], coeff);
  }

  // construct the variable constraints
  for (int i{0}; i < ncols; ++i) {
    mpq_init(lower[i]);
    mpq_init(upper[i]);
    if (i < nB.size()) {
      mpq_set(lower[i], mpq_ILL_MINDOUBLE);
      mpq_set(upper[i], mpq_ILL_MAXDOUBLE);
    } else if (i < ncols - 1) {
      mpreal_to_mpq(lower[i], dBounds[i - nB.size()].first);
      mpreal_to_mpq(upper[i], dBounds[i - nB.size()].second);
    } else {
      mpq_set(lower[i], mpq_ILL_MINDOUBLE);
      mpq_set(upper[i], mpq_ILL_MAXDOUBLE);
    }
  }

  // set the sense of the inequalities: "<="
  for (int i{0}; i < nrowsT; ++i) {
    sense[i] = 'L';
  }

  for (int i{0}; i < ncols; ++i) {
    cmatcnt[i] = nrowsT;
    cmatbeg[i] = i * nrowsT;
    for (int j{0}; j < nrowsT; ++j)
      cmatind[i * nrowsT + j] = j;
  }

  p = mpq_QSload_prob("diffcorr", ncols, nrowsT, cmatcnt, cmatbeg, cmatind,
                      cmatval, QS_MIN, obj, rhs, sense, lower, upper, nullptr,
                      nullptr);

  if (p == nullptr) {
    std::cerr << "Unable to load the problem...\n";

    delete[] cmatcnt;
    delete[] cmatbeg;
    delete[] cmatind;
    for (int i{0}; i < ncols * nrowsT; ++i)
      mpq_clear(cmatval[i]);
    delete[] cmatval;

    return false;
  }

  rval = QSexact_solver(p, nullptr, nullptr, nullptr, DUAL_SIMPLEX, &status);

  if (rval)
    std::cerr << "QSexact_solver failed\n";
  if (status != QS_LP_OPTIMAL) {
    std::cerr << "Did not find optimal solution.\n";
    std::cerr << "Status code: " << status << std::endl;
  }

  mpq_t objval;
  mpq_init(objval);
  rval = mpq_QSget_objval(p, &objval);
  if (rval) {
    std::cerr << "Could not get objective value, error code: " << rval << "\n;";
  } else {
    delta = mpfr::mpreal(objval);
  }

  mpq_t *xs = new mpq_t[ncols];
  for (int i{0}; i < ncols; ++i)
    mpq_init(xs[i]);
  rval = mpq_QSget_x_array(p, xs);
  if (rval) {
    std::cerr << "Could not get variable vector, error code: " << rval << "\n";
    success = false;
  } else {
    num.clear();
    den.clear();
    for (int i{0}; i < nB.size(); ++i)
      num.emplace_back(xs[i]);
    for (int i{0}; i < dB.size(); ++i)
      den.push_back(xs[nB.size() + i]);
  }

  mpq_QSfree_prob(p);

  delete[] cmatcnt;
  delete[] cmatbeg;
  delete[] cmatind;
  for (int i{0}; i < ncols; ++i)
    mpq_clear(xs[i]);
  delete[] xs;
  for (int i{0}; i < ncols * nrowsT; ++i)
    mpq_clear(cmatval[i]);
  delete[] cmatval;
  mpq_clear(objval);
  return success;
}

// TODO: this is a general routine. Need to also code a version where
// at least one of the denominator coefficients has an active bounding
// constraint (in a guaranteed way)
bool dckernel(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
              std::vector<mpfr::mpreal> &den,
              std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nB,
              std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dB,
              std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dBounds,
              std::function<mpfr::mpreal(mpfr::mpreal)> &nFix,
              std::function<mpfr::mpreal(mpfr::mpreal)> &dFix,
              std::vector<mpfr::mpreal> const &x,
              std::function<mpfr::mpreal(mpfr::mpreal)> &qk,
              std::function<mpfr::mpreal(mpfr::mpreal)> &f,
              std::function<mpfr::mpreal(mpfr::mpreal)> &w,
              mpfr::mpreal const &deltak) {
  int rval{0}, status{0}, ncols{(int)nB.size() + (int)dB.size() + 1},
      nrows{2 * (int)x.size()};
  bool success{true};

  int *cmatcnt = new int[ncols];
  int *cmatbeg = new int[ncols];
  int *cmatind = new int[ncols * nrows];
  char *sense = new char[nrows];

  mpq_t *cmatval = new mpq_t[ncols * nrows];
  mpq_t *obj = new mpq_t[ncols];
  mpq_t *rhs = new mpq_t[nrows];
  mpq_t *lower = new mpq_t[ncols];
  mpq_t *upper = new mpq_t[ncols];

  mpq_QSprob p;
  // construct the objective function
  for (int i{0}; i < ncols; ++i) {
    mpq_init(obj[i]);
    if (i < ncols - 1)
      mpq_set_ui(obj[i], 0u, 1u);
    else
      mpq_set_ui(obj[i], 1u, 1u);
  }

  mpfr::mpreal coeff;
  // construct the constraint matrix
  for (int i{(int)nB.size()}; i < ncols - 1; ++i)
    for (int j{0}; j < nrows; j += 2) {
      mpq_init(cmatval[i * nrows + j]);
      mpq_init(cmatval[i * nrows + j + 1]);

      coeff =
          (f(x[j / 2]) * w(x[j / 2]) - deltak) * dB[i - nB.size()](x[j / 2]);
      mpreal_to_mpq(cmatval[i * nrows + j], coeff);
      coeff =
          (-f(x[j / 2]) * w(x[j / 2]) - deltak) * dB[i - nB.size()](x[j / 2]);
      mpreal_to_mpq(cmatval[i * nrows + j + 1], coeff);
    }

  for (int i{0}; i < nB.size(); ++i)
    for (int j{0}; j < nrows; j += 2) {
      mpq_init(cmatval[i * nrows + j]);
      mpq_init(cmatval[i * nrows + j + 1]);

      coeff = -nB[i](x[j / 2]) * w(x[j / 2]);
      mpreal_to_mpq(cmatval[i * nrows + j], coeff);
      coeff = -coeff;
      mpreal_to_mpq(cmatval[i * nrows + j + 1], coeff);
    }

  for (int i{0}; i < nrows; i += 2) {
    mpq_init(cmatval[(ncols - 1) * nrows + i]);
    mpq_init(cmatval[(ncols - 1) * nrows + i + 1]);

    coeff = -qk(x[i / 2]);
    mpreal_to_mpq(cmatval[(ncols - 1) * nrows + i], coeff);
    mpq_set(cmatval[(ncols - 1) * nrows + i + 1],
            cmatval[(ncols - 1) * nrows + i]);

    // initialise the right hand side of the constraint matrix
    mpq_init(rhs[i]);
    mpq_init(rhs[i + 1]);
    coeff = 0;
    coeff = (deltak - f(x[i / 2]) * w(x[i / 2])) * dFix(x[i / 2]) +
            nFix(x[i / 2]) * w(x[i / 2]);
    mpreal_to_mpq(rhs[i], coeff);
    coeff = 0;
    coeff = (deltak + f(x[i / 2]) * w(x[i / 2])) * dFix(x[i / 2]) -
            nFix(x[i / 2]) * w(x[i / 2]);
    mpreal_to_mpq(rhs[i + 1], coeff);
  }

  // construct the variable constraints
  for (int i{0}; i < ncols; ++i) {
    mpq_init(lower[i]);
    mpq_init(upper[i]);
    if (i < nB.size()) {
      mpq_set(lower[i], mpq_ILL_MINDOUBLE);
      mpq_set(upper[i], mpq_ILL_MAXDOUBLE);
    } else if (i < ncols - 1) {
      mpreal_to_mpq(lower[i], dBounds[i - nB.size()].first);
      mpreal_to_mpq(upper[i], dBounds[i - nB.size()].second);
    } else {
      mpq_set(lower[i], mpq_ILL_MINDOUBLE);
      mpq_set(upper[i], mpq_ILL_MAXDOUBLE);
    }
  }

  // set the sense of the inequalities: "<="
  for (int i{0}; i < nrows; ++i) {
    sense[i] = 'L';
  }

  for (int i{0}; i < ncols; ++i) {
    cmatcnt[i] = nrows;
    cmatbeg[i] = i * nrows;
    for (int j{0}; j < nrows; ++j)
      cmatind[i * nrows + j] = j;
  }

  p = mpq_QSload_prob("diffcorr", ncols, nrows, cmatcnt, cmatbeg, cmatind,
                      cmatval, QS_MIN, obj, rhs, sense, lower, upper, nullptr,
                      nullptr);

  if (p == nullptr) {
    std::cerr << "Unable to load the problem...\n";

    delete[] cmatcnt;
    delete[] cmatbeg;
    delete[] cmatind;
    for (int i{0}; i < ncols * nrows; ++i)
      mpq_clear(cmatval[i]);
    delete[] cmatval;

    return false;
  }

  rval = QSexact_solver(p, nullptr, nullptr, nullptr, DUAL_SIMPLEX, &status);

  if (rval)
    std::cerr << "QSexact_solver failed\n";
  if (status != QS_LP_OPTIMAL) {
    std::cerr << "Did not find optimal solution.\n";
    std::cerr << "Status code: " << status << std::endl;
  }

  mpq_t objval;
  mpq_init(objval);
  rval = mpq_QSget_objval(p, &objval);
  if (rval) {
    std::cerr << "Could not get objective value, error code: " << rval << "\n;";
  } else {
    delta = mpfr::mpreal(objval);
  }

  mpq_t *xs = new mpq_t[ncols];
  for (int i{0}; i < ncols; ++i)
    mpq_init(xs[i]);
  rval = mpq_QSget_x_array(p, xs);
  if (rval) {
    std::cerr << "Could not get variable vector, error code: " << rval << "\n";
    success = false;
  } else {
    num.clear();
    den.clear();
    for (int i{0}; i < nB.size(); ++i)
      num.emplace_back(xs[i]);
    for (int i{0}; i < dB.size(); ++i)
      den.push_back(xs[nB.size() + i]);
  }

  mpq_QSfree_prob(p);

  delete[] cmatcnt;
  delete[] cmatbeg;
  delete[] cmatind;
  for (int i{0}; i < ncols; ++i)
    mpq_clear(xs[i]);
  delete[] xs;
  for (int i{0}; i < ncols * nrows; ++i)
    mpq_clear(cmatval[i]);
  delete[] cmatval;
  mpq_clear(objval);
  return success;
}

bool dckernelV2(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
                std::vector<mpfr::mpreal> &den,
                std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nB,
                std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dB,
                std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dBounds,
                std::function<mpfr::mpreal(mpfr::mpreal)> &nFix,
                std::function<mpfr::mpreal(mpfr::mpreal)> &dFix,
                std::vector<mpfr::mpreal> const &x,
                std::function<mpfr::mpreal(mpfr::mpreal)> &qk,
                std::function<mpfr::mpreal(mpfr::mpreal)> &f,
                std::function<mpfr::mpreal(mpfr::mpreal)> &w,
                mpfr::mpreal const &deltak) {
  int rval{0}, status{0}, ncols{(int)nB.size() + (int)dB.size() + 1},
      nrows{2 * (int)x.size()}, nrowsT{3 * (int)x.size()};
  bool success{true};

  int *cmatcnt = new int[ncols];
  int *cmatbeg = new int[ncols];
  int *cmatind = new int[ncols * nrowsT];
  char *sense = new char[nrowsT];

  mpq_t *cmatval = new mpq_t[ncols * nrowsT];
  mpq_t *obj = new mpq_t[ncols];
  mpq_t *rhs = new mpq_t[nrowsT];
  mpq_t *lower = new mpq_t[ncols];
  mpq_t *upper = new mpq_t[ncols];

  mpq_QSprob p;
  // construct the objective function
  for (int i{0}; i < ncols; ++i) {
    mpq_init(obj[i]);
    if (i < ncols - 1)
      mpq_set_ui(obj[i], 0u, 1u);
    else
      mpq_set_ui(obj[i], 1u, 1u);
  }

  mpfr::mpreal coeff;
  // construct the constraint matrix
  for (int i{(int)nB.size()}; i < ncols - 1; ++i) {
    for (int j{0}; j < nrows; j += 2) {
      mpq_init(cmatval[i * nrowsT + j]);
      mpq_init(cmatval[i * nrowsT + j + 1]);

      coeff =
          (f(x[j / 2]) * w(x[j / 2]) - deltak) * dB[i - nB.size()](x[j / 2]);
      mpreal_to_mpq(cmatval[i * nrowsT + j], coeff);
      coeff =
          (-f(x[j / 2]) * w(x[j / 2]) - deltak) * dB[i - nB.size()](x[j / 2]);
      mpreal_to_mpq(cmatval[i * nrowsT + j + 1], coeff);
    }
    for (int j{nrows}; j < nrowsT; ++j) {
      mpq_init(cmatval[i * nrowsT + j]);
      coeff = -dB[i - nB.size()](x[j - nrows]);
      mpreal_to_mpq(cmatval[i * nrowsT + j], coeff);
    }
  }

  for (int i{0}; i < nB.size(); ++i) {
    for (int j{0}; j < nrows; j += 2) {
      mpq_init(cmatval[i * nrowsT + j]);
      mpq_init(cmatval[i * nrowsT + j + 1]);

      coeff = -nB[i](x[j / 2]) * w(x[j / 2]);
      mpreal_to_mpq(cmatval[i * nrowsT + j], coeff);
      coeff = -coeff;
      mpreal_to_mpq(cmatval[i * nrowsT + j + 1], coeff);
    }
    for (int j{nrows}; j < nrowsT; ++j) {
      mpq_init(cmatval[i * nrowsT + j]);
      mpq_set_d(cmatval[i * nrowsT + j], 0.0);
    }
  }

  for (int i{0}; i < nrows; i += 2) {
    mpq_init(cmatval[(ncols - 1) * nrowsT + i]);
    mpq_init(cmatval[(ncols - 1) * nrowsT + i + 1]);

    coeff = -qk(x[i / 2]);
    mpreal_to_mpq(cmatval[(ncols - 1) * nrowsT + i], coeff);
    mpq_set(cmatval[(ncols - 1) * nrowsT + i + 1],
            cmatval[(ncols - 1) * nrowsT + i]);

    // initialise the right hand side of the constraint matrix
    mpq_init(rhs[i]);
    mpq_init(rhs[i + 1]);
    coeff = 0;
    coeff = (deltak - f(x[i / 2]) * w(x[i / 2])) * dFix(x[i / 2]) +
            nFix(x[i / 2]) * w(x[i / 2]);
    mpreal_to_mpq(rhs[i], coeff);
    coeff = 0;
    coeff = (deltak + f(x[i / 2]) * w(x[i / 2])) * dFix(x[i / 2]) -
            nFix(x[i / 2]) * w(x[i / 2]);
    mpreal_to_mpq(rhs[i + 1], coeff);
  }

  for (int i{nrows}; i < nrowsT; ++i) {
    mpq_init(cmatval[(ncols - 1) * nrowsT + i]);

    mpreal_to_mpq(cmatval[(ncols - 1) * nrowsT + i], coeff);
    mpq_set_d(cmatval[(ncols - 1) * nrowsT + i], 0.0);

    // initialise the right hand side of the constraint matrix
    mpq_init(rhs[i]);
    coeff = -1e-20;
    mpreal_to_mpq(rhs[i], coeff);
  }

  // construct the variable constraints
  for (int i{0}; i < ncols; ++i) {
    mpq_init(lower[i]);
    mpq_init(upper[i]);
    if (i < nB.size()) {
      mpq_set(lower[i], mpq_ILL_MINDOUBLE);
      mpq_set(upper[i], mpq_ILL_MAXDOUBLE);
    } else if (i < ncols - 1) {
      mpreal_to_mpq(lower[i], dBounds[i - nB.size()].first);
      mpreal_to_mpq(upper[i], dBounds[i - nB.size()].second);
    } else {
      mpq_set(lower[i], mpq_ILL_MINDOUBLE);
      mpq_set(upper[i], mpq_ILL_MAXDOUBLE);
    }
  }

  // set the sense of the inequalities: "<="
  for (int i{0}; i < nrowsT; ++i) {
    sense[i] = 'L';
  }

  for (int i{0}; i < ncols; ++i) {
    cmatcnt[i] = nrowsT;
    cmatbeg[i] = i * nrowsT;
    for (int j{0}; j < nrowsT; ++j)
      cmatind[i * nrowsT + j] = j;
  }

  p = mpq_QSload_prob("diffcorr", ncols, nrowsT, cmatcnt, cmatbeg, cmatind,
                      cmatval, QS_MIN, obj, rhs, sense, lower, upper, nullptr,
                      nullptr);

  if (p == nullptr) {
    std::cerr << "Unable to load the problem...\n";

    delete[] cmatcnt;
    delete[] cmatbeg;
    delete[] cmatind;
    for (int i{0}; i < ncols * nrowsT; ++i)
      mpq_clear(cmatval[i]);
    delete[] cmatval;

    return false;
  }

  rval = QSexact_solver(p, nullptr, nullptr, nullptr, DUAL_SIMPLEX, &status);

  if (rval)
    std::cerr << "QSexact_solver failed\n";
  if (status != QS_LP_OPTIMAL) {
    std::cerr << "Did not find optimal solution.\n";
    std::cerr << "Status code: " << status << std::endl;
  }

  mpq_t objval;
  mpq_init(objval);
  rval = mpq_QSget_objval(p, &objval);
  if (rval) {
    std::cerr << "Could not get objective value, error code: " << rval << "\n;";
  } else {
    delta = mpfr::mpreal(objval);
  }

  mpq_t *xs = new mpq_t[ncols];
  for (int i{0}; i < ncols; ++i)
    mpq_init(xs[i]);
  rval = mpq_QSget_x_array(p, xs);
  if (rval) {
    std::cerr << "Could not get variable vector, error code: " << rval << "\n";
    success = false;
  } else {
    num.clear();
    den.clear();
    for (int i{0}; i < nB.size(); ++i)
      num.emplace_back(xs[i]);
    for (int i{0}; i < dB.size(); ++i)
      den.push_back(xs[nB.size() + i]);
  }

  mpq_QSfree_prob(p);

  delete[] cmatcnt;
  delete[] cmatbeg;
  delete[] cmatind;
  for (int i{0}; i < ncols; ++i)
    mpq_clear(xs[i]);
  delete[] xs;
  for (int i{0}; i < ncols * nrowsT; ++i)
    mpq_clear(cmatval[i]);
  delete[] cmatval;
  mpq_clear(objval);
  return success;
}

bool diffcorr(std::vector<mpfr::mpreal> &num, std::vector<mpfr::mpreal> &den,
              mpfr::mpreal &errDC,
              std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nB,
              std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dB,
              std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dBounds,
              std::function<mpfr::mpreal(mpfr::mpreal)> &nFix,
              std::function<mpfr::mpreal(mpfr::mpreal)> &dFix,
              std::vector<mpfr::mpreal> const &x,
              std::function<mpfr::mpreal(mpfr::mpreal)> &f,
              std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log) {
  if (log)
    std::cout << "================DIFFCORR START=================\n";
  bool success{true};

  mpfr::mpreal del;
  dcinit(del, num, den, nB, dB, dBounds, nFix, dFix, x, f, w);

  std::function<mpfr::mpreal(mpfr::mpreal)> qk;
  std::function<mpfr::mpreal(mpfr::mpreal)> pk;

  qk = [den, dB, dFix](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = dFix(var);
    for (std::size_t i{0u}; i < den.size(); ++i)
      res += den[i] * dB[i](var);
    return res;
  };
  pk = [num, nB, nFix](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = nFix(var);
    for (std::size_t i{0u}; i < num.size(); ++i)
      res += num[i] * nB[i](var);
    return res;
  };

  mpfr::mpreal dk = mpfr::abs(w(x[0]) * (f(x[0]) - pk(x[0]) / qk(x[0])));
  for (std::size_t i{0}; i < x.size(); ++i)
    if (mpfr::abs(w(x[i]) * (f(x[i]) - pk(x[i]) / qk(x[i]))) > dk)
      dk = mpfr::abs(w(x[i]) * (f(x[i]) - pk(x[i]) / qk(x[i])));

  mpfr::mpreal minDen = qk(x[0]);
  for (auto &it : x)
    if (qk(it) < minDen)
      minDen = qk(it);

  mpfr::mpreal delta;
  success =
      dckernelV2(delta, num, den, nB, dB, dBounds, nFix, dFix, x, qk, f, w, dk);

  qk = [den, dB, dFix](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = dFix(var);
    for (std::size_t i{0u}; i < den.size(); ++i)
      res += den[i] * dB[i](var);
    return res;
  };
  pk = [num, nB, nFix](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = nFix(var);
    for (std::size_t i{0u}; i < num.size(); ++i)
      res += num[i] * nB[i](var);
    return res;
  };

  mpfr::mpreal ndk = mpfr::abs(w(x[0]) * (f(x[0]) - pk(x[0]) / qk(x[0])));
  for (std::size_t i{0}; i < x.size(); ++i)
    if (mpfr::abs(w(x[i]) * (f(x[i]) - pk(x[i]) / qk(x[i]))) > ndk)
      ndk = mpfr::abs(w(x[i]) * (f(x[i]) - pk(x[i]) / qk(x[i])));

  while (mpfr::abs(dk - ndk) / mpfr::abs(ndk) > 1e-18 && success) {
    if (log) {
      std::cout << "Convergence order = "
                << mpfr::abs(dk - ndk) / mpfr::abs(ndk) << std::endl;
      std::cout << "Differential correction delta = " << ndk << std::endl;
    }
    mpfr::mpreal minDen = qk(x[0]);
    for (auto &it : x)
      if (qk(it) < minDen)
        minDen = qk(it);
    dk = ndk;
    num.clear();
    den.clear();
    success = dckernelV2(delta, num, den, nB, dB, dBounds, nFix, dFix, x, qk, f,
                         w, dk);

    qk = [den, dB, dFix](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = dFix(var);
      for (std::size_t i{0u}; i < den.size(); ++i)
        res += den[i] * dB[i](var);
      return res;
    };
    pk = [num, nB, nFix](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = nFix(var);
      for (std::size_t i{0u}; i < num.size(); ++i)
        res += num[i] * nB[i](var);
      return res;
    };

    ndk = mpfr::abs(w(x[0]) * (f(x[0]) - pk(x[0]) / qk(x[0])));
    for (std::size_t i{0}; i < x.size(); ++i)
      if (mpfr::abs(w(x[i]) * (f(x[i]) - pk(x[i]) / qk(x[i]))) > ndk)
        ndk = mpfr::abs(w(x[i]) * (f(x[i]) - pk(x[i]) / qk(x[i])));
  }
  errDC = ndk;
  if (log)
    std::cout << "Differential correction delta = " << ndk << std::endl;
  minDen = qk(x[0]);
  for (auto &it : x)
    if (qk(it) < minDen)
      minDen = qk(it);
  if (log)
    std::cout << "================DIFFCORR FINISH=================\n";
  return success;
}

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
               std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log) {
  if (log)
    std::cout << "================ADIFFCORR START=================\n";
  bool success{true};

  mpfr::mpreal errDC0 = 0.0;
  mpfr::mpreal errDC1, errDCk;
  std::vector<mpfr::mpreal> x1 = xsmall;
  std::vector<std::size_t> ax1;

  if (x1.empty()) {
    // need to construct the small active subset
    std::size_t step =
        (std::size_t)ceil((double)x.size() / (nb.size() + db.size() + 1));
    for (std::size_t i{0u}; i < x.size(); i += step) {
      x1.push_back(x[i]);
      ax1.push_back(i);
    }
    x1.push_back(x[x.size() - 1]);
    ax1.push_back(x.size() - 1);
    remove_duplicate(x1);
    remove_duplicate(ax1);

    success = diffcorr(num, den, errDC1, nb, db, dbounds, nfix, dfix, x1, f, w);
  } else {
    // need to sort the active subset and determine the
    // correct indices in the large subset
    for (std::size_t i{0u}; i < x1.size(); ++i) {
      bool found{false};
      for (std::size_t j{0u}; j < x.size() && !found; ++j) {
        if (x1[i] == x[j]) {
          found = true;
          ax1.push_back(j);
        }
      }
    }
    errDC1 = errDC;
  }

  if (log)
    std::cout << "ADC delta = " << errDC1 << std::endl;

  std::function<mpfr::mpreal(mpfr::mpreal)> qk;
  std::function<mpfr::mpreal(mpfr::mpreal)> pk;
  std::function<mpfr::mpreal(mpfr::mpreal)> ek;

  qk = [den, db, dfix](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = dfix(var);
    for (std::size_t i{0u}; i < den.size(); ++i)
      res += den[i] * db[i](var);
    return res;
  };
  pk = [num, nb, nfix](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = nfix(var);
    for (std::size_t i{0u}; i < num.size(); ++i)
      res += num[i] * nb[i](var);
    return res;
  };
  // TODO: handle the case where a denominator value is
  // below a certain threshold value
  ek = [pk, qk, f, w](mpfr::mpreal var) -> mpfr::mpreal {
    return w(var) * (f(var) - pk(var) / qk(var));
  };
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> pts;
  for (auto &it : x1) {
    pts.push_back(std::make_pair(it, ek(it)));
  }

  std::vector<std::size_t> tr1;
  for (std::size_t i{0u}; i < x1.size(); ++i) {
    if (mpfr::abs(ek(x1[i])) >= errDC1 * (1.0 - 1e-8)) {
      tr1.push_back(ax1[i]);
    }
  }
  if (log)
    std::cout << "trk size = " << tr1.size() << std::endl;

  std::vector<std::size_t> a1;
  // perform the "uphill" search of the current point
  for (std::size_t i{0u}; i < tr1.size(); ++i) {
    int dir;
    if (tr1[i] > 0u && tr1[i] < x.size() - 1u)
      dir = (mpfr::abs(ek(x[tr1[i] - 1u])) > mpfr::abs(ek(x[tr1[i] + 1u]))) ? -1
                                                                            : 1;
    else if (tr1[i] == 0u)
      dir = 1;
    else // trk[i] == x.size() - 1u
      dir = -1;
    std::size_t toAdd{tr1[i]};
    while ((toAdd + dir < x.size()) && (toAdd + dir >= 0u) &&
           (mpfr::abs(ek(x[toAdd + dir])) > mpfr::abs(ek(x[toAdd])))) {
      toAdd = toAdd + dir;
    }
    if (mpfr::abs(ek(x[toAdd])) > errDC * (1.0 + 1e-8))
      a1.push_back(toAdd);
  }
  if (log)
    std::cout << "ak size = " << a1.size() << std::endl;

  std::vector<std::size_t> axk{tr1};
  for (auto &it : a1)
    axk.push_back(it);
  remove_duplicate(axk);
  std::sort(begin(axk), end(axk));
  std::vector<mpfr::mpreal> xk;
  for (auto &it : axk)
    xk.push_back(x[it]);

  if (log)
    std::cout << "xk size = " << xk.size() << std::endl;

  // TODO: test the fact that the current denominator takes
  // a positive value for all of the elements in x1 (if not,
  // we need to reinitialize the coefficient values)

  bool stop{false};
  std::vector<std::size_t> sk;
  std::vector<std::size_t> ak;
  std::vector<std::size_t> trk;

  while (!stop && success) {
    if (log) {
      std::cout << "newxk size = " << xk.size() << std::endl;
      success = diffcorr(num, den, errDCk, nb, db, dbounds, nfix, dfix, xk, f,
                         w, log);
      std::cout << "ADC delta = " << errDCk << std::endl;
    } else {
      success = diffcorr(num, den, errDCk, nb, db, dbounds, nfix, dfix, xk, f,
                         w, log);
    }

    stop = true;
    qk = [den, db, dfix](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = dfix(var);
      for (std::size_t i{0u}; i < den.size(); ++i)
        res += den[i] * db[i](var);
      return res;
    };
    pk = [num, nb, nfix](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = nfix(var);
      for (std::size_t i{0u}; i < num.size(); ++i)
        res += num[i] * nb[i](var);
      return res;
    };
    // TODO: handle the case where a denominator value is
    // below a certain threshold value
    ek = [pk, qk, f, w](mpfr::mpreal var) -> mpfr::mpreal {
      return w(var) * (f(var) - pk(var) / qk(var));
    };

    pts.clear();
    for (auto &it : xk) {
      pts.push_back(std::make_pair(it, ek(it)));
    }

    trk.clear();
    for (std::size_t i{0u}; i < xk.size(); ++i) {
      if (mpfr::abs(ek(xk[i])) >= errDCk * (1.0 - 1e-8))
        trk.push_back(axk[i]);
    }
    if (log)
      std::cout << "trk size = " << trk.size() << std::endl;

    ak.clear();
    // perform the "uphill" search of the current point
    for (std::size_t i{0u}; i < trk.size(); ++i) {
      int dir;
      if (trk[i] > 0u && trk[i] < x.size() - 1u)
        dir = (mpfr::abs(ek(x[trk[i] - 1u])) > mpfr::abs(ek(x[trk[i] + 1u])))
                  ? -1
                  : 1;
      else if (trk[i] == 0u)
        dir = 1;
      else // trk[i] == x.size() - 1u
        dir = -1;
      std::size_t toAdd{trk[i]};
      while ((toAdd + dir < x.size()) && (toAdd + dir >= 0u) &&
             (mpfr::abs(ek(x[toAdd + dir])) > mpfr::abs(ek(x[toAdd])))) {
        toAdd = toAdd + dir;
      }
      if (mpfr::abs(ek(x[toAdd])) > errDCk * (1.0 + 1e-8)) {
        stop = false;
        ak.push_back(toAdd);
      }
    }
    if (log)
      std::cout << "ak = " << ak.size() << std::endl;

    if (errDCk <= mpfr::max(errDC1, errDC0) * (1.0 + 1e-4))
      sk = axk;
    else
      sk = trk;

    std::vector<std::size_t> auxk{sk};
    for (auto &it : ak)
      auxk.push_back(it);
    if (errDCk < errDC1)
      for (auto &it : ax1)
        auxk.push_back(it);
    remove_duplicate(auxk);
    std::sort(begin(auxk), end(auxk));

    errDC0 = errDC1;
    errDC1 = errDCk;
    ax1 = axk;
    axk = auxk;
    xk.clear();
    for (auto &it : axk)
      xk.push_back(x[it]);
    if (log)
      std::cout << "xk size = " << xk.size() << std::endl;
  }
  errDC = errDCk;
  xsmall = xk;
  if (log)
    std::cout << "================ADIFFCORR FINISH=================\n";
  return success;
}

bool diffcorr(std::vector<mpfr::mpreal> &num, std::vector<mpfr::mpreal> &den,
              mpfr::mpreal &errDC,
              std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nB,
              std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dB,
              std::vector<mpfr::mpreal> const &x,
              std::function<mpfr::mpreal(mpfr::mpreal)> &f,
              std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log) {
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> dBounds(dB.size());
  for (size_t i{0u}; i < dBounds.size(); ++i) {
    dBounds[i].first = -1;
    dBounds[i].second = 1;
  }
  std::function<mpfr::mpreal(mpfr::mpreal)> nFix =
      [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0.0); };
  std::function<mpfr::mpreal(mpfr::mpreal)> dFix =
      [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0.0); };
  return diffcorr(num, den, errDC, nB, dB, dBounds, nFix, dFix, x, f, w, log);
}

bool diffcorr(std::vector<mpfr::mpreal> &num, std::vector<mpfr::mpreal> &den,
              mpfr::mpreal &errDC, std::pair<size_t, size_t> &type,
              std::vector<mpfr::mpreal> const &x,
              std::function<mpfr::mpreal(mpfr::mpreal)> &f,
              std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log) {
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> dBounds(type.second + 1u);
  for (size_t i{0u}; i < dBounds.size(); ++i) {
    dBounds[i].first = -1;
    dBounds[i].second = 1;
  }

  std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> nB(type.first + 1u);
  std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> dB(type.second + 1u);
  for (size_t i{0u}; i <= type.first; ++i)
    nB[i] = [i](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::pow(x, i); };
  for (size_t i{0u}; i <= type.second; ++i)
    dB[i] = [i](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::pow(x, i); };
  std::function<mpfr::mpreal(mpfr::mpreal)> nFix =
      [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0.0); };
  std::function<mpfr::mpreal(mpfr::mpreal)> dFix =
      [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0.0); };
  return diffcorr(num, den, errDC, nB, dB, dBounds, nFix, dFix, x, f, w, log);
}