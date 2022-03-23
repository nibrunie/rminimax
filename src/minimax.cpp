#include "minimax.h"
#include "bary.h"
#include "cheby.h"
#include "diffcorr.h"
#include "plotting.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <set>

using mpfr::mpreal;
using std::size_t;

typedef Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> MatrixXm;
typedef Eigen::Matrix<std::complex<mpfr::mpreal>, Eigen::Dynamic, 1> VectorXcm;

static const std::size_t maxit{100u};

template <typename T> void remove_duplicate(std::vector<T> &v) {
  std::set<T> s(begin(v), end(v));
  v.assign(begin(s), end(s));
}

void chgvar(std::vector<mpfr::mpreal> &out, std::vector<mpfr::mpreal> const &in,
            std::pair<mpfr::mpreal, mpfr::mpreal> const &type) {
  out = in;
  for (size_t i{0u}; i < in.size(); ++i)
    out[i] = fma((type.second - type.first) / 2, in[i],
                 (type.second + type.first) / 2);
}

void split(std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &doms,
           std::pair<mpfr::mpreal, mpfr::mpreal> const &dom, size_t n) {
  std::vector<mpfr::mpreal> s;
  for (size_t i{0u}; i <= n; ++i)
    s.emplace_back(dom.first + (dom.second - dom.first) * mpfr::mpreal(i) / n);
  for (size_t i{0u}; i < n; ++i)
    doms.emplace_back(std::make_pair(s[i], s[i + 1u]));
}

void split(std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &doms,
           std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
           std::vector<mpfr::mpreal> const &x) {
  auto s = x;
  std::sort(begin(s), end(s),
            [](mpfr::mpreal x, mpfr::mpreal y) -> bool { return x < y; });
  remove_duplicate(s);

  doms.clear();
  if (s[0] != dom.first)
    doms.emplace_back(std::make_pair(dom.first, s[0]));
  for (size_t i{0u}; i < s.size() - 1u; ++i)
    doms.emplace_back(std::make_pair(s[i], s[i + 1]));
  if (s[s.size() - 1u] != dom.second)
    doms.emplace_back(std::make_pair(s[s.size() - 1], dom.second));
}

void infnorm(std::pair<mpfr::mpreal, mpfr::mpreal> &norm,
             std::function<mpfr::mpreal(mpfr::mpreal)> &f,
             std::pair<mpfr::mpreal, mpfr::mpreal> const &dom) {
  mpfr::mpreal ia = -1;
  mpfr::mpreal ib = 1;
  size_t mdeg = 8u;
  size_t N = 512;
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> doms;
  split(doms, dom, N);
  norm.first = dom.first;
  norm.second = mpfr::abs(f(dom.first));

  std::vector<mpfr::mpreal> x;
  chebpts(x, mdeg + 1u);

  for (size_t i{0u}; i < doms.size(); ++i) {
    std::vector<mpfr::mpreal> nx;
    chgvar(nx, x, doms[i]);
    std::vector<mpfr::mpreal> fx(mdeg + 1u);
    for (size_t j{0u}; j < fx.size(); ++j)
      fx[j] = f(nx[j]);

    std::vector<mpfr::mpreal> c(mdeg + 1u);
    chebcoeffs(c, fx);
    std::vector<mpfr::mpreal> dc(mdeg);
    diffcoeffs(dc, c, ChebyshevKind::FIRST);
    std::vector<mpfr::mpreal> rs;
    roots(rs, dc, std::make_pair(ia, ib));
    chgvar(rs, rs, doms[i]);

    mpfr::mpreal candMax;
    for (const auto &r : rs) {
      candMax = mpfr::abs(f(r));
      if (candMax > norm.second) {
        norm.first = r;
        norm.second = candMax;
      }
    }
  }
}

void infnorm(std::pair<mpfr::mpreal, mpfr::mpreal> &norm,
             std::function<mpfr::mpreal(mpfr::mpreal)> &f,
             std::vector<mpfr::mpreal> &sp,
             std::pair<mpfr::mpreal, mpfr::mpreal> const &dom) {
  mpfr::mpreal ia = -1;
  mpfr::mpreal ib = 1;
  size_t mdeg = 16u;
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> doms;
  split(doms, dom, sp);
  norm.first = dom.first;
  norm.second = mpfr::abs(f(dom.first));

  std::vector<mpfr::mpreal> x;
  chebpts(x, mdeg + 1u);

  for (size_t i{0u}; i < doms.size(); ++i) {
    std::vector<mpfr::mpreal> nx;
    chgvar(nx, x, doms[i]);
    std::vector<mpfr::mpreal> fx(mdeg + 1u);
    for (size_t j{0u}; j < fx.size(); ++j)
      fx[j] = f(nx[j]);

    std::vector<mpfr::mpreal> c(mdeg + 1u);
    chebcoeffs(c, fx);
    std::vector<mpfr::mpreal> dc(mdeg);
    diffcoeffs(dc, c, ChebyshevKind::FIRST);
    std::vector<mpfr::mpreal> rs;
    roots(rs, dc, std::make_pair(ia, ib));
    chgvar(rs, rs, doms[i]);

    mpfr::mpreal candMax;
    for (const auto &r : rs) {
      candMax = mpfr::abs(f(r));
      if (candMax > norm.second) {
        norm.first = r;
        norm.second = candMax;
      }
    }
  }
}

void infnorm(std::pair<mpfr::mpreal, mpfr::mpreal> &norm,
             std::vector<mpfr::mpreal> &local, mpfr::mpreal const &errThr,
             std::function<mpfr::mpreal(mpfr::mpreal)> &f,
             std::vector<mpfr::mpreal> &sp,
             std::pair<mpfr::mpreal, mpfr::mpreal> const &dom) {
  mpfr::mpreal ia = -1;
  mpfr::mpreal ib = 1;
  size_t mdeg = 8u;
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> doms;
  split(doms, dom, sp);
  norm.first = dom.first;
  norm.second = mpfr::abs(f(dom.first));

  std::vector<mpfr::mpreal> x;
  chebpts(x, mdeg + 1u);

  for (size_t i{0u}; i < doms.size(); ++i) {
    std::vector<mpfr::mpreal> nx;
    chgvar(nx, x, doms[i]);
    std::vector<mpfr::mpreal> fx(mdeg + 1u);
    for (size_t j{0u}; j < fx.size(); ++j)
      fx[j] = f(nx[j]);

    std::vector<mpfr::mpreal> c(mdeg + 1u);
    chebcoeffs(c, fx);
    std::vector<mpfr::mpreal> dc(mdeg);
    diffcoeffs(dc, c, ChebyshevKind::FIRST);
    std::vector<mpfr::mpreal> rs;
    roots(rs, dc, std::make_pair(ia, ib));
    chgvar(rs, rs, doms[i]);

    mpfr::mpreal candMax;
    for (const auto &r : rs) {
      candMax = mpfr::abs(f(r));
      if (candMax > errThr)
        local.push_back(r);
      if (candMax > norm.second) {
        norm.first = r;
        norm.second = candMax;
      }
    }
  }
  remove_duplicate(local);
}

bool remez_1(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
             std::vector<mpfr::mpreal> &den,
             std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
             std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nbasis,
             std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dbasis,
             std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dbounds,
             std::function<mpfr::mpreal(mpfr::mpreal)> &nfix,
             std::function<mpfr::mpreal(mpfr::mpreal)> &dfix,
             std::vector<mpfr::mpreal> &x,
             std::function<mpfr::mpreal(mpfr::mpreal)> &f,
             std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log,
             mp_prec_t prec) {
  mp_prec_t prev = mpfr::mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  bool success{true};
  std::vector<mpfr::mpreal> xsmall;
  std::vector<mpfr::mpreal> local;

  mpfr::mpreal errDC;
  if (log)
    std::cout << "==============ITERATION 1 ==============\n\n";

  success = adiffcorr(num, den, errDC, nbasis, dbasis, dbounds, nfix, dfix, x,
                      xsmall, f, w, log);

  std::pair<mpfr::mpreal, mpfr::mpreal> cnorm;
  std::function<mpfr::mpreal(mpfr::mpreal)> qk;
  std::function<mpfr::mpreal(mpfr::mpreal)> pk;

  qk = [den, dbasis, dfix](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = dfix(var);
    for (size_t i{0u}; i < dbasis.size(); ++i)
      res += den[i] * dbasis[i](var);
    return res;
  };
  pk = [num, nbasis, nfix](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = nfix(var);
    for (size_t i{0u}; i < nbasis.size(); ++i)
      res += num[i] * nbasis[i](var);
    return res;
  };

  std::function<mpfr::mpreal(mpfr::mpreal)> err;
  err = [pk, qk, f, w](mpfr::mpreal var) -> mpfr::mpreal {
    return w(var) * (f(var) - pk(var) / qk(var));
  };

  infnorm(cnorm, local, errDC, err, xsmall, dom);

  if (log) {
    std::cout << "\nLocation\tMaximum error\n";
    std::cout << cnorm.first << "\t" << cnorm.second << std::endl;
  }

  mpfr::mpreal preverr = cnorm.second * 2;
  size_t it{2u};
  if (log) {
    std::cout << "Error norm = " << cnorm.second << std::endl;
    std::cout << "DC    norm = " << errDC << std::endl << std::endl;
  }
  while (mpfr::abs(cnorm.second - errDC) / cnorm.second > 1e-4 && success &&
         it < maxit) {
    if (log) {
      std::cout << "==============ITERATION " << it << " ==============\n";
    }
    preverr = cnorm.second;
    for (auto &it : local)
      x.emplace_back(it);
    remove_duplicate(x);
    std::sort(begin(x), end(x));
    if (log)
      std::cout << "Discrete set size = " << x.size() << std::endl;

    success = adiffcorr(num, den, errDC, nbasis, dbasis, dbounds, nfix, dfix, x,
                        xsmall, f, w, log);

    qk = [den, dbasis, dfix](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = dfix(var);
      for (size_t i{0u}; i < dbasis.size(); ++i)
        res += den[i] * dbasis[i](var);
      return res;
    };
    pk = [num, nbasis, nfix](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = nfix(var);
      for (size_t i{0u}; i < nbasis.size(); ++i)
        res += num[i] * nbasis[i](var);
      return res;
    };
    err = [pk, qk, f, w](mpfr::mpreal var) -> mpfr::mpreal {
      return w(var) * (f(var) - pk(var) / qk(var));
    };

    local.clear();
    infnorm(cnorm, local, errDC, err, xsmall, dom);

    ++it;
    if (log) {
      std::cout << "\nLocation\tMaximum error\n";
      std::cout << cnorm.first << "\t" << cnorm.second << std::endl;
      std::cout << "Error norm = " << cnorm.second << std::endl;
      std::cout << "DC    norm = " << errDC << std::endl << std::endl;
    }
  }

  delta = cnorm.second;
  x = xsmall;
  if (success) {
    std::string name = "error_R1";
    plotFunc(name, err, x, prec);
  }

  if (it == maxit) {
    if (log) {
      std::cout << "Warning! Max number of iterations = " << maxit
                << " reached in minimax" << std::endl
                << " approximation routine. Results might be inaccurate."
                << std::endl;
    }
  }

  mpreal::set_default_prec(prev);
  return success;
}

bool remez_1V2(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
               std::vector<mpfr::mpreal> &den,
               std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
               std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nbasis,
               std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dbasis,
               std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dbounds,
               std::function<mpfr::mpreal(mpfr::mpreal)> &nfix,
               std::function<mpfr::mpreal(mpfr::mpreal)> &dfix,
               std::vector<mpfr::mpreal> &x,
               std::function<mpfr::mpreal(mpfr::mpreal)> &f,
               std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log,
               mp_prec_t prec) {
  mp_prec_t prev = mpfr::mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  bool success{true};

  mpfr::mpreal errDC;
  if (log)
    std::cout << "==============ITERATION 1 ==============\n\n";

  success = diffcorr(num, den, errDC, nbasis, dbasis, dbounds, nfix, dfix, x, f,
                     w, log);

  std::pair<mpfr::mpreal, mpfr::mpreal> cnorm;
  std::function<mpfr::mpreal(mpfr::mpreal)> qk;
  std::function<mpfr::mpreal(mpfr::mpreal)> pk;

  qk = [den, dbasis, dfix](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = dfix(var);
    for (size_t i{0u}; i < dbasis.size(); ++i)
      res += den[i] * dbasis[i](var);
    return res;
  };
  pk = [num, nbasis, nfix](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = nfix(var);
    for (size_t i{0u}; i < nbasis.size(); ++i)
      res += num[i] * nbasis[i](var);
    return res;
  };

  std::function<mpfr::mpreal(mpfr::mpreal)> err;
  err = [pk, qk, f, w](mpfr::mpreal var) -> mpfr::mpreal {
    return w(var) * (f(var) - pk(var) / qk(var));
  };

  infnorm(cnorm, err, x, dom);

  if (log) {
    std::cout << "\nLocation\tMaximum error\n";
    std::cout << cnorm.first << "\t" << cnorm.second << std::endl;
  }

  mpfr::mpreal preverr = cnorm.second * 2;
  size_t it{2u};
  if (log) {
    std::cout << "Error norm = " << cnorm.second << std::endl;
    std::cout << "DC    norm = " << errDC << std::endl << std::endl;
  }
  while (mpfr::abs(cnorm.second - errDC) / cnorm.second > 1e-4 && success &&
         it < maxit) {
    if (log) {
      std::cout << "==============ITERATION " << it << " ==============\n";
    }
    preverr = cnorm.second;

    x.emplace_back(cnorm.first);
    remove_duplicate(x);
    std::sort(begin(x), end(x));

    success =
        diffcorr(num, den, errDC, nbasis, dbasis, dbounds, nfix, dfix, x, f, w);

    qk = [den, dbasis, dfix](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = dfix(var);
      for (size_t i{0u}; i < dbasis.size(); ++i)
        res += den[i] * dbasis[i](var);
      return res;
    };
    pk = [num, nbasis, nfix](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = nfix(var);
      for (size_t i{0u}; i < nbasis.size(); ++i)
        res += num[i] * nbasis[i](var);
      return res;
    };
    err = [pk, qk, f, w](mpfr::mpreal var) -> mpfr::mpreal {
      return w(var) * (f(var) - pk(var) / qk(var));
    };

    infnorm(cnorm, err, x, dom);

    ++it;
    if (log) {
      std::cout << "\nLocation\tMaximum error\n";
      std::cout << cnorm.first << "\t" << cnorm.second << std::endl;
      std::cout << "Error norm = " << cnorm.second << std::endl;
      std::cout << "DC    norm = " << errDC << std::endl << std::endl;
    }
  }

  delta = cnorm.second;

  if (success) {
    std::string name = "error_R1";
    plotFunc(name, err, x, prec);
  }

  if (it == maxit) {
    if (log) {
      std::cout << "Warning! Max number of iterations = " << maxit
                << " reached in minimax" << std::endl
                << " approximation routine. Results might be inaccurate."
                << std::endl;
    }
  }

  mpreal::set_default_prec(prev);
  return success;
}

std::pair<bool, std::vector<mpfr::mpreal>> minimax(
    mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
    std::vector<mpfr::mpreal> &den,
    std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nbasis,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dbasis,
    std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dbounds,
    std::function<mpfr::mpreal(mpfr::mpreal)> &nfix,
    std::function<mpfr::mpreal(mpfr::mpreal)> &dfix,
    std::vector<mpfr::mpreal> &x, std::function<mpfr::mpreal(mpfr::mpreal)> &f,
    std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log, mp_prec_t prec) {
  bool success = remez_1(delta, num, den, dom, nbasis, dbasis, dbounds, nfix,
                         dfix, x, f, w, log, prec);
  return std::make_pair(success, x);
}

std::pair<bool, std::vector<mpfr::mpreal>> minimax(
    mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
    std::vector<mpfr::mpreal> &den,
    std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nbasis,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dbasis,
    std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> &dbounds,
    std::vector<mpfr::mpreal> &x, std::function<mpfr::mpreal(mpfr::mpreal)> &f,
    std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log, mp_prec_t prec) {
  mp_prec_t prev = mpfr::mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  std::function<mpfr::mpreal(mpfr::mpreal)> nfix =
      [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0.0); };

  std::function<mpfr::mpreal(mpfr::mpreal)> dfix =
      [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0.0); };

  mpreal::set_default_prec(prev);

  bool success = remez_1(delta, num, den, dom, nbasis, dbasis, dbounds, nfix,
                         dfix, x, f, w, log, prec);
  return std::make_pair(success, x);
}

std::pair<bool, std::vector<mpfr::mpreal>>
minimax(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
        std::vector<mpfr::mpreal> &den,
        std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
        std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &nbasis,
        std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &dbasis,
        std::function<mpfr::mpreal(mpfr::mpreal)> &f,
        std::function<mpfr::mpreal(mpfr::mpreal)> &w, bool log,
        mp_prec_t prec) {
  mp_prec_t prev = mpfr::mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  std::vector<mpfr::mpreal> x;
  x.resize(100u * (nbasis.size() + dbasis.size()) + 1u);
  for (size_t i{0}; i < x.size(); ++i)
    x[i] = dom.first + (dom.second - dom.first) * i / (x.size() - 1);
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> dbounds;
  dbounds.resize(dbasis.size());
  for (size_t i{0u}; i < dbasis.size(); ++i) {
    dbounds[i].first = -1;
    dbounds[i].second = 1;
  }

  mpreal::set_default_prec(prev);

  return minimax(delta, num, den, dom, nbasis, dbasis, dbounds, x, f, w, log,
                 prec);
}

bool trialapprox(mpfr::mpreal &h, std::vector<mpfr::mpreal> &num,
                 std::vector<mpfr::mpreal> &den, std::vector<mpfr::mpreal> &fx,
                 std::vector<mpfr::mpreal> &wx, std::vector<mpfr::mpreal> &x,
                 std::pair<size_t, size_t> &type, mp_prec_t prec) {
  mp_prec_t prev = mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  size_t N = type.first + type.second;
  MatrixXm Vm(N + 2u, N + 2u);

  for (size_t i{0u}; i < N + 2u; ++i)
    Vm(i, 0) = 1.0;

  for (size_t i{0u}; i < N + 2u; ++i)
    for (size_t j{1u}; j < N + 2u; ++j)
      Vm(i, j) = Vm(i, j - 1u) * x[i];

  Eigen::HouseholderQR<MatrixXm> qr(Vm);
  MatrixXm C = qr.householderQ();

  MatrixXm ZL = C.block(0u, type.first + 1u, N + 2u, type.second + 1u);
  ZL.transposeInPlace();
  for (size_t i{0u}; i < ZL.rows(); ++i)
    for (size_t j{0u}; j < ZL.cols(); ++j)
      ZL(i, j) *= fx[j];
  ZL = ZL * C.block(0, 0, N + 2u, type.second + 1u);

  std::vector<mpfr::mpreal> wsigma(fx.size());
  for (size_t i{0u}; i < wsigma.size(); i += 2u)
    wsigma[i] = mpfr::mpreal(1.0) / wx[i];
  for (size_t i{1u}; i < wsigma.size(); i += 2u)
    wsigma[i] = mpfr::mpreal(-1.0) / wx[i];

  MatrixXm ZR = C.block(0u, type.first + 1u, N + 2u, type.second + 1u);
  ZR.transposeInPlace();
  for (size_t i{0u}; i < ZR.rows(); ++i)
    for (size_t j{0u}; j < ZR.cols(); ++j)
      ZR(i, j) *= wsigma[j];
  ZR = ZR * C.block(0, 0, N + 2u, type.second + 1u);

  Eigen::GeneralizedEigenSolver<MatrixXm> ges(ZL, ZR, true);
  MatrixXm d = ges.eigenvalues().real();
  MatrixXm v = ges.eigenvectors();

  MatrixXm qAll = C.block(0u, 0u, N + 2u, type.second + 1u) * v;

  bool valid{false};
  size_t vidx{0u};
  for (size_t i{0u}; i < qAll.cols(); ++i) {
    int sgnSum{0};
    for (size_t j{0u}; j < qAll.rows(); ++j)
      sgnSum += mpfr::sgn(qAll(j, i));
    if ((size_t)abs(sgnSum) == (size_t)qAll.rows()) {
      if (!valid) {
        valid = true;
        vidx = i;
      } else {
        std::cerr << "WARNING Remez: More than one candidate!\n";
        if (mpfr::abs(d(i)) < mpfr::abs(d(vidx)))
          vidx = i;
      }
    }
  }

  if (!valid) {
    std::cerr << "ERROR Remez: No candidate for current iteration!\n";

    mpreal::set_default_prec(prev);
    return false;
  }

  h = d(vidx);
  // retrieve the denominator coefficients in the
  // monomial basis
  MatrixXm bDen = Vm.block(0u, 0u, N + 2, type.second + 1u)
                      .fullPivHouseholderQr()
                      .solve(qAll.col(vidx));
  // ensure that the denominator takes positive values
  if (qAll(0, vidx) < 0)
    bDen = -bDen;
  den.clear();
  for (size_t i{0u}; i < bDen.rows(); ++i)
    den.push_back(bDen(i, 0));
  // construct the corresponding numerator
  num.clear();
  MatrixXm b = Vm.block(0u, 0u, N + 2u, type.second + 1u) * bDen;
  for (size_t i{0u}; i < b.rows(); i += 2)
    b(i, 0) *= (wx[i] * fx[i] - h);
  for (size_t i{1u}; i < b.rows(); i += 2)
    b(i, 0) *= (wx[i] * fx[i] + h);
  MatrixXm A = Vm.block(0u, 0u, N + 2u, type.first + 1u);
  for (size_t i{0u}; i < A.rows(); ++i)
    for (size_t j{0u}; j < A.cols(); ++j)
      A(i, j) *= wx[i];
  // solve the linear system for determining the numerator
  MatrixXm solNum = A.fullPivHouseholderQr().solve(b);
  for (size_t i{0u}; i < solNum.rows(); ++i)
    num.push_back(solNum(i, 0));
  // renormalization of the coefficients
  mpfr::mpreal renorm = den[0];
  for (size_t i{1u}; i < den.size(); ++i)
    if (mpfr::abs(den[i]) > mpfr::abs(renorm))
      renorm = den[i];

  for (size_t i{0u}; i < num.size(); ++i)
    num[i] /= renorm;
  for (size_t i{0u}; i < den.size(); ++i)
    den[i] /= renorm;

  mpreal::set_default_prec(prev);
  return true;
}

bool find_extrema(mpfr::mpreal &minimaxerr, std::vector<mpfr::mpreal> &num,
                  std::vector<mpfr::mpreal> &den,
                  mpfr::mpreal &convergenceOrder,
                  std::vector<mpfr::mpreal> &newRef,
                  std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
                  std::function<mpfr::mpreal(mpfr::mpreal)> &f,
                  std::function<mpfr::mpreal(mpfr::mpreal)> &w,
                  std::vector<mpfr::mpreal> &oldRef,
                  std::pair<size_t, size_t> &type, size_t Nmax, bool log,
                  mp_prec_t prec) {
  mp_prec_t prev = mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> intervals;
  mpfr::mpreal a = -1;
  mpfr::mpreal b = 1;
  split(intervals, dom, oldRef);

  mpfr::mpreal h;
  std::vector<mpfr::mpreal> oldFx(oldRef.size());
  std::vector<mpfr::mpreal> oldWx(oldRef.size());
  for (size_t i{0u}; i < oldFx.size(); ++i) {
    oldFx[i] = f(oldRef[i]);
    oldWx[i] = w(oldRef[i]);
  }

  if (!trialapprox(h, num, den, oldFx, oldWx, oldRef, type, prec)) {
    mpreal::set_default_prec(prec);
    return false;
  }
  std::vector<mpfr::mpreal> interpF(oldRef.size());
  for (size_t i{0u}; i < interpF.size(); i += 2u)
    interpF[i] = oldFx[i] - h / oldWx[i];
  for (size_t i{1u}; i < interpF.size(); i += 2u)
    interpF[i] = oldFx[i] + h / oldWx[i];

  std::vector<mpfr::mpreal> wx{oldRef};
  std::vector<mpfr::mpreal> wfx{interpF};
  size_t wSize = oldFx.size() - 1u;

  mpfr::mpreal buffer = wx[wSize];
  wx[wSize] = wx[wSize - 1u];
  wx[wSize - 1u] = buffer;
  wx.resize(wSize);

  buffer = wfx[wSize];
  wfx[wSize] = wfx[wSize - 1u];
  wfx[wSize - 1u] = buffer;
  wfx.resize(wSize);

  std::vector<mpfr::mpreal> xz, xnz, fxnz;

  for (size_t i{0u}; i < wfx.size(); ++i)
    if (wfx[i] == 0)
      xz.push_back(wx[i]);
    else {
      xnz.push_back(wx[i]);
      fxnz.push_back(wfx[i]);
    }

  std::vector<mpfr::mpreal> bw;
  std::function<mpfr::mpreal(mpfr::mpreal)> func;
  std::function<mpfr::mpreal(mpfr::mpreal)> err;
  size_t mu{xz.size()};
  std::pair<size_t, size_t> ntype{type};
  if (type.first - mu >= type.second) {
    ntype.first -= mu;
    bw.resize(xnz.size());
    baryweights(bw, xnz, fxnz, ntype, prec);
    func = [&](mpfr::mpreal val) -> mpfr::mpreal {
      mpfr::mpreal res;
      baryeval(res, val, xnz, fxnz, bw, prec);
      for (size_t j{0u}; j < xz.size(); ++j)
        res *= (res - xz[j]);
      return res;
    };

    err = [&](mpfr::mpreal val) -> mpfr::mpreal {
      mpfr::mpreal res;
      baryeval(res, val, xnz, fxnz, bw, prec);
      for (size_t j{0u}; j < xz.size(); ++j)
        res *= (res - xz[j]);
      return w(val) * (f(val) - res);
    };
  } else {
    ntype.first = type.second;
    ntype.second = type.first - mu;
    for (size_t i{0u}; i < fxnz.size(); ++i)
      fxnz[i] = mpfr::mpreal(1) / fxnz[i];
    bw.resize(xnz.size());
    baryweights(bw, xnz, fxnz, ntype, prec);
    func = [&](mpfr::mpreal val) -> mpfr::mpreal {
      mpfr::mpreal res;
      baryeval(res, val, xnz, fxnz, bw, prec);
      res = mpfr::mpreal(1) / res;
      for (size_t j{0u}; j < xz.size(); ++j)
        res *= (res - xz[j]);
      return res;
    };
    err = [&](mpfr::mpreal val) -> mpfr::mpreal {
      mpfr::mpreal res;
      baryeval(res, val, xnz, fxnz, bw, prec);
      res = mpfr::mpreal(1) / res;
      for (size_t j{0u}; j < xz.size(); ++j)
        res *= (res - xz[j]);
      return w(val) * (f(val) - res);
    };
  }

  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> pextrema;
  std::vector<mpfr::mpreal> chebx;
  chebpts(chebx, Nmax + 1u);

  mpfr::mpreal extrErrValL, extrErrValR, extrErrVal;
  std::vector<std::vector<mpfr::mpreal>> pExs(intervals.size());

  for (size_t i{0u}; i < intervals.size(); ++i) {
    std::vector<mpfr::mpreal> sichebx;
    chgvar(sichebx, chebx, intervals[i]);

    std::vector<mpfr::mpreal> fx(Nmax + 1u);
    for (size_t j{0u}; j < fx.size(); ++j)
      fx[j] = err(sichebx[j]);

    std::vector<mpfr::mpreal> c(Nmax + 1u);
    chebcoeffs(c, fx);
    std::vector<mpfr::mpreal> dc(Nmax);
    diffcoeffs(dc, c, ChebyshevKind::SECOND);
    std::vector<mpfr::mpreal> rs;
    roots(rs, dc, std::make_pair(a, b), ChebyshevKind::SECOND);
    chgvar(rs, rs, intervals[i]);

    for (size_t j{0u}; j < rs.size(); ++j)
      pExs[i].push_back(rs[j]);
    pExs[i].push_back(intervals[i].first);
    pExs[i].push_back(intervals[i].second);
  }

  for (size_t i{0u}; i < pExs.size(); ++i)
    for (size_t j{0u}; j < pExs[i].size(); ++j)
      pextrema.push_back(std::make_pair(pExs[i][j], err(pExs[i][j])));

  std::sort(begin(pextrema), end(pextrema),
            [](const std::pair<mpfr::mpreal, mpfr::mpreal> &lhs,
               const std::pair<mpfr::mpreal, mpfr::mpreal> &rhs) {
              return lhs.first < rhs.first;
            });

  newRef.clear();
  size_t eitx{0u};
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> alternatingExtrema;
  mpfr::mpreal minErr = INT_MAX;
  mpfr::mpreal maxErr = INT_MIN;
  mpfr::mpreal absErr;

  while (eitx < pextrema.size()) {
    std::pair<mpfr::mpreal, mpfr::mpreal> maxErrPoint = pextrema[eitx];
    while (eitx < pextrema.size() - 1u &&
           mpfr::sgn(maxErrPoint.second) *
                   mpfr::sgn(pextrema[eitx + 1u].second) >
               0) {
      ++eitx;
      if (mpfr::abs(maxErrPoint.second) < mpfr::abs(pextrema[eitx].second))
        maxErrPoint = pextrema[eitx];
    }
    alternatingExtrema.push_back(maxErrPoint);
    ++eitx;
  }

  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> buffExtrema;
  if (log) {
    std::cout << "Alternating exterma: " << alternatingExtrema.size() << "/"
              << oldRef.size() << std::endl;
  }

  if (alternatingExtrema.size() < oldRef.size()) {
    std::cerr << "ERROR Remez: Not enough alternating extrema!\n"
              << "POSSIBLE CAUSE: Nmax too small\n"
              << "Size: " << alternatingExtrema.size() << std::endl;
    convergenceOrder = 2.0;
    mpreal::set_default_prec(prev);
    return false;
  } else if (alternatingExtrema.size() > oldRef.size()) {
    size_t superflu = alternatingExtrema.size() - oldRef.size();
    if (superflu % 2 != 0) {
      mpfr::mpreal abs1 = mpfr::abs(alternatingExtrema[0u].second);
      mpfr::mpreal abs2 =
          mpfr::abs(alternatingExtrema[alternatingExtrema.size() - 1u].second);
      size_t sidx{0u};
      if (abs1 < abs2)
        sidx = 1u;
      for (size_t i{sidx}; i < alternatingExtrema.size() + sidx - 1u; ++i)
        buffExtrema.push_back(alternatingExtrema[i]);
      alternatingExtrema = buffExtrema;
      buffExtrema.clear();
    }

    while (alternatingExtrema.size() > oldRef.size()) {
      size_t remidx{0u};
      mpfr::mpreal minValRem =
          mpfr::max(mpfr::abs(alternatingExtrema[0].second),
                    mpfr::abs(alternatingExtrema[1].second));
      mpfr::mpreal remBuff;
      for (size_t i{1u}; i < alternatingExtrema.size() - 1u; ++i) {
        remBuff = mpfr::max(mpfr::abs(alternatingExtrema[i].second),
                            mpfr::abs(alternatingExtrema[i + 1u].second));
        if (remBuff < minValRem) {
          minValRem = remBuff;
          remidx = i;
        }
      }
      for (size_t i{0u}; i < remidx; ++i)
        buffExtrema.push_back(alternatingExtrema[i]);
      for (size_t i{remidx + 2u}; i < alternatingExtrema.size(); ++i)
        buffExtrema.push_back(alternatingExtrema[i]);
      alternatingExtrema = buffExtrema;
      buffExtrema.clear();
    }
  }

  if (alternatingExtrema.size() < oldRef.size()) {
    std::cerr << "ERROR Remez: Building new reference set failed\n";
    mpreal::set_default_prec(prev);
    return false;
  }

  newRef.clear();
  for (auto &it : alternatingExtrema) {
    newRef.push_back(it.first);
    absErr = mpfr::abs(it.second);
    minErr = mpfr::min(minErr, absErr);
    maxErr = mpfr::max(maxErr, absErr);
  }

  if (log) {
    std::cout << "Min error = " << minErr << std::endl;
    std::cout << "Max error = " << maxErr << std::endl;
  }
  convergenceOrder = (maxErr - minErr) / maxErr;
  if (log) {
    std::cout << "Convergence order = " << convergenceOrder << std::endl
              << std::endl;
  }
  minimaxerr = minErr;

  std::string name = "error_R2";
  plotFunc(name, err, dom.first, dom.second);

  mpreal::set_default_prec(prev);
  return true;
}

bool remez_2(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
             std::vector<mpfr::mpreal> &den,
             std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
             std::vector<mpfr::mpreal> &x,
             std::function<mpfr::mpreal(mpfr::mpreal)> &f,
             std::function<mpfr::mpreal(mpfr::mpreal)> &w,
             std::pair<size_t, size_t> &type, size_t Nmax,
             mpfr::mpreal const &convThreshold, bool log, mp_prec_t prec) {
  mp_prec_t prev = mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  mpfr::mpreal convOrder = 1;
  std::vector<mpfr::mpreal> nx;
  size_t itx{1u};
  do {
    if (log) {
      std::cout << "==============ITERATION " << itx << " ==============\n\n";
    }
    bool succ = find_extrema(delta, num, den, convOrder, nx, dom, f, w, x, type,
                             Nmax, log, prec);
    if (!succ) {
      mpreal::set_default_prec(prev);
      return false;
    }
    x = nx;
    std::cout << std::endl;
    ++itx;
  } while (convOrder > convThreshold && itx < maxit);
  x = nx;

  if (itx == maxit) {
    if (log) {
      std::cout << "Warning! Max number of iterations = " << maxit
                << " reached in minimax" << std::endl
                << " approximation routine. Results might be inaccurate."
                << std::endl;
    }
  }

  mpreal::set_default_prec(prev);
  return true;
}

std::pair<bool, std::vector<mpfr::mpreal>>
minimax(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
        std::vector<mpfr::mpreal> &den,
        std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
        std::function<mpfr::mpreal(mpfr::mpreal)> &f,
        std::function<mpfr::mpreal(mpfr::mpreal)> &w,
        std::pair<size_t, size_t> &type, AlgorithmType atype, bool log,
        mp_prec_t prec) {
  bool success{true};
  std::vector<mpfr::mpreal> x;
  switch (atype) {
  case AlgorithmType::REMEZ_FIRST: {
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> nBasis(type.first +
                                                                  1u);
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> dBasis(type.second +
                                                                  1u);
    nBasis[0] = [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(1); };
    dBasis[0] = [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(1); };
    for (size_t i{1u}; i <= type.first; ++i)
      nBasis[i] = [i](mpfr::mpreal val) -> mpfr::mpreal {
        return mpfr::pow(val, i);
      };
    for (size_t i{1u}; i <= type.second; ++i)
      dBasis[i] = [i](mpfr::mpreal val) -> mpfr::mpreal {
        return mpfr::pow(val, i);
      };
    std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> dBounds(dBasis.size());
    for (size_t i{0u}; i < dBounds.size(); ++i) {
      dBounds[i].first = -1;
      dBounds[i].second = 1;
    }
    std::function<mpfr::mpreal(mpfr::mpreal)> nFix =
        [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0); };

    std::function<mpfr::mpreal(mpfr::mpreal)> dFix =
        [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0); };

    x.resize(1000u * (type.first + type.second + 2u) + 1u);
    chebpts(x, x.size());
    chgvar(x, x, dom);
    success = remez_1(delta, num, den, dom, nBasis, dBasis, dBounds, nFix, dFix,
                      x, f, w, log, prec);
  }

  break;

  default: {
    x.resize(type.first + type.second + 2u);
    chebpts(x, x.size());
    chgvar(x, x, dom);

    success = remez_2(delta, num, den, dom, x, f, w, type, 8u,
                      mpfr::mpreal(1e-3), log, prec * 4);
    if (!success) {
      if (log)
        std::cout << "Attempting DC-based algorithm...\n";
      std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> nBasis(type.first +
                                                                    1u);
      std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> dBasis(
          type.second + 1u);
      nBasis[0] = [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(1); };
      dBasis[0] = [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(1); };
      for (size_t i{1u}; i <= type.first; ++i)
        nBasis[i] = [i](mpfr::mpreal val) -> mpfr::mpreal {
          return mpfr::pow(val, i);
        };
      for (size_t i{1u}; i <= type.second; ++i)
        dBasis[i] = [i](mpfr::mpreal val) -> mpfr::mpreal {
          return mpfr::pow(val, i);
        };
      std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> dBounds(dBasis.size());
      for (size_t i{0u}; i < dBounds.size(); ++i) {
        dBounds[i].first = -1;
        dBounds[i].second = 1;
      }
      std::function<mpfr::mpreal(mpfr::mpreal)> nFix =
          [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0); };

      std::function<mpfr::mpreal(mpfr::mpreal)> dFix =
          [](mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0); };

      x.resize(1000u * (type.first + type.second + 2u) + 1u);
      chebpts(x, x.size());
      chgvar(x, x, dom);

      success = remez_1(delta, num, den, dom, nBasis, dBasis, dBounds, nFix,
                        dFix, x, f, w, log, prec);
    }

  } break;
  }
  return std::make_pair(success, x);
}

std::pair<bool, std::vector<mpfr::mpreal>>
minimax(mpfr::mpreal &delta, std::vector<mpfr::mpreal> &num,
        std::vector<mpfr::mpreal> &den,
        std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
        std::function<mpfr::mpreal(mpfr::mpreal)> &f,
        std::function<mpfr::mpreal(mpfr::mpreal)> &w,
        std::pair<size_t, size_t> &type, bool log, mp_prec_t prec) {
  std::vector<mpfr::mpreal> x(type.first + type.second + 2u);
  chebpts(x, x.size());
  chgvar(x, x, dom);
  if (log) {
    std::cout << "Attempting Remez exchange algorithm...\n\n";
    if (!remez_2(delta, num, den, dom, x, f, w, type, 32u, mpfr::mpreal(1e-3),
                 log, prec)) {
      std::cout << "The Remez exchange algorithm implementation failed.\n\n"
                << "Attempting the general minimax rational approximation "
                   "code...\n\n";
      return minimax(delta, num, den, dom, f, w, type,
                     AlgorithmType::REMEZ_FIRST, log, prec);
    }
  } else {
    if (!remez_2(delta, num, den, dom, x, f, w, type, 32u, mpfr::mpreal(1e-3),
                 log, prec))
      return minimax(delta, num, den, dom, f, w, type,
                     AlgorithmType::REMEZ_FIRST, log, prec);
  }
  return std::make_pair(true, x);
}