#include "fpminimax.h"
#include "cheby.h"
#include "minimax.h"
#include <acb.h>
#include <arb.h>
#include <arb_fmpz_poly.h>
#include <cstdlib>
#include <flint/arith.h>
#include <fplll.h>
#include <fstream>
#include <sstream>

std::pair<mpz_class, mp_exp_t> mpfr_decomp(mpfr::mpreal const &val) {
  if (val == 0)
    return std::make_pair(mpz_class(0l), 0);
  mpz_t buffer;
  mpz_init(buffer);
  mp_exp_t exp;
  // compute normalized decomposition of the data
  // (in regard to the precision of the internal data
  // representation)
  exp = mpfr_get_z_2exp(buffer, val.mpfr_srcptr());
  while (mpz_divisible_ui_p(buffer, 2ul)) {
    mpz_divexact_ui(buffer, buffer, 2ul);
    ++exp;
  }
  mpz_class signif(buffer);
  mpz_clear(buffer);
  return std::make_pair(signif, exp);
}

void mav_vec_norm(mpz_t *norm, fplll::ZZ_mat<mpz_t> &basis) {

  mpz_t maxval;
  mpz_t iter;
  mpz_init(maxval);
  mpz_init(iter);
  mpz_set_ui(maxval, 0);

  for (int i = 0; i < basis.get_rows(); ++i) {
    mpz_set_ui(iter, 0);
    for (int j = 0; j < basis.get_cols(); ++j) {
      mpz_addmul(iter, basis(i, j).get_data(), basis(i, j).get_data());
    }
    if (mpz_cmp(iter, maxval) > 0)
      mpz_set(maxval, iter);
  }
  mpfr_t fmaxval;
  mpfr_init_set_z(fmaxval, maxval, GMP_RNDN);
  mpfr_sqrt(fmaxval, fmaxval, GMP_RNDN);
  mpfr_get_z(*norm, fmaxval, GMP_RNDD);

  mpz_clear(iter);
  mpz_clear(maxval);
  mpfr_clear(fmaxval);
}

void generate_cvp(fplll::ZZ_mat<mpz_t> &B, std::vector<mpz_class> &t,
                  std::vector<mpfr::mpreal> const &x,
                  std::vector<mpfr::mpreal> const &fx,
                  std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &basis,
                  mp_prec_t prec) {
  using mpfr::mpreal;
  mp_prec_t prev = mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  // store the min value for an exponent used to represent
  // in mantissa-exponent form the lattice basis and the
  // vector T to be approximated using the LLL approach
  mp_exp_t scale_exp = 0;
  // determine the scaling factor for the basis vectors and
  // the unknown vector to be approximated by T
  std::pair<mpz_class, mp_exp_t> decomp;

  for (std::size_t i{0u}; i < x.size(); ++i) {
    decomp = mpfr_decomp(fx[i]);
    if (decomp.second < scale_exp)
      scale_exp = decomp.second;
    for (std::size_t j{0u}; j < basis.size(); ++j) {
      decomp = mpfr_decomp(basis[j](x[i]));
      if (decomp.second < scale_exp)
        scale_exp = decomp.second;
    }
  }

  // scale the basis and vector t
  B.resize(basis.size(), x.size());
  mpz_t ibuff;
  mpz_init(ibuff);

  for (std::size_t i{0u}; i < basis.size(); ++i)
    for (std::size_t j{0u}; j < x.size(); ++j) {
      decomp = mpfr_decomp(basis[i](x[j]));
      decomp.second -= scale_exp;
      mpz_ui_pow_ui(ibuff, 2u, (unsigned int)decomp.second);
      mpz_mul(ibuff, ibuff, decomp.first.get_mpz_t());
      mpz_set(B(i, j).get_data(), ibuff);
    }
  t.clear();
  for (std::size_t i{0u}; i < fx.size(); ++i) {
    decomp = mpfr_decomp(fx[i]);
    decomp.second -= scale_exp;
    mpz_ui_pow_ui(ibuff, 2u, (unsigned int)decomp.second);
    mpz_mul(ibuff, ibuff, decomp.first.get_mpz_t());
    t.push_back(mpz_class(ibuff));
  }

  // clean-up
  mpz_clear(ibuff);

  mpreal::set_default_prec(prev);
}

void generate_cvp(fplll::ZZ_mat<mpz_t> &B, std::vector<mpz_class> &t,
                  std::vector<mpfr::mpreal> const &x,
                  std::vector<mpfr::mpreal> const &fx,
                  std::vector<mpfr::mpreal> const &wx,
                  std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &basis,
                  mp_prec_t prec) {
  using mpfr::mpreal;
  mp_prec_t prev = mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  // store the min value for an exponent used to represent
  // in mantissa-exponent form the lattice basis and the
  // vector T to be approximated using the LLL approach
  mp_exp_t scale_exp = 0;
  // determine the scaling factor for the basis vectors and
  // the unknown vector to be approximated by T
  std::pair<mpz_class, mp_exp_t> decomp;

  for (std::size_t i{0u}; i < x.size(); ++i) {
    decomp = mpfr_decomp(fx[i]);
    if (decomp.second < scale_exp)
      scale_exp = decomp.second;
    for (std::size_t j{0u}; j < basis.size(); ++j) {
      decomp = mpfr_decomp(wx[i] * basis[j](x[i]));
      if (decomp.second < scale_exp)
        scale_exp = decomp.second;
    }
  }

  // scale the basis and vector t
  B.resize(basis.size(), x.size());
  mpz_t ibuff;
  mpz_init(ibuff);

  for (std::size_t i{0u}; i < basis.size(); ++i)
    for (std::size_t j{0u}; j < x.size(); ++j) {
      decomp = mpfr_decomp(wx[j] * basis[i](x[j]));
      decomp.second -= scale_exp;
      mpz_ui_pow_ui(ibuff, 2u, (unsigned int)decomp.second);
      mpz_mul(ibuff, ibuff, decomp.first.get_mpz_t());
      mpz_set(B(i, j).get_data(), ibuff);
    }
  t.clear();
  for (std::size_t i{0u}; i < fx.size(); ++i) {
    decomp = mpfr_decomp(fx[i]);
    decomp.second -= scale_exp;
    mpz_ui_pow_ui(ibuff, 2u, (unsigned int)decomp.second);
    mpz_mul(ibuff, ibuff, decomp.first.get_mpz_t());
    t.push_back(mpz_class(ibuff));
  }

  // clean-up
  mpz_clear(ibuff);

  mpreal::set_default_prec(prev);
}

void kannan_embedding(fplll::ZZ_mat<mpz_t> &B, std::vector<mpz_class> &t) {

  int cols = B.get_cols() + 1;
  int rows = B.get_rows() + 1;

  mpz_t w; // the CVP t vector weight
  mpz_init(w);
  mav_vec_norm(&w, B);
  B.resize(rows, cols);

  for (int i = 0; i < rows - 1; ++i)
    mpz_set_ui(B(i, cols - 1).get_data(), 0);
  for (int j = 0; j < cols - 1; ++j)
    mpz_set(B(rows - 1, j).get_data(), t[j].get_mpz_t());
  mpz_set(B(rows - 1, cols - 1).get_data(), w);
  mpz_clear(w);
}

void acb_get_mpreal_real_imag(mpfr::mpreal &re, mpfr::mpreal &im, acb_t z) {
  acb_t mid;
  arb_t real, imag;
  acb_init(mid);
  arb_init(real);
  arb_init(imag);
  mpfr_t r, i;
  mpfr_init(r);
  mpfr_init(i);

  acb_get_mid(mid, z);
  acb_get_real(real, mid);
  arb_get_interval_mpfr(r, r, real);
  acb_get_imag(imag, mid);
  arb_get_interval_mpfr(i, i, imag);

  re = mpfr::mpreal(r);
  im = mpfr::mpreal(i);

  acb_clear(mid);
  arb_clear(real);
  arb_clear(imag);
  mpfr_clear(r);
  mpfr_clear(i);
}

void factorize(
    mpfr::mpreal &scale,
    std::vector<std::pair<std::vector<mpfr::mpreal>, std::size_t>> &factors,
    std::vector<mpfr::mpreal> const &coeffs, mp_prec_t prec) {
  using mpfr::mpreal;
  mp_prec_t prev = mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  factors.clear();
  mp_exp_t scale_exp = 0;
  std::pair<mpz_class, mp_exp_t> decomp;
  for (auto &coeff : coeffs) {
    decomp = mpfr_decomp(coeff);
    if (decomp.second < scale_exp)
      scale_exp = decomp.second;
  }

  fmpz_poly_t p;
  fmpz_poly_factor_t fac;
  acb_ptr roots;
  mpz_t ibuff;

  fmpz_poly_init(p);
  fmpz_poly_factor_init(fac);
  mpz_init(ibuff);

  for (std::size_t i{0u}; i < coeffs.size(); ++i) {
    decomp = mpfr_decomp(coeffs[i]);
    decomp.second -= scale_exp;
    mpz_ui_pow_ui(ibuff, 2u, (unsigned int)decomp.second);
    mpz_mul(ibuff, ibuff, decomp.first.get_mpz_t());
    fmpz_poly_set_coeff_mpz(p, i, ibuff);
  }
  scale = mpreal(ibuff);
  scale <<= scale_exp;

  fmpz_poly_factor_squarefree(fac, p);
  mpfr::mpreal re, im;

  for (int i{0}; i < fac->num; i++) {
    auto deg = fmpz_poly_degree(fac->p + i);
    roots = _acb_vec_init(deg);
    arb_fmpz_poly_complex_roots(roots, fac->p + i, 0, prec * 3);

    std::size_t idx{0u};
    while (idx < deg) {
      acb_get_mpreal_real_imag(re, im, roots + idx);
      std::vector<mpfr::mpreal> factor;
      if (im == 0) {
        // real root
        factor.push_back(-re);
        factor.push_back(mpfr::mpreal(1));
        ++idx;
      } else {
        // complex root; based on the description of the
        // arb_fmpz_poly_complex_roots function, there should
        // follow complex conjugate pairs for the roots
        factor.push_back(re * re + im * im);
        factor.push_back(-re * 2.0);
        factor.push_back(mpfr::mpreal(1));
        idx += 2u;
      }
      factors.push_back(std::make_pair(factor, fac->exp[i]));
    }

    _acb_vec_clear(roots, deg);
  }

  fmpz_poly_factor_clear(fac);
  fmpz_poly_clear(p);
  mpz_clear(ibuff);
  mpreal::set_default_prec(prev);
}

void fpirls(
    std::vector<mpfr::mpreal> &fpnum, std::vector<mpfr::mpreal> &fpden,
    std::function<mpfr::mpreal(mpfr::mpreal)> const &f,
    std::function<mpfr::mpreal(mpfr::mpreal)> const &r,
    std::function<mpfr::mpreal(mpfr::mpreal)> const &w,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> const &nbasis,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> const &dbasis,
    std::vector<mpfr::mpreal> const &num, std::vector<mpfr::mpreal> const &den,
    std::vector<mp_prec_t> const &nump, std::vector<mp_prec_t> const &denp,
    std::pair<mpfr::mpreal, mpfr::mpreal> const &dom, std::size_t idx,
    std::size_t itcount, mp_prec_t prec) {
  using mpfr::mpreal;
  mp_prec_t prev = mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  std::vector<mp_exp_t> numexp;
  std::vector<mp_exp_t> denexp;

  numexp.resize(nump.size());
  denexp.resize(denp.size());

  for (size_t i{0u}; i < num.size(); ++i) {
    mpfr::mpreal buff = num[i];
    buff.set_prec(nump[i], GMP_RNDN);
    auto decomp = mpfr_decomp(buff);
    numexp[i] = decomp.second;
  }

  for (size_t i{0u}; i < den.size(); ++i) {
    mpfr::mpreal buff = den[i];
    buff.set_prec(denp[i], GMP_RNDN);
    auto decomp = mpfr_decomp(buff);
    denexp[i] = decomp.second;
  }

  std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> bs;
  for (size_t i{0u}; i < nbasis.size(); ++i)
    bs.push_back([nbasis, i, w, numexp](mpfr::mpreal x) -> mpfr::mpreal {
      return (w(x) * nbasis[i](x)) << numexp[i];
    });

  for (size_t i{0u}; i < idx; ++i)
    bs.push_back([r, dbasis, i, w, denexp](mpfr::mpreal x) -> mpfr::mpreal {
      return (-r(x) * w(x) * dbasis[i](x)) << denexp[i];
    });

  for (size_t i{idx + 1u}; i < den.size(); ++i)
    bs.push_back([r, dbasis, i, w, denexp](mpfr::mpreal x) -> mpfr::mpreal {
      return (-r(x) * w(x) * dbasis[i](x)) << denexp[i];
    });

  std::function<mpfr::mpreal(mpfr::mpreal)> wr =
      [r, w, dbasis, den, idx](mpfr::mpreal x) -> mpfr::mpreal {
    return w(x) * den[idx] * r(x) * dbasis[idx](x);
  };

  std::vector<mpfr::mpreal> lllcoeffs;
  std::vector<std::vector<mpfr::mpreal>> svpcoeffs;

  std::size_t discsize = 10 * (nbasis.size() + dbasis.size());
  std::vector<mpfr::mpreal> disc(discsize);
  chebpts(disc, discsize);
  chgvar(disc, disc, dom);

  fpminimax_kernel(lllcoeffs, svpcoeffs, disc, wr, bs, prec);

  fpnum.resize(num.size());
  fpden.resize(den.size());
  for (size_t i{0u}; i < num.size(); ++i) {
    fpnum[i] = (lllcoeffs[i] << numexp[i]);
  }

  fpden[idx] = den[idx];
  for (size_t i{0u}; i < idx; ++i) {
    fpden[i] = (lllcoeffs[num.size() + i] << denexp[i]);
  }
  for (size_t i{idx + 1u}; i < den.size(); ++i) {
    fpden[i] = (lllcoeffs[num.size() + i - 1]) << denexp[i];
  }

  std::vector<mpfr::mpreal> wd(disc.size());
  for (auto &it : wd)
    it = 1.0;

  for (std::size_t itx{0u}; itx < itcount; ++itx) {
    std::function<mpfr::mpreal(mpfr::mpreal)> qfp;
    std::function<mpfr::mpreal(mpfr::mpreal)> pfp;

    qfp = [dbasis, fpden](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = 0;
      for (size_t i{0u}; i < fpden.size(); ++i)
        res += fpden[i] * dbasis[i](var);
      return res;
    };

    pfp = [nbasis, fpnum](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = 0;
      for (size_t i{0u}; i < fpnum.size(); ++i)
        res += fpnum[i] * nbasis[i](var);
      return res;
    };
    std::function<mpfr::mpreal(mpfr::mpreal)> errfp;
    errfp = [f, w, pfp, qfp](mpfr::mpreal var) -> mpfr::mpreal {
      return w(var) * (f(var) - pfp(var) / qfp(var));
    };

    std::pair<mpfr::mpreal, mpfr::mpreal> errfp_norm;
    infnorm(errfp_norm, errfp, dom);
    std::cout << "fpminimax approximation error                    = "
              << errfp_norm.second << std::endl;

    for (std::size_t i{0u}; i < wd.size(); ++i) {
      wd[i] *= mpfr::abs(errfp(disc[i]));
    }

    bs.clear();
    for (size_t i{0u}; i < nbasis.size(); ++i)
      bs.push_back([nbasis, i, w, numexp](mpfr::mpreal x) -> mpfr::mpreal {
        return (w(x) * nbasis[i](x)) << numexp[i];
      });

    for (size_t i{0u}; i < idx; ++i)
      bs.push_back([f, dbasis, i, w, denexp](mpfr::mpreal x) -> mpfr::mpreal {
        return (-f(x) * w(x) * dbasis[i](x)) << denexp[i];
      });

    for (size_t i{idx + 1u}; i < den.size(); ++i)
      bs.push_back([f, dbasis, i, w, denexp](mpfr::mpreal x) -> mpfr::mpreal {
        return (-f(x) * w(x) * dbasis[i](x)) << denexp[i];
      });

    std::function<mpfr::mpreal(mpfr::mpreal)> wf =
        [f, w, dbasis, den, idx](mpfr::mpreal x) -> mpfr::mpreal {
      return w(x) * den[idx] * f(x) * dbasis[idx](x);
    };

    fpirls_kernel(lllcoeffs, svpcoeffs, disc, wd, wf, bs, prec);

    fpnum.resize(num.size());
    fpden.resize(den.size());
    for (size_t i{0u}; i < num.size(); ++i) {
      fpnum[i] = (lllcoeffs[i] << numexp[i]);
    }

    fpden[idx] = den[idx];
    for (size_t i{0u}; i < idx; ++i) {
      fpden[i] = (lllcoeffs[num.size() + i] << denexp[i]);
    }
    for (size_t i{idx + 1u}; i < den.size(); ++i) {
      fpden[i] = (lllcoeffs[num.size() + i - 1]) << denexp[i];
    }
  }
}

void fpirls_kernel(
    std::vector<mpfr::mpreal> &lllcoeffs,
    std::vector<std::vector<mpfr::mpreal>> &svpcoeffs,
    std::vector<mpfr::mpreal> const &x, std::vector<mpfr::mpreal> const &wx,
    std::function<mpfr::mpreal(mpfr::mpreal)> const &target,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &basis,
    mp_prec_t prec) {
  using mpfr::mpreal;
  mp_prec_t prev = mpreal::get_default_prec();
  mpreal::set_default_prec(prec);
  std::size_t n = basis.size();

  std::vector<mpfr::mpreal> fx;
  fx.resize(x.size());

  for (std::size_t i{0u}; i < x.size(); ++i) {
    fx[i] = wx[i] * target(x[i]);
  }

  fplll::ZZ_mat<mpz_t> B;
  std::vector<mpz_class> t;
  generate_cvp(B, t, x, fx, wx, basis, prec);
  kannan_embedding(B, t);

  fplll::ZZ_mat<mpz_t> U(B.get_rows(), B.get_cols());
  fplll::lll_reduction(B, U, 0.99, 0.51);
  int xdp1 = (int)mpz_get_si(U(U.get_rows() - 1, U.get_cols() - 1).get_data());

  std::vector<mpz_class> illlcoeffs;
  std::vector<std::vector<mpz_class>> isvpcoeffs(n);
  xdp1 = (int)mpz_get_si(U(U.get_rows() - 1, U.get_cols() - 1).get_data());
  mpz_t coeffaux;
  mpz_init(coeffaux);
  switch (xdp1) {
  case 1:
    for (int i{0}; i < U.get_cols() - 1; ++i) {
      mpz_neg(coeffaux, U(U.get_rows() - 1, i).get_data());
      illlcoeffs.push_back(mpz_class(coeffaux));
      for (int j{0}; j < (int)n; ++j) {
        mpz_set(coeffaux, U(j, i).get_data());
        isvpcoeffs[j].push_back(mpz_class(coeffaux));
      }
    }
    break;
  case -1:
    for (int i = 0; i < U.get_cols() - 1; ++i) {
      illlcoeffs.push_back(mpz_class(U(U.get_rows() - 1, i).get_data()));
      for (int j = 0; j < (int)n; ++j) {
        mpz_set(coeffaux, U(j, i).get_data());
        isvpcoeffs[j].push_back(mpz_class(coeffaux));
      }
    }
    break;
  default:
    std::cout << "Failed to generate the approximation\n";
    exit(EXIT_FAILURE);
  }
  mpz_clear(coeffaux);

  svpcoeffs.resize(isvpcoeffs.size());
  lllcoeffs.resize(illlcoeffs.size());

  for (std::size_t i{0u}; i < illlcoeffs.size(); ++i) {

    mpreal newcoeff;
    lllcoeffs[i] = illlcoeffs[i].get_str();
    for (std::size_t j{0u}; j < svpcoeffs.size(); ++j) {
      mpreal buffer = isvpcoeffs[j][i].get_str();
      svpcoeffs[j].push_back(buffer);
    }
  }

  mpreal::set_default_prec(prev);
}

void fpminimax(
    std::vector<mpfr::mpreal> &fpnum, std::vector<mpfr::mpreal> &fpden,
    std::function<mpfr::mpreal(mpfr::mpreal)> const &r,
    std::function<mpfr::mpreal(mpfr::mpreal)> const &w,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> const &nbasis,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> const &dbasis,
    std::vector<mpfr::mpreal> const &num, std::vector<mpfr::mpreal> const &den,
    std::vector<mp_prec_t> const &nump, std::vector<mp_prec_t> const &denp,
    std::pair<mpfr::mpreal, mpfr::mpreal> const &dom, std::size_t idx,
    mp_prec_t prec) {
  using mpfr::mpreal;
  mp_prec_t prev = mpreal::get_default_prec();
  mpreal::set_default_prec(prec);

  std::vector<mp_exp_t> numexp;
  std::vector<mp_exp_t> denexp;

  numexp.resize(nump.size());
  denexp.resize(denp.size());

  for (size_t i{0u}; i < num.size(); ++i) {
    mpfr::mpreal buff = num[i];
    buff.set_prec(nump[i], GMP_RNDN);
    auto decomp = mpfr_decomp(buff);
    numexp[i] = decomp.second;
  }

  for (size_t i{0u}; i < den.size(); ++i) {
    mpfr::mpreal buff = den[i];
    buff.set_prec(denp[i], GMP_RNDN);
    auto decomp = mpfr_decomp(buff);
    denexp[i] = decomp.second;
  }

  std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> bs;
  for (size_t i{0u}; i < nbasis.size(); ++i)
    bs.push_back([nbasis, i, w, numexp](mpfr::mpreal x) -> mpfr::mpreal {
      return (w(x) * nbasis[i](x)) << numexp[i];
    });

  for (size_t i{0u}; i < idx; ++i)
    bs.push_back([r, dbasis, i, w, denexp](mpfr::mpreal x) -> mpfr::mpreal {
      return (-r(x) * w(x) * dbasis[i](x)) << denexp[i];
    });

  for (size_t i{idx + 1u}; i < den.size(); ++i)
    bs.push_back([r, dbasis, i, w, denexp](mpfr::mpreal x) -> mpfr::mpreal {
      return (-r(x) * w(x) * dbasis[i](x)) << denexp[i];
    });

  std::function<mpfr::mpreal(mpfr::mpreal)> wr =
      [r, w, dbasis, den, idx](mpfr::mpreal x) -> mpfr::mpreal {
    return w(x) * den[idx] * r(x) * dbasis[idx](x);
  };

  std::vector<mpfr::mpreal> lllcoeffs;
  std::vector<std::vector<mpfr::mpreal>> svpcoeffs;

  std::size_t discsize = 10 * (nbasis.size() + dbasis.size());
  std::vector<mpfr::mpreal> disc(discsize);
  chebpts(disc, discsize);
  chgvar(disc, disc, dom);

  fpminimax_kernel(lllcoeffs, svpcoeffs, disc, wr, bs, prec);

  fpnum.resize(num.size());
  fpden.resize(den.size());
  for (size_t i{0u}; i < num.size(); ++i) {
    fpnum[i] = (lllcoeffs[i] << numexp[i]);
  }

  fpden[idx] = den[idx];
  for (size_t i{0u}; i < idx; ++i) {
    fpden[i] = (lllcoeffs[num.size() + i] << denexp[i]);
  }
  for (size_t i{idx + 1u}; i < den.size(); ++i) {
    fpden[i] = (lllcoeffs[num.size() + i - 1]) << denexp[i];
  }

  mpreal::set_default_prec(prev);
}

void fpminimax_kernel(
    std::vector<mpfr::mpreal> &lllcoeffs,
    std::vector<std::vector<mpfr::mpreal>> &svpcoeffs,
    std::vector<mpfr::mpreal> const &x,
    std::function<mpfr::mpreal(mpfr::mpreal)> &target,
    std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> &basis,
    mp_prec_t prec) {
  using mpfr::mpreal;
  mp_prec_t prev = mpreal::get_default_prec();
  mpreal::set_default_prec(prec);
  std::size_t n = basis.size();

  std::vector<mpfr::mpreal> fx;
  fx.resize(x.size());

  for (std::size_t i{0u}; i < x.size(); ++i) {
    fx[i] = target(x[i]);
  }

  fplll::ZZ_mat<mpz_t> B;
  std::vector<mpz_class> t;
  generate_cvp(B, t, x, fx, basis, prec);
  kannan_embedding(B, t);

  fplll::ZZ_mat<mpz_t> U(B.get_rows(), B.get_cols());
  fplll::lll_reduction(B, U, 0.99, 0.51);
  int xdp1 = (int)mpz_get_si(U(U.get_rows() - 1, U.get_cols() - 1).get_data());

  std::vector<mpz_class> illlcoeffs;
  std::vector<std::vector<mpz_class>> isvpcoeffs(n);
  xdp1 = (int)mpz_get_si(U(U.get_rows() - 1, U.get_cols() - 1).get_data());
  mpz_t coeffaux;
  mpz_init(coeffaux);
  switch (xdp1) {
  case 1:
    for (int i{0}; i < U.get_cols() - 1; ++i) {
      mpz_neg(coeffaux, U(U.get_rows() - 1, i).get_data());
      illlcoeffs.push_back(mpz_class(coeffaux));
      for (int j{0}; j < (int)n; ++j) {
        mpz_set(coeffaux, U(j, i).get_data());
        isvpcoeffs[j].push_back(mpz_class(coeffaux));
      }
    }
    break;
  case -1:
    for (int i = 0; i < U.get_cols() - 1; ++i) {
      illlcoeffs.push_back(mpz_class(U(U.get_rows() - 1, i).get_data()));
      for (int j = 0; j < (int)n; ++j) {
        mpz_set(coeffaux, U(j, i).get_data());
        isvpcoeffs[j].push_back(mpz_class(coeffaux));
      }
    }
    break;
  default:
    std::cout << "Failed to generate the approximation\n";
    exit(EXIT_FAILURE);
  }
  mpz_clear(coeffaux);

  svpcoeffs.resize(isvpcoeffs.size());
  lllcoeffs.resize(illlcoeffs.size());

  for (std::size_t i{0u}; i < illlcoeffs.size(); ++i) {

    mpreal newcoeff;
    lllcoeffs[i] = illlcoeffs[i].get_str();
    for (std::size_t j{0u}; j < svpcoeffs.size(); ++j) {
      mpreal buffer = isvpcoeffs[j][i].get_str();
      svpcoeffs[j].push_back(buffer);
    }
  }

  mpreal::set_default_prec(prev);
}