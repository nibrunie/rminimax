#include "cheby.h"
#include "fpminimax.h"
#include "minimax.h"
#include "plotting.h"
#include "shuntingyard.h"
#include "tokenizer.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gmpxx.h>
#include <iostream>
#include <limits>
#include <mpfr.h>

extern "C" {
#include <qsopt_ex/QSopt_ex.h>
}

bool getCmdParameter(char *argv, const char *parameter, char *&value) {
  if (std::strstr(argv, parameter)) {
    value = argv + std::strlen(parameter);
    return true;
  } else {
    return false;
  }
}

void printShortHelp() {
  std::cout
      << "usage: ratapprox [OPTIONS]" << std::endl
      << std::endl
      << "---------------------------------------------------------------------"
         "-----------\n"
      << "General options:\n"
      << "Option:                       Meaning:\n"
      << "--help                        Prints this help.\n"
      << "--log                         Prints log information during "
         "execution.\n"
      << "--scalingSearch               Optionally searches for a better "
         "scaling factor.\n"
      << "--function=[string]           The function to approximate.\n"
      << "--weight=[string]             Weight function. By default this is\n"
      << "                              the reciprocal of the function.\n"
      << "--prec=uint                   Precision for internal MPFR "
         "floating-point\n"
      << "                              computations. It should usually be "
         "set\n"
      << "                              to a few hundred bits (default is "
         "200).\n"
      << "--num=[even|odd]              Specifies if the numerator should "
         "only\n"
      << "                              contain even/odd powered monomials.\n"
      << "--num=[string,string,...]     Specifies if the numerator should use "
         "a\n"
      << "                              custom basis where the functions are\n"
      << "                              given as parsable strings.\n"
      << "--den=[even|odd]              Specifies if the denominator should\n"
      << "                              contain only even/odd monomials.\n"
      << "--den=[string,string,...]     Specifies if the denominator should "
         "use\n"
      << "                              a custom basis where the functions "
         "are\n"
      << "                              given as parsable strings.\n"
      << "--denF=[int|string,...]       Specifies the list of floating point\n"
      << "                              formats for the den. coefficients.\n"
      << "                              They can either be numeric values or\n"
      << "                              strings, following Sollya notation:\n"
      << "                                  * HP - halfprecision\n"
      << "                                  * SG - singleprecision\n"
      << "                                  * D  - doubleprecision\n"
      << "                                  * DE - double extended\n"
      << "                                  * DD - double double\n"
      << "                                  * TD - triple double\n"
      << "                                  * QD - quad precision\n"
      << "                              In case it contains less elements "
         "than\n"
      << "                              the number of den. basis functions,\n"
      << "                              the last format is used for all\n"
      << "                              remaining coefficients. If it "
         "contains\n"
      << "                              more elements, then the extra formats\n"
      << "                              are discarded. By default, elements\n"
      << "                              will be optimized as double (D) "
         "values.\n"
      << "--numF=[int|string,...]       Similar to denF, only applied to the\n"
      << "                              numerator.\n"
      << "--type=[int,int]              Type of the rational approximation.\n"
      << "                              This option needs to be specified if\n"
      << "                              the monomial basis is used for the\n"
      << "                              numerator and denominator, or "
         "even/odd\n"
      << "                              alternatives.\n"
      << "--factorize                   Factorizes both the expressions in "
         "the\n"
      << "                              numerator and denominator into "
         "irreducible\n"
      << "                              degree one and degree two factors.\n"
      << "                              WARNING: results will be irrelevant if "
         "the\n"
      << "                              numerator and denominator are not "
         "polynomials\n"
      << "                              with a full set of basis functions.\n"
      << "--factorF=[string]            Specifies the floating point format "
         "for\n"
      << "                              the factorization coefficients. Can "
         "take\n"
      << "                              the same values as in numF and denF.\n"
      << "--dom=[double,double]         Approximation interval.\n"
      << "--output=[string]             Path to the output file.\n"
      << "                              By default string=./coeffs.sollya\n"
      << "--dispCoeff=[bin,dec,hex]     Display format for the approximation\n"
      << "                              coefficients in the output file, if\n"
      << "                              specified. 'hex' by default."
      << std::endl
      << std::endl;
}

bool is_uint(const std::string &s) {
  if (s.empty() || std::isspace(s[0]) || s[0] == '-')
    return false;
  char *p;
  strtol(s.c_str(), &p, 10);
  return (*p == 0);
}

void example(mp_prec_t prec) {
  // function to approximate
  std::function<mpfr::mpreal(mpfr::mpreal)> f =
      [](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::acos(x); };

  // weight function (here I just do unweighted approximation)
  // you can use 1/f(x) if you want relative error approximations
  std::function<mpfr::mpreal(mpfr::mpreal)> w =
      [](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::mpreal(1.0); };

  // the approximation interval (in this case [-1, 1])
  auto dom = std::make_pair(mpfr::mpreal(-1), mpfr::mpreal(1));
  // the type of the approximation (in this case a degree 5 polynomial)
  // the type is of the form (m, n) for a rational function, where m is the
  // degree of the numerator and n is the degree of the denominator
  auto type = std::make_pair(5u, 5u);

  // the basis functions for the numerator (1, x, x^2,...)
  std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> nbasis;
  // the basis function for the denominator (1, x, x^2,...)
  std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> dbasis;
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> dbounds;

  // construct the actual basis functions in the numerator and denominator
  for (size_t i{0u}; i <= type.first; ++i)
    nbasis.push_back(
        [i](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::pow(x, i); });
  for (size_t i{0u}; i <= type.second; ++i)
    dbasis.push_back(
        [i](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::pow(x, i); });

  // the general minimax approximation algorithm that computes the
  // best approximation with real (i.e., multiprecision) coefficients
  // (minimax is the analogous to Remez in Sollya)

  mpfr::mpreal delta;
  std::vector<mpfr::mpreal> num;
  std::vector<mpfr::mpreal> den;
  auto success = minimax(delta, num, den, dom, nbasis, dbasis, f, w, prec);

  // parameters for the floating-point coefficient approximation
  // (here single precision, float32 values)
  std::vector<mp_prec_t> numPrec;
  std::vector<mp_prec_t> denPrec;
  std::vector<mp_exp_t> numExp;
  std::vector<mp_exp_t> denExp;
  std::vector<mpfr::mpreal> fpnum;
  std::vector<mpfr::mpreal> fpden;
  for (size_t i{0u}; i < nbasis.size(); ++i)
    numPrec.push_back(24ul); // 24 bits of mantissa (including
                             // the implicit leading 1)
  for (size_t i{0u}; i < dbasis.size(); ++i)
    denPrec.push_back(24ul);

  // function handles for the approximation r = p/q
  std::function<mpfr::mpreal(mpfr::mpreal)> q;
  std::function<mpfr::mpreal(mpfr::mpreal)> p;

  q = [dbasis, den](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = 0;
    for (size_t i{0u}; i < den.size(); ++i)
      res += den[i] * dbasis[i](var);
    return res;
  };

  p = [nbasis, num](mpfr::mpreal var) -> mpfr::mpreal {
    mpfr::mpreal res = 0;
    for (size_t i{0u}; i < num.size(); ++i)
      res += num[i] * nbasis[i](var);
    return res;
  };

  std::function<mpfr::mpreal(mpfr::mpreal)> err;
  err = [f, w, p, q](mpfr::mpreal var) -> mpfr::mpreal {
    return w(var) * (f(var) - p(var) / q(var));
  };

  // the index of the normalizing coefficient needed for constructing
  // the CVP lattice problem solved inside fpminimax
  size_t idx{0u};
  for (size_t i{0u}; i < den.size(); ++i) {
    if (den[i] == 1 || den[i] == -1)
      idx = i;
  }

  std::function<mpfr::mpreal(mpfr::mpreal)> r;

  r = [p, q](mpfr::mpreal var) -> mpfr::mpreal { return p(var) / q(var); };

  // the fpminimax algorithm that will compute the approximation
  // with floating point coefficients
  fpminimax(fpnum, fpden, r, w, nbasis, dbasis, num, den, numPrec, denPrec, dom,
            idx, prec);

  // compute error of the fpminimax-type result
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
  std::string fpName = "fpminimax";
  plotFunc(fpName, errfp, dom.first, dom.second, prec);
  std::pair<mpfr::mpreal, mpfr::mpreal> errfp_norm;
  infnorm(errfp_norm, errfp, dom);

  std::cout << "fpminimax error = " << errfp_norm.second << std::endl;
}

void rminimax(int argc, char *argv[]) {
  mp_prec_t prec = 200;
  mpfr::mpreal::set_default_prec(prec);
  std::pair<mpfr::mpreal, mpfr::mpreal> dom =
      std::make_pair(mpfr::mpreal(-1), mpfr::mpreal(1));
  std::pair<size_t, size_t> type = std::make_pair(5u, 5u);
  bool useDomArg{false}, useTypeArg{false}, useWeightArg{false},
      useCustomNumBasis{false}, EONumBasis{false}, EODenBasis{false},
      numBasisEven{false}, useCustomDenBasis{false}, denBasisEven{false},
      useLog{false}, allowScaling{false}, factorizationEnabled{false},
      useCustomDenFormat{false}, useCustomNumFormat{false},
      useCustomInternalPrec{false};

  std::stringstream domstream;
  std::stringstream typestream;
  std::stringstream numbasisstream;
  std::stringstream denbasisstream;
  std::stringstream numformatstream;
  std::stringstream denformatstream;
  std::string facformatstr = "D";
  std::string functionstr = "exp(x)";
  std::string weightstr = "1";
  std::string outfilename = "coeffs.sollya";
  std::string coeffDisplay = "%Ra";
  std::string precstr;

  shuntingyard sh;
  std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> nbasis;
  std::vector<std::vector<std::string>> nbToks;
  std::vector<std::function<mpfr::mpreal(mpfr::mpreal)>> dbasis;
  std::vector<std::vector<std::string>> dbToks;
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> dbounds;

  for (int i{1}; i < argc; ++i) {
    char *value;
    if (getCmdParameter(argv[i], "--help", value)) {
      printShortHelp();
      QSexactClear();
      exit(-1);
    } else if (getCmdParameter(argv[i], "--scalingSearch", value)) {
      allowScaling = true;
    } else if (getCmdParameter(argv[i], "--log", value)) {
      useLog = true;
    } else if (getCmdParameter(argv[i], "--dom=", value)) {
      useDomArg = true;
      domstream << value;
    } else if (getCmdParameter(argv[i], "--prec=", value)) {
      useCustomInternalPrec = true;
      precstr = value;
    } else if (getCmdParameter(argv[i], "--type=", value)) {
      useTypeArg = true;
      typestream << value;
    } else if (getCmdParameter(argv[i], "--output=", value)) {
      outfilename = value;
    } else if (getCmdParameter(argv[i], "--function=", value)) {
      functionstr = value;
    } else if (getCmdParameter(argv[i], "--weight=", value)) {
      useWeightArg = true;
      weightstr = value;
    } else if (getCmdParameter(argv[i], "--num=", value)) {
      useCustomNumBasis = true;
      numbasisstream << value;
    } else if (getCmdParameter(argv[i], "--den=", value)) {
      useCustomDenBasis = true;
      denbasisstream << value;
    } else if (getCmdParameter(argv[i], "--denF=", value)) {
      useCustomDenFormat = true;
      denformatstream << value;
    } else if (getCmdParameter(argv[i], "--numF=", value)) {
      useCustomNumFormat = true;
      numformatstream << value;
    } else if (getCmdParameter(argv[i], "--factorize", value)) {
      factorizationEnabled = true;
    } else if (getCmdParameter(argv[i], "--factorF=", value)) {
      facformatstr = value;
    } else if (getCmdParameter(argv[i], "--dispCoeff=", value)) {
      if (strcmp("dec", value) == 0) {
        coeffDisplay = "%Re";
      } else if (strcmp("bin", value) == 0) {
        coeffDisplay = "%Rb";
      } else {
        coeffDisplay = "%Ra";
      }
    } else {
      std::cout << "Error: Illegal option: " << argv[i] << std::endl;
      printShortHelp();
      QSexactClear();
      exit(-1);
    }
  }

  if (useCustomInternalPrec) {
    mp_prec_t nprec = strtoul(precstr.c_str(), nullptr, 10);
    if (nprec == 0ul) {
      std::cout
          << "Error while parsing value in prec: Incompatible input format\n";
      printShortHelp();
      QSexactClear();
      exit(-1);
    }
    prec = nprec;
  }

  mpf_QSset_precision(prec);
  mpfr::mpreal::set_default_prec(prec);

  if (useDomArg) {
    double left, right;
    int parsed = sscanf(domstream.str().c_str(), "[%la,%la]", &left, &right);
    dom.first = left;
    dom.second = right;
    if (parsed != 2) {
      std::cout
          << "Error while parsing values in dom: Incompatible input format\n";
      printShortHelp();
      QSexactClear();
      exit(-1);
    }
  }

  if (useTypeArg) {
    char leftBracket = typestream.get();
    if (leftBracket != '[') {
      std::cout
          << "Error while parsing '[' in type: Incompatible input format\n";
      printShortHelp();
      QSexactClear();
      exit(-1);
    }
    char boundaryDelimiter;
    typestream >> type.first >> boundaryDelimiter >> type.second;
    if (boundaryDelimiter != ',') {
      std::cout
          << "Error while parsing ',' in type: Incompatible input format\n";
      printShortHelp();
      QSexactClear();
      exit(-1);
    }
    char rightBracket = typestream.get();
    if (rightBracket != ']') {
      std::cout
          << "Error while parsing ']' in type: Incompatible input format\n";
      printShortHelp();
      QSexactClear();
      exit(-1);
    }
  }

  std::vector<std::string> numFuncs;

  if (useCustomNumBasis) {
    std::string numbasisstr = numbasisstream.str();
    if (numbasisstr.size() <= 2u) {
      std::cout << "Error while parsing num: value too short\n";
      printShortHelp();
      QSexactClear();
      exit(-1);
    }
    char basisDelim = numbasisstream.peek();
    if (basisDelim == '[') {
      if (numbasisstr.back() != ']') {
        std::cout
            << "Error while parsing ']' in num: no corresponding ']' found\n";
        printShortHelp();
        QSexactClear();
        exit(-1);
      }
      numbasisstr = numbasisstr.substr(1, numbasisstr.size() - 2u);
      numbasisstream.str(numbasisstr);
      while (std::getline(numbasisstream, numbasisstr, ','))
        numFuncs.push_back(numbasisstr);

      nbasis.resize(numFuncs.size());
      nbToks.resize(numFuncs.size());
      for (size_t i{0u}; i < numFuncs.size(); ++i) {
        nbToks[i] = tokenizer(numFuncs[i]).getTokens();
        nbasis[i] = [sh, i, &nbToks](mpfr::mpreal x) -> mpfr::mpreal {
          return sh.evaluate(nbToks[i], x);
        };
      }
    } else {
      if (numbasisstr == "even") {
        numBasisEven = true;
        EONumBasis = true;
      } else if (numbasisstr == "odd") {
        numBasisEven = false;
        EONumBasis = true;
      } else {
        std::cout << "Error while parsing num: invalid string\n";
        printShortHelp();
        QSexactClear();
        exit(-1);
      }
    }
  }

  std::vector<std::string> denFuncs;
  if (useCustomDenBasis) {
    std::string denbasisstr = denbasisstream.str();

    if (denbasisstr.size() <= 2u) {
      std::cout << "Error while parsing den: value too short\n";
      printShortHelp();
      QSexactClear();
      exit(-1);
    }
    char basisDelim = denbasisstream.peek();
    if (basisDelim == '[') {
      if (denbasisstr.back() != ']') {
        std::cout
            << "Error while parsing ']' in den: no corresponding ']' found\n";
        printShortHelp();
        QSexactClear();
        exit(-1);
      }
      denbasisstr = denbasisstr.substr(1, denbasisstr.size() - 2u);
      denbasisstream.str(denbasisstr);
      while (std::getline(denbasisstream, denbasisstr, ','))
        denFuncs.push_back(denbasisstr);

      dbasis.resize(denFuncs.size());
      dbToks.resize(denFuncs.size());
      for (size_t i{0u}; i < denFuncs.size(); ++i) {
        dbToks[i] = tokenizer(denFuncs[i]).getTokens();
        dbasis[i] = [sh, i, &dbToks](mpfr::mpreal x) -> mpfr::mpreal {
          return sh.evaluate(dbToks[i], x);
        };
      }
    } else {
      if (denbasisstr == "even") {
        denBasisEven = true;
        EODenBasis = true;
      } else if (denbasisstr == "odd") {
        denBasisEven = false;
        EODenBasis = true;
      } else {
        std::cout << "Error while parsing den: invalid string\n";
        printShortHelp();
        QSexactClear();
        exit(-1);
      }
    }
  }

  if (EODenBasis) {
    size_t start{0u};
    if (!denBasisEven)
      start = 1u;
    for (size_t i{start}; i <= type.second; i += 2u)
      dbasis.push_back(
          [i](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::pow(x, i); });

    if (!EONumBasis && !useCustomNumBasis)
      for (size_t i{0u}; i <= type.first; ++i)
        nbasis.push_back(
            [i](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::pow(x, i); });
  }

  if (EONumBasis) {
    size_t start{0u};
    if (!numBasisEven)
      start = 1u;
    for (size_t i{start}; i <= type.first; i += 2u)
      nbasis.push_back(
          [i](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::pow(x, i); });

    if (!EODenBasis && !useCustomDenBasis) {
      for (size_t i{0u}; i <= type.second; ++i)
        dbasis.push_back(
            [i](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::pow(x, i); });
    }
  }

  if (useCustomNumBasis || useCustomDenBasis) {
    dbounds.resize(dbasis.size());
    for (size_t i{0u}; i < dbasis.size(); ++i) {
      dbounds[i].first = -1;
      dbounds[i].second = 1;
    }
  }

  if (useCustomNumBasis ^ useCustomDenBasis) {
    if (!EONumBasis && !EODenBasis) {
      std::cout
          << "Error while parsing num & den: incomplete basis information\n";
      printShortHelp();
      QSexactClear();
      exit(-1);
    }
  }

  auto ftok = tokenizer(functionstr).getTokens();
  auto wtok = tokenizer(weightstr).getTokens();
  std::function<mpfr::mpreal(mpfr::mpreal)> f =
      [sh, &ftok](mpfr::mpreal x) -> mpfr::mpreal {
    return sh.evaluate(ftok, x);
  };

  std::function<mpfr::mpreal(mpfr::mpreal)> w;
  if (useWeightArg) {
    w = [sh, &wtok](mpfr::mpreal x) -> mpfr::mpreal {
      return sh.evaluate(wtok, x);
    };
  } else {
    w = [f](mpfr::mpreal x) -> mpfr::mpreal {
      return mpfr::mpreal(1.0) / f(x);
    };
  }

  std::vector<mpfr::mpreal> num;
  std::vector<mpfr::mpreal> den;
  mpfr::mpreal delta;

  std::vector<mpfr::mpreal> x(type.first + type.second + 2u);
  chebpts(x, x.size());
  chgvar(x, x, dom);

  std::pair<bool, std::vector<mpfr::mpreal>> success = std::make_pair(true, x);
  if (useCustomNumBasis || useCustomDenBasis) {
    x.resize(100u * (nbasis.size() + dbasis.size()) + 1u);
    for (size_t i{0}; i < x.size(); ++i)
      x[i] = dom.first + (dom.second - dom.first) * i / (x.size() - 1);
    if (useLog)
      std::cout << "Attempting general minimax approximation algorithm...\n\n";

    success = minimax(delta, num, den, dom, nbasis, dbasis, dbounds, x, f, w,
                      useLog, prec);
  } else {
    // prepare the basis function vectors
    nbasis.resize(type.first + 1u);
    dbasis.resize(type.second + 1u);
    for (int i{0}; i <= type.first; ++i)
      nbasis[i] = [i](mpfr::mpreal x) -> mpfr::mpreal {
        return mpfr::pow(x, i);
      };
    for (int i{0}; i <= type.second; ++i)
      dbasis[i] = [i](mpfr::mpreal x) -> mpfr::mpreal {
        return mpfr::pow(x, i);
      };
    success = minimax(delta, num, den, dom, f, w, type, useLog, prec);
  }

  if (success.first) {
    size_t idx{0u};
    for (size_t i{0u}; i < den.size(); ++i) {
      if (den[i] == 1 || den[i] == -1)
        idx = i;
    }

    std::function<mpfr::mpreal(mpfr::mpreal)> q;
    std::function<mpfr::mpreal(mpfr::mpreal)> p;

    if (useLog) {
      std::cout << "Final coefficients:\n";
      std::cout << "Numerator:\n";
      for (auto &it : num)
        std::cout << it << " ";
      std::cout << "\nDenominator:\n";
      for (auto &it : den)
        std::cout << it << " ";
      std::cout << std::endl;
    }

    q = [dbasis, den](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = 0;
      for (size_t i{0u}; i < den.size(); ++i)
        res += den[i] * dbasis[i](var);
      return res;
    };

    p = [nbasis, num](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = 0;
      for (size_t i{0u}; i < num.size(); ++i)
        res += num[i] * nbasis[i](var);
      return res;
    };

    std::function<mpfr::mpreal(mpfr::mpreal)> err;
    err = [f, w, p, q](mpfr::mpreal var) -> mpfr::mpreal {
      return w(var) * (f(var) - p(var) / q(var));
    };

    std::string initName = "init";
    plotFunc(initName, err, dom.first, dom.second, prec);

    mpfr::mpreal numscale;
    std::vector<std::pair<std::vector<mpfr::mpreal>, std::size_t>> numfactors;
    mpfr::mpreal denscale;
    std::vector<std::pair<std::vector<mpfr::mpreal>, std::size_t>> denfactors;
    std::function<mpfr::mpreal(mpfr::mpreal)> fp, fq, ferr;
    std::pair<mpfr::mpreal, mpfr::mpreal> ferr_norm;

    mpfr::mpreal rnumscale;
    std::vector<std::pair<std::vector<mpfr::mpreal>, std::size_t>> rnumfactors;
    mpfr::mpreal rdenscale;
    std::vector<std::pair<std::vector<mpfr::mpreal>, std::size_t>> rdenfactors;
    std::function<mpfr::mpreal(mpfr::mpreal)> rfp, rfq, rferr;
    std::pair<mpfr::mpreal, mpfr::mpreal> rferr_norm;

    if (factorizationEnabled) {
      if (useLog)
        std::cout << "Starting factorization...\n";
      factorize(numscale, numfactors, num, prec);
      factorize(denscale, denfactors, den, prec);
      if (useLog)
        std::cout << "Finished factorization...\n";

      rnumscale = numscale;
      rdenscale = denscale;
      rnumfactors = numfactors;
      rdenfactors = denfactors;

      fp = [numscale, numfactors](mpfr::mpreal var) -> mpfr::mpreal {
        mpfr::mpreal res = numscale;
        for (std::size_t i{0u}; i < numfactors.size(); ++i) {
          mpfr::mpreal factor = 0ul;
          for (std::size_t j{0u}; j < numfactors[i].first.size(); ++j)
            factor += numfactors[i].first[j] * mpfr::pow(var, j);
          for (std::size_t j{0u}; j < numfactors[i].second; ++j)
            res *= factor;
        }
        return res;
      };

      fq = [denscale, denfactors](mpfr::mpreal var) -> mpfr::mpreal {
        mpfr::mpreal res = denscale;
        for (std::size_t i{0u}; i < denfactors.size(); ++i) {
          mpfr::mpreal factor = 0ul;
          for (std::size_t j{0u}; j < denfactors[i].first.size(); ++j)
            factor += denfactors[i].first[j] * mpfr::pow(var, j);
          for (std::size_t j{0u}; j < denfactors[i].second; ++j)
            res *= factor;
        }
        return res;
      };

      ferr = [f, w, fp, fq](mpfr::mpreal var) -> mpfr::mpreal {
        return w(var) * (f(var) - fp(var) / fq(var));
      };

      std::string factorName = "factorize";
      plotFunc(factorName, ferr, dom.first, dom.second, prec);
      infnorm(ferr_norm, ferr, dom);
    }

    std::vector<mp_prec_t> numPrec;
    std::vector<mp_prec_t> denPrec;
    std::vector<mp_exp_t> numExp;
    std::vector<mp_exp_t> denExp;
    std::vector<mpfr::mpreal> fpnum, rdnum;
    std::vector<mpfr::mpreal> fpden, rdden;

    if (useCustomDenFormat) {
      std::string denformatstr = denformatstream.str();
      if (denformatstr.size() <= 2u) {
        std::cout << "Error while parsing denF: value too short\n";
        printShortHelp();
        QSexactClear();
        exit(-1);
      }
      char formatDelim = denformatstream.peek();
      if (formatDelim == '[') {
        if (denformatstr.back() != ']') {
          std::cout << "Error while parsing ']' in denF: no corresponding ']' "
                       "found\n";
          printShortHelp();
          QSexactClear();
          exit(-1);
        }
        denformatstr = denformatstr.substr(1, denformatstr.size() - 2u);
        denformatstream.str(denformatstr);
        while (std::getline(denformatstream, denformatstr, ',')) {
          if (is_uint(denformatstr)) {
            mpfr::mpreal buff = denformatstr;
            denPrec.push_back(buff.toULong());
          } else if (denformatstr == "HP") {
            denPrec.push_back(11ul);
          } else if (denformatstr == "SG") {
            denPrec.push_back(24ul);
          } else if (denformatstr == "D") {
            denPrec.push_back(53ul);
          } else if (denformatstr == "DE") {
            denPrec.push_back(64ul);
          } else if (denformatstr == "DD") {
            denPrec.push_back(107ul);
          } else if (denformatstr == "TD") {
            denPrec.push_back(113ul);
          } else if (denformatstr == "QD") {
            denPrec.push_back(161ul);
          } else {
            std::cout << "Error while parsing denF: invalid string\n";
            printShortHelp();
            QSexactClear();
            exit(-1);
          }
        }
        if (denPrec.size() > den.size())
          denPrec.resize(den.size());
        else if (denPrec.size() < den.size()) {
          mp_prec_t buff = denPrec[denPrec.size() - 1];
          std::size_t oldsize = denPrec.size();
          for (int i{0}; i < den.size() - oldsize; ++i)
            denPrec.push_back(buff);
        }
      } else {
        std::cout << "Error while parsing denF: invalid string\n";
        printShortHelp();
        QSexactClear();
        exit(-1);
      }
    } else {
      for (size_t i{0u}; i < den.size(); ++i)
        denPrec.push_back(53ul);
    }

    if (useCustomNumFormat) {
      std::string numformatstr = numformatstream.str();
      if (numformatstr.size() <= 2u) {
        std::cout << "Error while parsing numF: value too short\n";
        printShortHelp();
        QSexactClear();
        exit(-1);
      }
      char formatDelim = numformatstream.peek();
      if (formatDelim == '[') {
        if (numformatstr.back() != ']') {
          std::cout << "Error while parsing ']' in numF: no corresponding ']' "
                       "found\n";
          printShortHelp();
          QSexactClear();
          exit(-1);
        }
        numformatstr = numformatstr.substr(1, numformatstr.size() - 2u);
        numformatstream.str(numformatstr);
        while (std::getline(numformatstream, numformatstr, ',')) {
          if (is_uint(numformatstr)) {
            mpfr::mpreal buff = numformatstr;
            numPrec.push_back(buff.toULong());
          } else if (numformatstr == "HP") {
            numPrec.push_back(11ul);
          } else if (numformatstr == "SG") {
            numPrec.push_back(24ul);
          } else if (numformatstr == "D") {
            numPrec.push_back(53ul);
          } else if (numformatstr == "DE") {
            numPrec.push_back(64ul);
          } else if (numformatstr == "DD") {
            numPrec.push_back(107ul);
          } else if (numformatstr == "TD") {
            numPrec.push_back(113ul);
          } else if (numformatstr == "QD") {
            numPrec.push_back(161ul);
          } else {
            std::cout << "Error while parsing numF: invalid string\n";
            printShortHelp();
            QSexactClear();
            exit(-1);
          }
        }
        if (numPrec.size() > num.size())
          numPrec.resize(num.size());
        else if (numPrec.size() < num.size()) {
          mp_prec_t buff = numPrec[numPrec.size() - 1];
          std::size_t oldsize = numPrec.size();
          for (int i{0}; i < num.size() - oldsize; ++i)
            numPrec.push_back(buff);
        }
      } else {
        std::cout << "Error while parsing numF: invalid string\n";
        printShortHelp();
        QSexactClear();
        exit(-1);
      }
    } else {
      for (size_t i{0u}; i < num.size(); ++i)
        numPrec.push_back(53ul);
    }

    if (factorizationEnabled) {
      if (facformatstr == "SG") {
        rnumscale = rnumscale.toFloat();
        rdenscale = rdenscale.toFloat();
        for (std::size_t i{0u}; i < rnumfactors.size(); ++i)
          for (std::size_t j{0u}; j < rnumfactors[i].first.size(); ++j)
            rnumfactors[i].first[j] = rnumfactors[i].first[j].toFloat();
        for (std::size_t i{0u}; i < rdenfactors.size(); ++i)
          for (std::size_t j{0u}; j < rdenfactors[i].first.size(); ++j)
            rdenfactors[i].first[j] = rdenfactors[i].first[j].toFloat();
      } else if (facformatstr == "D") {
        rnumscale = rnumscale.toDouble();
        rdenscale = rdenscale.toDouble();
        for (std::size_t i{0u}; i < rnumfactors.size(); ++i)
          for (std::size_t j{0u}; j < rnumfactors[i].first.size(); ++j)
            rnumfactors[i].first[j] = rnumfactors[i].first[j].toDouble();
        for (std::size_t i{0u}; i < rdenfactors.size(); ++i)
          for (std::size_t j{0u}; j < rdenfactors[i].first.size(); ++j)
            rdenfactors[i].first[j] = rdenfactors[i].first[j].toDouble();
      } else if (facformatstr == "DE") {
        rnumscale = rnumscale.toLDouble();
        rdenscale = rdenscale.toLDouble();
        for (std::size_t i{0u}; i < rnumfactors.size(); ++i)
          for (std::size_t j{0u}; j < rnumfactors[i].first.size(); ++j)
            rnumfactors[i].first[j] = rnumfactors[i].first[j].toLDouble();
        for (std::size_t i{0u}; i < rdenfactors.size(); ++i)
          for (std::size_t j{0u}; j < rdenfactors[i].first.size(); ++j)
            rdenfactors[i].first[j] = rdenfactors[i].first[j].toLDouble();
      } else {
        std::cout << "Error while parsing numF: invalid string\n";
        printShortHelp();
        QSexactClear();
        exit(-1);
      }

      rfp = [rnumscale, rnumfactors](mpfr::mpreal var) -> mpfr::mpreal {
        mpfr::mpreal res = rnumscale;
        for (std::size_t i{0u}; i < rnumfactors.size(); ++i) {
          mpfr::mpreal factor = 0ul;
          for (std::size_t j{0u}; j < rnumfactors[i].first.size(); ++j)
            factor += rnumfactors[i].first[j] * mpfr::pow(var, j);
          for (std::size_t j{0u}; j < rnumfactors[i].second; ++j)
            res *= factor;
        }
        return res;
      };

      rfq = [rdenscale, rdenfactors](mpfr::mpreal var) -> mpfr::mpreal {
        mpfr::mpreal res = rdenscale;
        for (std::size_t i{0u}; i < rdenfactors.size(); ++i) {
          mpfr::mpreal factor = 0ul;
          for (std::size_t j{0u}; j < rdenfactors[i].first.size(); ++j)
            factor += rdenfactors[i].first[j] * mpfr::pow(var, j);
          for (std::size_t j{0u}; j < rdenfactors[i].second; ++j)
            res *= factor;
        }
        return res;
      };

      rferr = [f, w, rfp, rfq](mpfr::mpreal var) -> mpfr::mpreal {
        return w(var) * (f(var) - rfp(var) / rfq(var));
      };

      std::string rfactorName = "factorize_rounded";
      plotFunc(rfactorName, rferr, dom.first, dom.second, prec);
      infnorm(rferr_norm, rferr, dom);
    }

    // compute information regarding the naive rounding result

    rdnum.resize(num.size());
    rdden.resize(den.size());

    for (size_t i{0u}; i < num.size(); ++i) {
      mpfr::mpreal buff = num[i];
      buff.set_prec(numPrec[i], GMP_RNDN);
      rdnum[i] = buff;
    }

    for (size_t i{0u}; i < den.size(); ++i) {
      mpfr::mpreal buff = den[i];
      buff.set_prec(denPrec[i], GMP_RNDN);
      rdden[i] = buff;
    }

    std::function<mpfr::mpreal(mpfr::mpreal)> qr;
    std::function<mpfr::mpreal(mpfr::mpreal)> pr;

    qr = [dbasis, rdden](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = 0;
      for (size_t i{0u}; i < rdden.size(); ++i)
        res += rdden[i] * dbasis[i](var);
      return res;
    };

    pr = [nbasis, rdnum](mpfr::mpreal var) -> mpfr::mpreal {
      mpfr::mpreal res = 0;
      for (size_t i{0u}; i < rdnum.size(); ++i)
        res += rdnum[i] * nbasis[i](var);
      return res;
    };

    std::function<mpfr::mpreal(mpfr::mpreal)> errr;
    errr = [f, w, pr, qr](mpfr::mpreal var) -> mpfr::mpreal {
      return w(var) * (f(var) - pr(var) / qr(var));
    };

    std::pair<mpfr::mpreal, mpfr::mpreal> errr_norm;
    std::string roundName = "round";
    plotFunc(roundName, errr, dom.first, dom.second, prec);
    infnorm(errr_norm, errr, dom);

    std::function<mpfr::mpreal(mpfr::mpreal)> r;

    r = [p, q](mpfr::mpreal var) -> mpfr::mpreal { return p(var) / q(var); };

    if (!allowScaling) {
      if (useLog)
        std::cout << "Starting fpminimax...\n";
      fpminimax(fpnum, fpden, r, w, nbasis, dbasis, num, den, numPrec, denPrec,
                dom, idx, prec);
      if (useLog)
        std::cout << "Finished fpminimax...\n";
    } else {
      if (useLog)
        std::cout << "Starting scaling factor exploration...\n";
      mpfr::mpreal scaleFactor = 1.0;
      std::size_t vals = 1000;
      std::vector<mpfr::mpreal> qsnum, snum;
      std::vector<mpfr::mpreal> qsden, sden;
      snum.resize(num.size());
      sden.resize(den.size());
      mpfr::mpreal minerr;
      for (int itx{0}; itx < vals - 1; ++itx) {
        std::function<mpfr::mpreal(mpfr::mpreal)> sp, sq;
        for (std::size_t i{0u}; i < num.size(); ++i) {
          snum[i] = num[i] * scaleFactor;
        }
        for (std::size_t i{0}; i < den.size(); ++i) {
          sden[i] = den[i] * scaleFactor;
        }

        sq = [dbasis, sden](mpfr::mpreal var) -> mpfr::mpreal {
          mpfr::mpreal res = 0;
          for (std::size_t i{0u}; i < sden.size(); ++i)
            res += sden[i] * dbasis[i](var);
          return res;
        };
        sp = [nbasis, snum](mpfr::mpreal var) -> mpfr::mpreal {
          mpfr::mpreal res = 0;
          for (std::size_t i{0u}; i < snum.size(); ++i)
            res += snum[i] * nbasis[i](var);
          return res;
        };

        // the index of the normalizing coefficient needed for constructing
        // the CVP lattice problem solved inside fpminimax
        std::size_t idx{0u};
        for (size_t i{0u}; i < sden.size(); ++i) {
          if (sden[i] == scaleFactor || sden[i] == -scaleFactor)
            idx = i;
        }

        std::function<mpfr::mpreal(mpfr::mpreal)> r;
        r = [sp, sq](mpfr::mpreal var) -> mpfr::mpreal {
          return sp(var) / sq(var);
        };

        fpminimax(qsnum, qsden, r, w, nbasis, dbasis, snum, sden, numPrec,
                  denPrec, dom, idx, prec);

        // compute error of the fpminimax-type result
        std::function<mpfr::mpreal(mpfr::mpreal)> qsp;
        std::function<mpfr::mpreal(mpfr::mpreal)> qsq;

        qsq = [dbasis, qsden](mpfr::mpreal var) -> mpfr::mpreal {
          mpfr::mpreal res = 0;
          for (size_t i{0u}; i < qsden.size(); ++i)
            res += qsden[i] * dbasis[i](var);
          return res;
        };

        qsp = [nbasis, qsnum](mpfr::mpreal var) -> mpfr::mpreal {
          mpfr::mpreal res = 0;
          for (size_t i{0u}; i < qsnum.size(); ++i)
            res += qsnum[i] * nbasis[i](var);
          return res;
        };

        std::function<mpfr::mpreal(mpfr::mpreal)> errfp;
        errfp = [f, w, qsp, qsq](mpfr::mpreal var) -> mpfr::mpreal {
          return w(var) * (f(var) - qsp(var) / qsq(var));
        };

        std::pair<mpfr::mpreal, mpfr::mpreal> qerr;
        infnorm(qerr, errfp, dom);
        if (itx == 0) {
          minerr = qerr.second;
          fpnum = qsnum;
          fpden = qsden;
        } else {
          if (qerr.second < minerr) {
            minerr = qerr.second;
            fpnum = qsnum;
            fpden = qsden;
          }
        }

        scaleFactor = 1.0 + (double)(itx + 1) / (double)vals;
      }
      if (useLog)
        std::cout << "Finished scaling factor exploration...\n";
    }

    // compute information regarding the fpminimax result
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
    std::string fpName = "fpminimax";
    plotFunc(fpName, errfp, dom.first, dom.second, prec);
    std::pair<mpfr::mpreal, mpfr::mpreal> errfp_norm;
    infnorm(errfp_norm, errfp, dom);

    if (useLog) {
      std::cout << "-----------------------INFORMATION--------------------\n";
      std::cout << "Function = " << functionstr << std::endl;
      if (useWeightArg) {
        std::cout << "Weight = " << weightstr << std::endl;
      } else {
        std::cout << "Weight = reciprocal of " << functionstr << std::endl;
      }
      std::cout << "Domain: [" << dom.first << ", " << dom.second << "]\n";
      if (useCustomDenBasis || useCustomNumBasis) {
        if (EONumBasis) {
          if (numBasisEven) {
            std::cout << "Numerator basis = [1";
            for (size_t i{2u}; i <= type.first; i += 2u)
              std::cout << ", x^" << i;
            std::cout << "]\n";
          } else {
            std::cout << "Numerator basis = [x";
            for (size_t i{3u}; i <= type.first; i += 2u)
              std::cout << ", x^" << i;
            std::cout << "]\n";
          }
        } else if (useCustomNumBasis) {
          std::cout << "Numerator basis = " << numbasisstream.str()
                    << std::endl;
        } else {
          std::cout << "Numerator basis = [1";
          for (size_t i{1u}; i <= type.first; ++i)
            std::cout << ", x^" << i;
          std::cout << "]\n";
        }

        if (EODenBasis) {
          if (denBasisEven) {
            std::cout << "Denominator basis = [1";
            for (size_t i{2u}; i <= type.second; i += 2u)
              std::cout << ", x^" << i;
            std::cout << "]\n";
          } else {
            std::cout << "Denominator basis = [x";
            for (size_t i{3u}; i <= type.second; i += 2u)
              std::cout << ", x^" << i;
            std::cout << "]\n";
          }
        } else if (useCustomDenBasis) {
          std::cout << "Denominator basis = " << denbasisstream.str()
                    << std::endl;
        } else {
          std::cout << "Denominator basis = [1";
          for (size_t i{1u}; i <= type.second; ++i)
            std::cout << ", x^" << i;
          std::cout << "]\n";
        }
      } else {
        std::cout << "Type: (" << type.first << ", " << type.second << ")\n";
      }
      std::cout << "Approximation error                              = "
                << delta << std::endl;
      std::cout << "Naive rounding approximation error               = "
                << errr_norm.second << std::endl;
      std::cout << "fpminimax approximation error                    = "
                << errfp_norm.second << std::endl;
      if (factorizationEnabled) {
        std::cout << "factorization approximation error                = "
                  << ferr_norm.second << std::endl;
        std::cout << "Naive rounding factorization approximation error = "
                  << rferr_norm.second << std::endl;
      }

      std::cout << "------------------------------------------------------\n";

      std::cout << "Writing results to the file " << outfilename << std::endl
                << std::endl;
    }

    std::ofstream coeffFile;
    coeffFile.open(outfilename);

    coeffFile << "Numerator = [|" << endl;
    for (std::size_t i{0u}; i < num.size() - 1; ++i)
      coeffFile << fpnum[i].toString(coeffDisplay) << "," << std::endl;
    coeffFile << fpnum[num.size() - 1].toString(coeffDisplay) << "|];"
              << std::endl;

    coeffFile << "Denominator = [|" << endl;
    for (std::size_t i{0u}; i < den.size() - 1; ++i)
      coeffFile << fpden[i].toString(coeffDisplay) << "," << std::endl;
    coeffFile << fpden[den.size() - 1].toString(coeffDisplay) << "|];"
              << std::endl;

    coeffFile.close();
  }
}

int main(int argc, char *argv[]) {
  QSexactStart();

  rminimax(argc, argv);

  QSexactClear();
}