//
// Created by Silviu-Ioan Filip on 08/06/2018.
//

#ifndef RMINIMAX_SHUNTINGYARD
#define RMINIMAX_SHUNTINGYARD
#include "tokenutils.h"
#include <functional>
#include <map>
#include <mpreal.h>
#include <string>
#include <utility>
#include <vector>

using namespace std;

class shuntingyard {
public:
  shuntingyard();
  mpfr::mpreal evaluate(std::vector<std::string> &tokens, mpfr::mpreal x) const;

private:
  std::map<std::string, std::pair<int, int>> operatorMap;
  std::map<std::string, std::function<mpfr::mpreal(mpfr::mpreal)>> funcMap;
  std::map<std::string, std::function<mpfr::mpreal(mpfr::mpreal, mpfr::mpreal)>>
      opMap;

  bool isOperator(std::string const &token) const;
  bool isFunction(std::string const &token) const;
  bool isOperand(std::string const &token) const;
  bool isAssociative(std::string const &token, AssociativityType type) const;
  int cmpPrecedence(std::string const &token1, std::string const &token2) const;
};

#endif
