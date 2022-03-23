#include "shuntingyard.h"
#include <iostream>
#include <queue>
#include <stack>

shuntingyard::shuntingyard() {
  operatorMap["+"] = std::pair<int, AssociativityType>(0, LEFT_ASSOC);
  operatorMap["-"] = std::pair<int, AssociativityType>(0, LEFT_ASSOC);
  operatorMap["*"] = std::pair<int, AssociativityType>(4, LEFT_ASSOC);
  operatorMap["/"] = std::pair<int, AssociativityType>(4, LEFT_ASSOC);
  operatorMap["^"] = std::pair<int, AssociativityType>(6, RIGHT_ASSOC);

  funcMap["sin"] = [](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::sin(x); };
  funcMap["cos"] = [](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::cos(x); };
  funcMap["tan"] = [](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::tan(x); };
  funcMap["asin"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::asin(x);
  };
  funcMap["acos"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::acos(x);
  };
  funcMap["atan"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::atan(x);
  };
  funcMap["sinh"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::sinh(x);
  };
  funcMap["cosh"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::cosh(x);
  };
  funcMap["tanh"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::tanh(x);
  };
  funcMap["asinh"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::asinh(x);
  };
  funcMap["acosh"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::acosh(x);
  };
  funcMap["atanh"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::atanh(x);
  };
  funcMap["log"] = [](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::log(x); };
  funcMap["log2"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::log2(x);
  };
  funcMap["log10"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::log10(x);
  };
  funcMap["exp"] = [](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::exp(x); };
  funcMap["exp2"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::exp2(x);
  };
  funcMap["exp10"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::exp10(x);
  };
  funcMap["sqrt"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::sqrt(x);
  };
  funcMap["abs"] = [](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::abs(x); };
  funcMap["gamma"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::gamma(x);
  };
  funcMap["lngamma"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::lngamma(x);
  };
  funcMap["zeta"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::zeta(x);
  };
  funcMap["erf"] = [](mpfr::mpreal x) -> mpfr::mpreal { return mpfr::erf(x); };
  funcMap["erfc"] = [](mpfr::mpreal x) -> mpfr::mpreal {
    return mpfr::erfc(x);
  };

  opMap["+"] = [](mpfr::mpreal x, mpfr::mpreal y) -> mpfr::mpreal {
    return x + y;
  };
  opMap["-"] = [](mpfr::mpreal x, mpfr::mpreal y) -> mpfr::mpreal {
    return x - y;
  };
  opMap["*"] = [](mpfr::mpreal x, mpfr::mpreal y) -> mpfr::mpreal {
    return x * y;
  };
  opMap["/"] = [](mpfr::mpreal x, mpfr::mpreal y) -> mpfr::mpreal {
    return x / y;
  };
  opMap["^"] = [](mpfr::mpreal x, mpfr::mpreal y) -> mpfr::mpreal {
    return mpfr::pow(x, y);
  };
}

mpfr::mpreal shuntingyard::evaluate(std::vector<std::string> &tokens,
                                    mpfr::mpreal x) const {

  std::stack<std::pair<std::string, int>> operatorStack;
  std::stack<mpfr::mpreal> operandStack;

  for (std::size_t i{0u}; i < tokens.size(); ++i) {
    tokens[i] = trim(tokens[i]);

    if (isOperator(tokens[i])) {
      if (operandStack.empty()) {
        operatorStack.push(make_pair(tokens[i], 1)); // push a unary operator
      } else {
        if (i > 0 && (isOperator(tokens[i - 1]) || tokens[i - 1] == "(")) {
          // prefix unary operators
          while (!operatorStack.empty() &&
                 isOperator(operatorStack.top().first) &&
                 operatorStack.top().second == 1) {
            if (operatorStack.top().first == "+")
              operatorStack.pop();
            else if (operatorStack.top().first == "-") {
              operatorStack.pop();
              mpfr::mpreal result = -operandStack.top();
              operandStack.pop();
              operandStack.push(result);
            } else {
              std::cerr << "Badly formatted function string: processing error!";
              exit(EXIT_FAILURE);
            }
          }
          operatorStack.push(make_pair(tokens[i], 1));
        } else if (i > 0 &&
                   (isOperand(tokens[i - 1]) || tokens[i - 1] == ")")) {
          // binary operators or postfix unary operators
          while (!operatorStack.empty() &&
                 isOperator(operatorStack.top().first)) {
            if ((isAssociative(tokens[i], LEFT_ASSOC) &&
                 cmpPrecedence(tokens[i], operatorStack.top().first) <= 0) ||
                (cmpPrecedence(tokens[i], operatorStack.top().first) < 0)) {
              if (operatorStack.top().second == 2) {
                mpfr::mpreal op2 = operandStack.top();
                operandStack.pop();
                mpfr::mpreal op1 = operandStack.top();
                operandStack.pop();
                auto f = opMap.find(operatorStack.top().first);
                operandStack.push(f->second(op1, op2));
                operatorStack.pop();
              } else {
                if (operatorStack.top().first == "+")
                  operatorStack.pop();
                else if (operatorStack.top().first == "-") {
                  operatorStack.pop();
                  mpfr::mpreal result = -operandStack.top();
                  operandStack.pop();
                  operandStack.push(result);
                } else {
                  std::cerr
                      << "Badly formatted function string: processing error!";
                  exit(EXIT_FAILURE);
                }
              }
              continue;
            }
            break;
          }
          operatorStack.push(make_pair(tokens[i], 2));
        }
      }
    } else if (isFunction(tokens[i]) || tokens[i] == "(") {
      operatorStack.push(make_pair(tokens[i], 3));
    } else if (tokens[i] == ")") {
      while (operatorStack.top().first != "(") {
        if (isFunction(operatorStack.top().first)) {
          mpfr::mpreal op = operandStack.top();
          operandStack.pop();
          auto f = funcMap.find(operatorStack.top().first);
          operandStack.push(f->second(op));
        } else if (operatorStack.top().second == 1) {
          if (operatorStack.top().first == "+") {
          } else if (operatorStack.top().first == "-") {
            mpfr::mpreal result = -operandStack.top();
            operandStack.pop();
            operandStack.push(result);
          } else {
            std::cerr << "Badly formatted function string: processing error!";
            exit(EXIT_FAILURE);
          }
        } else {
          mpfr::mpreal op2 = operandStack.top();
          operandStack.pop();
          mpfr::mpreal op1 = operandStack.top();
          operandStack.pop();
          auto f = opMap.find(operatorStack.top().first);
          operandStack.push(f->second(op1, op2));
        }
        operatorStack.pop();
      }
      operatorStack.pop();
      if (!operatorStack.empty() && isFunction(operatorStack.top().first)) {
        mpfr::mpreal op = operandStack.top();
        operandStack.pop();
        auto f = funcMap.find(operatorStack.top().first);
        operandStack.push(f->second(op));
        operatorStack.pop();
      }
    } else {
      // numeric token
      if (tokens[i] == "pi") // the constant Pi
        operandStack.push(mpfr::const_pi());
      else if (tokens[i] == "x") // the unknown argument for the expression
        operandStack.push(x);
      else // assume a numeric constant otherwise,
           // parsable by the mprf/mpreal constructors
      {
        operandStack.push(mpfr::mpreal(tokens[i]));
      }
    }
  }

  while (!operatorStack.empty()) {
    if (isFunction(operatorStack.top().first)) {
      mpfr::mpreal op = operandStack.top();
      operandStack.pop();
      auto f = funcMap.find(operatorStack.top().first);
      operandStack.push(f->second(op));
    } else if (operatorStack.top().second == 1) {
      if (operatorStack.top().first == "+") {
      } else if (operatorStack.top().first == "-") {
        mpfr::mpreal result = -operandStack.top();
        operandStack.pop();
        operandStack.push(result);
      } else {
        std::cerr << "Badly formatted input: processing error!";
        exit(EXIT_FAILURE);
      }
    } else {
      mpfr::mpreal op2 = operandStack.top();
      operandStack.pop();
      mpfr::mpreal op1 = operandStack.top();
      operandStack.pop();
      auto f = opMap.find(operatorStack.top().first);
      operandStack.push(f->second(op1, op2));
    }
    operatorStack.pop();
  }

  if (operandStack.size() != 1) {
    std::cerr << "Badly formatted function string: processing error!\n";
    exit(EXIT_FAILURE);
  }

  mpfr::mpreal result = operandStack.top();
  return result;
}

bool shuntingyard::isOperator(std::string const &token) const {
  return token == "+" || token == "-" || token == "*" || token == "/" ||
         token == "^";
}

bool shuntingyard::isFunction(std::string const &token) const {
  return token == "sin" || token == "cos" || token == "tan" || token == "log" ||
         token == "log2" || token == "log10" || token == "asin" ||
         token == "atan" || token == "sqrt" || token == "exp" ||
         token == "exp2" || token == "exp10" || token == "acos" ||
         token == "abs" || token == "gamma" || token == "zeta" ||
         token == "lngamma" || token == "erf" || token == "erfc" ||
         token == "sinh" || token == "cosh" || token == "tanh" ||
         token == "asinh" || token == "acosh" || token == "atanh";
}

bool shuntingyard::isOperand(std::string const &token) const {
  return !isOperator(token) && !isFunction(token) && token != "(" &&
         token != ")";
}

bool shuntingyard::isAssociative(std::string const &token,
                                 AssociativityType type) const {
  std::pair<int, int> p = operatorMap.find(token)->second;
  if (p.second == type) {
    return true;
  }
  return false;
}

int shuntingyard::cmpPrecedence(std::string const &token1,
                                std::string const &token2) const {
  std::pair<int, int> p1 = operatorMap.find(token1)->second;
  std::pair<int, int> p2 = operatorMap.find(token2)->second;

  return p1.first - p2.first;
}
