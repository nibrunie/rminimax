#include "tokenutils.h"
#include <algorithm>
#include <cctype>
#include <functional>
#include <locale>
#include <utility>

std::string &ltrim(std::string &s) {
  const char *t = " \t\n\r\f\v";
  s.erase(0, s.find_first_not_of(t));
  return s;
}

std::string &rtrim(std::string &s) {
  const char *t = " \t\n\r\f\v";
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

std::string &trim(std::string &s) { return ltrim(rtrim(s)); }
