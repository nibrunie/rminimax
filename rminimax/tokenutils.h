#ifndef RMINIMAX_TOKENUTILS_H
#define RMINIMAX_TOKENUTILS_H
#include <string>

enum TokenType { OPERATOR, OPERAND, UNKNOWN, FUNCTION, PAREN };

enum AssociativityType { LEFT_ASSOC, RIGHT_ASSOC };

// trim from the start
std::string &ltrim(std::string &s);

// trim from the end
std::string &rtrim(std::string &s);

// trim from both ends
std::string &trim(std::string &s);

#endif
