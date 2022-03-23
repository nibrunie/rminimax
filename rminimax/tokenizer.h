#ifndef RMINIMAX_TOKENIZER
#define RMINIMAX_TOKENIZER
#include <string>
#include <vector>

class tokenizer {
public:
  tokenizer();
  explicit tokenizer(std::string mData);
  ~tokenizer();
  std::vector<std::string> &getTokens();
  void setData(std::string newData);

private:
  void tokenize();

  std::string data;
  std::vector<std::string> tokens;
};

#endif
