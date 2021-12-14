#pragma once

namespace graph {

struct EdgeImpl {
  long weight = 0;
  int mark = 0;
  bool isTried = false;

  EdgeImpl() = default;

  EdgeImpl(int weight) : weight(weight) {}
};

}
