#ifndef BIODYNAMO_NODE_INDEX_H
#define BIODYNAMO_NODE_INDEX_H

#include <tuple>

namespace bdm {
namespace spatial_organization {

using std::tuple;

struct NodeIndex {

  explicit NodeIndex(tuple<uint, uint, uint>);

  NodeIndex(uint64_t code);

  uint Level();
  NodeIndex Parent();

  tuple<uint, uint, uint> Position();

  union {
    uint64_t code;
    struct {
      uint64_t x : 1;
      uint64_t y : 1;
      uint64_t z : 1;
    };
  };
};


NodeIndex::NodeIndex(uint64_t code_) : code(code_) {}


NodeIndex NodeIndex::Parent() {
  return NodeIndex(this->code >> 3);
}
uint NodeIndex::Level() {
  uint level = 0;
  while (code >> 3*(level+1)) {
    level += 1;
  }
  return level;
}

tuple<uint, uint, uint> NodeIndex::Position() {
  uint x_ = 0, y_ = 0, z_ = 0;
  for (int i = Level() - 1; i >= 0; i--) {
    auto code_at_level = code >> 3 * i;
    x_ = (x_ << 1) |  code_at_level & 1;
    y_ = (y_ << 1) | (code_at_level >> 1) & 1;
    z_ = (z_ << 1) | (code_at_level >> 2) & 1;
  }
  return std::make_tuple(x_, y_, z_);
};


}  // namespace spatial_organization
}  // namespace bdm

#endif //BIODYNAMO_NODE_INDEX_H_H
