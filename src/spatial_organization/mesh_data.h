#ifndef BIODYNAMO_MESH_DATA_H
#define BIODYNAMO_MESH_DATA_H

#include <unordered_map>
#include <vector>



namespace bdm {
namespace spatial_organization {

using std::vector;
using std::array;
using std::unordered_map;


template <typename T>
//struct MeshData<T>::Node {
struct Node {
  uint_fast64_t location_code;
  T value;
  bool is_leaf;

  inline NodeIndex Index() { return NodeIndex(location_code); }
};

template <typename T>
class MeshData {
 public:


  MeshData(const Bound &bnd = {0, 0, 0, 1, 1, 1},
           T empty_value = 0, int max_level = 10);

  void Put(const Point &p, T obj);
  T At(const Point &p);

 protected:

  NodeIndex NodeIndexAt(const Point &p_) const;
  NodeIndex NodeIndexAt(const Point &p_, int level) const;

  Bound GetNodeBound(NodeIndex index) const;
  Point NodeSizeAtLevel(uint level) const;

  void Split(Node<T> *node);
  void Split(NodeIndex node_index, NodeIndex need_depth_node_index);

  Node<T>* CreateNode(NodeIndex);
  Node<T>* GetNode(NodeIndex);
  Node<T>* FindLeaf(NodeIndex);
  Node<T>* BSLeaf(NodeIndex min, NodeIndex max);

  void DestroyNode(NodeIndex);
  static constexpr auto root_index = 0b0001;

  using HashTable = unordered_map<uint_fast64_t, Node<T>>;
  vector<HashTable> levels;

  Bound bound_;
  uint max_level_;
  T empty_val_;

 private:


};

template <typename T>
MeshData<T>::MeshData(const Bound &bnd, T value, int max_level) :
        bound_(bnd), empty_val_(value),
        max_level_(max_level), levels(max_level) {

  /// Root Node initialize
  Node<T> n = {root_index, value};
  n.is_leaf = true;

  for (int i = 0; i < max_level; i++) {
    levels.push_back({});
    levels[i].reserve(1 << i);
  }

  levels[0][root_index] = n;

}

template <typename T>
void MeshData<T>::Put(const Point &p, T obj) {
  auto max_depth_node_index = NodeIndexAt(p);
  auto leaf = FindLeaf(max_depth_node_index);

  if (leaf->value != obj) {
    if (leaf->Index().Level() < max_level_) {
      Split(leaf->Index(), max_depth_node_index);
      leaf = GetNode(max_depth_node_index);
    }
    leaf->value = obj;
  }
}

template <typename T>
T MeshData<T>::At(const Point &p) {
  auto max_depth_node_index = NodeIndexAt(p);
  return FindLeaf(max_depth_node_index)->value;
}

template <typename T>
NodeIndex MeshData<T>::NodeIndexAt(const Point &p_) const {
  return NodeIndexAt(p_, max_level_);
}

template <typename T>
NodeIndex MeshData<T>::NodeIndexAt(const Point &p_, int level) const {
  auto node_size = NodeSizeAtLevel(level);
  auto p = p_ + bound_.far_left_bottom_point_ * (-1);
  uint x = p.x_ / node_size.x_;
  uint y = p.y_ / node_size.y_;
  uint z = p.z_ / node_size.z_;
  auto position = std::make_tuple(x, y, z);
  return NodeIndex(position);
}

template <typename T>
Point MeshData<T>::NodeSizeAtLevel(uint level) const {
  auto size =
          bound_.near_right_top_point_
          + bound_.far_left_bottom_point_
            * (-1);
  return size * (1.0 / (1 << level));
}

template <typename T>
Bound MeshData<T>::GetNodeBound(NodeIndex index) const {
  auto node_size = NodeSizeAtLevel(index.Level());
  uint x, y, z;
  std::tie(x, y, z) = index.Position();

  Point local_position = {node_size.x_*x,
                          node_size.y_*y,
                          node_size.z_*z};

  auto center = bound_.far_left_bottom_point_ + local_position;
  return {center + node_size*(-1), center + node_size};
}

template <typename T>
void MeshData<T>::Split(Node<T> *node) {
  if (!node->is_leaf)
    return;

  node->is_leaf = false;
  for (uint i = 0; i < 8; i++) {
    auto child_index = node->location_code << 3 | i;
    auto child = CreateNode(NodeIndex(child_index));
    child->location_code = child_index;
    child->is_leaf = true;
    child->value = node->value;
  }
}

template <typename T>
void MeshData<T>::Split(NodeIndex node_index, NodeIndex need_depth_node_index) {
  while (node_index.Level() < need_depth_node_index.Level()) {
    auto node = GetNode(node_index);
    Split(node);
    auto level_offset =
            need_depth_node_index.Level() - node_index.Level() - 1;
    node_index = NodeIndex(need_depth_node_index.code >> 3 * (level_offset));
  }
}

template <typename T>
Node<T>* MeshData<T>::CreateNode(NodeIndex node_index) {
  auto &nodes = levels[node_index.Level()];
  nodes[node_index.code] = {node_index.code, empty_val_, true};
  auto node = &nodes[node_index.code];
  return node;
}

template <typename T>
Node<T>* MeshData<T>::GetNode(NodeIndex node_index) {
//  BSLeaf
  auto &nodes = levels[node_index.Level()];
  auto node_iter = nodes.find(node_index.code);
  if (node_iter == nodes.end()) {
    return nullptr;
  } else {
    return &node_iter->second;
  }
}

template <typename T>
Node<T>* MeshData<T>::FindLeaf(NodeIndex nodeIndex) {
  return BSLeaf(NodeIndex(root_index), nodeIndex);
}
template <typename T>
Node<T>* MeshData<T>::BSLeaf(NodeIndex min, NodeIndex max) {
  auto offset = (max.Level() - min.Level()+1)/2;

  NodeIndex middle(max.code >> 3*(offset));
//  NodeIndex middle = middle_child.Parent();

  auto node = GetNode(middle);

  if (node == nullptr) {
//    if (min.code == max.code)
//      return nullptr;
//    else
      return BSLeaf(min, middle);
  }
  if (node->is_leaf)
    return node;
  else
    return BSLeaf(max.code >> 3*(offset-1), max);
}

template <typename T>
void MeshData<T>::DestroyNode(NodeIndex index) {
  auto &nodes = levels[index.Level()];
  auto node_iter = nodes.find(index.code);
  if (node_iter != nodes.end()) {
    nodes.erase(node_iter);
  }
}

}  // namespace spatial_organization
}  // namespace bdm


#endif //BIODYNAMO_MESH_DATA_H
