#ifndef BIODYNAMO_MESH_DATA_H
#define BIODYNAMO_MESH_DATA_H

#include <unordered_map>
#include <vector>


#include "spatial_organization/bound.h"
#include "spatial_organization/node_index.h"

namespace bdm {
namespace spatial_organization {

using std::vector;
using std::array;
using std::unordered_map;

template <typename T> class MeshData;

/// Describes a hashed-octree node
template <typename T>
struct Node {
  uint_fast64_t location_code;
  T value;
  bool is_leaf;

  inline NodeIndex Index() {
    return NodeIndex(location_code);
  };
};


/// Implements hashed-octree.
/// Represent (a geometric object) as a set of finite elements
/// @tparam T - type of the data to be stored in the elements
template <typename T>
class MeshData {
 public:

  /// Empty constructor, initializes with bounds {0, 0, 0, 1, 1, 1},
  /// maximum depth 10,
  /// value of root with 0
  /// @tparam T - type of the data to be stored in the mesh
  MeshData();

  /// Constructor
  /// @tparam T  - type of the data to be stored in the mesh
  /// @param bnd - Bound of the mesh
  /// @param init_value - value of root node of the tree
  /// @param max_level - maximum possible depth of the tree. After reaching that
  /// point, nodes won't split. Maximum depth is 21
  MeshData(const Bound &bnd, T init_value, int max_level);

  void Put(const Point &p, T val);

  T At(const Point &p);

 protected:

  /// Returns node index on maximum level where point p is located
  /// @param p - point in space
  NodeIndex NodeIndexAt(const Point &p_) const;

  /// Returns node index on level l where point p is located
  /// @param p - point in space
  /// @param l - level of node index
  NodeIndex NodeIndexAt(const Point &p_, int l) const;

  /// Returns bound of node
  /// @param index - node index
  Bound GetNodeBound(NodeIndex index) const;

  /// Returns size of nodes on level l
  /// @param l - level of nodes
  Point NodeSizeAtLevel(uint l) const;

  /// Splits node to 8 equal subspaces
  /// @tparam T - type of mesh value
  /// @node - pointer to node which will split
  virtual void Split(Node<T> *node);

  /// Splits
  void Split(NodeIndex node_index, NodeIndex need_depth_node_index);


  /// Creates a node by the index where the node is located
  /// @param index - index where the node is located
  /// @return pointer to created node
  Node<T>* CreateNode(NodeIndex index);

  ///
  Node<T>* GetNode(NodeIndex);

  /// Destroy Node
  void DestroyNode(NodeIndex);

  ///
  Node<T>* FindLeaf(NodeIndex);

  /// Binary search leaf level
  Node<T>* BSLeaf(int min_l, int max_l, NodeIndex max_depth_index);

  /// Linear search leaf level
  Node<T>* LSLeaf(NodeIndex max_depth_index);


  using HashTable = unordered_map<uint_fast64_t, Node<T>>;

  /// Tree levels
  vector<HashTable> levels;

  /// Mesh bound
  Bound bound_;

  /// Maximum level of the mesh
  uint max_level_;

  /// Size of root node
  Point root_size_;

 private:


};

template <typename T>
MeshData<T>::MeshData()
        : MeshData({0, 0, 0, 1, 1, 1}, 0, 10) {

}

template <typename T>
MeshData<T>::MeshData(const Bound &bnd, T value, int max_level)
        : levels(max_level), bound_(bnd), max_level_(max_level) {

  /// Reserve memory for each level
  for (int i = 0; i < max_level; i++) {
    levels.push_back({});
    levels[i].reserve(1 << i);
  }

  /// Root Node initialize
  Node<T> n = {0b001, value};
  n.is_leaf = true;
  levels[0][0b001] = n;

  /// Calculate root node size
  root_size_ = bound_.near_right_top_point_
               + bound_.far_left_bottom_point_
                 * (-1);

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

  return NodeIndex(x, y, z, level);
}

template <typename T>
Point MeshData<T>::NodeSizeAtLevel(uint level) const {
  return root_size_ * (1.0 / (1l << level));
}

template <typename T>
Bound MeshData<T>::GetNodeBound(NodeIndex index) const {
  auto node_size = NodeSizeAtLevel(index.Level());
  uint x = index.X();
  uint y = index.Y();
  uint z = index.Z();

  Point local_position = {node_size.x_*x,
                          node_size.y_*y,
                          node_size.z_*z};

  auto far_left_bottom_point = bound_.far_left_bottom_point_ + local_position;
  return {far_left_bottom_point, far_left_bottom_point + node_size};
}

template <typename T>
void MeshData<T>::Split(Node<T> *node) {
  if (!node->is_leaf)
    return;

  node->is_leaf = false;
  for (uint i = 0; i < 8; i++) {
    auto child_index = (node->location_code << 3) | i;
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
    this->Split(node);
    auto level_offset =
            need_depth_node_index.Level() - node_index.Level() - 1;
    node_index = NodeIndex(need_depth_node_index.code >> (3*level_offset));
  }
}

template <typename T>
Node<T>* MeshData<T>::CreateNode(NodeIndex node_index) {
  auto &nodes = levels[node_index.Level()];
  nodes[node_index.code] = {node_index.code, 0, true};
  auto node = &nodes[node_index.code];
  return node;
}

template <typename T>
Node<T>* MeshData<T>::GetNode(NodeIndex node_index) {
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
//  return BSLeaf(0, max_level_, nodeIndex);
  return LSLeaf(nodeIndex);
}

//TODO(Sotnem13): Binary search have bug, needs to be fixed
template <typename T>
Node<T>* MeshData<T>::BSLeaf(int min_l, int max_l, NodeIndex max_depth_index) {
  auto middle_l = (max_l - min_l)/2;

  NodeIndex middle_index(max_depth_index.code >> 3*(max_level_-middle_l));

  auto node = GetNode(middle_index);

  if (node == nullptr) {
      return BSLeaf(min_l, middle_l, max_depth_index);
  }
  if (node->is_leaf)
    return node;
  else
    return BSLeaf(middle_l, max_l, max_depth_index);
}

template <typename T>
Node<T>* MeshData<T>::LSLeaf(NodeIndex max_depth_index) {
  auto node = GetNode(max_depth_index);
  int i = 0;
  while(node == nullptr) {
    i++;
    node = GetNode(NodeIndex(max_depth_index.code >> 3 * i));
  }
  return node;
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
