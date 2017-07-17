#ifndef BIODYNAMO_ADAPTIVE_MESH_H
#define BIODYNAMO_ADAPTIVE_MESH_H

#include <vector>
#include <stack>
#include <unordered_map>
#include <functional>

#include "spatial_organization/node_index.h"
#include "spatial_organization/mesh_data.h"

namespace bdm {
namespace spatial_organization {

using std::vector;
using std::unordered_map;
using std::function;


template <typename T>
class AdaptiveMesh : public MeshData<T> {
 public:

  /// Empty constructor, initializes with bounds {0, 0, 0, 1, 1, 1},
  /// maximum depth 10,
  /// value of root with 0
  /// @tparam T - type of the data to be stored in the mesh
  AdaptiveMesh();

  /// Constructor
  /// @tparam T  - type of the data to be stored in the mesh
  /// @param bnd - Bound of the mesh
  /// @param init_value - value of root node of the tree
  /// @param max_level - maximum possible depth of the tree. After reaching that
  /// point, nodes won't split. Maximum depth is 21
  AdaptiveMesh(const Bound &bnd, T empty_val, int max_level);

  /// Destructor
  ~AdaptiveMesh();

  /// Something like access policy
  class Element;

  void Refine(function<void(Element&)> refine_func,
          bool mode = false);

//private:
protected:
  ///
  void Split(Node<T> *node);
  ///
  void Coarse(Node<T> *node);



};


/// Octree node wrapper
template <typename T>
class AdaptiveMesh<T>::Element {
public:

  /// Is value changed
  inline bool Changed();

  ///
  inline T Value();

  ///
  inline void SetValue(T new_value);

  /// Return mesh element bound
  inline Bound GetBound();

  /// Return vector of neighbor value and delta position
  vector<pair<T, Point>> Neighbors();

private:
  Element(Node<T>* node, vector<NodeIndex>* neighbors, AdaptiveMesh<T>* parent)
          : node_(node), neighbors_(neighbors),
            parent_(parent), value(node->value) {}

  Node<T> *node_;
  vector<NodeIndex> * neighbors_;
  AdaptiveMesh<T>* parent_;
  T value;
  friend AdaptiveMesh;
};

template <typename T>
AdaptiveMesh<T>::AdaptiveMesh()
        : AdaptiveMesh(Bound{0, 0, 0, 1, 1, 1}, 0, 10) {}

template <typename T>
AdaptiveMesh<T>::AdaptiveMesh(const Bound &bnd, T empty_val, int max_depth)
        : MeshData<T>(bnd, empty_val, max_depth) {
}

template <typename T>
AdaptiveMesh<T>::~AdaptiveMesh() {}

template <typename T>
void AdaptiveMesh<T>::Split(Node<T> *node) {
  MeshData<T>::Split(node);
  auto node_index = node->Index();
  auto neighbors =
          NodeIndex::GetAdjacentIndex<true>(node_index, node_index.Level());

  for (auto neighbor_index : neighbors) {
    if (neighbor_index.Parent().code != node->Index().Parent().code) {
      auto neighbor = this->FindLeaf(neighbor_index);
      if (neighbor->Index().Level() < node->Index().Level()) {
        Split(neighbor);
      }
    }
  }
}
template <typename T>
void AdaptiveMesh<T>::Coarse(Node<T>* node) {
  if (node == nullptr || node->is_leaf) return;

  auto node_index = node->Index();
  auto first_child = this->GetNode(NodeIndex(node_index.code << 3));
  auto need_coarse = true;

  for (int i = 0; i < 8 && need_coarse; i++) {
    auto child_index = NodeIndex(node->location_code << 3 | i);
    auto child = this->GetNode(child_index);
    if (child && !child->is_leaf) {
      Coarse(child);
    }
    need_coarse = child->is_leaf &&
            first_child->value == child->value;
  }

  if (need_coarse) {
    auto neighbors = NodeIndex::GetAdjacentIndex<true>(node_index,
                                                       node_index.Level());
    for (int i = 0; i < neighbors.size() && need_coarse; i++) {
      auto neighbor = this->GetNode(neighbors[i]);
      if (neighbor == nullptr) {
        neighbor = this->GetNode(node_index.Parent());
      }
      need_coarse = neighbor && neighbor->is_leaf &&
              neighbor->value == first_child->value;
    }
    if (need_coarse) {
      node->value = first_child->value;
      node->is_leaf = true;
      for (int i = 0; i <8; i++) {
        auto child_index = NodeIndex(node_index.code << 3 | i);
        this->DestroyNode(child_index);
      }
    }
  }
}


template <typename T>
void AdaptiveMesh<T>::Refine(function<void(Element&)> refine_func,
                             bool mode) {

  auto &leafs = this->levels[this->max_level_];

  auto new_values = vector<pair<Node<T>*, T>>();
  std::stack<NodeIndex> stack;
  new_values.reserve(leafs.size());

  auto neighbor_func = [](bool m){
    if (m)
      return NodeIndex::GetAdjacentIndex<true>;
    else
      return NodeIndex::GetAdjacentIndex<false>;
  }(mode);

  for (auto &it : leafs) {
    auto &node = it.second;
    auto neighbors = neighbor_func(node.Index(), this->max_level_);
    Element elem(&node, &neighbors, this);
    for (auto &neighbor_index : neighbors) {
      auto neighbor = this->LSLeaf(neighbor_index);
      if (neighbor->Index().Level() != this->max_level_) {
        stack.push(neighbor_index);
      }
    }

    refine_func(elem);

    if (elem.Changed()) {
      new_values.push_back({&node, elem.Value()});
    }
  }

  new_values.reserve(stack.size());

  while (!stack.empty()) {
    auto node_index = stack.top();
    Split(this->GetNode(node_index.Parent()));

    auto node = this->GetNode(node_index);
    auto neighbors = neighbor_func(node->Index(), this->max_level_);
    auto elem = Element(node, &neighbors, this);

    refine_func(elem);

    if (elem.Changed()) {
      new_values.push_back({node, elem.Value()});
    }
    stack.pop();
  }

  for (auto &n : new_values) {
    auto node = n.first;
    node->value = n.second;
  }

  for (auto &n : new_values) {
    auto node = n.first;
    Coarse(this->GetNode(node->Index().Parent()));
  }
}


template <typename T>
inline bool AdaptiveMesh<T>::Element::Changed() {
  return Value() != node_->value;
}

template <typename T>
inline T AdaptiveMesh<T>::Element::Value() {
  return value;
}

template <typename T>
inline void AdaptiveMesh<T>::Element::SetValue(T new_value) {
  value = new_value;
}

template <typename T>
inline Bound AdaptiveMesh<T>::Element::GetBound() {
  return parent_->GetNodeBound(node_->Index());
}

template <typename T>
vector<pair<T, Point>> AdaptiveMesh<T>::Element::Neighbors() {
  vector<pair<T, Point>> result;
  result.reserve(neighbors_->size());
  auto node_index = node_->Index();
  auto node_size = parent_->NodeSizeAtLevel(node_index.Level());

  uint x = node_index.X();
  uint y = node_index.Y();
  uint z = node_index.Z();

  for (auto neighbor_index : *neighbors_) {

    uint x_ = neighbor_index.X();
    uint y_ = neighbor_index.Y();
    uint z_ = neighbor_index.Z();

    auto neighbor_value = parent_->FindLeaf(neighbor_index)->value;
    auto delta = Point(x - x_, y - y_, z - z_);
    delta.x_*=node_size.x_;
    delta.y_*=node_size.y_;
    delta.z_*=node_size.z_;
    result.push_back(make_pair(neighbor_value, delta));
  }
  return result;
}


}  // namespace spatial_organization
}  // namespace bdm



#endif //BIODYNAMO_ADAPTIVE_MESH_H
