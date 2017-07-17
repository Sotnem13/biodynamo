#ifndef BIODYNAMO_ADAPTIVE_MESH_H
#define BIODYNAMO_ADAPTIVE_MESH_H

#include <iostream>
#include <cstdint>
#include <vector>
#include <stack>
#include <array>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <bitset>
#include <functional>


#include "spatial_organization/octree_node.h"

#include "spatial_organization/mesh_tree_node.h"
#include "spatial_organization/voxel.h"
#include "spatial_organization/node_index.h"
#include "spatial_organization/mesh_data.h"

namespace bdm {
namespace spatial_organization {

using namespace std;
using std::vector;
using std::array;
using std::unordered_map;
using std::function;


template <typename T>
class AdaptiveMesh : public MeshData<T> {
 public:

  class Element;

  AdaptiveMesh();
  AdaptiveMesh(const Bound &bnd, T empty_val, int max_level);

  ~AdaptiveMesh();

  void Refine(function<void(Element&)> refine_func,
          bool mode = false);


  class Element {
   public:
    bool Changed() {
      return Value() != node_->value;
    }
    T Value() {
      return value;
    }
    void SetValue(T new_value) {
      value = new_value;
    }

    Point GetPosition() {
      return parent_->GetNodeBound(node_->Index());
    }

    vector<pair<T, Point>> Neighbors() {
      vector<pair<T, Point>> result;
      result.reserve(neighbors_->size());
      auto node_index = node_->Index();
      auto node_size = parent_->NodeSizeAtLevel(node_index.Level());
      int64_t x, y, z;
      tie(x, y, z) = node_index.Position();

      for (auto neighbor_index : *neighbors_) {
//        auto neighbor_value = parent_->GetNode(neighbor_index);
//        if (neighbor == nullptr) {
//          parent_
//          parent_->GetNode(neighbor_index.Parent());
//        }
        int64_t x_, y_, z_;
        tie(x_, y_, z_) = neighbor_index.Position();
        auto neighbor_value = parent_->At(Point(x_, y_, z_));
        auto delta = Point(x - x_, y - y_, z - z_);
        delta.x_*=node_size.x_;
        delta.y_*=node_size.y_;
        delta.z_*=node_size.z_;
        result.push_back(make_pair(neighbor_value, delta));
      }
      return result;
    }

  private:
    Element(Node<T>* node, vector<NodeIndex>* neighbors, AdaptiveMesh<T>* parent) :
            node_(node), neighbors_(neighbors), parent_(parent), value(node->value){
    }

    Node<T> *node_;
    vector<NodeIndex> * neighbors_;
    AdaptiveMesh<T>* parent_;
    T value;
    friend AdaptiveMesh;
  };

//private:
protected:
  void Split(Node<T> *node);
  void Coarse(Node<T> *node);

//  static vector<NodeIndex> GetAdjacentIndex(NodeIndex index, int at_level);
//  template <typename>
//  friend AccessPolicy;

//  AccessPolicy<T>

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

  auto neighbors = NodeIndex::GetAdjacentIndex<true>(node->Index(), node->Index().Level());

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

  auto first_child = this->GetNode(NodeIndex(node->location_code << 3));
  auto need_coarse = true;

  for (int i = 0; i < 8 && need_coarse; i++) {
    auto child_index = NodeIndex(node->location_code << 3 | i);
    auto child = this->GetNode(child_index);
    if(child && !child->is_leaf) {
      Coarse(child);
    }
    need_coarse = child->is_leaf && first_child->value == child->value;
  }

  if (need_coarse) {
    auto neighbors = NodeIndex::GetAdjacentIndex<true>(node->Index(), node->Index().Level());
    for (int i = 0; i < neighbors.size() && need_coarse; i++) {
      auto neighbor = this->GetNode(neighbors[i]);
      if (!neighbor) {
        neighbor = this->GetNode(neighbor->Index().Parent());
      }
      need_coarse = neighbor && neighbor->is_leaf &&
              neighbor->value == first_child->value;
    }
    if (need_coarse) {
      node->value = first_child->value;
      node->is_leaf = true;
      for (int i = 0; i <8; i++) {
        auto child_index = NodeIndex(node->location_code << 3 | i);
        this->DestroyNode(child_index);
      }
    }
  }
}


template <typename T>
void AdaptiveMesh<T>::Refine(
        function<void(Element&)> refine_func,
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
      auto neighbor_level = neighbor_index.Level();
//      auto neighbor = this->BSLeaf(neighbor_level - 1, neighbor_level, neighbor_index);
      auto neighbor = this->LSLeaf(neighbor_index);
      if (neighbor->Index().Level() != this->max_level_ + 1) {
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
    cout << node->location_code  << n.second << endl;
    node->value = n.second;
  }

  for (auto &n : new_values) {
    auto node = n.first;
    Coarse(this->GetNode(node->Index().Parent()));
  }
}




}  // namespace spatial_organization
}  // namespace bdm



#endif //BIODYNAMO_ADAPTIVE_MESH_H
