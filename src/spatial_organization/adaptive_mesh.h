#ifndef BIODYNAMO_ADAPTIVE_MESH_H
#define BIODYNAMO_ADAPTIVE_MESH_H

#include <vector>
#include <stack>
#include <unordered_map>
#include <functional>
#include <bitset>

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
  class Face;

  void Refine(function<void(Element&)> refine_func,
          bool mode = false);

  void prepareChildList();

  void prepareFaces();

  void PrepareFacesXY(const Node<T>& node);
  void PrepareFacesXZ(const Node<T>& node);
  void PrepareFacesYZ(const Node<T>& node);

//  void Refine(const function<void(Face&)> &refine_func) {
//
////    prepareChildList();
////    prepareFaces();
//
//    for (auto &face : faces) {
//      refine_func(face);
//    }
//
//    changes.clear();
//    for (const auto &face : faces) {
//      changes[face.front] += face.flux;
//      changes[face.back] -= face.flux;
//    }
//
////    coarse.clear();
//    for (const auto &change : changes) {
//      const auto node_index = NodeIndex(change.first);
//      const auto dvalue = change.second;
//      auto node = FindOrCreateNode(node_index);
//      node->value += dvalue;
//      coarse[node_index.Parent().code]++;
//
//    }
//
//    Coarse();
//
//  }

//  void Coarse() {
//
//    for (auto c : coarse) {
//      if (c.second == 8) {
//        auto node = GetNode(NodeIndex(c.first));
//        Coarse(node);
//      }
//    }
//
//  }

  void Refine(const function<void(Face&)> &refine_func) {

//    prepareChildList();
    prepareFaces();

    changes.clear();
    // #pragma omp parallel for
    for (size_t i = 0; i < faces.size(); i++) {
      const auto &face = faces[i];
      const auto n1 = FindOrCreateNode(face.first);
      const auto n2 = FindOrCreateNode(face.second);

//      node_faces.push
      Face f(n1, n2);
      refine_func(f);
      changes[face.first.code] += f.flux;
      changes[face.second.code] -= f.flux;
    }

//    changes.clear();
////
//    for (const auto &face : faces) {
//      changes[face.front] += face.flux;
//      changes[face.back] -= face.flux;
//    }

//    #pragma omp parallel for
//    for (size_t i = 0; i < changes.bucket_count(); i++) {
//
//      for (auto change = changes.begin(i); change != changes.end(i); ++change) {

//      }
      for (auto &change : changes) {
        const auto node_index = NodeIndex(change.first);
        const auto dvalue = change.second;
        auto node = this->GetNode(node_index);
        node->value += dvalue;
      }

//    }

  }



  Node<T>* FindOrCreateNode(NodeIndex index) {
    auto node = this->GetNode(index);
    if (node == nullptr) {
      auto parent = this->GetNode(index.Parent());
      this->Split(parent);
      node = this->GetNode(index);



//      auto leaf = FindLeaf(index);
//      Split(leaf->Index(), index);
//      node = GetNode(index);
    }
    return node;
  }

//  void Put(const Point &p, T obj) {
//    auto node_index = NodeIndexAt(p);
//    auto leaf = FindLeaf(node_index);
//
//    if (leaf->value != obj) {
//      if (leaf->Index().Level() < max_level_) {
//        this->Split(leaf->Index(), node_index);
//        leaf = GetNode(node_index);
//      }
//      leaf->value = obj;
////      auto parent  = GetNode(node_index.Parent());
//
//
////      leaf->Faces();
////      NodeFace face(leaf->location_code, {0,0,0});
////      NodeFace face(leaf->location_code, {0,0,1});
////      NodeFace face(leaf->location_code, {0,1,0});
////      NodeFace face(leaf->location_code, {1,0,0});
//
////      faces.push_back(
//
//    }
//  }

//private:
protected:
  ///
  void Split(Node<T> *node);
  ///
  void Coarse(Node<T> *node);

  vector<pair<NodeIndex, NodeIndex>> faces;
  vector<Face> node_faces;
  vector<uint_fast64_t> childs_for_refine;
  unordered_map<uint_fast64_t, T> changes;

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
class AdaptiveMesh<T>::Face {
 public:

  inline T FrontValue() {
    return fvalue;
  };

  inline T BackValue() {
    return bvalue;
  };

  inline void SetFlux(T dvalue) {
    flux = dvalue;
  };

  Face(const Node<T> *n1, const Node<T> *n2) : flux(0), fvalue(n1->value), bvalue(n2->value) {
    front = n1->location_code;
    back  = n2->location_code;
  }


  Face() : flux(0), fvalue(0), bvalue(0), front(0), back(0) {

  }
 private:

  T flux;
  T fvalue;
  T bvalue;
  uint_fast64_t front;
  uint_fast64_t back;

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

  for (int i = 1; i < 8 && need_coarse; i++) {
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

  auto &leafs = this->levels.back();

  auto new_values = vector<pair<Node<T>*, T>>();
  std::stack<NodeIndex> stack;
  new_values.reserve(leafs.size());

//  buffer

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
void AdaptiveMesh<T>::PrepareFacesYZ(const Node<T> &node) {
  const auto max_coord = (1 << this->max_level_);
  const auto index = node.Index();
  const auto x = index.X();

  if (x > 0) {
    const auto neighbor = index.AdjacentX(-1);
    const auto np = this->GetNode(neighbor);
    if (np == nullptr) {
      faces.emplace_back(index, neighbor);
    }
  }

  if (x + 1 < max_coord) {
    const auto neighbor = index.AdjacentX(1);
    const auto np = this->FindLeaf(neighbor);
    if (np && np->value != node.value) {
      faces.emplace_back(index, neighbor);
    }
  }
}

template <typename T>
void AdaptiveMesh<T>::PrepareFacesXY(const Node<T> &node) {
  const auto max_coord = (1 << this->max_level_);
  const auto index = node.Index();
  const auto z = index.Z();

  if (z > 0) {
    const auto neighbor = index.AdjacentZ(-1);
    const auto np = this->GetNode(neighbor);
    if (np == nullptr) {
      faces.emplace_back(index, neighbor);
    }
  }

  if (z + 1 < max_coord) {
    const auto neighbor = index.AdjacentZ(1);
    const auto np = this->FindLeaf(neighbor);
    if (np && np->value != node.value) {
      faces.emplace_back(index, neighbor);
    }
  }
}

template <typename T>
void AdaptiveMesh<T>::PrepareFacesXZ(const Node<T> &node) {
  const auto max_coord = (1 << this->max_level_);
  const auto index = node.Index();
  const auto y = index.Y();


  if (y > 0) {
    const auto neighbor = index.AdjacentY(-1);
    const auto np = this->GetNode(neighbor);
    if (np == nullptr) {
      faces.emplace_back(index, neighbor);
    }
  }

  if (y + 1 < max_coord) {
    const auto neighbor = index.AdjacentY(1);
    const auto np = this->FindLeaf(neighbor);
    if (np && np->value != node.value) {
      faces.emplace_back(index, neighbor);
    }
  }
}
//template <typename T>
//void AdaptiveMesh<T>::prepareChildList() {



//  static const array<uint, 4> child_codes = {0b000, 0b011, 0b101, 0b110};
//  const auto &leaf_parents = this->levels[this->max_level_ - 1];



//  childs_for_refine.resize(0);

//  #pragma omp parallel
//  {
//    std::vector<int> vec_private;
//
//    #pragma omp for
//    for (size_t i = 0; i < leaf_parents.bucket_count(); i++) {
//
//      for (auto parent = leaf_parents.begin(i); parent != leaf_parents.end(i); ++parent) {
//        if (parent->second.is_leaf) continue;
//        const auto lc = parent->first << 3;
//
//        for (const auto &child_code : child_codes) {
//          vec_private.emplace_back(lc | child_code);
//        }
//      }
//
//    }
//    #pragma omp critical
//    childs_for_refine.insert(childs_for_refine.end(), vec_private.begin(), vec_private.end());
//
//  };




//  for (const auto &parent : leaf_parents) {
//    if (parent.second.is_leaf) continue;
//
//
//
//    const auto lc = parent.first << 3;
//
//    for (const auto &child_code : child_codes) {
//      childs_for_refine.emplace_back(lc | child_code);
//    }
//  }
//}

template <typename T>
void AdaptiveMesh<T>::prepareFaces() {

  const auto &leafs = this->levels[this->max_level_];
  const auto max_coord = (1 << this->max_level_);

  faces.clear();

  for (const auto &leaf : leafs) {
    const auto node = leaf.second;
//    const auto index = node.Index();

    PrepareFacesXY(node);
    PrepareFacesYZ(node);
    PrepareFacesXZ(node);

  }

//  const auto index = node.Index();
//  const auto x = index.X();

//  if (x > 0) {
//    const auto neighborX = index.AdjacentX(-1);
//    const auto np = GetNode(neighborX);
//    if (np == nullptr) {
//      faces.emplace_back(index, neighborX);
//    }
//  }
//
//  if (x + 1 < max_coord) {
//    const auto neighborX = index.AdjacentX(1);
//    const auto np = GetNode(neighborX);
//    if (np && np->value != node.value) {
//      faces.emplace_back(index, neighborX);
//    }
//  }


////  faces[]
//
//  faces.resize(0);
////  #pragma omp parallel for
////  for (size_t j = 0; j < childs_for_refine.size(); j++) {
////    const auto &child_index = childs_for_refine[j];
////  }
//
//
//  for (const auto &child_index : childs_for_refine) {
//    const auto child_neighbors = NodeIndex::GetAdjacentIndex(NodeIndex{child_index});
//    const auto child = this->GetNode(NodeIndex{child_index});
//
//    for (const auto &neighbor_index : child_neighbors) {
//      const auto neighbor = FindOrCreateNode(neighbor_index);
//      if (neighbor->value != child->value)
//        faces.emplace_back(child, neighbor);
////        faces.push_back(Face{child, neighbor});
//    }
//  }
//
//
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
