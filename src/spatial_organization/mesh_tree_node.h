#ifndef BIODYNAMO_MESHTREENODE_H
#define BIODYNAMO_MESHTREENODE_H

#include <vector>
#include <memory>

#include "spatial_organization/spatial_tree_node.h"
#include "spatial_organization/bound.h"
#include "spatial_organization/adaptive_mesh.h"
#include "spatial_organization/voxel.h"

namespace bdm {
namespace spatial_organization {


using std::vector;
using std::array;
using std::pair;

using std::unique_ptr;

//template <typename T> class AdaptiveMesh;

template <typename T>
class MeshTreeNode {
public:

  MeshTreeNode();
  MeshTreeNode(const Bound &bnd, T empty_value, int level);

  ~MeshTreeNode();

  T At(const Point &p) const;
  void Put(const Point &p, T obj);  // override;

  bool IsLeaf() const;  // override;

  /// Splits node to 8 equal subspaces
  /// @tparam T - type of objects
  void Split();

 private:
  using Node = MeshTreeNode<T>;

  bool is_leaf_;
  int level_;

  Bound bound_;

  array<unique_ptr<Node>, 8> childrens_;
  T value_;

  void Put(const Point &p, T obj, vector<Node*> neighbors);


  /// returns in what octa subspace point p is located
  /// @tparam T - type of the object
  /// @param p - point in space
  /// @return number of subspace
  int GetChildID(Point const &p) const;

  //
  bool Coarse();

  //
//  size_t GetChildrenSize() const override;

  //
//  const vector<pair<Point, T>> &GetObjects() const override;
//  SpatialTreeNode<T> **GetChildrenNodes() const override ;

  template <typename T1>
  friend class AdaptiveMesh;
};

template <typename T>
MeshTreeNode<T>::MeshTreeNode()
        : MeshTreeNode(Bound(0, 0, 0, 1, 1, 1), T(0), 0) {}

template <typename T>
MeshTreeNode<T>::MeshTreeNode(const Bound &bnd, T init_val, int level)
        : bound_(bnd), is_leaf_(true),
          level_(level), value_{init_val} {  //, SpatialTreeNode<T>(bnd) {

  this->voxel_ = unique_ptr<VoxelType>(
          new VoxelType(make_pair(this, value_))
  );
}

template <typename T>
MeshTreeNode<T>::~MeshTreeNode() {

}

template <typename T>
int MeshTreeNode<T>::GetChildID(const Point &p) const {
  auto center = bound_.Center();
  int x = center.x_ > p.x_;
  int y = center.y_ > p.y_;
  int z = center.z_ > p.z_;
  return (z << 2) | (y << 1) | x;
}

template <typename T>
void MeshTreeNode<T>::Put(Point const &p, T obj) {
  Put(p, obj, {});
//  if (IsLeaf()) {
//    if (value_ != obj) {
//      if (level_ > 0) {
//        Split();
//        Put(p, obj);
//      } else {
//        value_ = obj;
//      }
//    }
//  } else {
//    auto child_index = GetChildID(p);
//    auto neighbors = vector<Node*>();
//
//    neighbors.reserve(7);
//    for (int i = 0; i < childrens_.size(); ++i) {
//      if (i != child_index) {
//        auto child_ptr = childrens_[i].get();
//        neighbors.push_back(child_ptr);
//      }
//    }
//
//    childrens_[child_index]->Put(p, obj, std::move(neighbors));
//  }
}

template <typename T>
void MeshTreeNode<T>::Put(const Point &point, T obj, vector<Node*> neighbors) {
  if (IsLeaf()) {
    if (voxel_->value.second != obj) {
      if (level_ > 0) {
        Split();
        Put(point, obj, neighbors);
      } else {
        value_ = obj;
      }
    }
  } else  {
    auto center = bound_.Center();
    int x = center.x_ > point.x_;
    int y = center.y_ > point.y_;
    int z = center.z_ > point.z_;

    auto child_index = (z << 2) | (y << 1) | x;

    vector<Node *> neighbors_;
    neighbors_.reserve(26);  // 3 * 3 * 3 - 1

    for (int i = 0; i < childrens_.size(); ++i) {
      if (i != child_index) {
        auto child_ptr = childrens_[i].get();
        neighbors_.push_back(child_ptr);
      }
    }


    for (auto neighbor : neighbors) {
      auto n_bnd = neighbor->bound_;
      auto n_center = n_bnd.Center();

//      n_center.coords.map( []() )
//      for (auto coor : n_center.coords) {
//        coor
//      }

      int x_ = center.x_ > n_center.x_;
      int y_ = center.y_ > n_center.y_;
      int z_ = center.z_ > n_center.z_;

      int xx = x == x_ || n_center.x_ == center.x_;
      int yy = y == y_ || n_center.y_ == center.y_;
      int zz = z == z_ || n_center.z_ == center.z_;

      if (xx && yy & zz) {
        if (neighbor->level_ - this->level_ > 0) {
          if (neighbor->is_leaf_) {
            neighbor->Split();
          }
          auto n_child_index = neighbor->GetChildID(point);
          auto n_child = neighbor->childrens_[n_child_index].get();
          neighbors_.push_back(n_child);
        } else {
          neighbors_.push_back(neighbor);
        }
      }
    }

    childrens_[child_index]->Put(point, obj, neighbors_);
  }
}


template <typename T>
T MeshTreeNode<T>::At(const Point &p) const {
  if (IsLeaf()) {
    return value_;
  } else {
    auto index = GetChildID(p);
    return childrens_[index]->At(p);
  }
}


template <typename T>
inline bool MeshTreeNode<T>::IsLeaf() const {
  return is_leaf_ || level_ == 0;
}

template <typename T>
void MeshTreeNode<T>::Split() {
  if (IsLeaf() && level_ > 0) {
    auto center = this->bound_.Center();
    auto far_left_bottom_point = this->bound_.far_left_bottom_point_;
    auto near_right_top_point  = this->bound_.near_right_top_point_;

    auto children_size =
            (near_right_top_point + far_left_bottom_point*(-1))*0.5;

    //
    for (auto face : this->voxel_->faces) {
      auto large_voxel = face->LargeVoxel();
      if (large_voxel != this->voxel_.get()) {
        auto node = large_voxel->GetValue().first;
        if ((node->level_ - this->level_) > 0) {
          node->Split();
        }
      }
    }

    //      auto smaller_voxels = face->SmallerVoxels();

    for (int i = 0; i < childrens_.size(); ++i) {
      int x = i & Coord::X ? 1 : 0;
      int y = i & Coord::Y ? 1 : 0;
      int z = i & Coord::Z ? 1 : 0;

      auto offset = Point(children_size.x_ * x,
                          children_size.y_ * y,
                          children_size.z_ * z);

      auto child_bound = Bound(far_left_bottom_point + offset, center + offset);

      auto child = new Node(child_bound, value_, level_ - 1);
      childrens_[i] = unique_ptr<Node>(child);
    }



    this->is_leaf_ = false;
//    this->voxel_ = nullptr;
  }
}



//      auto x1_face = this->voxel_->faces[x + 0];
//      auto y1_face = this->voxel_->faces[y + 2];
//      auto z1_face = this->voxel_->faces[z + 4];

//      childrens_[( z << 2) | ( y << 1) | !x];
//      childrens_[( z << 2) | (!y << 1) |  x];
//      childrens_[(!z << 2) | ( y << 1) |  x];

//      child->voxel->faces[face1_index] = face1;
//      child->voxel->faces[face2_index] = face2;



//      for (int j = 0; j < 3; ++j) {
//        int  coord = (i >> j) & 1;
//
//        auto face1_index =  coord + j*2;
//        auto face2_index = !coord + j*2;

//        auto opposite_child_index = ( z << 2) | ( y << 1) | !x

//        auto face_1 = this->voxel_->faces[face1_index];
//
//        if (childrens_[i]) {
//
//        }
//
//        auto x_face = this->voxel_->faces[x];
//      }



//      for (auto face : this->voxel_->faces) {
//        auto larger_voxel = face->LargeVoxel();
//        auto smaller_voxels = face->SmallerVoxels();
//
//        if (smaller_voxels.size() == 1) {
//          if (larger_voxel == this->voxel_.get()) {
//          } else {
//            larger_voxel.face
//          }
//        } else {
//
//        }
//
//          if (smaller_voxels.size() > 1) {
//
//          } else {
//
//          }
//        } else {
//
//        }
//        if (larger_voxel != this->voxel_.get()) {
////          auto node = larger_voxel->GetValue().first;
////          if ((node->level_ - this->level_) > 0) {
////            node->Split();
////          }
//        }
//      }

//      auto cvoxel = child->voxel_;
//      child->level_;

//    }





    //        auto face->SmallVoxels();
//      auto opposite_voxel = face->GetOpposite(this->voxel_.get());
//      auto opposite_nodes = face->OppositeNodes(this);
//      if (opposite_nodes.size() = 1) {
//        opposite_nodes[0]->Split();
//      } else {

//      }
//    }

//      auto faces = voxel_->GetFaces(x, y, z);
//      child->voxel_->faces[0] =
//      auto face_x = voxel_->faces[x];
//      auto face_x = voxel_->FrontFace();

//      if()


//    }

//      auto value = voxel_.value;
//      Voxel<T> voxel;
//      voxel.value = value;

//      voxel.SetValue()

//      auto face_x = &voxel_.faces[x + 0]; // 0 1
//      auto face_y = &voxel_.faces[y + 2]; // 2 3
//      auto face_z = &voxel_.faces[z + 4]; // 4 5

//      voxel.faces[x] = face_x;
//      voxel.faces[y + 2] = face_y;
//      voxel.faces[z + 4] = face_z;


//      auto child_x = childrens_[!x | ( y << 1) | ( z << 2)];
//      auto child_y = childrens_[ x | (!y << 1) | ( z << 2)];
//      auto child_z = childrens_[ x | ( y << 1) | (!z << 2)];

//      for (int j = 0; j < 3; j++ ) {
//        int coord = (i >> j) & 1;
//        int face_index = coord + j*2;

//        auto f1 = voxel_.faces[face_index];
//        auto f2 = fa

//        auto child = childrens_[!x | ( y << 1) | ( z << 2)];
//      }
//      if (child_x) {
//        auto child_voxel = &(child_x->voxel_);
//        voxel.faces[!x] = child_voxel->faces[x];
//      }
//      if (child_y) {
//        auto child_voxel = &(child_y->voxel_);
//        voxel.faces[!y + 2] = child_voxel->faces[x];
//      }
//      if (child_z) {
//        auto child_voxel = &(child_z->voxel_);
//        voxel.faces[!z + 4] = child_voxel->faces[x];
//      }

//      voxel.faces[!x] =




//    }





}  //  namespace spatial_organization
}  //  namespace bdm


#endif //BIODYNAMO_MESHTREENODE_H
