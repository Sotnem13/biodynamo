#ifndef BIODYNAMO_VOXEL_H
#define BIODYNAMO_VOXEL_H

#include <array>
#include <vector>

namespace bdm {
namespace spatial_organization {
//namespace geometric_shapes {

using std::array;


//class Ghost<T> {
//
//};


//Node


template <typename T> class Face;

template <typename T>
class Voxel {
 public:

  Voxel() = default;

  explicit Voxel(T init_val) : value(init_val) {
//    for (auto &face : faces) {

//    }
  }

  ~Voxel() = default;


  T& GetValue() {
    return value;
  }

  void SetValue(T new_value) {
    value = new_value;
  }

  array<Face<T>*, 6>& GetFaces() {
    return faces;
  }



 private:
  Point position;
  array<Face<T>*, 6> faces;
  T value;
  template <typename T1>
  friend class MeshTreeNode;
};



template <typename T>
class Face {
 public:

  Voxel<T>* LargeVoxel() {
    return large;
  }

  vector<Voxel<T>*> SmallerVoxels() const {
    return smaller;
  }

private:
  Voxel<T> *large;
  vector<Voxel<T>*> smaller;
};


//}  //  namespace geometric_shapes
}  //  namespace spatial_organization
}  //  namespace bdm
#endif //BIODYNAMO_VOXEL_H
