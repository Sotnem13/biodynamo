#ifndef BIODYNAMO_DIFFUSION_OP_H
#define BIODYNAMO_DIFFUSION_OP_H

#include <vector>

namespace bdm {
using spatial_organization::SpatialTreeNode;
using spatial_organization::Point;
using spatial_organization::Voxel;
using spatial_organization::AdaptiveMesh;


using std::vector;
using std::array;

// TODO(Sotnem13):
template <typename T>
struct ConcentrationChange {
  Voxel<T>* top;
  Voxel<T>* bottom;
  double delta_value;
};


class DiffusionOp {

 public:
  DiffusionOp() {}
  explicit DiffusionOp(int iteration_count = 1, double diffusion_constant = 0.5)
      : iteration_count_(iteration_count), diffusion_constant_(diffusion_constant){}
  ~DiffusionOp() {}
  DiffusionOp(const DiffusionOp&) = delete;
  DiffusionOp& operator=(const DiffusionOp&) = delete;

  template <typename T>
  void Compute(AdaptiveMesh<T>* mesh) const {
    for (int i = 0; i < iteration_count_; i++) {

      mesh->Refine([D = diffusion_constant_]
                           (typename AdaptiveMesh<T>::Element &mesh_element) {
        auto neighbors = mesh_element.Neighbors();
        auto delta_value = 0; ;
        for (auto neighbor: neighbors) {
          delta_value += (neighbor.first - mesh_element.Value());
          std::cout << 111 << std::endl;
        }
        auto new_value = mesh_element.Value();
        new_value += delta_value*D/6.0;
        mesh_element.SetValue(new_value);
      }, false);
    }
  }

 private:
  int iteration_count_;
  double diffusion_constant_;
};
}  // namespace bdm

#endif //BIODYNAMO_DIFFUSION_OP_H
