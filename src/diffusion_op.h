#ifndef BIODYNAMO_DIFFUSION_OP_H
#define BIODYNAMO_DIFFUSION_OP_H

#include "spatial_organization/adaptive_mesh.h"

namespace bdm {

using spatial_organization::AdaptiveMesh;


class DiffusionOp {

 public:
  explicit DiffusionOp(int iteration_count = 1, double diffusion_constant = 0.5)
      : iteration_count_(iteration_count), diffusion_constant_(diffusion_constant){}
  ~DiffusionOp() {}
  DiffusionOp(const DiffusionOp&) = delete;
  DiffusionOp& operator=(const DiffusionOp&) = delete;

  template <typename T>
  void Compute1(AdaptiveMesh<T>* mesh) const {
    const double D = diffusion_constant_;
    for (int i = 0; i < iteration_count_; i++) {
      mesh->Refine([D](typename AdaptiveMesh<T>::Element &mesh_element) {
        auto neighbors = mesh_element.Neighbors();
        auto delta_value = 0.0;
        for (auto neighbor : neighbors) {
          delta_value += (neighbor.first - mesh_element.Value());
        }
        auto new_value = mesh_element.Value();
        new_value += delta_value*D/6.0;
        mesh_element.SetValue(new_value);
      }, false);
    }
  }

  template <typename T>
  void Compute(AdaptiveMesh<T>* mesh) const {
    const double D = diffusion_constant_;

    const auto &diffusion_flow = [D](typename AdaptiveMesh<T>::Face &face) {
      const auto concentration_difference = face.BackValue() - face.FrontValue();
      const auto flux = concentration_difference*D/6;
      face.SetFlux(flux);
    };

    for (int i = 0; i < iteration_count_; i++) {
      mesh->Refine(diffusion_flow);
    }
  }

 private:
  int iteration_count_;
  double diffusion_constant_;
};
}  // namespace bdm

#endif //BIODYNAMO_DIFFUSION_OP_H
