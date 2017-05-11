#ifndef BIODYNAMO_DIFFUSION_OP_H
#define BIODYNAMO_DIFFUSION_OP_H

namespace bdm {
using spatial_organization::SpatialTreeNode;
using spatial_organization::Point;

using std::array;

// TODO(Sotnem13): move to separate file
class Voxel {
public:
  Point position_;
  double concentration = 0;
};

// TODO(Sotnem13):
struct G {
  Voxel* top;
  Voxel* bottom;
  double dvalue;
};

template <
        template < typename T1 > class T >
class DiffusionOp {
 public:
  DiffusionOp() {}
  explicit DiffusionOp(int iterationCount) : iterationCount_(iterationCount) {}
  ~DiffusionOp(){}
  DiffusionOp(const DiffusionOp&) = delete;
  DiffusionOp& operator=(const DiffusionOp&) = delete;

  template <typename TContainer>
  void Compute(TContainer* voxels) const {

    double search_radius = 1;
    NeighborOp op(search_radius);

    for (int i = 0; i < iterationCount_; i++) {
//      G{voxels};

      for (auto &g : gg) {
        auto t = g.top;
        auto b = g.bottom;

//        if (!t) {
//          auto pos = b->getPosition() + g.dpos;
//          t = createAndAddVoxelTo(voxels, pos);
//        } else if (!b) {
//          auto pos = t->getPosition() - g.dpos;
//          b = createAndAddVoxelTo(voxels, pos);
//        }
        t->concentration += g.dvalue;
        b->concentration -= g.dvalue;
      }
    }
  }

 private:
    int iterationCount_ = 100;
};
}  // namespace bdm

#endif //BIODYNAMO_DIFFUSION_OP_H
