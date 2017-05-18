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

  array<Ghost*, 6> ghosts;
  array<Voxel*, 6> neighbor;

};


// TODO(Sotnem13):
struct Ghost {
  Voxel* top;
  Voxel* bottom;
  double delta_value;
  Point delta_pos;  //

};


class DiffusionOp {

    const auto D = ;

 public:
  DiffusionOp() {}
  explicit DiffusionOp(int iteration_count)
      : iteration_count_(iteration_count) {}
  ~DiffusionOp() {}
  DiffusionOp(const DiffusionOp&) = delete;
  DiffusionOp& operator=(const DiffusionOp&) = delete;

  template <typename TVoxel>
  inline Ghost makeGhost(int ghost_index, TVoxel top, TVoxel bottom) {

  }

  template <typename TContainer>
  void Compute(TContainer* voxels) const {

// TODO(Sotnem13): use refinement octree

    double search_radius = 1;
    NeighborOp op(search_radius);
    op.Compute(voxels);

    sortNeighbors(voxels);
    std::vector<Ghost> ghosts;

    for (int j = 0; j < iterationCount_; j++) {

      for (auto voxel : voxels) {
        for (int i = 0; i < maxNeighborsCount; i++) {
          auto& ghost = voxel.ghosts[i];
          if (!ghost) {
            auto& neighbor = voxel.neighbor[i];

            // Computing concentration change
            auto delta_value = -voxel.concentration;
            if (neighbor) {
              delta_value += neighbor->concentration;
            }
            delta_value *= D / maxNeighborsCount;

            //
            ghost = makeGhost(i, voxel, neighbor);
            ghost->delta_value = delta_value;
            ghosts.push_back(ghost);
          }
        }
      }

      for (auto &ghost : ghosts) {
        auto top = *ghost.top;
        auto bottom = *ghost.bottom;

        auto delta_value = ghost.delta_value;
        if( delta_value > param::kMinimalDifferenceConcentrationForExtracacellularDiffusion) {
          //
          if (!bottom) {
            auto pos = top->getPosition() - ghost.delta_pos;
            bottom = createAndAddVoxelTo(voxels, pos);
          }
          // Update concentration
          top->concentration    += delta_value;
          bottom->concentration -= delta_value;
        }
      }
    }
  }

 private:
    int iteration_count_ = 100;
};
}  // namespace bdm

#endif //BIODYNAMO_DIFFUSION_OP_H
