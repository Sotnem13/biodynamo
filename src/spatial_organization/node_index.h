#ifndef BIODYNAMO_NODE_INDEX_H
#define BIODYNAMO_NODE_INDEX_H

#include <vector>
namespace bdm {
namespace spatial_organization {

using std::vector;
/// Implement Morton keys
/// The key k(n) of a node n can be generated recursively from the
/// octree hierarchy: the key of the root is 1, and the key of the child
/// of n is the concatenation of k(n) with the 3 bits of the child octant.
/// The key k(n) can also be generated from the depth l of n and
/// the position (x, y, z) of its center: assuming that
/// the root is the unit cube [0, 1]^3 , the key k(n) is computed
/// by interleaving the bits of x, y and z:
/// k(n) = 1  xl yl zl  xl−1 yl−1 zl−1 ... x1 y1 z1
/// Interleaving was accelerated by integer dilation.
struct NodeIndex {

  /// Empty constructor, initializes key with 1 (root index)
  NodeIndex();

  /// Constructor
  /// Generate key from the depth l of node and
  /// the position (x, y, z)
  /// @param x - node coordinate on axis x
  /// @param y - node coordinate on axis y
  /// @param z - node coordinate on axis z
  /// @param level - depth of node
  NodeIndex(uint x, uint y, uint z, uint level);

  /// Constructor
  /// @param code - key of the node
  explicit NodeIndex(uint_fast64_t code);

  /// The depth of node is floor(log2(k)/3)
  uint Level();

  /// The key of the parent of node is obtained
  /// by truncating the 3 least significant bits of k(n)
  NodeIndex Parent();

  /// Restore the coordinate from key
  uint X();
  uint Y();
  uint Z();

  /// Integer dilation, for interleaving the bits of coordinate
  /// For more information, see
  /// https://www.researchgate.net/publication/3649954_Integer_dilation_and_contraction_for_quadtrees_and_octrees
  /// @param coord - coordinate of the node
  static inline uint_fast64_t oct_dilate(uint coord);

  /// Reverse operation of integer dilation
  /// @param code - key of the node
  static inline uint oct_contract(uint_fast64_t code);

  /// For get adjacent index need restore the position (x, y, z) of node and
  /// just generate indexes by shifting the coordinates
  /// like a (x+-1, y, z) (x, y+-1, z) (x, y, z+-1).
  /// if template argument is true then result will include angular node indexes
  /// like a (x-1, y-1, z-1) (x+1, y-1, z-1)
  /// @tparam T - type of the return adjacent Index
  template <bool>
  static vector<NodeIndex> GetAdjacentIndex(NodeIndex index, int at_level);


  union {
    /// Key k(n)
    uint_fast64_t code;
    struct {
      /// Last 3 bits of key k(n)
      uint_fast64_t x : 1;
      uint_fast64_t y : 1;
      uint_fast64_t z : 1;
    };
  };

};

NodeIndex::NodeIndex() : code(1) {}

NodeIndex::NodeIndex(uint_fast64_t code_) : code(code_) {}

NodeIndex::NodeIndex(uint x, uint y, uint z, uint level) {
  code = oct_dilate(x)    |
       oct_dilate(y) << 1 |
       oct_dilate(z) << 2;
  code ^= 1l << (3*level);
}

NodeIndex NodeIndex::Parent() {
  return NodeIndex(this->code >> 3);
}

uint NodeIndex::Level() {
  static const int tab64[64] = {
          63,  0, 58,  1, 59, 47, 53,  2,
          60, 39, 48, 27, 54, 33, 42,  3,
          61, 51, 37, 40, 49, 18, 28, 20,
          55, 30, 34, 11, 43, 14, 22,  4,
          62, 57, 46, 52, 38, 26, 32, 41,
          50, 36, 17, 19, 29, 10, 13, 21,
          56, 45, 25, 31, 35, 16,  9, 12,
          44, 24, 15,  8, 23,  7,  6,  5
  };

  uint_fast64_t level = code;
  level |= level >> 1;
  level |= level >> 2;
  level |= level >> 4;
  level |= level >> 8;
  level |= level >> 16;
  level |= level >> 32;

  return tab64[((uint64_t)((level - (level >> 1))*0x07EDD5E59A4E28C2)) >> 58]/3;
}

uint NodeIndex::X() {
  return oct_contract(code ^ (1l << 3*Level()) );
}
uint NodeIndex::Y() {
  return oct_contract(code >> 1);
}
uint NodeIndex::Z() {
  return oct_contract(code >> 2);
}

inline uint_fast64_t NodeIndex::oct_dilate(uint coord) {
  uint_fast64_t code = coord;
  code = (code | code << 32) & 0x001F00000000FFFF;
  code = (code | code << 16) & 0x001F0000FF0000FF;
  code = (code | code <<  8) & 0x100F00F00F00F00F;
  code = (code | code <<  4) & 0x10C30C30C30C30C3;
  code = (code | code <<  2) & 0x1249249249249249;
  return code;
}
inline uint NodeIndex::oct_contract(uint_fast64_t code) {
  code = code & 0x1249249249249249;
  code = (code | code >>  2) & 0x10C30C30C30C30C3;
  code = (code | code >>  4) & 0x100F00F00F00F00F;
  code = (code | code >>  8) & 0x001F0000FF0000FF;
  code = (code | code >> 16) & 0x001F00000000FFFF;
  code = (code | code >> 32) & 0x00000000001FFFFF;
  return code;
}

/// Didn't invent of anything better than this =(
template <>
vector<NodeIndex> NodeIndex::GetAdjacentIndex<true>(NodeIndex index, int at_level) {
  int node_level = index.Level();
  vector<NodeIndex> result;

  uint x = index.X();
  uint y = index.Y();
  uint z = index.Z();

  at_level = std::max(at_level,node_level);
  int level_offset = at_level - node_level;

  int nodes_count = 1 << level_offset;

  int ey = ((y + 1 < 1 << node_level) ? 1 : 0) + nodes_count;
  int ez = ((z + 1 < 1 << node_level) ? 1 : 0) + nodes_count;

  int by = y > 0 ? -1 : 0;
  int bz = z > 0 ? -1 : 0;


  if (x > 0) {
    auto nx = (x << level_offset) - 1;

    for (int dy = by; dy < ey; dy++) {
      for (int dz = bz; dz < ez; dz++) {
        int ny = (y << level_offset) + dy;
        int nz = (z << level_offset) + dz;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }

  if (x + 1 < 1 << node_level) {
    auto nx = (x + 1) << level_offset;
    for (int dy = by; dy < ey; dy++) {
      for (int dz = bz; dz < ez; dz++) {
        int ny = (y << level_offset) + dy;
        int nz = (z << level_offset) + dz;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }

  if (y > 0) {
    auto ny = (y << level_offset) - 1;
    for (int dx = 0; dx < nodes_count; dx++) {
      for (int dz = bz; dz < ez; dz++) {
        int nx = (x << level_offset) + dx;
        int nz = (z << level_offset) + dz;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }

  if (y + 1 < 1 << node_level) {
    auto ny = (y + 1) << level_offset;
    for (int dx = 0; dx < nodes_count; dx++) {
      for (int dz = bz; dz < ez; dz++) {
        int nx = (x << level_offset) + dx;
        int nz = (z << level_offset) + dz;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }

  if (z > 0) {
    auto nz = (z << level_offset) - 1;
    for (int dy = 0; dy < nodes_count; dy++) {
      for (int dx = 0; dx < nodes_count; dx++) {
        int ny = (y << level_offset) + dy;
        int nx = (x << level_offset) + dx;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }

  if (z + 1 < 1 << node_level) {
    auto nz = (z + 1) << level_offset;
    for (int dy = 0; dy < nodes_count; dy++) {
      for (int dx = 0; dx < nodes_count; dx++) {
        int ny = (y << level_offset) + dy;
        int nx = (x << level_offset) + dx;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }
  return result;
}

template <>
vector<NodeIndex> NodeIndex::GetAdjacentIndex<false>(NodeIndex index, int at_level) {
  int node_level = index.Level();
  vector<NodeIndex> result;

  uint x = index.X();
  uint y = index.Y();
  uint z = index.Z();

  at_level = at_level > node_level ? at_level : node_level;
  int level_offset = at_level - node_level;

  int nodes_count = 1 << level_offset;
  if (x > 0) {
    auto nx = (x << level_offset) - 1;
    for (int dy = 0; dy < nodes_count; dy++) {
      for (int dz = 0; dz < nodes_count; dz++) {
        int ny = (y << level_offset) + dy;
        int nz = (z << level_offset) + dz;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }

  if (x + 1 < 1 << node_level) {
    auto nx = (x + 1) << level_offset;
    for (int dy = 0; dy < nodes_count; dy++) {
      for (int dz = 0; dz < nodes_count; dz++) {
        int ny = (y << level_offset) + dy;
        int nz = (z << level_offset) + dz;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }

  if (y > 0) {
    auto ny = (y << level_offset) - 1;
    for (int dx = 0; dx < nodes_count; dx++) {
      for (int dz = 0; dz < nodes_count; dz++) {
        int nx = (x << level_offset) + dx;
        int nz = (z << level_offset) + dz;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }

  if (y + 1 < 1 << node_level) {
    auto ny = (y + 1) << level_offset;
    for (int dx = 0; dx < nodes_count; dx++) {
      for (int dz = 0; dz < nodes_count; dz++) {
        int nx = (x << level_offset) + dx;
        int nz = (z << level_offset) + dz;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }

  if (z > 0) {
    auto nz = (z << level_offset) - 1;
    for (int dy = 0; dy < nodes_count; dy++) {
      for (int dx = 0; dx < nodes_count; dx++) {
        int ny = (y << level_offset) + dy;
        int nx = (x << level_offset) + dx;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }

  if (z + 1 < 1 << node_level) {
    auto nz = (z + 1) << level_offset;
    for (int dy = 0; dy < nodes_count; dy++) {
      for (int dx = 0; dx < nodes_count; dx++) {
        int ny = (y << level_offset) + dy;
        int nx = (x << level_offset) + dx;
        result.push_back(NodeIndex(nx, ny, nz, at_level));
      }
    }
  }
  return result;
}

}  // namespace spatial_organization
}  // namespace bdm

#endif //BIODYNAMO_NODE_INDEX_H_H
