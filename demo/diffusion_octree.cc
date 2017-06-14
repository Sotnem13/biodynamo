#include <iostream>
#include <array>

#include "spatial_organization/octree_node.h"
#include "spatial_organization/adaptive_mesh.h"
#include "diffusion_op.h"

using namespace bdm::spatial_organization;
using bdm::DiffusionOp;
using namespace std;

template <int level>
struct m {
  bitset<3*level> m;
};


int main(int args, char** argv) {

//  NodeIndex i = 1;
//  cout
//  bitset<level>

//
////  MeshNode<double> m;
////  MeshData<double> m(b, 10);
////  auto s = m.size_;
//
//
////  std::cout << s.x_ * s.y_ * s.z_;

  auto x = 149.;
  auto y = 50.;
  auto z = 128.;
  auto p = Point(x, y, z);
  Bound b = {10, 10, 10, 210, 210, 210};

//  MeshData<double> m(b, 0, 10);
  AdaptiveMesh<double> m(b, 0, 10);

  m.Put(p,1);
  DiffusionOp diffusionOp(1);

  diffusionOp.Compute(&m);
  cout << m.At({x, y, z}) << endl;
  cout << m.At({x-1, y, z}) << endl;
  cout << m.At({x+1, y, z}) << endl;
  cout << m.At({x, y-1, z}) << endl;
//  auto index = m.get_node_index(p);
//  uint x_, y_, z_;
//  tie(x_, y_, z_) = index.Coord();
//
//
////  cout << x_ << " " << y_ << " " << z_ << endl;
////  cout << bitset<30>(index.index) << endl;
//
//  auto size = m.get_node_size_at_level(index.level());
//  double dis = 100;
//  double d   = dis / size.x_;
//  auto index2 = NodeIndex(make_tuple(x_-d, y_, z_));
//  tie(x_, y_, z_) = index2.Coord();
//
//  auto center = m.get_node_bound(index2).Center();
//
//
//  cout << bitset<30>(index2.index) << endl;
//  m.Put(p, 20);
//  cout  << m.At(p) << endl;
//  cout << d << endl;
//  cout << size.x_ << " " << size.y_ <<  " " << size.z_ << endl;
//  cout << center.x_ << " " << center.y_ <<  " " << center.z_ << endl;

//  auto bs = bitset<4>(8);
//  cout << bs.size() << endl;

//  m.Put(p, 50);

//  MeshData<double>::NodeIndex index = 1001001101;

//  for
//  auto p1 = Point{ 99,  99,  99};
//  auto p2 = Point{101, 101, 101};
//  auto p3 = Point{100, 100, 100};
//
//  m.Put(p1, 55);
//  m.Put(p2, 45);
//  m.Put(p3, 35);
//
//  cout << m.At(p1) << endl;
//  cout << m.At(p2) << endl;
//  cout << m.At(p3) << endl;

//  auto bs = bitset<4>(8);
//  cout << bs.count() << endl;

//  auto index = m.get_node_index(p);
//  auto bound = m.get_node_bound(index);
//  auto center = bound.Center();
//  auto size = m.get_node_size_at_level(index.level());

//  p1 = bound.far_left_bottom_point_;
//  p2 = bound.near_right_top_point_;
//
//  cout << p.x_ << " " << p.y_ <<  " " << p.z_ << endl;
//  cout << center.x_ << " " << center.y_ <<  " " << center.z_ << endl;
//  cout << p1.x_ << " " << p1.y_ <<  " " << p1.z_ << endl;
//  cout << p2.x_ << " " << p2.y_ <<  " " << p2.z_ << endl;

//  cout << bound. << endl;

//  auto bs = bitset<4>(8);
//  cout << bs.count() << endl;

//  cout << bitset<15>(index.index) << endl;

//  auto neighbor_nodes_index = m.get_adjacent_neighbors_index(index);


//  cout << neighbor_nodes_index.size() << endl ;
//
//  cout << bitset<15>(index.index) << endl << endl;
//
//  for (auto index : neighbor_nodes_index)
//    cout << bitset<15>(index.index) << endl;

//  auto index = m.get_node_index(p1);

//  index.index = 0b1;
//  auto size = m.get_node_size(index);
//  cout << size.x_ << " " << size.y_ <<  " " << size.z_ << endl;





//  auto coor = 0b001;
//  uint_fast64_t l = -1;
//  size_t l = 0;
//  while (index >> l*3) {
//    l += 1;
//  }
//  for (long i = index.level(); i >= 0; i--) {
//    auto xyz = (index.index >> i*3) & 0b111;
//    for (auto j = 0; j < 3; j++) {
//      cout << ((xyz >> (2 - j)) & 1);
//    }
//    cout << " ";
//  }
//  while (index) {
//    cout << " ";

//    cout << ": ";
//    index >> 1;
//  }



//  cout << m.get_child_code(p) << endl;

//  cout << m.root_index << endl;
//  cout << m.At(p) << endl;
//  cout << m.nodes[1].value << endl;
//  cout << m.nodes.size() << endl;
//
//
//
//
////  AdaptiveMesh<double> m(b, 10);
////  MeshTreeNode<double> m(b, 0, 10);
//  MeshTreeNode<double> m(b, {}, 10);
//
////  m.SetRefinementCriterion();
//  m.Split();
//
//  auto concentration = 5;
//
//  auto voxel = Voxel<double>();
//  voxel.SetValue(concentration);
//
//  m.Put(point, concentration);
//
////  DiffusionOp diff_op;

//  union {
//    int i = 0;
//    struct {
//      int : 5;
//      int x : 1, y : 1, z : 1;
//    };
//  } index;
//  };
//  int x = center.x_ > p.x_;
//  int y = center.y_ > p.y_;
//  int z = center.z_ > p.z_;
//  index.x = 1;
//  index.y = 0;
//  index.z = 1;
//  index.j
//  std::cout << index.i << std::endl;
//  std::cout << m.At(point) << std::endl;
//  std::cout << (5 & 0b100 ? 1 : 0)<< std::endl;

//  diff_op.Compute(m);



//    SpatialTreeNode<size_t>* tree = new OctreeNode<size_t>(
//            Bound(-10000.0, -10000.0, -10000.0, 10000.0, 10000.0, 10000.0), 100,
//            10);

}