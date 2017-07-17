#include <iostream>
#include <array>
#include <bitset>


//#include "spatial_organization/octree_node.h"
#include "spatial_organization/adaptive_mesh.h"
#include "diffusion_op.h"
#include <spatial_organization/mesh_data.h>
#include <spatial_organization/bound.h>

using namespace bdm::spatial_organization;
using bdm::DiffusionOp;
using namespace std;






int main(int args, char** argv) {

//  NodeIndex index(1,1,1,2);

//  cout << index.Level() << endl;
//  cout << bitset<10>(index.code) << endl;
//  int level  = 2;
//  int at_level = 2;
//
//  int x = 1, y = 1, z = 1;

//  auto result = NodeIndex::GetAdjacentIndex<true>(index, 3);

//  cout << result.size() << endl;
//  for (auto res : result)
//    cout << bitset<12>(res.code) << endl;
//

  auto x = 149.;
  auto y = 50.;
  auto z = 128.;
  auto p = Point(x, y, z);
//
  Bound b = {0, 0, 0, 100, 100, 100};
////
//  MeshData<double> m(b, 0, 10);
//
//  NodeIndex n(1207701440);
////
//  cout << get<0>(n.Position()) << endl;



//  m.Put({100,100,100},5);
//  cout << m.At({100,100,100}) << endl;
//  cout << m.At({100,100,101}) << endl;
  AdaptiveMesh<double> m(b, 0, 10);
//
  m.Put(p, 1);
  DiffusionOp diffusionOp(1);

  diffusionOp.Compute(&m);

  cout << m.At({x, y, z}) << endl;
  cout << m.At({x-1, y, z}) << endl;
  cout << m.At({x+1, y, z}) << endl;
  cout << m.At({x, y-1, z}) << endl;




}