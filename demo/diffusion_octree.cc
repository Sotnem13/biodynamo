#include <iostream>

#include "spatial_organization/adaptive_mesh.h"
#include "diffusion_op.h"

using namespace bdm::spatial_organization;
using bdm::DiffusionOp;
using namespace std;


int main(int args, char** argv) {


  Bound b = {0, 0, 0, 128, 128, 128};
//  MeshData<double> m(b, 0, 21);
  AdaptiveMesh<double> m(b, 0, 7);


  auto x = 64.0;
  auto y = 64.0;
  auto z = 64.0;

  Point p{x, y, z};

  DiffusionOp op(1, 0.3);

  for (int i = 0; i < 10; i++) {
    m.Put(p, 1);
    op.Compute(&m);
  }

//  ofstream out("output");

  for (int x_ = -11; x_ < 12; x_++) {
    for (int y_ = -11; y_ < 12; y_++) {
      cout << m.At({ x + x_,y + y_,z})* 1e10 << " ";
    }
    cout << endl;
  }

//  out.close();

}