#include <iostream>
#include <chrono>

#include "spatial_organization/adaptive_mesh.h"
#include "diffusion_op.h"



using namespace bdm::spatial_organization;
using bdm::DiffusionOp;
using namespace std;
using namespace chrono;



int main(int args, char** argv) {



  Bound b = {0, 0, 0, 128, 128, 128};
  AdaptiveMesh<double> m(b, 0, 7);
  AdaptiveMesh<double> m1(b, 0, 7);
  DiffusionOp op(1, 0.3);




  auto p = Point{64, 64, 64};


  m.Put(p, 1);
  m1.Put(p, 1);

  auto count = 50;

  high_resolution_clock c;

  auto bt = c.now();
  for (int i = 0; i < count; i++) {
    op.Compute(&m);
  }

  auto nt = duration_cast<duration<double>>(c.now() - bt).count();
  cout << "NT: " << nt << endl;

  bt = c.now();
  for (int i = 0; i < count; i++) {
    op.Compute1(&m1);
  }


  auto ot = duration_cast<duration<double>>(c.now() - bt).count();
  cout << "OT: " << ot << endl;
  cout << "Speed Up: " << ot / nt << endl;

  cout << m.At(p) << endl;
  cout << m1.At(p) << endl;

//  auto &l = m1.levels.back();
//  cout << l.bucket_count() << endl;


//  for (int i = 0; i < l.bucket_count(); i++) {
//    cout << l.bucket_size(i) << endl;
//  }

//  cout << l.bucket_count() << endl;


//
//
//  auto x = 64.0;
//  auto y = 64.0;
//  auto z = 64.0;
//
//  Point p{x, y, z};
//
//
//  for (int i = 0; i < 10; i++) {
//  }
//
////  ofstream out("output");
//
//  for (int x_ = -11; x_ < 12; x_++) {
//    for (int y_ = -11; y_ < 12; y_++) {
//      cout << m.At({ x + x_,y + y_,z})* 1e10 << " ";
//    }
//    cout << endl;
//  }
//
////  out.close();

}