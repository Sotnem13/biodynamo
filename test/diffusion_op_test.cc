#include <gtest/gtest.h>

#include <spatial_organization/point.h>
#include "diffusion_op.h"

namespace bdm {
namespace diffusion_op_test_internal {

using spatial_organization::Point;
void runDiffusionStep(double**** Conc, double**** tempConc, int L, double D) {
  // computes the changes in substance concentrations due to diffusion
  int i1,i2,i3, subInd;

  for (i1 = 0; i1 < L; i1++) {
    for (i2 = 0; i2 < L; i2++) {
      for (i3 = 0; i3 < L; i3++) {
        tempConc[0][i1][i2][i3] = Conc[0][i1][i2][i3];
        tempConc[1][i1][i2][i3] = Conc[1][i1][i2][i3];
      }
    }
  }

  int xUp, xDown, yUp, yDown, zUp, zDown;

  for (i1 = 0; i1 < L; i1++) {
    for (i2 = 0; i2 < L; i2++) {
      for (i3 = 0; i3 < L; i3++) {
        xUp = (i1+1);
        xDown = (i1-1);
        yUp = (i2+1);
        yDown = (i2-1);
        zUp = (i3+1);
        zDown = (i3-1);

        for (subInd = 0; subInd < 2; subInd++) {
          if (xUp<L) {
            Conc[subInd][i1][i2][i3] += (tempConc[subInd][xUp][i2][i3]-tempConc[subInd][i1][i2][i3])*D/6;
          }
          if (xDown>=0) {
            Conc[subInd][i1][i2][i3] += (tempConc[subInd][xDown][i2][i3]-tempConc[subInd][i1][i2][i3])*D/6;
          }
          if (yUp<L) {
            Conc[subInd][i1][i2][i3] += (tempConc[subInd][i1][yUp][i3]-tempConc[subInd][i1][i2][i3])*D/6;
          }
          if (yDown>=0) {
            Conc[subInd][i1][i2][i3] += (tempConc[subInd][i1][yDown][i3]-tempConc[subInd][i1][i2][i3])*D/6;
          }
          if (zUp<L) {
            Conc[subInd][i1][i2][i3] += (tempConc[subInd][i1][i2][zUp]-tempConc[subInd][i1][i2][i3])*D/6;
          }
          if (zDown>=0) {
            Conc[subInd][i1][i2][i3] += (tempConc[subInd][i1][i2][zDown]-tempConc[subInd][i1][i2][i3])*D/6;
          }
        }
      }
    }
  }
}

TEST(DiffusionOp, CompareWithNewCastleDiffusionTest) {

  double diffusion_constant = 0.1;
  int iteration_count = 15;
  DiffusionOp op(iteration_count, diffusion_constant);

  int L = 128;
  AdaptiveMesh<float> m({0, 0, 0, 1.*L, 1.*L, 1.*L}, 0, 7);

  // create 3D concentration matrix
  double**** Conc;
  double**** tempConc;
  Conc = new double***[L];
  tempConc = new double***[L];
  for (int i1 = 0; i1 < 2; i1++) {
    Conc[i1] = new double**[L];
    tempConc[i1] = new double**[L];
    for (int i2 = 0; i2 < L; i2++) {
      Conc[i1][i2] = new double*[L];
      tempConc[i1][i2] = new double*[L];
      for (int i3 = 0; i3 < L; i3++) {
        Conc[i1][i2][i3] = new double[L];
        tempConc[i1][i2][i3] = new double[L];
        for (int i4 = 0; i4 < L; i4++) {
          Conc[i1][i2][i3][i4] = 0;
          tempConc[i1][i2][i3][i4] = 0;
        }
      }
    }
  }

  Point p = {L/2., L/2., L/2.};
  m.Put(p, 1);

  Conc[0][L/2][L/2][L/2] = 1;

  op.Compute(&m);

  for (int i = 0; i < iteration_count; i++) {
    runDiffusionStep(Conc, tempConc, L, diffusion_constant);
  }

  for (int x = 0; x < L; x++) {
    for (int y = 0; y < L; y++) {
      for (int z = 0; z < L; z++) {
        auto diff = Conc[0][x][y][z] - m.At({x*1., y*1., z*1.});

//        ASSERT_TRUE( std::abs(diff) < Param::kEpsilon);
        ASSERT_EQ(Conc[0][x][y][z], m.At({x*1., y*1., z*1.}));
      }
    }
  }


  for (int i1 = 0; i1 < 2; i1++) {
    for (int i2 = 0; i2 < L; i2++) {
      for (int i3 = 0; i3 < L; i3++) {
        delete[] Conc[i1][i2][i3];
        delete[] tempConc[i1][i2][i3];
      }
      delete[] Conc[i1][i2];
      delete[] tempConc[i1][i2];
    }
    delete[] Conc[i1];
    delete[] tempConc[i1];
  }
  delete[] Conc;
  delete[] tempConc;

}


TEST(DiffusionOp, LittleManualTest) {

  double diffusion_constant = 0.1;
  DiffusionOp op(1, diffusion_constant);

  AdaptiveMesh<double> m({0, 0, 0, 1024, 1024, 1024}, 0, 10);

  Point p = {512, 512, 512};
  m.Put(p, 1);

  op.Compute(&m);
  ASSERT_EQ(m.At({512, 512, 512}), 1 - diffusion_constant);

  ASSERT_EQ(m.At({511, 512, 512}), diffusion_constant/6);
  ASSERT_EQ(m.At({513, 512, 512}), diffusion_constant/6);
  ASSERT_EQ(m.At({512, 511, 512}), diffusion_constant/6);
  ASSERT_EQ(m.At({512, 513, 512}), diffusion_constant/6);
  ASSERT_EQ(m.At({512, 512, 511}), diffusion_constant/6);
  ASSERT_EQ(m.At({512, 512, 513}), diffusion_constant/6);

}

} // namespace diffusion_op_test_internal
}  // namespace bdm