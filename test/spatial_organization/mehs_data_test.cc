#include <gtest/gtest.h>

#include "spatial_organization/mesh_data.h"
#include "spatial_organization/adaptive_mesh.h"


namespace bdm {
namespace spatial_organization {

TEST(MeshDataTest, DefaultValueTest) {
  MeshData<double> m1({0, 0, 0, 100, 100, 100}, 0, 10);

  ASSERT_EQ(m1.At({50, 50, 50}), 0);

  MeshData<double> m2({0, 0, 0, 100, 100, 100}, 5, 10);

  ASSERT_EQ(m2.At({50, 50, 50}), 5);
}

TEST(MeshDataTest, InsertAndFoundTest) {
  MeshData<double> m({0, 0, 0, 100, 100, 100}, 0, 10);

  m.Put({50, 50, 50},5);
  ASSERT_EQ(m.At({50, 50, 50}), 5);

  m.Put({50, 50, 50},6);

  ASSERT_EQ(m.At({50, 51, 50}), 0);
  ASSERT_EQ(m.At({50, 50, 51}), 0);
  ASSERT_EQ(m.At({51, 50, 50}), 0);
  ASSERT_EQ(m.At({50, 50, 50}), 6);

}


}  // namespace spatial_organization
}  // namespace bdm
