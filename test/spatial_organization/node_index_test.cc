#include <gtest/gtest.h>
#include <ctime>
#include <random>
#include <algorithm>

#include "spatial_organization/node_index.h"

namespace bdm {
namespace spatial_organization {

TEST(NodeIndexTest, AdjacentIndexAtSameLevelTest2) {

  NodeIndex index(5, 5, 5, 5);

  std::vector<NodeIndex> n;

  n.push_back(NodeIndex(4, 5, 5, 5));
  n.push_back(NodeIndex(6, 5, 5, 5));

  n.push_back(NodeIndex(5, 4, 5, 5));
  n.push_back(NodeIndex(5, 6, 5, 5));

  n.push_back(NodeIndex(5, 5, 4, 5));
  n.push_back(NodeIndex(5, 5, 6, 5));


  auto neighbors = NodeIndex::GetAdjacentIndex<false>(index, 5);

  auto comparator = [](NodeIndex a, NodeIndex b) {
    return a.code < b.code;
  };
  std::sort(n.begin(), n.end(), comparator);
  std::sort(neighbors.begin(), neighbors.end(), comparator);

  for (int i = 0; i < neighbors.size(); i++) {
    ASSERT_EQ(n[i].code,neighbors[i].code);
  }

  neighbors = NodeIndex::GetAdjacentIndex<true>(index, 5);

  n.push_back(NodeIndex(4, 4, 4, 5));
  n.push_back(NodeIndex(4, 4, 5, 5));
  n.push_back(NodeIndex(4, 4, 6, 5));

  n.push_back(NodeIndex(4, 5, 4, 5));
  n.push_back(NodeIndex(4, 5, 6, 5));

  n.push_back(NodeIndex(4, 6, 4, 5));
  n.push_back(NodeIndex(4, 6, 5, 5));
  n.push_back(NodeIndex(4, 6, 6, 5));

  n.push_back(NodeIndex(5, 4, 4, 5));
  n.push_back(NodeIndex(5, 4, 6, 5));

  n.push_back(NodeIndex(5, 6, 4, 5));
  n.push_back(NodeIndex(5, 6, 6, 5));

  n.push_back(NodeIndex(6, 4, 4, 5));
  n.push_back(NodeIndex(6, 4, 5, 5));
  n.push_back(NodeIndex(6, 4, 6, 5));

  n.push_back(NodeIndex(6, 5, 4, 5));
  n.push_back(NodeIndex(6, 5, 6, 5));

  n.push_back(NodeIndex(6, 6, 4, 5));
  n.push_back(NodeIndex(6, 6, 5, 5));
  n.push_back(NodeIndex(6, 6, 6, 5));

  std::sort(n.begin(), n.end(), comparator);
  std::sort(neighbors.begin(), neighbors.end(), comparator);

  for (int i = 0; i < neighbors.size(); i++) {
    ASSERT_EQ(n[i].code,neighbors[i].code);
  }


}
TEST(NodeIndexTest, AdjacentIndexTest1) {
  NodeIndex index(1, 1, 1, 4);


  ASSERT_EQ(index.AdjacentX(-1).code, NodeIndex(0, 1, 1, 4).code);
  ASSERT_EQ(index.AdjacentX( 1).code, NodeIndex(2, 1, 1, 4).code);
  ASSERT_EQ(index.AdjacentY(-1).code, NodeIndex(1, 0, 1, 4).code);
  ASSERT_EQ(index.AdjacentY( 1).code, NodeIndex(1, 2, 1, 4).code);
  ASSERT_EQ(index.AdjacentZ(-1).code, NodeIndex(1, 1, 0, 4).code);
  ASSERT_EQ(index.AdjacentZ( 1).code, NodeIndex(1, 1, 2, 4).code);

  ASSERT_EQ(index.AdjacentZ( 2).code, NodeIndex(1, 1, 3, 4).code);

//  NodeIndex n1(1, 1, 0, 4);
//
//
//
//
//
//  NodeIndex n2(0, 1, 0, 4);
//  NodeIndex n3(1, 0, 0, 4);
//
//  std::vector<NodeIndex> n;
//
//  n.push_back(n1);
//  n.push_back(n2);
//  n.push_back(n3);
//
//  auto neighbors = NodeIndex::GetAdjacentIndex<false>(index, 3);
//
//  auto comparator = [](NodeIndex a, NodeIndex b) {
//    return a.code < b.code;
//  };
//  std::sort(n.begin(), n.end(), comparator);
//  std::sort(neighbors.begin(), neighbors.end(), comparator);
//
//  for (int i = 0; i < neighbors.size(); i++) {
//    ASSERT_EQ(n[i].code,neighbors[i].code);
//  }
//
//
//  neighbors = NodeIndex::GetAdjacentIndex<true>(index, 3);
//
//  n.push_back(NodeIndex(0, 1, 1, 3));
//  n.push_back(NodeIndex(1, 0, 1, 3));
//  n.push_back(NodeIndex(1, 1, 0, 3));
//  n.push_back(NodeIndex(1, 1, 1, 3));
//
//  std::sort(n.begin(), n.end(), comparator);
//  std::sort(neighbors.begin(), neighbors.end(), comparator);
//
//  for (int i = 0; i < neighbors.size(); i++) {
//    ASSERT_EQ(n[i].code, neighbors[i].code);
//  }

}


TEST(NodeIndexTest, AdjacentIndexAtSameLevelTest1) {
  NodeIndex index(0, 0, 0, 3);

  NodeIndex n1(0, 0, 1, 3);
  NodeIndex n2(0, 1, 0, 3);
  NodeIndex n3(1, 0, 0, 3);

  std::vector<NodeIndex> n;

  n.push_back(n1);
  n.push_back(n2);
  n.push_back(n3);

  auto neighbors = NodeIndex::GetAdjacentIndex<false>(index, 3);

  auto comparator = [](NodeIndex a, NodeIndex b) {
    return a.code < b.code;
  };
  std::sort(n.begin(), n.end(), comparator);
  std::sort(neighbors.begin(), neighbors.end(), comparator);

  for (int i = 0; i < neighbors.size(); i++) {
    ASSERT_EQ(n[i].code,neighbors[i].code);
  }


  neighbors = NodeIndex::GetAdjacentIndex<true>(index, 3);

  n.push_back(NodeIndex(0, 1, 1, 3));
  n.push_back(NodeIndex(1, 0, 1, 3));
  n.push_back(NodeIndex(1, 1, 0, 3));
  n.push_back(NodeIndex(1, 1, 1, 3));

  std::sort(n.begin(), n.end(), comparator);
  std::sort(neighbors.begin(), neighbors.end(), comparator);

  for (int i = 0; i < neighbors.size(); i++) {
    ASSERT_EQ(n[i].code, neighbors[i].code);
  }

}



TEST(NodeIndexTest, LevelTest) {
  NodeIndex index1(0, 0, 0, 3);
  NodeIndex index2(0, 1, 0, 5);
  NodeIndex index3(0, 0, 4, 20);

  ASSERT_EQ(0,  NodeIndex().Level());
  ASSERT_EQ(3,  index1.Level());
  ASSERT_EQ(5,  index2.Level());
  ASSERT_EQ(20, index3.Level());
}

TEST(NodeIndexTest, DilateTest) {
  std::srand(time(nullptr));

  auto coord = std::rand()% (1<<20);
  auto code = NodeIndex::oct_dilate(coord);

  for (int i = 0; i < 21; i++) {
    auto coord_bit = coord >> i & 1;
    auto code_bit  = (code >> (3l*i)) & 1;
    ASSERT_EQ(coord_bit, code_bit);
  }
}

TEST(NodeIndexTest, ContractTest) {
  std::srand(time(nullptr));

  auto code = std::rand() % (1l<<63);
  auto x = NodeIndex::oct_contract(code);
  auto y = NodeIndex::oct_contract(code >> 1);
  auto z = NodeIndex::oct_contract(code >> 2);

  for (int i = 0; i < 21; i++) {
    auto x_bit  = (code >> (3*i    )) & 1;
    auto y_bit  = (code >> (3*i + 1)) & 1;
    auto z_bit  = (code >> (3*i + 2)) & 1;

    ASSERT_EQ(x_bit, (x >> i) & 1);
    ASSERT_EQ(y_bit, (y >> i) & 1);
    ASSERT_EQ(z_bit, (z >> i) & 1);
  }
}

}  // namespace spatial_organization
}  // namespace bdm
