#include "gtest/gtest.h"
#include "backend.h"
#include "simulation_object.h"

namespace bdm {

// The following tests check if code insertion in new classes works as intended
// Therefore BdmSimObject is extended in two stages: first by CellExt and
// then by NeuronExt

template <typename Base = BdmSimObject<>>
class CellExt : public Base {
  BDM_CLASS_HEADER(CellExt, CellExt<>,
                   CellExt<typename Base::template Self<Backend>>, position_,
                   diameter_);

 public:
  explicit CellExt(const std::array<real_v, 3>& pos) : position_{{pos}} {}

  CellExt() : position_{{1, 2, 3}} {}

  const std::array<real_v, 3>& GetPosition() const { return position_[idx_]; }
  const real_v& GetDiameter() const { return diameter_[idx_]; }

  void SetDiameter(const real_v& diameter) { diameter_[idx_] = diameter; }

 protected:
  BDM_PROTECTED_MEMBER(Container<std::array<real_v COMMA() 3>>, position_);
  BDM_PROTECTED_MEMBER(Container<real_v>, diameter_) = {real_v(6.28)};
};

// -----------------------------------------------------------------------------
// libraries for specific specialities add functionality - e.g. Neuroscience
class Neurite {
 public:
  Neurite() : id_(0) {}
  explicit Neurite(std::size_t id) : id_(id) {}
  std::size_t id_;
};

// add Neurites to BaseCell
template <typename Base = CellExt<>>
class NeuronExt : public Base {
  BDM_CLASS_HEADER(NeuronExt, NeuronExt<>,
                   NeuronExt<typename Base::template Self<Backend>>, neurites_);

 public:
  template <class... A>
  explicit NeuronExt(const SimdArray<std::vector<Neurite>>& neurites,
                     const A&... a)
      : Base(a...) {
    neurites_[idx_] = neurites;
  }

  NeuronExt() = default;

  const SimdArray<std::vector<Neurite>>& GetNeurites() const {
    return neurites_[idx_];
  }

 private:
  BDM_PRIVATE_MEMBER(Container<SimdArray<std::vector<Neurite>>>,
                     neurites_) = {{}};

  FRIEND_TEST(SimulationObjectUtilTest, SoaBackend_clear);
  FRIEND_TEST(SimulationObjectUtilTest, SoaBackend_reserve);
  FRIEND_TEST(SimulationObjectUtilTest,
              SoaBackend_push_backScalarOnNonEmptySoa);
};

// define easy to use templated type alias
template <typename Backend = VcBackend>
using Neuron = NeuronExt<CellExt<BdmSimObject<SelectAllMembers, Backend>>>;

TEST(SimulationObjectUtilTest, DefaultConstructor) {
  // are data members in all extensions correctly initialized?
  Neuron<VcBackend> neuron;

  EXPECT_TRUE((VcBackend::real_v(6.28) == neuron.GetDiameter()).isFull());
  auto positions = neuron.GetPosition();
  EXPECT_TRUE((VcBackend::real_v(1) == positions[0]).isFull());
  EXPECT_TRUE((VcBackend::real_v(2) == positions[1]).isFull());
  EXPECT_TRUE((VcBackend::real_v(3) == positions[2]).isFull());

  auto neurites_array = neuron.GetNeurites();
  for (auto& neurites : neurites_array) {
    EXPECT_EQ(0u, neurites.size());
  }
}

TEST(SimulationObjectUtilTest, NonDefaultConstructor) {
  // are data members in all extensions correctly initialized?
  using real_v = VcBackend::real_v;

  std::vector<Neurite> neurites;
  neurites.push_back(Neurite(2));
  neurites.push_back(Neurite(3));
  VcBackend::SimdArray<std::vector<Neurite>> neurite_v;
  for (std::size_t i = 0; i < neurite_v.size(); i++) neurite_v[i] = neurites;

  Neuron<VcBackend> neuron(
      neurite_v, std::array<real_v, 3>{real_v(4), real_v(5), real_v(6)});

  EXPECT_TRUE((VcBackend::real_v(6.28) == neuron.GetDiameter()).isFull());
  auto positions = neuron.GetPosition();
  EXPECT_TRUE((VcBackend::real_v(4) == positions[0]).isFull());
  EXPECT_TRUE((VcBackend::real_v(5) == positions[1]).isFull());
  EXPECT_TRUE((VcBackend::real_v(6) == positions[2]).isFull());

  auto& neurites_array = neuron.GetNeurites();
  for (auto& neurites : neurites_array) {
    EXPECT_EQ(2u, neurites.size());
  }
}

TEST(SimulationObjectUtilTest, NewScalar) {
  using real_v = ScalarBackend::real_v;
  auto neuron = Neuron<VcBackend>::NewScalar();

  EXPECT_TRUE((real_v(6.28) == neuron.GetDiameter()).isFull());
  auto positions = neuron.GetPosition();
  EXPECT_TRUE((real_v(1) == positions[0]).isFull());
  EXPECT_TRUE((real_v(2) == positions[1]).isFull());
  EXPECT_TRUE((real_v(3) == positions[2]).isFull());

  auto neurites_array = neuron.GetNeurites();
  for (auto& neurites : neurites_array) {
    EXPECT_EQ(0u, neurites.size());
  }
}

TEST(SimulationObjectUtilTest, NewEmptySoa) {
  auto neurons = Neuron<>::NewEmptySoa();
  neurons.size();
  EXPECT_EQ(0u, neurons.size());
  EXPECT_EQ(0u, neurons.vectors());
  EXPECT_EQ(0u, neurons.elements());
}

TEST(SimulationObjectUtilTest, GetSoaRef) {
  using real_v = VcBackend::real_v;
  Neuron<VcSoaBackend> neurons;
  auto neurons_ref = neurons.GetSoaRef();
  neurons_ref.SetDiameter(real_v(12.34));
  EXPECT_TRUE((real_v(12.34) == neurons.GetDiameter()).isFull());
}

TEST(SimulationObjectUtilTest,
     SoaBackend_push_backVector_AndSubscriptOperator) {
  using real_v = VcBackend::real_v;

  std::vector<Neurite> neurites;
  neurites.push_back(Neurite(2));
  neurites.push_back(Neurite(3));
  VcBackend::SimdArray<std::vector<Neurite>> neurite_v;
  for (std::size_t i = 0; i < neurite_v.size(); i++) neurite_v[i] = neurites;

  Neuron<VcBackend> neuron_v1(
      neurite_v, std::array<real_v, 3>{real_v(4), real_v(5), real_v(6)});

  neurites.push_back(Neurite(4));
  for (std::size_t i = 0; i < neurite_v.size(); i++) neurite_v[i] = neurites;
  Neuron<VcBackend> neuron_v2(
      neurite_v, std::array<real_v, 3>{real_v(9), real_v(8), real_v(7)});

  auto neurons = Neuron<>::NewEmptySoa();
  neurons.push_back(neuron_v1);
  neurons.push_back(neuron_v2);

  // switch to neuron_v2
  auto& result1 = neurons[1];
  std::cout << typeid(result1).name() << std::endl;

  // operator[] returns reference to *this
  EXPECT_TRUE(&result1 == &neurons);

  EXPECT_TRUE((VcBackend::real_v(6.28) == result1.GetDiameter()).isFull());
  auto& positions = result1.GetPosition();
  EXPECT_TRUE((VcBackend::real_v(9) == positions[0]).isFull());
  EXPECT_TRUE((VcBackend::real_v(8) == positions[1]).isFull());
  EXPECT_TRUE((VcBackend::real_v(7) == positions[2]).isFull());

  auto neurite_v_actual = result1.GetNeurites();
  for (auto& neurites_s : neurite_v_actual) {
    EXPECT_EQ(3u, neurites_s.size());
  }
}

TEST(SimulationObjectUtilTest, SoaBackend_push_backScalarOnEmptySoa) {
  auto neurons = Neuron<>::NewEmptySoa();
  EXPECT_EQ(0u, neurons.size());

  using real_v = ScalarBackend::real_v;
  std::vector<Neurite> neurites;
  neurites.push_back(Neurite(22));
  neurites.push_back(Neurite(33));
  Neuron<ScalarBackend> single_neuron(
      {neurites}, std::array<real_v, 3>{real_v(6), real_v(3), real_v(9)});
  single_neuron.SetDiameter({1.2345});

  neurons.push_back(single_neuron);
  EXPECT_EQ(1u, neurons.size());
  EXPECT_EQ(1u, neurons.vectors());
  EXPECT_EQ(1u, neurons.elements());

  // check if scalar data members have been copied correctly
  EXPECT_EQ(1.2345, neurons.GetDiameter()[0]);
  auto& position = neurons.GetPosition();
  EXPECT_EQ(6, position[0][0]);
  EXPECT_EQ(3, position[1][0]);
  EXPECT_EQ(9, position[2][0]);
  auto& neurites_actual = neurons.GetNeurites()[0];
  EXPECT_EQ(2u, neurites_actual.size());
  EXPECT_EQ(22u, neurites_actual[0].id_);
  EXPECT_EQ(33u, neurites_actual[1].id_);
}

TEST(SimulationObjectUtilTest, SoaBackend_push_backScalarOnNonEmptySoa) {
  Neuron<VcSoaBackend> neurons;  // stores one vector neuron with default values
  EXPECT_EQ(1u, neurons.size());
  EXPECT_EQ(1u, neurons.vectors());
  auto expected_elements = VcBackend::kVecLen;
  EXPECT_EQ(expected_elements, neurons.elements());

  // simulate that the first vector only holds one scalar
  neurons.size_last_vector_ = 1;

  using real_v = ScalarBackend::real_v;
  std::vector<Neurite> neurites;
  neurites.push_back(Neurite(22));
  neurites.push_back(Neurite(33));
  Neuron<ScalarBackend> single_neuron(
      {neurites}, std::array<real_v, 3>{real_v(6), real_v(3), real_v(9)});
  single_neuron.SetDiameter({1.2345});

  neurons.push_back(single_neuron);
  EXPECT_EQ(1u, neurons.size());
  EXPECT_EQ(1u, neurons.vectors());
  EXPECT_EQ(2u, neurons.elements());

  // check if scalar data members have been copied correctly
  EXPECT_EQ(1.2345, neurons.GetDiameter()[1]);
  auto& position = neurons.GetPosition();
  EXPECT_EQ(6, position[0][1]);
  EXPECT_EQ(3, position[1][1]);
  EXPECT_EQ(9, position[2][1]);
  auto& neurites_actual = neurons.GetNeurites()[1];
  EXPECT_EQ(2u, neurites_actual.size());
  EXPECT_EQ(22u, neurites_actual[0].id_);
  EXPECT_EQ(33u, neurites_actual[1].id_);
}

TEST(SimulationObjectUtilTest, SoaBackend_clear) {
  Neuron<VcSoaBackend> neurons;
  EXPECT_EQ(1u, neurons.size());
  neurons.clear();
  EXPECT_EQ(0u, neurons.size());
  EXPECT_EQ(0u, neurons.neurites_.size());
  EXPECT_EQ(0u, neurons.diameter_.size());
  EXPECT_EQ(0u, neurons.position_.size());
}

TEST(SimulationObjectUtilTest, SoaBackend_reserve) {
  Neuron<VcSoaBackend> neurons;
  neurons.reserve(10);
  EXPECT_EQ(10u, neurons.neurites_.capacity());
  EXPECT_EQ(10u, neurons.diameter_.capacity());
  EXPECT_EQ(10u, neurons.position_.capacity());
}

TEST(SimulationObjectUtilTest, SoaBackend_Gather) {
  auto objects = Neuron<>::NewEmptySoa();

  // create objects
  for (size_t i = 0; i < 10; i++) {
    using real_v = ScalarBackend::real_v;
    std::vector<Neurite> neurites;
    neurites.push_back(Neurite(i));
    Neuron<ScalarBackend> scalar(
        {neurites}, std::array<real_v, 3>{real_v(i), real_v(i), real_v(i)});
    scalar.SetDiameter(i);
    objects.push_back(scalar);
  }

  aosoa<Neuron<VcBackend>, VcBackend> gathered;
  InlineVector<int, 8> indexes;
  indexes.push_back(5);
  indexes.push_back(3);
  indexes.push_back(9);
  indexes.push_back(2);
  objects.Gather(indexes, &gathered);
  // check if it returns the correct objects
  size_t target_n_vectors =
      4 / VcBackend::kVecLen + (4 % VcBackend::kVecLen ? 1 : 0);
  EXPECT_EQ(target_n_vectors, gathered.vectors());
  size_t counter = 0;
  for (size_t i = 0; i < gathered.vectors(); i++) {
    for (size_t j = 0; j < VcBackend::kVecLen; j++) {
      EXPECT_EQ(indexes[counter], gathered[i].GetDiameter()[j]);
      EXPECT_EQ(indexes[counter], gathered[i].GetPosition()[0][j]);
      EXPECT_EQ(indexes[counter], gathered[i].GetPosition()[1][j]);
      EXPECT_EQ(indexes[counter], gathered[i].GetPosition()[2][j]);
      EXPECT_EQ((unsigned)indexes[counter],
                gathered[i].GetNeurites()[j][0].id_);
      counter++;
    }
  }
}

TEST(SimulationObjectUtilTest, VectorBackend_push_backScalar) {
  Neuron<VcBackend> neurons;  // stores one vector neuron with default values
  EXPECT_EQ(4u, neurons.size());  // replace with VcBackend::kVecLen
  EXPECT_EQ(1u, neurons.vectors());
  auto expected_elements = VcBackend::kVecLen;
  EXPECT_EQ(expected_elements, neurons.elements());

  // simulate that the vector only holds one scalar - remaining slots are free
  neurons.SetSize(1);

  using real_v = ScalarBackend::real_v;
  std::vector<Neurite> neurites;
  neurites.push_back(Neurite(22));
  neurites.push_back(Neurite(33));
  Neuron<ScalarBackend> single_neuron(
      {neurites}, std::array<real_v, 3>{real_v(6), real_v(3), real_v(9)});
  single_neuron.SetDiameter({1.2345});

  neurons.push_back(single_neuron);
  EXPECT_EQ(2u, neurons.size());
  EXPECT_EQ(1u, neurons.vectors());
  EXPECT_EQ(2u, neurons.elements());

  // check if scalar data members have been copied correctly
  EXPECT_EQ(1.2345, neurons.GetDiameter()[1]);
  auto& position = neurons.GetPosition();
  EXPECT_EQ(6, position[0][1]);
  EXPECT_EQ(3, position[1][1]);
  EXPECT_EQ(9, position[2][1]);
  auto& neurites_actual = neurons.GetNeurites()[1];
  EXPECT_EQ(2u, neurites_actual.size());
  EXPECT_EQ(22u, neurites_actual[0].id_);
  EXPECT_EQ(33u, neurites_actual[1].id_);
}

// TODO(lukas) test assignment operator

}  // namespace bdm