#include "geometry.h"

#include "io.h"
#include "mesh.h"
#include "tester.h"

using namespace helix;

UT_TEST_SUITE(geometry_test_suite)

UT_TEST_CASE(test_egads) {
  std::string filename = "../tools/csmextract/data/gallery/Hypersonic.egads";
  filename = "tire.egads";
  Geometry geometry(GeometryKernel::kEgads, filename);
  UT_ASSERT_EQUALS(geometry.n_bodies(), 1);
  UT_ASSERT_EQUALS(geometry.n_entities(), 123);

  Operation op;
  op.type = OperationType::kSubtraction;
  op.centering = CenteringMethod::kModel;
  op.scaling = LengthScaling::kRelative;
  op.shape = OperationShape::kSphere;
  op.set_radius(2.0);
  geometry.operate(op);

  Mesh mesh(3);
  geometry.tessellate(mesh);
  meshb::write(mesh, "tire.meshb");
}
UT_TEST_CASE_END(test_egads)

UT_TEST_SUITE_END(geometry_test_suite)