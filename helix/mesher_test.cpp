#include "mesher.h"

#include "geometry.h"
#include "io.h"
#include "tester.h"

using namespace helix;

UT_TEST_SUITE(mesher_test_suite)

UT_TEST_CASE(test1) {
#if 0
  Mesh mesh(3);
  meshb::read("bunny.meshb", mesh);
  meshb::read("tire.meshb", mesh);
  TetGenMesher mesher(mesh);
  MeshingParameters params{.constrained = true};
  params.quality = 1.1;
  mesher.generate(params);
  meshb::write(mesh, "results/tire-vol.meshb");
#else
  std::string filename = "../tools/csmextract/data/bottle.egads";
  // filename = "DPW6/DPW6_CRM_wbnpt_ih+0_v09_2016-01-28_cf.stp";
  Geometry geometry(GeometryKernel::kEgads, filename);
  Operation op;
  op.type = OperationType::kSubtraction;
  op.centering = CenteringMethod::kModel;
  op.scaling = LengthScaling::kRelative;
  op.shape = OperationShape::kSphere;
  op.set_radius(1.0);
  // geometry.operate(op);

  Mesh mesh(3);
  geometry.tessellate(mesh);
  meshb::write(mesh, "results/test-surface.meshb");

  TetGenMesher mesher(mesh);
  MeshingParameters params{.constrained = true};
  params.verbose = false;
  // params.quality = 1.6;
  mesher.generate(params);
  meshb::write(mesh, "results/test.meshb");
#endif
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(mesher_test_suite)