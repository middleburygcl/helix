#include "library.h"

#include "io.h"
#include "tester.h"

using namespace helix;

UT_TEST_SUITE(library_test_suite)

UT_TEST_CASE(grid_triangle_test) {}
UT_TEST_CASE_END(grid_triangle_test)

UT_TEST_CASE(polygon_test) {
  int n = 10;
  Grid<Polygon> mesh({n, n});
  meshb::write(mesh, "results/polygrid.meshb");
}
UT_TEST_CASE_END(polygon_test)

UT_TEST_CASE(polyhedron_test) {
  int n = 10;
  Grid<Polyhedron> mesh({n, n, n});
  meshb::write(mesh, "results/polyhedra.meshb");
}
UT_TEST_CASE_END(polyhedron_test)

UT_TEST_CASE(sphere_test) {
  Sphere mesh(2);
  mesh.add_prism_layers(1.5, 5, 2.0);
  meshb::write(mesh, "results/sphere.meshb");
}
UT_TEST_CASE_END(sphere_test)

UT_TEST_SUITE_END(library_test_suite)
