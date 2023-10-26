#include "mesher.h"

#include <string>
#include <unordered_set>

#include "geometry.h"
#include "log.h"
#include "mesh.h"
#include "stlext.h"

#define TETLIBRARY
#define REAL double
#include "tetgen.h"

namespace helix {

TetGenMesher::TetGenMesher(Mesh& surface) : mesh_(surface) {}

void TetGenMesher::generate(const MeshingParameters& params) {
  mesh_.tetrahedra().clear(false);

  // -p tetrahedralizes a piecewise linear complex (PLC).
  // -Y preserves the input surface mesh (does not modify it).
  // -V verbose
  // -Q quiet
  // -q mesh quality (maximum radius-edge ratio)/(minimum dihedral angle)
  // -a maximum volume constraint
  // -f provides the interior + boundary triangular faces
  // -nn to get tet neighbors for each triangular face
  // -m applies a mesh sizing function
  // -C checks the consistency of the final mesh.
  std::string options = "z";
  if (params.constrained) options += "pY";
  options += (params.verbose) ? "V" : "Q";
  if (params.quality > 0 && params.constrained) {
    options += "q" + std::to_string(params.quality).substr(0, 3);
  } else {
    options += "O0";
  }
  options += "fnn";
  LOG << "running tetgen with options: " << options;

  // all memory in the tetgenio objects will be freed when they are
  // destructed (see the "clean_memory" function in tetgen.h)
  tetgenio input, output;
  tetgenio::facet* f;
  tetgenio::polygon* p;

  input.numberofpoints = mesh_.vertices().n();
  input.pointlist = new REAL[3 * input.numberofpoints];
  input.pointmarkerlist = new int[input.numberofpoints];

  for (int k = 0; k < mesh_.vertices().n(); k++) {
    for (int d = 0; d < 3; d++)
      input.pointlist[3 * k + d] = mesh_.vertices()[k][d];
    input.pointmarkerlist[k] = mesh_.vertices().group(k);
  }

  // set the constraints
  if (params.constrained) {
    input.numberoffacets = mesh_.triangles().n();
    input.facetlist = new tetgenio::facet[input.numberoffacets];
    for (index_t k = 0; k < mesh_.triangles().n(); k++) {
      f = &input.facetlist[k];
      f->numberofpolygons = 1;
      f->polygonlist = new tetgenio::polygon[1];
      f->numberofholes = 0;
      f->holelist = nullptr;

      p = &f->polygonlist[0];
      p->numberofvertices = 3;
      p->vertexlist = new int[3];
      for (int i = 0; i < 3; i++)
        p->vertexlist[i] = mesh_.triangles()(k, i);
    }
  } else {
    input.numberoffacets = 0;
  }

  // tetrahedralize
  LOG << "tetrahedralizing...";
  tetrahedralize((char*)options.c_str(), &input, &output);
  LOG << "done";

  // add new points if necessary
  if (output.numberofpoints != mesh_.vertices().n()) {
    ASSERT(output.numberofpoints > mesh_.vertices().n())
        << output.numberofpoints;
    mesh_.vertices().reserve(output.numberofpoints);
    for (index_t i = mesh_.vertices().n(); i < output.numberofpoints; i++)
      mesh_.vertices().add(output.pointlist + 3 * i);
  }

  // create a dictionary of all the triangles
  LOG << "creating dictionary of surface triangles";
  std::unordered_map<std::array<int, 3>, index_t> triangle_map;
  triangle_map.reserve(mesh_.triangles().n());
  std::array<int, 3> t;
  for (int k = 0; k < mesh_.triangles().n(); k++) {
    for (int d = 0; d < 3; d++)
      t[d] = mesh_.triangles()[k][d];
    triangle_map.insert({t, k});
  }
  auto is_surface_triangle = [&triangle_map](int f0, int f1, int f2) {
    if (triangle_map.find({f0, f1, f2}) != triangle_map.end()) return 1;
    if (triangle_map.find({f2, f0, f1}) != triangle_map.end()) return 2;
    if (triangle_map.find({f1, f2, f0}) != triangle_map.end()) return 3;
    if (triangle_map.find({f1, f0, f2}) != triangle_map.end()) return -1;
    if (triangle_map.find({f2, f1, f0}) != triangle_map.end()) return -2;
    if (triangle_map.find({f0, f2, f1}) != triangle_map.end()) return -3;
    return 0;
  };

  // find one tet that is on the boundary
  LOG << "determining which tetrahedra to keep";
  int* neighbors = output.neighborlist;
  int tbnd = -1;
  for (int k = 0; k < output.numberoftetrahedra; k++) {
    for (int j = 0; j < 4; j++) {
      if (neighbors[4 * k + j] < 0) {
        tbnd = k;
        break;
      }
    }
    if (tbnd >= 0) break;
  }
  ASSERT(tbnd >= 0);

  std::unordered_set<int> visited;
  visited.reserve(output.numberoftetrahedra);
  std::vector<int> stack;
  stack.reserve(64);
  stack.push_back(tbnd);
  while (!stack.empty()) {
    int t = stack.back();
    visited.insert(t);
    stack.pop_back();
    for (int j = 0; j < 4; j++) {
      int n = neighbors[4 * t + j];
      if (n < 0) continue;  // exterior boundary
      if (visited.find(n) != visited.end()) continue;
      // skip faces that are part of the surface triangles
      // TODO: determine if this is an interior surface in which
      // we want tetrahedra on both sides
      int f = output.tet2facelist[4 * t + j];
      int face[3];
      for (int d = 0; d < 3; d++)
        face[d] = output.trifacelist[3 * f + d];
      if (is_surface_triangle(face[0], face[1], face[2]) != 0) continue;
      stack.push_back(n);
    }
  }

  // create tetrahedra
  index_t tet[4];
  mesh_.tetrahedra().reserve(visited.size());
  for (int k = 0; k < output.numberoftetrahedra; k++) {
    if (visited.find(k) == visited.end()) continue;
    for (int j = 0; j < 4; j++)
      tet[j] = output.tetrahedronlist[4 * k + j];
    mesh_.tetrahedra().add(tet);
  }
  LOG << fmt::format(
      "created tetrahedralization with {} tetrahedra (discarded {})",
      mesh_.tetrahedra().n(), output.numberoftetrahedra - visited.size());
}

}  // namespace helix