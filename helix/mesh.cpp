#include "mesh.h"

#include <set>

#include "elements.h"
#include "geometry.h"

namespace helix {

template <typename T>
void get_element_edges(const Topology<T>& elems,
                       int k,
                       std::vector<Edge>& edges) {
  for (int j = 0; j < T::n_edges; j++) {
    auto e0 = elems(k, T::edges[2 * j]);
    auto e1 = elems(k, T::edges[2 * j + 1]);

    if (e0 > e1) std::swap(e0, e1);
    edges.push_back({e0, e1});
  }
}

template <>
void get_element_edges<Polygon>(const Topology<Polygon>& elems,
                                int k,
                                std::vector<Edge>& edges) {
  auto n_vertices = elems.length(k);
  for (int j = 0; j < n_vertices; j++) {
    auto e0 = elems(k, j);
    auto e1 = elems(k, (j + 1) % n_vertices);

    if (e0 > e1) std::swap(e0, e1);
    edges.push_back({e0, e1});
  }
}

template <typename T>
void Topology<T>::append_edges(std::vector<Edge>& edges) const {
  std::set<Edge> E;
  std::vector<Edge> element_edges;

  for (int k = 0; k < n(); k++) {
    element_edges.clear();
    get_element_edges(*this, k, element_edges);
    for (const auto& e : element_edges) {
      if (E.find(e) == E.end()) {
        E.insert(e);
        edges.push_back(e);
      }
    }
  }
}

template <>
void Topology<Triangle>::flip_orientation() {
  for (int k = 0; k < n(); k++) {
    index_t t1 = (*this)(k, 1);
    index_t t2 = (*this)(k, 2);
    (*this)(k, 2) = t1;
    (*this)(k, 1) = t2;
  }
}

void Mesh::get_edges(std::vector<Edge>& edges) const {
  edges.clear();
  triangles_.append_edges(edges);
  tetrahedra_.append_edges(edges);
  polygons_.append_edges(edges);
}

template <>
const Topology<Triangle>& Mesh::get<Triangle>() const {
  return triangles_;
}

template <>
Topology<Triangle>& Mesh::get<Triangle>() {
  return triangles_;
}

template <>
const Topology<Polygon>& Mesh::get<Polygon>() const {
  return polygons_;
}

template <>
Topology<Polygon>& Mesh::get<Polygon>() {
  return polygons_;
}

template <>
const Topology<Quad>& Mesh::get<Quad>() const {
  return quads_;
}

template <>
Topology<Quad>& Mesh::get<Quad>() {
  return quads_;
}

template <>
const Topology<Tet>& Mesh::get<Tet>() const {
  return tetrahedra_;
}

template <>
Topology<Tet>& Mesh::get<Tet>() {
  return tetrahedra_;
}

template <>
const Topology<Polyhedron>& Mesh::get<Polyhedron>() const {
  return polyhedra_;
}

template <>
Topology<Polyhedron>& Mesh::get<Polyhedron>() {
  return polyhedra_;
}

int Mesh::get_surface_connected_components(std::vector<int>& components) const {
  ASSERT(quads_.n() == 0);
  ASSERT(polygons_.n() == 0);
  NOT_IMPLEMENTED;
  return -1;
}

void Vertices::print() const {
  for (int k = 0; k < n(); k++) {
    std::cout << fmt::format("v[{}] = (", k);
    for (int d = 0; d < dim(); d++) {
      std::cout << (*this)[k][d];
      if (d + 1 < dim())
        std::cout << ", ";
      else
        std::cout << ") ";
    }
    if (entity_[k] != nullptr) {
      std::cout << (*entity_[k]);
    } else
      std::cout << "interior";
    std::cout << std::endl;
  }
}

}  // namespace helix
