#include "library.h"

#include <cmath>
#include <map>
#include <vector>

namespace helix {

template <typename type>
Grid<type>::Grid(const std::vector<int>& sizes, int dim)
    : Mesh((dim < 0) ? sizes.size() : dim), sizes_(sizes) {
  build();
}

template <>
void Grid<Triangle>::build() {
  int nx = sizes_[0];
  int ny = sizes_[1];

  double dx = 1. / double(nx);
  double dy = 1. / double(ny);

  vertices_.reserve((nx + 1) * (ny + 1));
  triangles_.reserve(2 * nx * ny);

  std::vector<double> x(vertices_.dim(), 0.0);
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      x[0] = i * dx;
      x[1] = j * dy;

      vertices_.add(x.data());
    }
  }

  index_t t[3];
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int i0 = j * (nx + 1) + i;
      int i1 = i0 + 1;
      int i2 = i1 + nx + 1;
      int i3 = i2 - 1;

      t[0] = i0;
      t[1] = i1;
      t[2] = i2;
      triangles_.add(t);

      t[0] = i0;
      t[1] = i2;
      t[2] = i3;
      triangles_.add(t);
    }
  }
}

template <>
void Grid<Quad>::build() {
  // these are actually quads
  int nx = sizes_[0];
  int ny = sizes_[1];

  double dx = 1. / double(nx);
  double dy = 1. / double(ny);

  vertices_.reserve((nx + 1) * (ny + 1));
  quads_.reserve(nx * ny);

  std::vector<double> x(vertices_.dim(), 0.0);
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      x[0] = i * dx;
      x[1] = j * dy;

      vertices_.add(x.data());
    }
  }

  index_t p[4];
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int i0 = j * (nx + 1) + i;
      int i1 = i0 + 1;
      int i2 = i1 + nx + 1;
      int i3 = i2 - 1;

      p[0] = i0;
      p[1] = i1;
      p[2] = i2;
      p[3] = i3;
      quads_.add(p);
    }
  }
}

template <>
void Grid<Polygon>::build() {
  // these are actually quads
  int nx = sizes_[0];
  int ny = sizes_[1];

  double dx = 1. / double(nx);
  double dy = 1. / double(ny);

  vertices_.reserve((nx + 1) * (ny + 1));
  polygons_.reserve(nx * ny);

  std::vector<double> x(vertices_.dim(), 0.0);
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      x[0] = i * dx;
      x[1] = j * dy;

      vertices_.add(x.data());
    }
  }

  index_t p[4];
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int i0 = j * (nx + 1) + i;
      int i1 = i0 + 1;
      int i2 = i1 + nx + 1;
      int i3 = i2 - 1;

      p[0] = i0;
      p[1] = i1;
      p[2] = i2;
      p[3] = i3;
      polygons_.add(p, 4);
    }
  }
}

template <>
void Grid<Tet>::build() {
  int nx = sizes_[0];
  int ny = sizes_[1];
  int nz = sizes_[2];

  double dx = 1. / double(nx);
  double dy = 1. / double(ny);
  double dz = 1. / double(nz);

  vertices_.reserve((nx + 1) * (ny + 1) * (nz + 1));
  tetrahedra_.reserve(6 * nx * ny * nz);

  // create the vertices
  std::vector<double> x(vertices_.dim(), 0.0);
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        x[0] = i * dx;
        x[1] = j * dy;
        x[2] = k * dz;
        vertices_.add(x.data());
      }
    }
  }

  const int hextets[6][4] = {{0, 1, 2, 4}, {1, 4, 3, 2}, {6, 2, 3, 4},

                             {1, 4, 5, 3}, {4, 6, 5, 3}, {7, 6, 3, 5}};

  const int joffset = (nx + 1);
  const int koffset = (nx + 1) * (ny + 1);

  unsigned long t[4];
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        const int n0 = k * koffset + j * joffset + i;

        // all the vertex indices that make up an individial hex
        const int hexnodes[8] = {n0 + 0,
                                 n0 + 1,
                                 n0 + joffset + 0,
                                 n0 + joffset + 1,

                                 n0 + koffset + 0,
                                 n0 + koffset + 1,
                                 n0 + koffset + joffset + 0,
                                 n0 + koffset + joffset + 1};

        // loop over all tets that make up a hex
        for (int tet = 0; tet < 6; tet++) {
          // map the nodes from the hex for each tet
          for (int n = 0; n < 4; n++) t[n] = hexnodes[hextets[tet][n]];

          // add the tet indices
          tetrahedra_.add(t);
        }
      }
    }
  }
}

void Honeycomb::build(int n) {
  int nx = n;
  double dx = 1.0 / (nx - 1.0);
  double dy = 0.5 * std::sqrt(3.0) * dx;
  int ny = 1.0 / dy + 1;

  vertices_.reserve(nx * ny);
  std::vector<coord_t> point(vertices_.dim(), 0.0);
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++) {
      double x = i * dx - 0.25 * dx;
      if (j % 2 == 1) x += 0.5 * dx;

      point[0] = x;
      point[1] = j * dy + 0.25 * dy;

      vertices_.add(point.data());
    }

  /*
  auto add_centroid = [&]( index_t* t ) {
    double x[3] = {0,0};
    for (int i = 0; i < 3; i++)
    for (int d = 0; d < 2; d++)
      x[d] += vertices_[t[i]][d] / 3.0;
    vertices_.add(x);
  };
  */

  Topology<Triangle> triangles;
  index_t t1[3], t2[3];
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int i0 = j * (nx + 1) + i;
      int i1 = i0 + 1;
      int i2 = i1 + nx + 1;
      int i3 = i2 - 1;

      if (i % 2 == 0) {
        t1[0] = i0;
        t1[1] = i1;
        t1[2] = i3;

        t2[0] = i1;
        t2[1] = i2;
        t2[2] = i3;
      } else {
        t1[0] = i0;
        t1[1] = i1;
        t1[2] = i2;

        t2[0] = i0;
        t2[1] = i2;
        t2[2] = i3;
      }
      triangles_.add(t1);
      triangles_.add(t2);
    }
  }
}

void Sphere::build(int n) {
  // icosahedron vertices
  coord_t t = (1.0 + std::sqrt(5.0)) / 2.0;

  coord_t coordinates[12][3] = {{t, 1, 0}, {-t, 1, 0}, {t, -1, 0}, {-t, -1, 0},
                                {1, 0, t}, {1, 0, -t}, {-1, 0, t}, {-1, 0, -t},
                                {0, t, 1}, {0, -t, 1}, {0, t, -1}, {0, -t, -1}};

  for (int i = 0; i < 12; i++) {
    for (int d = 0; d < 3; d++) coordinates[i][d] /= std::sqrt(1 + t * t);
    vertices_.add(coordinates[i]);
  }

  index_t triangles[20][3] = {
      {0, 8, 4},   // 0
      {0, 5, 10},  // 1
      {2, 4, 9},   // 2
      {2, 11, 5},  // 3
      {1, 6, 8},   // 4
      {1, 10, 7},  // 5
      {3, 9, 6},   // 6
      {3, 7, 11},  // 7
      {0, 10, 8},  // 8
      {1, 8, 10},  // 9
      {2, 9, 11},  // 10
      {3, 11, 9},  // 11
      {4, 2, 0},   // 12
      {5, 0, 2},   // 13
      {6, 1, 3},   // 14
      {7, 3, 1},   // 15
      {8, 6, 4},   // 16
      {9, 4, 6},   // 17
      {10, 5, 7},  // 18
      {11, 7, 5}   // 19
  };

  for (int i = 0; i < 20; i++) {
    triangles_.add(triangles[i]);
    triangles_.set_group(i, 0);
  }

  for (int i = 0; i < n; i++) subdivide();
}

void Sphere::subdivide() {
  std::map<Edge, index_t> edges;

  Topology<Triangle> triangles;
  for (int k = 0; k < triangles_.n(); k++) {
    int edge_indices[3];
    for (int j = 0; j < 3; j++) {
      index_t e0 = triangles_(k, j);
      index_t e1 = triangles_(k, (j + 1) % 3);

      // does this edge exist?
      index_t idx;
      auto itr = edges.find({e1, e0});
      if (itr == edges.end()) {
        // create a new point on this edge
        std::array<coord_t, 3> q;
        coord_t l = 0;
        for (int d = 0; d < 3; d++) {
          q[d] = 0.5 * (vertices_[e0][d] + vertices_[e1][d]);
          l += q[d] * q[d];
        }
        l = std::sqrt(l);
        // normalize to place point on unit sphere
        for (int d = 0; d < 3; d++) q[d] /= l;

        idx = vertices_.n();
        vertices_.add(q.data());
        edges.insert({{e0, e1}, idx});
      } else {
        idx = itr->second;
      }
      edge_indices[j] = idx;
    }

    // create the four new triangles from the subdivision
    index_t t0 = triangles_(k, 0);
    index_t t1 = triangles_(k, 1);
    index_t t2 = triangles_(k, 2);

    index_t e0 = edge_indices[0];
    index_t e1 = edge_indices[1];
    index_t e2 = edge_indices[2];

    index_t triangle0[3] = {t0, e0, e2};
    index_t triangle1[3] = {e0, t1, e1};
    index_t triangle2[3] = {e2, e1, t2};
    index_t triangle3[3] = {e0, e1, e2};

    triangles.add(triangle0);
    triangles.add(triangle1);
    triangles.add(triangle2);
    triangles.add(triangle3);
  }

  triangles_.clear(false);
  for (index_t k = 0; k < triangles.n(); k++) {
    triangles_.add(triangles[k]);
    triangles_.set_group(k, 0);
  }
}

void Sphere::add_prism_layers(double rf, int n_layers, double growth_rate) {
  using vec3d = std::array<coord_t, 3>;
  ASSERT(rf > 1.0);

  double r0 = 1.0;
  double dr = (rf - r0) / n_layers;

  int n = vertices_.n();
  int m = triangles_.n();  // surface triangles

  for (int k = 0; k < n_layers; k++) {
    for (int j = 0; j < n; j++) {
      vec3d p;
      coord_t l = 0.0;
      for (int d = 0; d < 3; d++) {
        p[d] = vertices_[j][d];
        l += p[d];
      }
      l = std::sqrt(l);
      vec3d q;
      for (int d = 0; d < 3; d++)
        q[d] = p[d] + std::pow((k + 1) * dr, growth_rate) * p[d] / l;
      vertices_.add(q.data());
    }
  }

  for (int k = 0; k < n_layers; k++) {
    for (int j = 0; j < m; j++) {
      index_t prism[6];
      for (int i = 0; i < 3; i++) {
        prism[i] = triangles_(j, i) + k * n;
        prism[3 + i] = triangles_(j, i) + (k + 1) * n;
      }

      prisms_.add(prism);
    }
  }

  // add a boundary group at the outer surface
  for (int j = 0; j < m; j++) {
    index_t triangle[3];
    for (int i = 0; i < 3; i++) triangle[i] = triangles_(j, i) + n_layers * n;

    triangles_.add(triangle);
    triangles_.set_group(m + j, 1);
  }
}

template class Grid<Triangle>;
template class Grid<Quad>;
template class Grid<Polygon>;
template class Grid<Tet>;
template class Grid<Polyhedron>;

}  // namespace helix
