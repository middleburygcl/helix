#include "util.h"

#include "mesh.h"
#include "morton-nd/mortonND_LUT.h"
#include "stlext.h"

using morton_t = uint64_t;
namespace mortonnd {
constexpr auto encoder64_2d = MortonNDLutEncoder<2, 32, 10>();
constexpr auto encoder64_3d = MortonNDLutEncoder<3, 21, 10>();
constexpr auto encoder64_4d = MortonNDLutEncoder<4, 16, 10>();

template <int dim>
struct Resolution;

template <>
struct Resolution<2> {
  // 32 bits for each coordinate
  static constexpr double value = 4294967296.0;
};

template <>
struct Resolution<3> {
  // 21 bits for each coordinate
  static constexpr double value = 2091752.0;
};

template <>
struct Resolution<4> {
  // 16 bits for each coordinate
  static constexpr double value = 65536.0;
};

template <int dim>
morton_t encode(const std::array<uint64_t, dim>& x);

template <>
morton_t encode<2>(const std::array<uint64_t, 2>& x) {
  return encoder64_2d.Encode(x[0], x[1]);
}

template <>
morton_t encode<3>(const std::array<uint64_t, 3>& x) {
  return encoder64_3d.Encode(x[0], x[1], x[2]);
}

template <>
morton_t encode<4>(const std::array<uint64_t, 4>& x) {
  return encoder64_4d.Encode(x[0], x[1], x[2], x[3]);
}

}  // namespace mortonnd

namespace helix {

template <typename T, int dim>
morton_t encode_morton(const T* values, const T* xmin, const T* xmax) {
  const auto n = mortonnd::Resolution<dim>::value;
  std::array<uint64_t, dim> u;
  for (int d = 0; d < dim; d++) {
    T x = (values[d] - xmin[d]) / (xmax[d] - xmin[d]);
    x = std::min(std::max(x * n, 0.0), n);
    u[d] = static_cast<uint64_t>(x);
  }
  return mortonnd::encode<dim>(u);
}

template <typename T>
morton_t encode_morton(const T* values, int dim, const T* xmin, const T* xmax) {
  if (dim == 2)
    return encode_morton<T, 2>(values, xmin, xmax);
  else if (dim == 3)
    return encode_morton<T, 3>(values, xmin, xmax);
  assert(dim == 4);
  return encode_morton<T, 4>(values, xmin, xmax);
}

template <typename coord_t>
void get_bounding_box(const coord_t* points,
                      int64_t n_points,
                      int8_t dim,
                      coord_t* xmin,
                      coord_t* xmax) {
  std::fill(xmin, xmin + dim, std::numeric_limits<coord_t>::max());
  std::fill(xmax, xmax + dim, std::numeric_limits<coord_t>::min());
  for (int64_t i = 0; i < n_points; i++) {
    for (int d = 0; d < dim; d++) {
      const auto& x = points[dim * i + d];
      if (x < xmin[d]) xmin[d] = x;
      if (x > xmax[d]) xmax[d] = x;
    }
  }
}

template <typename coord_t, typename index_t>
void sort_points_on_zcurve(const coord_t* points,
                           uint64_t n_points,
                           int8_t dim,
                           std::vector<index_t>& order) {
  std::vector<coord_t> xmin(dim), xmax(dim);
  get_bounding_box(points, n_points, dim, xmin.data(), xmax.data());
  std::vector<std::pair<int64_t, morton_t>> z(order.size());
  std::parafor_i(0, n_points, [&](int thread_id, int k) {
    z[k] = {k, encode_morton(&points[k * dim], dim, xmin.data(), xmax.data())};
  });
  std::parasort(z.begin(), z.end(), [](const auto& p, const auto& q) {
    return p.second < q.second;
  });
  for (uint64_t k = 0; k < n_points; k++)
    order[k] = z[k].first;
}

#define INSTANTIATE_ZCURVE(INDEX_T)                           \
  template void sort_points_on_zcurve<float, INDEX_T>(        \
      const float*, uint64_t, int8_t, std::vector<INDEX_T>&); \
  template void sort_points_on_zcurve<double, INDEX_T>(       \
      const double*, uint64_t, int8_t, std::vector<INDEX_T>&);

INSTANTIATE_ZCURVE(int32_t)
INSTANTIATE_ZCURVE(uint32_t)
INSTANTIATE_ZCURVE(int64_t)
INSTANTIATE_ZCURVE(uint64_t)
#undef INSTANTIATE_ZCURVE

// helper struct to represent a 3d vector
struct vec3d : std::array<coord_t, 3> {
  vec3d() {
    for (int d = 0; d < 3; d++)
      (*this)[d] = 0;
  }
  vec3d(const coord_t* x) {
    for (int d = 0; d < 3; d++)
      (*this)[d] = x[d];
  }
};

void sample_surface(const Mesh& mesh, Vertices& sites, index_t n) {
  for (index_t i = 0; i < n; i++) {
    // pick a triangle
    index_t t = floor(mesh.triangles().n() * double(rand()) / double(RAND_MAX));

    // generate a random point in a triangle
    double r1 = double(rand()) / double(RAND_MAX);
    double r2 = double(rand()) / double(RAND_MAX);
    double alpha = 1.0 - std::sqrt(r1);
    double beta = std::sqrt(r1) * (1.0 - r2);
    double gamma = 1.0 - alpha - beta;

    vec3d u(mesh.vertices()[mesh.triangles()(t, 0)]);
    vec3d v(mesh.vertices()[mesh.triangles()(t, 1)]);
    vec3d w(mesh.vertices()[mesh.triangles()(t, 2)]);

    vec3d p;
    for (int d = 0; d < 3; d++)
      p[d] = alpha * u[d] + beta * v[d] + gamma * w[d];
    sites.add(p.data());
  }
}

void sample_volume(const Mesh& mesh, Vertices& sites, index_t n) {
  for (index_t i = 0; i < n; i++) {
    // pick a tetrahedron
    index_t T =
        floor(mesh.tetrahedra().n() * double(rand()) / double(RAND_MAX));

    // generate a random point in a tetrahedron
    double s = double(rand()) / double(RAND_MAX);
    double t = double(rand()) / double(RAND_MAX);
    double u = double(rand()) / double(RAND_MAX);
    if (s + t > 1.0) {
      s = 1.0 - s;
      t = 1.0 - t;
    }
    if (t + u > 1.0) {
      double tmp = u;
      u = 1.0 - s - t;
      t = 1.0 - tmp;

    } else if (s + t + u > 1.0) {
      double tmp = u;
      u = s + t + u - 1.0;
      s = 1 - t - tmp;
    }
    double a = 1 - s - t - u;

    vec3d p0(mesh.vertices()[mesh.tetrahedra()(T, 0)]);
    vec3d p1(mesh.vertices()[mesh.tetrahedra()(T, 1)]);
    vec3d p2(mesh.vertices()[mesh.tetrahedra()(T, 2)]);
    vec3d p3(mesh.vertices()[mesh.tetrahedra()(T, 3)]);

    vec3d p;
    for (int d = 0; d < 3; d++)
      p[d] = a * p0[d] + s * p1[d] + t * p2[d] + u * p3[d];
    sites.add(p.data());
  }
}

}  // namespace helix
