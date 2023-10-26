#pragma once

#include <fmt/format.h>

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "defs.h"
#include "log.h"

namespace helix {

class Mesh;
enum EntityType : uint8_t {
  kNode,
  kEdge,
  kLoop,
  kFace,
  kShell,
  kBody,
  kVolume,
  kModel
};
static const std::map<EntityType, std::string> EntityName{
    {kNode, "node"}, {kEdge, "edge"}, {kFace, "face"}};

static const std::set<EntityType> TessellatableEntities{kNode, kEdge, kFace};

enum GeometryKernel : uint8_t { kEgads, kPlc };

class Entity {
 public:
  Entity(int dim, EntityType type, int id) : dim_(dim), type_(type), id_(id) {}
  virtual ~Entity() {}
  virtual void evaluate(const coord_t* u, coord_t* x) const = 0;
  virtual void project(const coord_t* x, coord_t* u,
                       coord_t* u_guess = nullptr) const = 0;
  virtual void get_normal(const coord_t* u, coord_t* n) const = 0;
  virtual void get_params(const Entity* e, const coord_t* u_lower,
                          coord_t* u_upper) const = 0;

  void set_parent(Entity* p) { parent_ = p; }
  const Entity* parent() const { return parent_; }
  auto& children() { return children_; }
  bool above(const Entity* e) const;
  bool below(const Entity* e) const;
  bool sibling(const Entity* e) const;
  Entity* common(Entity* e);

  const Entity* intersect(const Entity* e) const;
  const Entity* intersect(const Entity* e0, const Entity* e1) const;
  const Entity* intersect(const Entity* e0, const Entity* e1,
                          const Entity* e2) const;

  int dim() const { return dim_; }
  EntityType type() const { return type_; }
  int64_t id() const { return id_; }
  bool tessellatable() const {
    return TessellatableEntities.find(type_) != TessellatableEntities.end();
  }

  friend std::ostream& operator<<(std::ostream& out, const Entity& e) {
    out << fmt::format("{}: dim {}, id {}", EntityName.at(e.type()), e.dim(),
                       e.id());
    return out;
  }

 protected:
  int dim_{-1};
  EntityType type_;
  int64_t id_{-1};

  std::vector<Entity*> children_;
  Entity* parent_;
  bool interior_{false};
};

enum class OperationType : uint8_t {
  kNone,
  kUnion,
  kSubtraction,
  kIntersection,
};
enum class LengthScaling : uint8_t { kRelative, kAbsolute };
enum class OperationShape : uint8_t {
  kNone,
  kBox,
  kSphere,
  kHalfSpace,
  kSymmetry
};
enum class CenteringMethod : uint8_t {
  kNone,
  kUser,
  kModel,
  kBoxMaxX,
  kBoxMinX,
  kBoxMaxY,
  kBoxMinY,
  kBoxMaxZ,
  kBoxMinZ,
  kCenterX,
  kCenterY,
  kCenterZ
};

struct Operation {
  OperationType type{OperationType::kNone};
  CenteringMethod centering{CenteringMethod::kModel};
  LengthScaling scaling{LengthScaling::kRelative};
  OperationShape shape{OperationShape::kNone};
  bool existing_model_is_source{false};
  double center[3];
  double values[6];
  // kBox: values represent length of box (lx, ly, lz)
  // kSphere: first value is radius of sphere (r)
  // if the centering method is kSymmetry, values[3,4,5] are the normal
  void set_normal(double nx, double ny, double nz) {
    values[3] = nx;
    values[4] = ny;
    values[5] = nz;
  }
  double normal(int d) const { return values[d + 3]; }
  void set_radius(double r) { values[0] = r; }
  void set_box_lengths(double lx, double ly, double lz) {
    values[0] = lx;
    values[1] = ly;
    values[2] = lz;
  }
  double get_length(int d) const { return values[d]; }
  void set_center(double cx, double cy, double cz) {
    center[0] = cx;
    center[1] = cy;
    center[2] = cz;
  }
  // double center(int d) const { return center[d]; }
  double tol{0.0};
};

class Geometry;
class Context {
 public:
  Context(GeometryKernel t) : type_(t) {}
  virtual ~Context() {}
  virtual void load(const std::string& filename, Geometry& geometry) = 0;
  virtual void operate(const Operation& op) = 0;
  virtual void tessellate(Mesh& mesh) const = 0;
  virtual void get_control_mesh(Mesh& mesh) const = 0;
  virtual void save(const std::string& filename) const = 0;
  virtual Entity* create_body() = 0;
  virtual Entity* create_face() = 0;
  virtual Entity* create_edge() = 0;
  virtual Entity* create_node() = 0;

  GeometryKernel type() const { return type_; }

  size_t n_bodies() const { return bodies_.size(); }
  size_t n_entities() const { return entities_.size(); }

  const Entity& body(size_t k) const { return *bodies_[k]; }
  const Entity& entity(size_t k) const { return *entities_[k]; }

 protected:
  GeometryKernel type_;
  std::vector<std::unique_ptr<Entity>> entities_;
  std::vector<Entity*> bodies_;
};

class Geometry {
 public:
  Geometry(GeometryKernel type);
  Geometry(GeometryKernel type, const std::string& filename);

  void operate(const Operation& op) { context_->operate(op); }
  void tessellate(Mesh& mesh) { context_->tessellate(mesh); }
  void get_control_mesh(Mesh& mesh) { context_->get_control_mesh(mesh); }
  void save(const std::string& filename) const { context_->save(filename); }

  size_t n_bodies() const { return context_->n_bodies(); }
  size_t n_entities() const { return context_->n_entities(); }

  Entity* create_body() { return context_->create_body(); }
  Entity* create_face() { return context_->create_face(); }
  Entity* create_edge() { return context_->create_edge(); }
  Entity* create_node() { return context_->create_node(); }

 private:
  std::unique_ptr<Context> context_;
};

}  // namespace helix