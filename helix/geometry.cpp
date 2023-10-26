#include "geometry.h"

#include <unordered_map>

#include "mesh.h"

#define HELIX_WITH_EGADS 1
#if HELIX_WITH_EGADS
#include <egads.h>

#include "GeomConvert.hxx"
#include "GeomConvert_BSplineSurfaceToBezierSurface.hxx"
#include "egadsClasses.h"
#define EGADS_ENABLE(X) X
#else
struct egObject {};
typedef egObject* ego;
#define EGADS_ENABLE(X) NOT_CONFIGURED;
#endif

#include "log.h"

namespace helix {

static const std::map<int, std::string> EgadsErrorMessage = {
    {-37, "extrapol"}, {-36, "effctobj"}, {-35, "uvmap"},
    {-34, "sequerr"},  {-33, "cntxhrd"},  {-32, "readerr"},
    {-31, "tesstate"}, {-30, "exists"},   {-29, "attrerr"},
    {-28, "topocnt"},  {-27, "ocsegflt"}, {-26, "badscale"},
    {-25, "notortho"}, {-24, "degen"},    {-23, "consterr"},
    {-22, "topoerr"},  {-21, "geomerr"},  {-20, "notbody"},
    {-19, "writerr"},  {-18, "notmodel"}, {-17, "noload"},
    {-16, "rangerr"},  {-15, "notgeom"},  {-14, "nottess"},
    {-13, "empty"},    {-12, "nottopo"},  {-11, "referce"},
    {-10, "notxform"}, {-9, "notcntx"},   {-8, "mixcntx"},
    {-7, "nodata"},    {-6, "noname"},    {-5, "indexerr"},
    {-4, "malloc"},    {-3, "notobj"},    {-2, "nullobj"},
    {-1, "notfound"},  {0, "success"},    {1, "outside/not the same"}};

#define EG_CHECK(X)                        \
  {                                        \
    int status = (X);                      \
    if (status != EGADS_SUCCESS) {         \
      LOG << EgadsErrorMessage.at(status); \
    }                                      \
  }

#define EG_ASSERT(X)                                                 \
  {                                                                  \
    int status = (X);                                                \
    ASSERT(status == EGADS_SUCCESS) << EgadsErrorMessage.at(status); \
  }

bool Entity::above(const Entity* e) const {
  if (!e) return false;
  if (e->dim() >= dim_) return false;
  for (const auto* child : children_) {
    ASSERT(child) << (*this) << " " << children_.size();
    if (child == e) return true;
    bool result = child->above(e);
    if (result) return true;
  }
  return false;
}

bool Entity::below(const Entity* e) const {
  if (!e) return false;
  return !e->above(this);
}

bool Entity::sibling(const Entity* e) const {
  if (!e) return false;
  if (e == this) return true;
  if (e->parent() == parent_) return true;
  return false;
}

Entity* Entity::common(Entity* e) {
  if (!e) return nullptr;
  if (e == this) return e;
  if (parent_ == e) return e;
  if (e->parent() == this) return this;
  if (above(e)) return this;
  if (e->above(this)) return e;
  return nullptr;
}

namespace EGADS {
class Context : public helix::Context {
 public:
  Context(const helix::Context& ctx)
      : helix::Context(GeometryKernel::kEgads), owner_(false) {
    ASSERT(ctx.type() == GeometryKernel::kEgads);
    const auto& ectx = static_cast<const Context&>(ctx);
    context_ = ectx.context();
  }
  Context() : helix::Context(GeometryKernel::kEgads), owner_(true) {
    EGADS_ENABLE(EG_open(&context_));
  }
  ~Context() {
    if (owner_) {
      EGADS_ENABLE(EG_setOutLevel(context_, 0));
      EGADS_ENABLE(EG_close(context_));
    }
  }

  ego context() const { return context_; }

  void load(const std::string& filename, Geometry& geometry);
  void tessellate(Mesh& mesh) const;
  void get_control_mesh(Mesh& mesh) const;
  void save(const std::string& filename) const {
    EG_ASSERT(EG_saveModel(model_, filename.c_str()));
  }
  void operate(const Operation& op);
  void build();

  Entity* create_body() {
    NOT_IMPLEMENTED;
    return nullptr;
  }
  Entity* create_face() {
    NOT_IMPLEMENTED;
    return nullptr;
  }
  Entity* create_edge() {
    NOT_IMPLEMENTED;
    return nullptr;
  }
  Entity* create_node();

 private:
  ego context_;
  ego model_;
  bool owner_;
  std::unordered_map<helix::Entity*, ego> entity_map_;
  std::unordered_map<ego, helix::Entity*> ego_map_;
};

class Entity : public helix::Entity {
 public:
  Entity(int dim, EntityType type, ego object, int id)
      : helix::Entity(dim, type, id), object_(object) {}
  void evaluate(const coord_t* u, coord_t* x) const {
    double param[2];
    for (int d = 0; d < dim_; d++)
      param[d] = u[d];
    double xyz[18];
    EG_ASSERT(EG_evaluate(object_, param, xyz));
    for (int d = 0; d < 3; d++)
      x[d] = xyz[d];
  }
  void project(const coord_t* x, coord_t* u, coord_t* u_guess = nullptr) const {
    double xyz[3], params[2], result[3];
    for (int d = 0; d < 3; d++)
      xyz[d] = x[d];
    if (!u_guess) {
      EG_ASSERT(EG_invEvaluate(object_, xyz, params, result));
      for (int d = 0; d < dim_; d++)
        u[d] = params[d];
    } else
      NOT_IMPLEMENTED;
  }
  void get_normal(const coord_t* u, coord_t* n) const { NOT_IMPLEMENTED; }
  void get_params(const helix::Entity* e,
                  const coord_t* u_lower,
                  coord_t* u_upper) const {
    NOT_IMPLEMENTED;
  }

 private:
  ego object_;
};

helix::Entity* Context::create_node() {
  entities_.emplace_back(
      std::make_unique<EGADS::Entity>(0, EntityType::kNode, nullptr, 0));
  return entities_.back().get();
}

void Context::load(const std::string& filename, Geometry& geometry) {
#if HELIX_WITH_EGADS
  EG_loadModel(context_, 2, filename.c_str(), &model_);
#endif
  build();
}

void Context::build() {
#if HELIX_WITH_EGADS
  bodies_.clear();
  entities_.clear();
  ego_map_.clear();
  entity_map_.clear();

  // extract the bodies
  ego ref, *bodies;
  int oclass, mtype, nchild, *senses;
  double reals[4];
  EG_ASSERT(EG_getTopology(model_, &ref, &oclass, &mtype, reals, &nchild,
                           &bodies, &senses));

  entity_map_.reserve(nchild);
  ego_map_.reserve(nchild);

  // create bodies
  bodies_.resize(nchild);
  for (int k = 0; k < nchild; k++) {
    entities_.emplace_back(
        std::make_unique<EGADS::Entity>(3, EntityType::kBody, bodies[k], k));
    helix::Entity* body = entities_.back().get();
    bodies_[k] = body;
    ego_map_.insert({bodies[k], body});
    entity_map_.insert({body, bodies[k]});

    // get all the faces, edges and nodes
    int n_nodes, n_edges, n_faces;
    ego* faces;
    ego* edges;
    ego* nodes;
    EG_ASSERT(EG_getBodyTopos(bodies[k], nullptr, NODE, &n_nodes, &nodes));
    EG_ASSERT(EG_getBodyTopos(bodies[k], nullptr, EDGE, &n_edges, &edges));
    EG_ASSERT(EG_getBodyTopos(bodies[k], nullptr, FACE, &n_faces, &faces));
    // LOG << fmt::format("body has {} faces, {} edges, {} nodes", n_faces,
    //                    n_edges, n_nodes);

    size_t n_entities = ego_map_.max_load_factor() * ego_map_.bucket_count();
    entity_map_.reserve(n_nodes + n_edges + n_faces + n_entities);
    ego_map_.reserve(n_nodes + n_edges + n_faces + n_entities);

    for (int j = 0; j < n_nodes; j++) {
      int id = EG_indexBodyTopo(bodies[k], nodes[j]);
      entities_.emplace_back(
          std::make_unique<EGADS::Entity>(0, EntityType::kNode, nodes[j], id));
      helix::Entity* node = entities_.back().get();
      ego_map_.insert({nodes[j], node});
      entity_map_.insert({node, nodes[j]});
    }

    for (int j = 0; j < n_edges; j++) {
      int id = EG_indexBodyTopo(bodies[k], edges[j]);
      entities_.emplace_back(
          std::make_unique<EGADS::Entity>(1, EntityType::kEdge, edges[j], id));
      helix::Entity* edge = entities_.back().get();
      ego_map_.insert({edges[j], edge});
      entity_map_.insert({edge, edges[j]});

      int n_objects;
      ego* objects;
      EG_ASSERT(
          EG_getBodyTopos(bodies[k], edges[j], NODE, &n_objects, &objects));
      ASSERT(n_objects <= 2);

      // add child nodes
      edge->children().reserve(n_objects);
      for (int i = 0; i < n_objects; i++) {
        auto eit = ego_map_.find(objects[i]);
        ASSERT(eit != ego_map_.end());
        edge->children().push_back(eit->second);
        eit->second->set_parent(edge);
      }
      EG_free(objects);
    }

    for (int j = 0; j < n_faces; j++) {
      int id = EG_indexBodyTopo(bodies[k], faces[j]);
      entities_.emplace_back(
          std::make_unique<EGADS::Entity>(2, EntityType::kFace, faces[j], id));
      helix::Entity* face = entities_.back().get();
      ego_map_.insert({faces[j], face});
      entity_map_.insert({face, faces[j]});

      int n_objects;
      ego* objects;
      EG_ASSERT(
          EG_getBodyTopos(bodies[k], faces[j], EDGE, &n_objects, &objects));

      // add child edges
      face->children().reserve(n_objects);
      for (int i = 0; i < n_objects; i++) {
        auto eit = ego_map_.find(objects[i]);
        ASSERT(eit != ego_map_.end());
        face->children().push_back(eit->second);
        eit->second->set_parent(face);
      }
      EG_free(objects);
    }

    EG_free(nodes);
    EG_free(edges);
    EG_free(faces);
  }
#endif
}

void Context::tessellate(Mesh& mesh) const {
  std::map<ego, ego> body2tess;
  LOG << "tessellating " << bodies_.size() << " bodies";
  for (auto& body : bodies_) {
    ego ebody = entity_map_.at(body);
    ego tess;
    double box[6];
    EG_ASSERT(EG_getBoundingBox(ebody, box));
    double size = std::sqrt((box[0] - box[3]) * (box[0] - box[3]) +
                            (box[1] - box[4]) * (box[1] - box[4]) +
                            (box[2] - box[5]) * (box[2] - box[5]));
    ego* faces;
    int n_faces;
    EG_ASSERT(EG_getBodyTopos(ebody, nullptr, FACE, &n_faces, &faces));
#if 0
    for (int iface = 0; iface < n_faces; ++iface) {
      double fbox[6];
      EG_ASSERT(EG_getBoundingBox(faces[iface], fbox));
      double fsize = std::sqrt((fbox[0] - fbox[3]) * (fbox[0] - fbox[3]) +
                               (fbox[1] - fbox[4]) * (fbox[1] - fbox[4]) +
                               (fbox[2] - fbox[5]) * (fbox[2] - fbox[5]));
      if (fsize < size) size = fsize;
    }
#endif
    double parms[3];
    parms[0] = 0.025 * size;
    parms[1] = 0.001 * size;
    parms[2] = 30.0;
    EG_ASSERT(EG_makeTessBody(ebody, parms, &tess));
    body2tess.insert({ebody, tess});

    // add all the vertices from this body
    int nvert, state = 0;
    ego ref;
    EG_ASSERT(EG_statusTessBody(tess, &ref, &state, &nvert));
    ASSERT(ref == ebody);
    ASSERT(state == 1);
    int n = mesh.vertices().n();
    for (int j = 1; j <= nvert; j++) {
      double coord[3];
      int type, indx;
      EG_ASSERT(EG_getGlobal(tess, j, &type, &indx, coord));
      mesh.vertices().add(coord);
    }

    // label the nodes
    int n_nodes;
    ego* nodes;
    EG_ASSERT(EG_getBodyTopos(ebody, nullptr, NODE, &n_nodes, &nodes));
    for (int j = 0; j < n_nodes; j++) {
      int v;
      int id = EG_indexBodyTopo(ebody, nodes[j]);
      if (EG_localToGlobal(tess, 0, id, &v) != EGADS_SUCCESS) continue;
      helix::Entity* entity = ego_map_.at(nodes[j]);
      ASSERT(entity->dim() == 0);
      mesh.vertices().set_entity(n + v - 1, entity);
    }

    // add the edges
    int n_edges;
    ego* edges;
    EG_ASSERT(EG_getBodyTopos(ebody, nullptr, EDGE, &n_edges, &edges));
    for (int j = 0; j < n_edges; j++) {
      if (edges[j]->mtype == DEGENERATE) continue;
      int id = EG_indexBodyTopo(ebody, edges[j]);
      int nvert;
      const double *xyz, *t;
      EG_ASSERT(EG_getTessEdge(tess, id, &nvert, &xyz, &t));

      int type, indx;
      int edge[2];
      EG_ASSERT(EG_localToGlobal(tess, -(j + 1), 1, &edge[0]));
      EG_ASSERT(EG_getGlobal(tess, edge[0], &type, &indx, nullptr));
      ASSERT(type == 0);  // node

      edge[0] = edge[0] - 1 + n;
      for (int i = 1; i < nvert; i++) {
        EG_ASSERT(EG_localToGlobal(tess, -(j + 1), i + 1, &edge[1]));
        EG_ASSERT(EG_getGlobal(tess, edge[1], &type, &indx, nullptr));
        edge[1] = edge[1] - 1 + n;
        int i_edge = mesh.lines().n();
        mesh.lines().add(edge);
        mesh.lines().set_group(i_edge, id);
        edge[0] = edge[1];
        if (i + 1 < nvert) {
          ASSERT(type > 0) << type;  // edge
          mesh.vertices().set_entity(edge[1], ego_map_.at(edges[j]));
          mesh.vertices().set_param(edge[1], &t[i], 1);
        }
      }
    }

    // add the triangles
    for (int j = 0; j < n_faces; j++) {
      int id = EG_indexBodyTopo(ebody, faces[j]);
      int ntri, nvert;
      const double *xyz, *uv;
      const int *type, *indx, *tris, *tric;
      EG_ASSERT(EG_getTessFace(tess, id, &nvert, &xyz, &uv, &type, &indx, &ntri,
                               &tris, &tric));

      int tri[3];
      for (int i = 0; i < ntri; i++) {
        for (int d = 0; d < 3; d++) {
          EG_ASSERT(EG_localToGlobal(tess, j + 1, tris[3 * i + d], &tri[d]));
          int type, indx;
          EG_ASSERT(EG_getGlobal(tess, tri[d], &type, &indx, nullptr));
          tri[d] = tri[d] - 1 + n;
          if (type < 0) {  // face
            mesh.vertices().set_entity(tri[d], ego_map_.at(faces[j]));
            mesh.vertices().set_param(tri[d], &uv[tris[3 * i + d]], 2);
          }
        }
        int i_tri = mesh.triangles().n();
        mesh.triangles().add(tri);
        mesh.triangles().set_group(i_tri, id);
      }
    }
  }
  LOG << fmt::format("tessellation has {} triangles, {} edges and {} vertices",
                     mesh.triangles().n(), mesh.lines().n(),
                     mesh.vertices().n());
}

void Context::get_control_mesh(Mesh& mesh) const {
  bool use_bezier = true;
  auto add_grid = [&mesh](int npu, int npv, int n) {
    for (int j = 0; j < npv; j++) {
      for (int i = 0; i < npu - 1; i++) {
        int p = j * npu + i;
        int q = p + 1;
        int edge[2] = {n + p, n + q};
        mesh.lines().add(edge);
      }
    }

    for (int i = 0; i < npu; i++) {
      for (int j = 0; j < npv - 1; j++) {
        int p = j * npu + i;
        int q = p + npu;
        int edge[2] = {n + p, n + q};
        mesh.lines().add(edge);
      }
    }
  };

  for (auto& body : bodies_) {
    ego ebody = entity_map_.at(body);
    ego* faces;
    int n_faces;
    EG_ASSERT(EG_getBodyTopos(ebody, nullptr, FACE, &n_faces, &faces));
    for (int k = 0; k < n_faces; k++) {
      ego bspline;
      if (EG_convertToBSpline(faces[k], &bspline) != EGADS_SUCCESS) continue;

      if (use_bezier) {
#if 0
        auto* surface = static_cast<egadsSurface*>(bspline->blind);
        Handle(Geom_BSplineSurface) bspline_surface =
            GeomConvert::SurfaceToBSplineSurface(surface->handle);
        GeomConvert_BSplineSurfaceToBezierSurface bezier(bspline_surface);
        LOG << fmt::format("surface {} split into {} x {} bezier patches", k,
                           bezier.NbUPatches(), bezier.NbVPatches());
        for (int jj = 0; jj < bezier.NbVPatches(); jj++) {
          for (int ii = 0; ii < bezier.NbUPatches(); ii++) {
            const auto& patch = bezier.Patch(ii + 1, jj + 1);
            int npu = patch->NbUPoles();
            int npv = patch->NbVPoles();  // pole === control point
            int n = mesh.vertices().n();
            for (int j = 0; j < npv; j++) {
              for (int i = 0; i < npu; i++) {
                const auto& pnt = patch->Pole(i + 1, j + 1);
                double xyz[3] = {pnt.X(), pnt.Y(), pnt.Z()};
                mesh.vertices().add(xyz);
              }
            }
            add_grid(npu, npv, n);
          }
        }
#endif
      } else {
        int oclass, mtype;
        ego ref;
        int* ints;
        double* reals;
        EG_ASSERT(
            EG_getGeometry(bspline, &oclass, &mtype, &ref, &ints, &reals));
        ASSERT(mtype == BSPLINE);
        int npu = ints[2], npv = ints[5];
        int nku = ints[3], nkv = ints[6];
        const double* pts = reals + nku + nkv;
        int n = mesh.vertices().n();
        for (int j = 0; j < npv; j++) {
          for (int i = 0; i < npu; i++) {
            mesh.vertices().add(&pts[3 * (j * npu + i)]);
          }
        }
        add_grid(npu, npv, n);
        EG_free(ints);
        EG_free(reals);
      }
    }
  }
}

void Context::operate(const Operation& op) {
  // get the bounding box of the current model
  double box[6];
  EG_ASSERT(EG_getBoundingBox(model_, box));
  double length[3] = {box[3] - box[0], box[4] - box[1], box[5] - box[2]};
  double lmax = std::max(length[0], std::max(length[1], length[2]));
  LOG << fmt::format(
      "box = [{:3.2f}, {:3.2f}] x [{:3.2f}, {:3.2f}] x [{:3.2f} {:03.2f}]",
      box[0], box[3], box[1], box[4], box[2], box[5]);

  double center[3];
  for (int d = 0; d < 3; d++)
    center[d] = 0.5 * (box[d] + box[d + 3]);

  // create the body
  double scale = 1.0;
  if (op.scaling == LengthScaling::kRelative) scale = lmax;
  ego body;
  double data[8];
  int stype = BOX;
  if (op.shape == OperationShape::kSphere) {
    stype = SPHERE;
    for (int d = 0; d < 3; d++)
      data[d] = center[d];
    data[3] = scale * op.values[0];
    if (op.centering == CenteringMethod::kBoxMinX) data[0] = box[0];
    if (op.centering == CenteringMethod::kBoxMaxX) data[0] = box[3];
    if (op.centering == CenteringMethod::kBoxMinY) data[1] = box[1];
    if (op.centering == CenteringMethod::kBoxMaxY) data[1] = box[4];
    if (op.centering == CenteringMethod::kBoxMinZ) data[2] = box[2];
    if (op.centering == CenteringMethod::kBoxMaxZ) data[2] = box[5];
  } else if (op.shape == OperationShape::kBox) {
    stype = BOX;
    for (int d = 0; d < 3; d++) {
      double ld = op.get_length(d);
      if (op.scaling == LengthScaling::kRelative) ld *= (box[3 + d] - box[d]);
      data[d] = center[d] - 0.5 * ld;
      data[d + 3] = ld;
    }
    if (op.centering == CenteringMethod::kCenterX) data[0] = center[0];
    if (op.centering == CenteringMethod::kCenterY) data[1] = center[1];
    if (op.centering == CenteringMethod::kCenterZ) data[2] = center[2];
    if (op.centering == CenteringMethod::kBoxMinX) {
      data[0] = box[0];
      data[3] = data[3] - box[0];
    }
    if (op.centering == CenteringMethod::kBoxMaxX) data[3] = box[3] - data[0];
    if (op.centering == CenteringMethod::kBoxMinY) {
      data[1] = box[1];
      data[4] = data[4] - box[1];
    }
    if (op.centering == CenteringMethod::kBoxMaxY) data[4] = box[4] - data[1];
    if (op.centering == CenteringMethod::kBoxMinZ) {
      data[2] = box[2];
      data[5] = data[5] - box[2];
    }
    if (op.centering == CenteringMethod::kBoxMaxZ) data[5] = box[5] - data[2];
    LOG << fmt::format(
        "box: [x, y, z] = ({}, {}, {}), [dx, dy, dz] = ({}, {}, {})", data[0],
        data[1], data[2], data[3], data[4], data[5]);

  } else if (op.shape == OperationShape::kSymmetry) {
    stype = BOX;
    for (int d = 0; d < 3; d++) {
      data[d] = box[d] - 0.01 * length[d];
      data[d + 3] = lmax * 100;
    }
    if (op.centering == CenteringMethod::kCenterX) data[0] = center[0];
    if (op.centering == CenteringMethod::kCenterY) data[1] = center[1];
    if (op.centering == CenteringMethod::kCenterZ) data[2] = center[2];
    LOG << fmt::format(
        "box: [x, y, z] = ({}, {}, {}), [dx, dy, dz] = ({}, {}, {})", data[0],
        data[1], data[2], data[3], data[4], data[5]);
  }
  EG_ASSERT(EG_makeSolidBody(context_, stype, data, &body));

  // determine the order of the bodies in the operation
  ego src, tool;
  if (op.existing_model_is_source) {
    src = model_;
    tool = body;
  } else {
    src = body;
    tool = model_;
  }

  int oper = INTERSECTION;
  if (op.type == OperationType::kSubtraction)
    oper = SUBTRACTION;
  else {
    NOT_IMPLEMENTED;
  }
  EG_ASSERT(EG_generalBoolean(src, tool, oper, op.tol, &model_));
  build();
}

}  // namespace EGADS

namespace PLC {
class Context : public helix::Context {
 public:
  Context() : helix::Context(GeometryKernel::kPlc) {}
  ~Context() {}
  void load(const std::string& filename, Geometry& geometry) {
    NOT_IMPLEMENTED;
  }
  void tessellate(Mesh& mesh) const { NOT_IMPLEMENTED; }
  void get_control_mesh(Mesh& mesh) const { NOT_IMPLEMENTED; }
  void operate(const Operation& op) { NOT_IMPLEMENTED; }
  void save(const std::string& filename) const { NOT_IMPLEMENTED; }

  Entity* create_body() {
    NOT_IMPLEMENTED;
    return nullptr;
  }
  Entity* create_face() {
    NOT_IMPLEMENTED;
    return nullptr;
  }
  Entity* create_edge() {
    NOT_IMPLEMENTED;
    return nullptr;
  }
  Entity* create_node() {
    NOT_IMPLEMENTED;
    return nullptr;
  }
};

}  // namespace PLC

Geometry::Geometry(GeometryKernel type) {
  if (type == GeometryKernel::kEgads)
    context_ = std::make_unique<EGADS::Context>();
  else if (type == GeometryKernel::kPlc)
    context_ = std::make_unique<PLC::Context>();
  else
    NOT_POSSIBLE;
}

Geometry::Geometry(GeometryKernel type, const std::string& filename)
    : Geometry(type) {
  context_->load(filename, *this);
}

}  // namespace helix
