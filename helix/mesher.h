#pragma once

namespace helix {

class Mesh;
class Mesher {};

struct MeshingParameters {
  double quality{-1};
  bool verbose{false};
  bool constrained{true};
};

class TetGenMesher : public Mesher {
 public:
  TetGenMesher(Mesh& surface);
  void generate(const MeshingParameters& params);

 private:
  Mesh& mesh_;
};

}  // namespace helix