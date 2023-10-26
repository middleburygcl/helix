#include <argparse/argparse.hpp>
#include <csignal>
#include <filesystem>

#include "geometry.h"
#include "io.h"
#include "log.h"
#include "mesh.h"
#include "mesher.h"
// #include "subdivision.h"

void signal_handler(int signal) {
  LOG << "caught signal " << signal;
  std::exit(0);
}

namespace helix {

namespace {
std::pair<std::string, std::string> file_ext(const std::string& filename) {
  auto idx = filename.rfind('.');  // find the '.' in reverse order
  ASSERT(idx != std::string::npos);
  return {filename.substr(0, idx), filename.substr(idx + 1)};
}

bool is_geometry_file(const std::string& filename) {
  auto ext = file_ext(filename).second;
  std::transform(ext.begin(), ext.end(), ext.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return (ext == "step" || ext == "stp" || ext == "egads" || ext == "iges" ||
          ext == "igs");
}

}  // namespace

void run_cad(argparse::ArgumentParser& program) {
  auto input = program.get<std::string>("input");
  auto file_info = file_ext(input);
  LOG << "loading geometry [" << file_info.first << "] with ext "
      << file_info.second;
  std::string output;
  if (program.is_used("--output"))
    output = program.get<std::string>("--output");
  else
    output = file_info.first + ".output." + file_info.second;

  std::vector<double> axis;
  std::vector<std::string> lengths_str;
  std::string shape_str, radius_str, height_str, symmetry;
  bool enclose = false;
  auto center_info = program.get<std::string>("--center");
  if (program.is_used("--axis"))
    axis = program.get<std::vector<double>>("--axis");
  if (program.is_used("--radius"))
    radius_str = program.get<std::string>("--radius");
  if (program.is_used("--lengths"))
    lengths_str = program.get<std::vector<std::string>>("--lengths");
  if (program.is_used("--height"))
    height_str = program.get<std::string>("--height");
  if (program.is_used("--symmetry")) {
    symmetry = program.get<std::string>("--symmetry");
    ASSERT(symmetry.size() == 1)
        << "specify either x, y or z, received " << symmetry;
  }
  if (program.is_used("--enclose")) {
    shape_str = program.get<std::string>("--enclose");
    enclose = true;
  }

  // read the geometry
  Geometry geometry(GeometryKernel::kEgads, input);

  const std::map<std::string, CenteringMethod> center_map = {
      {"cx", CenteringMethod::kCenterX}, {"cy", CenteringMethod::kCenterY},
      {"cz", CenteringMethod::kCenterZ}, {"X", CenteringMethod::kBoxMaxX},
      {"x", CenteringMethod::kBoxMinX},  {"Y", CenteringMethod::kBoxMaxY},
      {"y", CenteringMethod::kBoxMinY},  {"Z", CenteringMethod::kBoxMaxZ},
      {"z", CenteringMethod::kBoxMinZ},  {"", CenteringMethod::kModel}};

  const std::map<std::string, OperationShape> shape_map = {
      {"sphere", OperationShape::kSphere}, {"box", OperationShape::kBox}};

  auto get_value = [](const std::string& s, Operation& op) {
    if (s.back() == 'a') {
      op.scaling = LengthScaling::kAbsolute;
      return std::atoi(s.substr(0, s.size() - 1).c_str());
    }
    op.scaling = LengthScaling::kRelative;
    return std::atoi(s.c_str());
  };

  // option to enclose the geometry within an outer shape (add a farfield)
  if (enclose) {
    Operation op;
    op.type = OperationType::kSubtraction;
    op.shape = shape_map.at(shape_str);
    if (!radius_str.empty()) {
      ASSERT(op.shape == OperationShape::kSphere);
      double radius = get_value(radius_str, op);
      op.set_radius(radius);
    }
    if (!lengths_str.empty()) {
      ASSERT(op.shape == OperationShape::kBox);
      ASSERT(lengths_str.size() == 3);
      double lengths[3];
      for (int d = 0; d < 3; d++) {
        auto s = op.scaling;
        lengths[d] = get_value(lengths_str[d], op);
        if (d > 0) ASSERT(op.scaling == s);
      }
      op.set_box_lengths(lengths[0], lengths[1], lengths[2]);
    }
    op.centering = center_map.at(center_info);
    op.existing_model_is_source = false;
    geometry.operate(op);
  }

  // option to add a symmetry plane
  if (!symmetry.empty()) {
    Operation sym;
    sym.type = OperationType::kSubtraction;
    sym.shape = OperationShape::kSymmetry;
    sym.centering = center_map.at("c" + symmetry);
    sym.existing_model_is_source = true;
    geometry.operate(sym);
  }

  // write the geometry
  LOG << "writing geometry [" << output << "]";
  std::filesystem::remove(output);
  geometry.save(output);
}

void run_mesher(const argparse::ArgumentParser& program) {
  auto input = program.get<std::string>("input");
  auto file_info = file_ext(input);
  bool is_geometry = is_geometry_file(input);
  Mesh mesh(3);
  if (is_geometry) {
    // load the geometry
    Geometry geometry(GeometryKernel::kEgads, input);

    // tessellate
    // TODO control tessellation options
    geometry.tessellate(mesh);

  } else {
    read_mesh(input, mesh);
  }

  std::string output;
  if (program.is_used("--output"))
    output = program.get<std::string>("--output");
  else
    output = file_info.first + ".output.meshb";

  auto surface_mode = program.get<std::string>("surface");
  LOG << "surface mode: " << surface_mode;  // TODO
  auto volume_mode = program.get<std::string>("volume");
  LOG << "volume mode: " << volume_mode;
  bool minimal = volume_mode == "minimal";

  TetGenMesher mesher(mesh);
  MeshingParameters params{.constrained = true};
  params.verbose = true;
  if (minimal && !program.is_used("--quality"))
    params.quality = -1;
  else {
    params.quality = program.get<double>("--quality");
  }
  mesher.generate(params);

  LOG << "writing mesh [" << output << "]";
  meshb::write(mesh, output);
}

}  // namespace helix

int main(int argc, const char** argv) {
  std::signal(SIGINT, signal_handler);
  argparse::ArgumentParser program("helix", "1.0");

  argparse::ArgumentParser cmd_cad("cad");
  cmd_cad.add_description("perform CAD operations on a geometry");
  cmd_cad.add_argument("input").help("input geometry");
  cmd_cad.add_argument("-o", "--output").help("output geometry file");
  cmd_cad.add_argument("--enclose")
      .default_value("")
      .help("enclose the input geometry within a farfield");
  cmd_cad.add_argument("--subtract")
      .default_value("")
      .help("subtract a geometry from the input");
  cmd_cad.add_argument("--radius")
      .help("radius")
      .default_value("1.0")
      .metavar("R");
  cmd_cad.add_argument("--lengths")
      .help("lengths")
      .nargs(3)
      .default_value(std::vector<std::string>{"1.0", "1.0", "1.0"});
  cmd_cad.add_argument("--height")
      .help("height")
      .default_value("1.0")
      .metavar("H");
  cmd_cad.add_argument("--center")
      .help("center on bounding box: [x, X], [y, Y], [z, Z], or cx, cy, cz")
      .default_value("")
      .metavar("C");
  cmd_cad.add_argument("--axis")
      .help("axis")
      .default_value(std::vector<double>{0.0, 0.0, 1.0})
      .metavar("A")
      .scan<'g', double>();
  cmd_cad.add_argument("--symmetry")
      .metavar("PLANE")
      .default_value("")
      .help("x, y, z");

  argparse::ArgumentParser cmd_mesh("mesh");
  cmd_mesh.add_description("create a mesh within a domain");
  cmd_mesh.add_argument("input").help("input mesh or geometry file");
  cmd_mesh.add_argument("--g", "--geometry")
      .help("geometry file attached to mesh");
  cmd_mesh.add_argument("--surface")
      .default_value("tess")
      .help("mode used to mesh the surface: tess/input (default), metric")
      .metavar("MODE");
  cmd_mesh.add_argument("--volume")
      .default_value("")
      .help("mode used to mesh the volume: quality (default), minimal")
      .metavar("MODE");
  cmd_mesh.add_argument("--quality")
      .default_value(1.8)
      .help("target mesh quality")
      .scan<'g', double>();
  cmd_mesh.add_argument("-o", "--output").help("output mesh file");

  program.add_subparser(cmd_cad);
  program.add_subparser(cmd_mesh);

  try {
    program.parse_args(argc, argv);
  } catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  if (program.is_subcommand_used("cad")) {
    helix::run_cad(program.at<argparse::ArgumentParser>("cad"));
  } else if (program.is_subcommand_used("mesh")) {
    helix::run_mesher(program.at<argparse::ArgumentParser>("mesh"));
  } else {
    std::cout << program.help().str();
  }

  return 0;
}