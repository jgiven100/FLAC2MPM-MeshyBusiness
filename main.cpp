#include "external/json.hpp"
using json = nlohmann::json;

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Check compatibility
void check_compatibility(int &nelem_x, int &nelem_y, const double xmin,
                         const double xmax, const double ymin,
                         const double ymax, const double he);

// Save mesh txt file
void write_mesh(const int nelem_x, const int nelem_y, const double xmin,
                const double ymin, const double he);

// Save entity set JSON file
void write_entity_set(const int nnode_x, const int nnode_y, const double he);

// Read FLAC file
void load_flac(std::string fname,
               std::map<int, std::pair<double, double>> &gridpoints,
               std::map<int, std::tuple<int, int, int, int>> &q_zones,
               std::map<int, std::tuple<int, int, int>> &t_zones,
               std::vector<std::vector<int>> &zone_groups,
               std::vector<std::string> &name_zone_groups,
               std::vector<std::tuple<int, int, int>> &faces,
               std::vector<std::vector<int>> &face_groups);

// Find particle
bool find_particle_quad(std::map<int, std::pair<double, double>> &gridpoints,
                        std::tuple<int, int, int, int> &nodes,
                        std::pair<double, double> &particle_coord);

// Find particle
bool find_particle_tri(std::map<int, std::pair<double, double>> &gridpoints,
                       std::tuple<int, int, int> &nodes,
                       std::pair<double, double> &particle_coord);

// Save particle coordinates txt file
void write_particle_coords(
    const int id,
    const std::vector<std::pair<double, double>> &particle_coords);

// Save particle cells txt file
void write_particle_cells(const std::vector<int> &particle_cells);

// Save particle stresses txt file
void write_particle_stresses();

// Save particle stresses beginning txt file
void write_particle_stresses_beginning();

// Save particle volumes txt file
void write_particle_volumes(
    const std::vector<std::pair<double, double>> &particle_coords,
    const double vol);

// Generate mesh and node sets
void generate_mesh_files(const double xmin, const double xmax,
                         const double ymin, const double ymax,
                         const double he) {
  // Number of elements in each direction
  int nelem_x, nelem_y;
  check_compatibility(nelem_x, nelem_y, xmin, xmax, ymin, ymax, he);

  // Number of nodes in each direction
  const int nnode_x = nelem_x + 1;
  const int nnode_y = nelem_y + 1;

  // Generate mesh file
  write_mesh(nelem_x, nelem_y, xmin, ymin, he);

  // Generate entity set file
  write_entity_set(nnode_x, nnode_y, he);
}

void check_compatibility(int &nelem_x, int &nelem_y, const double xmin,
                         const double xmax, const double ymin,
                         const double ymax, const double he) {
  // Check x-dir compatibility
  const double deltax = xmax - xmin;
  double remainder = std::fmod(deltax, he);
  if (remainder != 0.) {
    std::cout << "Bad combination of x-dir dimensions and he... quitting"
              << std::endl;
    exit(1);
  }

  // Check y-dir compatibility
  const double deltay = ymax - ymin;
  remainder = std::fmod(deltay, he);
  if (remainder != 0.) {
    std::cout << "Bad combination of y-dir dimensions and he... quitting"
              << std::endl;
    exit(1);
  }

  // Number of elements in each direction
  nelem_x = static_cast<int>(deltax / he);
  nelem_y = static_cast<int>(deltay / he);
}

void write_mesh(const int nelem_x, const int nelem_y, const double xmin,
                const double ymin, const double he) {
  // Number of nodes
  const int nnode_x = nelem_x + 1;
  const int nnode_y = nelem_y + 1;
  const int nnode = nnode_x * nnode_y;

  // Nodes
  std::vector<std::pair<double, double>> nodes;
  nodes.resize(nnode, std::make_pair(0., 0.));
  for (int i = 0; i < nnode_y; i++) {
    for (int j = 0; j < nnode_x; j++) {
      const int index = (i * nnode_x) + j;
      nodes[index] = std::make_pair(static_cast<double>(j) * he + xmin,
                                    static_cast<double>(i) * he + ymin);
    }
  }

  // Total number of elements
  const int nelem = nelem_x * nelem_y;

  // Elements
  std::vector<std::tuple<int, int, int, int>> elems;
  elems.resize(nelem, std::make_tuple(0, 0, 0, 0));
  for (int i = 0; i < nelem_y; i++) {
    for (int j = 0; j < nelem_x; j++) {
      const int index = (i * nelem_x) + j;
      elems[index] = std::make_tuple((i * nnode_x) + j, (i * nnode_x) + j + 1,
                                     (i * nnode_x) + j + 1 + nnode_x,
                                     (i * nnode_x) + j + nnode_x);
    }
  }

  // Write mesh file
  std::ofstream outfile("mesh.txt");
  if (outfile.is_open()) {
    // Write header
    outfile << nnode << " " << nelem << std::endl;
    outfile << std::fixed << std::setprecision(10);

    // Write nodes
    for (std::pair<double, double> node : nodes)
      outfile << node.first << " " << node.second << std::endl;

    // Write elements
    for (std::tuple<int, int, int, int> elem : elems)
      outfile << std::get<0>(elem) << " " << std::get<1>(elem) << " "
              << std::get<2>(elem) << " " << std::get<3>(elem) << std::endl;

    // Close
    outfile.close();
  } else {
    std::cout << "Unable to open mesh file for saving" << std::endl;
    exit(1);
  }
}

void write_entity_set(const int nnode_x, const int nnode_y, const double he) {
  // Left bondary
  std::vector<int> left;
  left.reserve(nnode_y);
  for (int i = 0; i < nnode_y; i++) {
    left.push_back(i * nnode_x);
  }

  // Right boundary
  std::vector<int> right;
  right.reserve(nnode_y);
  for (int i = 0; i < nnode_y; i++) {
    right.push_back((i * nnode_x) + (nnode_x - 1));
  }

  // Bottom boundary
  std::vector<int> bottom;
  bottom.reserve(nnode_x);
  for (int i = 0; i < nnode_x; i++) {
    bottom.push_back(i);
  }

  // Save
  std::vector<const std::vector<int> *> node_sets;
  node_sets.push_back(&left);
  node_sets.push_back(&right);
  node_sets.push_back(&bottom);

  // Write entity sets
  std::ofstream outfile("entity-sets.json");
  if (outfile.is_open()) {
    // Write header
    outfile << "{\n  \"node_sets\": [\n";

    // Save each node set
    int count = 0;
    for (const auto &node_set : node_sets) {

      outfile << "    {\n      \"id\": " << count << ",\n";
      outfile << "      \"set\": [\n";

      int ccount = 0;
      for (int node : *node_set) {
        ++ccount;

        outfile << "        " << node;
        if (count == 0 || count == 1) {
          if (ccount < nnode_y)
            outfile << ",";
        } else {
          if (ccount < nnode_x)
            outfile << ",";
        }
        outfile << "\n";
      }

      outfile << "      ]\n    }";
      if (count < 2)
        outfile << ",";
      outfile << "\n";

      ++count;
    }

    // Write footer
    outfile << "  ]\n}";

    // Close
    outfile.close();
  } else {
    std::cout << "Unable to open entity set file for saving" << std::endl;
    exit(1);
  }
}

void generate_particle_files(const double xmin, const double xmax,
                             const double ymin, const double ymax,
                             const double he, const int ppc,
                             const std::string fname) {
  // Number of elements in each direction
  int nelem_x, nelem_y;
  check_compatibility(nelem_x, nelem_y, xmin, xmax, ymin, ymax, he);

  // Set 1d particle coordinates [0, 1]
  std::vector<double> pts;
  for (int i = 0; i < ppc; i++) {
    pts.push_back((static_cast<double>(i) + 0.5) / static_cast<double>(ppc));
  }

  // Particle coordinates and cells
  std::map<int, std::pair<double, double>> dense_particle_coords;
  std::map<int, int> dense_particle_cells;

  // Particle and cell id
  int pid = 0;
  int cid = 0;

  // Initialize dense particle distribution in mesh
  for (int i = 0; i < nelem_y; i++) {
    for (int j = 0; j < nelem_x; j++) {
      for (auto pt_x : pts) {
        for (auto pt_y : pts) {
          // Global coordinates
          const double x = (static_cast<double>(j) * he) + (pt_x * he) + xmin;
          const double y = (static_cast<double>(i) * he) + (pt_y * he) + ymin;

          // Save particle
          dense_particle_coords[pid] = std::make_pair(x, y);
          dense_particle_cells[pid] = cid;

          ++pid;
        }
      }
      ++cid;
    }
  }

  // Load FLAC grippoints and zones
  std::map<int, std::pair<double, double>> gridpoints;
  std::map<int, std::tuple<int, int, int, int>> q_zones;
  std::map<int, std::tuple<int, int, int>> t_zones;
  std::vector<std::vector<int>> zone_groups;
  std::vector<std::string> name_zone_groups;
  std::vector<std::tuple<int, int, int>> faces;
  std::vector<std::vector<int>> face_groups;

  load_flac(fname, gridpoints, q_zones, t_zones, zone_groups, name_zone_groups,
            faces, face_groups);

  // Bounding box
  double xmin_flac = 1.E+6;
  double xmax_flac = -1.E+6;
  double ymin_flac = 1.E+6;
  double ymax_flac = -1.E+6;
  for (const auto &gridpoint : gridpoints) {
    if (gridpoint.second.first < xmin_flac)
      xmin_flac = gridpoint.second.first;
    if (gridpoint.second.first > xmax_flac)
      xmax_flac = gridpoint.second.first;
    if (gridpoint.second.second < ymin_flac)
      ymin_flac = gridpoint.second.second;
    if (gridpoint.second.second > ymax_flac)
      ymax_flac = gridpoint.second.second;
  }

  // Trim dense particles outside boudning box
  std::map<int, std::pair<double, double>> dense_particle_coords_copy =
      dense_particle_coords;
  for (auto pitr = dense_particle_coords_copy.begin();
       pitr != dense_particle_coords_copy.end(); ++pitr) {
    if (((pitr->second).first < xmin_flac) ||
        ((pitr->second).first > xmax_flac))
      dense_particle_coords.erase(pitr->first);

    if (((pitr->second).second < ymin_flac) ||
        ((pitr->second).second > ymax_flac))
      dense_particle_coords.erase(pitr->first);
  }

  // Loop each zone group
  for (unsigned i = 0; i < name_zone_groups.size(); i++) {
    std::cout << "particles-" << i << ".txt -> " << name_zone_groups[i]
              << std::endl;

    // Update copy of dense particles
    dense_particle_coords_copy = dense_particle_coords;

    // Get list of zones for current zone group
    const std::vector<int> zone_group = zone_groups[i];

    //
    std::vector<std::pair<double, double>> particle_coords_zone;
    particle_coords_zone.clear();

    // Loop each zone
    for (int zid : zone_group) {

      // Check if zone id is quad or triangle
      if (q_zones.find(zid) != q_zones.end()) {
        // Loop particles
        for (auto pitr = dense_particle_coords_copy.begin();
             pitr != dense_particle_coords_copy.end(); ++pitr) {
          if (find_particle_quad(gridpoints, q_zones[zid], pitr->second)) {
            particle_coords_zone.emplace_back(pitr->second);
            dense_particle_coords.erase(pitr->first);
          }
        }

      } else if (t_zones.find(zid) != t_zones.end()) {
        // Loop particles
        for (auto pitr = dense_particle_coords_copy.begin();
             pitr != dense_particle_coords_copy.end(); ++pitr) {
          if (find_particle_tri(gridpoints, t_zones[zid], pitr->second)) {
            particle_coords_zone.emplace_back(pitr->second);
            dense_particle_coords.erase(pitr->first);
          }
        }

      } else {
        std::cout << "Unable to find zone id" << std::endl;
        exit(1);
      }
    }
    write_particle_coords(i, particle_coords_zone);
  }

  // TODO REMOVE ME
  // for (const auto &gridpoint : gridpoints)
  //   std::cout << std::fixed << std::setprecision(10) << gridpoint.first << "
  //   " << gridpoint.second.first << " " << gridpoint.second.second <<
  //   std::endl;

  // for (const auto &q_zone : q_zones)
  //   std::cout << q_zone.first << " " << std::get<0>(q_zone.second) << " " <<
  //   std::get<1>(q_zone.second) << " " << std::get<2>(q_zone.second) << " " <<
  //   std::get<3>(q_zone.second) << std::endl;

  // for (const auto &t_zone : t_zones)
  //   std::cout << t_zone.first << " " << std::get<0>(t_zone.second) << " " <<
  //   std::get<1>(t_zone.second) << " " << std::get<2>(t_zone.second) <<
  //   std::endl;

  std::vector<std::tuple<double, double, double, double, double, double>>
      particle_stresses;
  std::vector<std::tuple<double, double, double, double, double, double>>
      particle_stresses_beginning;

  // write_particle_cells(particle_cells);

  // // write_particle_stresses();

  // // write_particle_stresses_beginning();

  // // Particle volume
  // const double vol = (he * he) / (ppc * ppc);
  // write_particle_volumes(particle_coords, vol);
}

void load_flac(std::string fname,
               std::map<int, std::pair<double, double>> &gridpoints,
               std::map<int, std::tuple<int, int, int, int>> &q_zones,
               std::map<int, std::tuple<int, int, int>> &t_zones,
               std::vector<std::vector<int>> &zone_groups,
               std::vector<std::string> &name_zone_groups,
               std::vector<std::tuple<int, int, int>> &faces,
               std::vector<std::vector<int>> &face_groups) {
  // Read FLAC
  std::ifstream infile(fname);
  if (infile.is_open()) {
    // Read line by line
    std::string line;
    while (std::getline(infile, line)) {
      // Ignore comments
      if ((line.find("*") == std::string::npos) && (line != "")) {
        //
        std::istringstream iss(line);

        if (line[0] == 'G') {
          // Grid points
          char c;
          int id;
          double x, y;
          iss >> c >> id >> x >> y;
          gridpoints[id] = std::make_pair(x, y);

        } else if ((line[0] == 'Z') and
                   (line.find("ZGROUP") == std::string::npos)) {
          // Zones
          char c;
          std::string type;
          int id, n0, n1, n2, n3;
          iss >> c >> type >> id >> n0 >> n1 >> n2 >> n3;

          // Quads versus triangles
          if (type == "Q4") {
            q_zones[id] = std::make_tuple(n0, n1, n2, n3);
          } else if (type == "T3") {
            t_zones[id] = std::make_tuple(n0, n1, n2);
          }

        } else if (line.find("ZGROUP") != std::string::npos) {
          // Zone groups
          std::string type, zone_name, slot, slot_name;

          iss >> type;
          // Handle zone group name (which may have spaces)
          zone_name.clear();
          char c;
          if (iss >> c && c == '"') {
            while (iss.get(c) && c != '"') {
              zone_name.push_back(c);
            }
          }
          iss >> slot;
          // Handle slot name (which may have spaces)
          slot_name.clear();
          if (iss >> c && c == '"') {
            while (iss.get(c) && c != '"') {
              slot_name.push_back(c);
            }
          }

          // Save "Construction" zones
          if (slot_name == "Construction") {

            std::vector<int> zone_group;

            bool get_next_line = true;
            while (get_next_line) {

              std::getline(infile, line);
              std::istringstream iss(line);

              int zid;
              while (iss >> zid) {
                zone_group.push_back(zid);
              }

              if ((infile.peek() == 'Z') || (infile.peek() == '*'))
                get_next_line = false;
            }

            // Save zone group names and ids
            name_zone_groups.emplace_back(zone_name);
            zone_groups.emplace_back(zone_group);
          }
        }
        // IGNORE: Faces
        // IGNORE: Face groups
      }
    }

    // Close
    infile.close();
  } else {
    std::cout << "Unable to open FLAC file" << std::endl;
    exit(1);
  }
}

bool find_particle_quad(std::map<int, std::pair<double, double>> &gridpoints,
                        std::tuple<int, int, int, int> &nodes,
                        std::pair<double, double> &particle_coord) {
  // Get nodes
  const int n0 = std::get<0>(nodes);
  const int n1 = std::get<1>(nodes);
  const int n2 = std::get<2>(nodes);
  const int n3 = std::get<3>(nodes);

  // Check triangle 0-1-2
  std::tuple<int, int, int> nodes_tmp = std::make_tuple(n0, n1, n2);
  const bool t012 = find_particle_tri(gridpoints, nodes_tmp, particle_coord);
  if (t012)
    return true;

  // Check triangle 1-3-2
  nodes_tmp = std::make_tuple(n1, n3, n2);
  const bool t132 = find_particle_tri(gridpoints, nodes_tmp, particle_coord);
  if (t132)
    return true;

  return false;
}

bool find_particle_tri(std::map<int, std::pair<double, double>> &gridpoints,
                       std::tuple<int, int, int> &nodes,
                       std::pair<double, double> &particle_coord) {
  // Get nodes
  const int n0 = std::get<0>(nodes);
  const int n1 = std::get<1>(nodes);
  const int n2 = std::get<2>(nodes);

  // Get points
  const std::pair<double, double> p0 = gridpoints[n0];
  const std::pair<double, double> p1 = gridpoints[n1];
  const std::pair<double, double> p2 = gridpoints[n2];

  const double area =
      0.5 * (-p1.second * p2.first + p0.second * (-p1.first + p2.first) +
             p0.first * (p1.second - p2.second) + p1.first * p2.second);
  const double s = 1 / (2 * area) *
                   (p0.second * p2.first - p0.first * p2.second +
                    (p2.second - p0.second) * particle_coord.first +
                    (p0.first - p2.first) * particle_coord.second);
  const double t = 1 / (2 * area) *
                   (p0.first * p1.second - p0.second * p1.first +
                    (p0.second - p1.second) * particle_coord.first +
                    (p1.first - p0.first) * particle_coord.second);

  return (s > 0 && t > 0 && (1 - s - t) > 0);
}

void write_particle_coords(
    const int id,
    const std::vector<std::pair<double, double>> &particle_coords) {
  // Number of particles
  const int nparticles = particle_coords.size();

  // Write particle coords
  std::string fname =
      std::string("particles-") + std::to_string(id) + std::string(".txt");
  std::ofstream outfile(fname);
  if (outfile.is_open()) {
    // Write header
    outfile << nparticles << std::endl;
    outfile << std::fixed << std::setprecision(10);

    // Write particle
    for (std::pair<double, double> particle_coord : particle_coords)
      outfile << particle_coord.first << " " << particle_coord.second
              << std::endl;

    // Close
    outfile.close();
  } else {
    std::cout << "Unable to open particle coordinate file for saving"
              << std::endl;
    exit(1);
  }
}

void write_particle_cells(const std::vector<int> &particle_cells) {
  // Write particle cells
  std::ofstream outfile("particle-cells.txt");
  if (outfile.is_open()) {
    // Write particle
    for (unsigned i = 0; i < particle_cells.size(); i++)
      outfile << i << " " << particle_cells[i] << std::endl;

    // Close
    outfile.close();
  } else {
    std::cout << "Unable to open particle cell file for saving" << std::endl;
    exit(1);
  }
}

void write_particle_volumes(
    const std::vector<std::pair<double, double>> &particle_coords,
    const double vol) {
  // Number of particles
  const int nparticles = particle_coords.size();

  // Write particle volumes
  std::ofstream outfile("particles-volumes.txt");
  if (outfile.is_open()) {
    // Write particle
    for (int i = 0; i < nparticles; i++)
      outfile << i << " " << std::fixed << std::setprecision(10) << vol
              << std::endl;

    // Close
    outfile.close();
  } else {
    std::cout << "Unable to open particle volume file for saving" << std::endl;
    exit(1);
  }
}

int main(int argc, char *argv[]) {

  // Check number of arguments
  if (argc != 3) {
    std::cout << "Invalid number of arguments... Usage: ./main -i input.json"
              << std::endl;
    return 1;
  }

  // Check flag
  if (std::string(argv[1]) != "-i") {
    std::cout << "Invalid flag... Usage: ./main -i input.json" << std::endl;
    return 1;
  }

  // Open JSON file
  std::ifstream file(argv[2]);
  if (!file.is_open()) {
    std::cout << "Input JSON not opened... Usage: ./main -i input.json"
              << std::endl;
    return 1;
  }

  // Save input JSON file as JSON object
  json input_json;
  try {
    file >> input_json;
  } catch (const std::exception &exception) {
    std::cout << "Failed to parse JSON: " << exception.what() << std::endl;
    return 1;
  }

  // Generate mesh with node sets
  try {
    // Get mesh parameters
    const double xmin = input_json.at("xmin").get<double>();
    const double xmax = input_json.at("xmax").get<double>();
    const double ymin = input_json.at("ymin").get<double>();
    const double ymax = input_json.at("ymax").get<double>();
    const double he = input_json.at("he").get<double>();

    // Generate mesh and node sets
    generate_mesh_files(xmin, xmax, ymin, ymax, he);
  } catch (const std::exception &exception) {
    std::cout << "Failed to find mesh settings: " << exception.what()
              << std::endl;
    return 1;
  }

  // Generate particle files
  try {
    // Get mesh parameters
    const double xmin = input_json.at("xmin").get<double>();
    const double xmax = input_json.at("xmax").get<double>();
    const double ymin = input_json.at("ymin").get<double>();
    const double ymax = input_json.at("ymax").get<double>();
    const double he = input_json.at("he").get<double>();
    const int ppc = input_json.at("ppc").get<double>();
    const std::string fname = input_json.at("file").get<std::string>();

    // Generate particle files
    generate_particle_files(xmin, xmax, ymin, ymax, he, ppc, fname);
  } catch (const std::exception &exception) {
    std::cout << "Failed to find particle settings: " << exception.what()
              << std::endl;
    return 1;
  }

  return 0;
}