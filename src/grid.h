#pragma once

#include <vector>

#include "matrix.h"

class IMolecule;

struct GridPoint {
  scfcpp::Vector3 pos;
  double weight;
};

std::vector<GridPoint> generate_molecular_grid(const IMolecule& mol,
                                               int radial_points,
                                               int angular_points);
