#include "elements.h"

namespace helix {

int Tet::edges[12] = {0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3};
int Tet::faces[12] = {1, 2, 3, 2, 0, 3, 0, 1, 3, 0, 2, 1};

int Quad::faces[8] = {0, 1, 1, 2, 2, 3, 3, 0};
int Triangle::faces[6] = {1, 2, 2, 0, 0, 1};
int Triangle::edges[6] = {0, 1, 1, 2, 2, 0};

}