#ifndef ELEMENT_H
#define ELEMENT_H

#include <string>
#include <vector>

#include "Node.h"

struct Element {
    const int id;
    // VSK_TYPE
    const int type;
    const std::vector<Node> nodes;

    Element(int id, int type, const std::vector<Node>& nodes) : id(id), type(type), nodes(nodes) {}
};

#endif  // ELEMENT_H
