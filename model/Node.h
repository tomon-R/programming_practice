#ifndef NODE_H
#define NODE_H

#include <string>
#include <vector>

struct Node {
    const int id;
    const std::vector<double> coordinates;

    Node(int id, const std::vector<double>& coordinates) : id(id), coordinates(coordinates) {}
};

#endif  // NODE_H
