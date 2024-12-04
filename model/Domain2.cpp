#include "Domain2.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// putting boundary.
void Domain2d::setBoundary(double e1Min, double e1Max, double e2Min, double e2Max, const std::function<double(std::vector<double>)>& g,
                           const std::function<double(std::vector<double>)>& q, std::array<std::string, 4> boundaryConditions) {
    this->e1Min = e1Min;
    this->e1Max = e1Max;
    this->e2Min = e2Min;
    this->e2Max = e2Max;
    this->g = g;
    this->q = q;
    this->boundaryConditions = boundaryConditions;
}

// setting nodes. distributes equally for now (it might be updated soon).
void Domain2d::setNodes(int N, int M) {
    this->nodes.clear();
    this->u.clear();
    this->nodes.reserve(N * M);
    this->u.resize(N * M, 0.0);

    double e1Step = (this->e1Max - this->e1Min) / (N - 1);
    double e2Step = (this->e2Max - this->e2Min) / (M - 1);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            int nodeId = i * M + j;
            double e1 = this->e1Min + j * e1Step;
            double e2 = this->e2Min + i * e2Step;
            std::vector<double> coordinates = {e1, e2};

            // conditional branch by borders
            // start with Neumann BC then Dirichlet BC
            if (e2 == e2Min) {
                if (this->boundaryConditions[0] == "neumann" && e1 < this->e1Max) {
                    const std::vector<int> integralLine = {nodeId, nodeId + 1};
                    this->neumannBoundaryNodes.emplace_back(integralLine);
                } else if (this->boundaryConditions[0] == "dirichlet") {
                    this->dirichletBoundaryNodes.emplace_back(nodeId);
                    this->u[nodeId] = this->g(coordinates);
                }
            }
            if (e1 == this->e1Max) {
                if (this->boundaryConditions[1] == "neumann" && e2 < this->e2Max) {
                    const std::vector<int> integralLine = {nodeId, nodeId + M};
                    this->neumannBoundaryNodes.emplace_back(integralLine);
                } else if (this->boundaryConditions[1] == "dirichlet") {
                    this->dirichletBoundaryNodes.emplace_back(nodeId);
                    this->u[nodeId] = this->g(coordinates);
                }
            }
            if (e2 == this->e2Max) {
                if (this->boundaryConditions[2] == "neumann" && e1 > this->e1Min) {
                    const std::vector<int> integralLine = {nodeId, nodeId - 1};
                    this->neumannBoundaryNodes.emplace_back(integralLine);
                } else if (this->boundaryConditions[2] == "dirichlet") {
                    this->dirichletBoundaryNodes.emplace_back(nodeId);
                    this->u[nodeId] = this->g(coordinates);
                }
            }
            if (e1 == this->e1Min) {
                if (this->boundaryConditions[3] == "neumann" && e2 > this->e2Min) {
                    const std::vector<int> integralLine = {nodeId, nodeId - M};
                    this->neumannBoundaryNodes.emplace_back(integralLine);
                } else if (this->boundaryConditions[3] == "dirichlet") {
                    this->dirichletBoundaryNodes.emplace_back(nodeId);
                    this->u[nodeId] = this->g(coordinates);
                }
            }

            this->nodes.emplace_back(Node(nodeId, coordinates));
        }
    }

    std::cout << std::endl;
}

Node Domain2d::findNode(int nodeId) const {
    for (const auto& node : this->nodes) {
        if (node.id == nodeId) {
            return node;
        }
    }
    throw std::invalid_argument("Node not found");
}

void Domain2d::setElements(int N, int M, std::string interpolationMethod) {
    this->elements.clear();
    this->elements.reserve((N - 1) * (M - 1) * 2);

    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < M - 1; ++j) {
            int elementId = j * (M - 1) + i;
            Node node1 = findNode(j * M + i);
            Node node2 = findNode(j * M + i + 1);
            Node node4 = findNode((j + 1) * M + i + 1);
            int elementType = 5;  // VTK_TRIANGLE

            const std::vector<Node>& upsideTriangle = {node1, node2, node4};

            this->elements.emplace_back(Element(elementId, elementType, upsideTriangle));
        }
    }

    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < M - 1; ++j) {
            int elementId = j * (M - 1) + i + (N - 1) * (M - 1);
            Node node1 = findNode(j * M + i);
            Node node3 = findNode((j + 1) * M + i);
            Node node4 = findNode((j + 1) * M + i + 1);
            int elementType = 5;  // VTK_TRIANGLE

            const std::vector<Node>& downsideTriangle = {node1, node4, node3};

            this->elements.emplace_back(Element(elementId, elementType, downsideTriangle));
        }
    }
    this->interpolationMethod = interpolationMethod;
}

void Domain2d::printStatus(std::string fileName) const {
    std::ofstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing: " << fileName << std::endl;
        return;
    }
    file << std::setprecision(4) << std::fixed << std::scientific;

    file << "Domain Status:" << std::endl;
    file << std::endl;
    file << "e1 range: [" << this->e1Min << ", " << this->e1Max << "]" << std::endl;
    file << "e2 range: [" << this->e2Min << ", " << this->e2Max << "]" << std::endl;
    file << std::endl;
    file << "Number of nodes: " << this->nodes.size() << std::endl;
    file << "Number of elements: " << this->elements.size() << std::endl;
    file << std::endl;

    file << "Dirichlet boundary nodes:" << std::endl;
    file << "[";
    if (!this->dirichletBoundaryNodes.empty()) {
        for (size_t i = 0; i < this->dirichletBoundaryNodes.size() - 1; ++i) {
            file << this->dirichletBoundaryNodes[i] << ", ";
        }
        file << this->dirichletBoundaryNodes.back();
    }
    file << "]" << std::endl;
    file << std::endl;

    file << "Neumann boundary nodes:" << std::endl;
    file << "[" << std::endl;
    for (size_t i = 0; i < this->neumannBoundaryNodes.size(); ++i) {
        const auto& nodeSegment = this->neumannBoundaryNodes[i];
        file << "    [";
        for (size_t j = 0; j < nodeSegment.size(); ++j) {
            file << nodeSegment[j];
            if (j < nodeSegment.size() - 1) {
                file << ", ";
            }
        }
        file << "]";
        if (i < this->neumannBoundaryNodes.size() - 1) {
            file << ",";  // Add comma here for elements within the array
        }
        file << std::endl;
    }
    file << "           ]" << std::endl;
}

void Domain2d::outputVTK(std::string fileName) const {
    std::ofstream file(fileName);

    file << std::setprecision(2) << std::fixed << std::scientific;

    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing: " << fileName << std::endl;
        return;
    }

    file << "# vtk DataFile Version 3.0\n";
    file << "VTK file example\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << this->nodes.size() << " float\n";

    for (const auto& node : this->nodes) {
        file << node.coordinates[0] << " " << node.coordinates[1] << " " << node.coordinates[2] << "\n";
    }

    file << "CELLS " << this->elements.size() << " " << this->elements.size() * 4 << "\n";  // 4 = 1 (num of nodes in element) + 3 (nodes' indices)
    for (const auto& element : this->elements) {
        file << "3 " << element.nodes[0].id << " " << element.nodes[1].id << " " << element.nodes[2].id << "\n";
    }

    file << "CELL_TYPES " << this->elements.size() << "\n";
    for (const auto& element : this->elements) {
        file << element.type << "\n";
    }

    file << "POINT_DATA " << this->nodes.size() << "\n";
    file << "SCALARS node_scalars float\n";
    file << "LOOKUP_TABLE default\n";
    for (const double& value : this->u) {
        file << value << "\n";
    }

    file.close();
}