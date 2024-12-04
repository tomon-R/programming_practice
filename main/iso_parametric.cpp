#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "../model/Domain2.h"
#include "../solver/Matrix.h"
#include "../src/outputMatrix.h"

#define N 30   // Number of numbers along x1
#define M 30   // Number of numbers along x2
#define G 1.0  // Dirichlet Boundary Conditions
#define Q 0.0  // Neumann Boundary Conditions
#define F 0.0  // External Force

namespace fs = std::filesystem;

// Simple functions for e1, e2, g, q
double e1Func(std::vector<double> x) { return x[0]; }
double e2Func(std::vector<double> x) { return x[1]; }
double gFunc(std::vector<double> x) {
    static const double pi = 3.141592653589793;
    if (x[0] == 1.0) {
        return G * sin(pi * x[1]);
        // return G;
    } else {
        return 0;
    }

    // return G;
}
double qFunc(std::vector<double> x) { return Q; }

double theoreticalValue(std::vector<double> x_vec) {
    static const double pi = 3.141592653589793;
    double x = x_vec[0];
    double y = x_vec[1];

    return sinh(pi * x) * sin(pi * y) / sinh(pi);
}

Matrix buildK(const std::vector<Node>& nodes, const std::vector<Element>& elements, const std::string& interpolationMethod,
              const std::vector<int>& dirichletBoundaryNodes);
std::vector<double> buildB(const std::vector<Node>& nodes, const std::vector<std::vector<int>>& neumannBoundaryNodes,
                           std::function<double(std::vector<double>)>& q, const std::vector<int>& dirichletBoundaryNodes,
                           std::function<double(std::vector<double>)>& g);
std::vector<double> buildF(const std::vector<Node>& nodes, const std::vector<Element>& elements, const std::string& interpolationMethod,
                           const std::vector<int>& dirichletBoundaryNodes);

int main() {
    std::function<double(std::vector<double>)> e1 = e1Func;
    std::function<double(std::vector<double>)> e2 = e2Func;
    std::function<double(std::vector<double>)> g = gFunc;
    std::function<double(std::vector<double>)> q = qFunc;

    Domain2d domain(e1, e2);
    std::array<std::string, 4> boundaryConditions = {"dirichlet", "dirichlet", "dirichlet", "dirichlet"};
    domain.setBoundary(0.0, 1.0, 0.0, 1.0, g, q, boundaryConditions);
    domain.setNodes(N, M);
    domain.setElements(N, M, "linear");

    Matrix K = buildK(domain.nodes, domain.elements, domain.interpolationMethod, domain.dirichletBoundaryNodes);
    std::vector<double> b = buildB(domain.nodes, domain.neumannBoundaryNodes, q, domain.dirichletBoundaryNodes, g);
    std::vector<double> f = buildF(domain.nodes, domain.elements, domain.interpolationMethod, domain.dirichletBoundaryNodes);

    std::vector<double> RHS(f.size());
    for (size_t i = 0; i < f.size(); ++i) {
        RHS[i] = f[i] + b[i];
    }

    std::vector<double> u = K.solve(RHS);

    domain.u = u;

    fs::path outputDir = OUTPUT_DIR;

    std::string kFile = (outputDir / "ISO/ISO_K.txt").string();
    std::string bFile = (outputDir / "ISO/ISO_b.txt").string();
    std::string rhsFile = (outputDir / "ISO/ISO_rhs.txt").string();
    std::string fFile = (outputDir / "ISO/ISO_f.txt").string();
    std::string uFile = (outputDir / "ISO/ISO_u.txt").string();

    outputMatrix(kFile, K);
    outputVector(bFile, b);
    outputVector(rhsFile, RHS);
    outputVector(fFile, f);
    outputVector(uFile, u);

    std::string statusFilePath = (outputDir / "ISO/ISO_domain_status.txt").string();
    std::string vtkFilePath = (outputDir / "ISO/ISO_domain.vtk").string();

    domain.printStatus(statusFilePath);
    domain.outputVTK(vtkFilePath);

    std::vector<double> theoretical(domain.nodes.size());
    for (size_t i = 0; i < theoretical.size(); i++) {
        Node n = domain.nodes[i];
        theoretical[i] = theoreticalValue(n.coordinates);
    }

    domain.u = theoretical;

    std::string theoreticalFile = (outputDir / "ISO/ISO_theoretical.vtk").string();
    domain.outputVTK(theoreticalFile);

    std::vector<double> err(domain.nodes.size());
    for (size_t i = 0; i < domain.nodes.size(); ++i) {
        err[i] = std::abs(u[i] - theoretical[i]);
    }

    domain.u = err;

    std::string errFile = (outputDir / "ISO/ISO_err.vtk").string();
    domain.outputVTK(errFile);

    return 0;
}

Matrix buildK(const std::vector<Node>& nodes, const std::vector<Element>& elements, const std::string& interpolationMethod,
              const std::vector<int>& dirichletBoundaryNodes) {
    Matrix K(nodes.size(), nodes.size());
    Matrix Ke(3, 3);

    for (const auto& element : elements) {
        std::cout << "element: " << element.id << std::endl;

        // VTK_TRIANGLE
        if (element.type == 5) {
            const Node& node0 = element.nodes[0];
            const double x0 = node0.coordinates[0];
            const double y0 = node0.coordinates[1];
            const Node& node1 = element.nodes[1];
            const double x1 = node1.coordinates[0];
            const double y1 = node1.coordinates[1];
            const Node& node2 = element.nodes[2];
            const double x2 = node2.coordinates[0];
            const double y2 = node2.coordinates[1];

            std::cout << "node0(" << node0.id << "): (" << x0 << ", " << y0 << ")" << std::endl;
            std::cout << "node1(" << node1.id << "): (" << x1 << ", " << y1 << ")" << std::endl;
            std::cout << "node2(" << node2.id << "): (" << x2 << ", " << y2 << ")" << std::endl;

            double hx = x1 - x0;
            double hy = y1 - y0;
            std::cout << "hx: " << hx << std::endl;
            std::cout << "hy: " << hx << std::endl;

            // LINEAR INTERPOLATION
            if (interpolationMethod == "linear") {
                double D = x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1);
                double a0 = (x1 * y2 - x2 * y1) / D;
                double a1 = (x2 * y0 - x0 * y2) / D;
                double a2 = (x0 * y1 - x1 * y0) / D;
                double b0 = (y1 - y2) / D;
                double b1 = (y2 - y0) / D;
                double b2 = (y0 - y1) / D;
                double c0 = (x2 - x1) / D;
                double c1 = (x0 - x2) / D;
                double c2 = (x1 - x0) / D;

                std::cout << "D: " << D << std::endl;

                std::cout << "b:[" << b0 << ", " << b1 << ", " << b2 << "]" << std::endl;
                std::cout << "c:[" << c0 << ", " << c1 << ", " << c2 << "]" << std::endl;

                Ke[0][0] = 2.0 / 2.0;
                Ke[0][1] = -1.0 / 2.0;
                Ke[0][2] = -1.0 / 2.0;
                Ke[1][0] = -1.0 / 2.0;
                Ke[1][1] = 1.0 / 2.0;
                Ke[1][2] = 0.0 / 2.0;
                Ke[2][0] = -1.0 / 2.0;
                Ke[2][1] = 0.0 / 2.0;
                Ke[2][2] = 1.0 / 2.0;

                fs::path outputDir = OUTPUT_DIR;
                std::string KeFile = (outputDir / "ISO/ISO_Ke.txt").string();
                outputMatrix(KeFile, Ke);

                K[node0.id][node0.id] += Ke[0][0];
                K[node0.id][node1.id] += Ke[0][1];
                K[node0.id][node2.id] += Ke[0][2];
                K[node1.id][node0.id] += Ke[1][0];
                K[node1.id][node1.id] += Ke[1][1];
                K[node1.id][node2.id] += Ke[1][2];
                K[node2.id][node0.id] += Ke[2][0];
                K[node2.id][node1.id] += Ke[2][1];
                K[node2.id][node2.id] += Ke[2][2];
            }
        }
    }

    for (int j = 0; j < K.getCols(); ++j) {
        auto isDirichlet = [j, &dirichletBoundaryNodes]() {
            for (const int& dirichletNode : dirichletBoundaryNodes) {
                if (j == dirichletNode) {
                    return true;
                }
            }
            return false;
        };

        if (isDirichlet()) {
            // for (int index = 0; index < K.getRows(); ++index) {
            //     K[index][j] = 0.0;
            // }
            for (int index = 0; index < K.getCols(); ++index) {
                K[j][index] = 0.0;
            }
            K[j][j] = 1.0;
        }
    }
    return K;
}

std::vector<double> buildB(const std::vector<Node>& nodes, const std::vector<std::vector<int>>& neumannBoundaryNodes,
                           std::function<double(std::vector<double>)>& q, const std::vector<int>& dirichletBoundaryNodes,
                           std::function<double(std::vector<double>)>& g) {
    std::vector<double> b(nodes.size(), 0.0);

    for (const auto& neumannBoundary : neumannBoundaryNodes) {
        int nodeNumber0 = neumannBoundary[0];
        Node node0 = nodes[nodeNumber0];
        double x0 = node0.coordinates[0];
        double y0 = node0.coordinates[1];
        int nodeNumber1 = neumannBoundary[1];
        Node node1 = nodes[nodeNumber1];
        double x1 = node1.coordinates[0];
        double y1 = node1.coordinates[1];
        double l = std::sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));

        double b0 = 1.0 / 2.0 * l * q(node0.coordinates);
        double b1 = 1.0 / 2.0 * l * q(node1.coordinates);

        b[nodeNumber0] += b0;
        b[nodeNumber1] += b1;

        std::cout << "node0(" << node0.id << "): (" << x0 << ", " << y0 << ")" << std::endl;
        std::cout << "node1(" << node1.id << "): (" << x1 << ", " << y1 << ")" << std::endl;

        std::cout << "b: [" << b0 << ", " << b1 << "]" << std::endl;
    }

    for (int j = 0; j < b.size(); ++j) {
        auto isDirichlet = [j, &dirichletBoundaryNodes]() {
            for (const int& dirichletNode : dirichletBoundaryNodes) {
                if (j == dirichletNode) {
                    return true;
                }
            }
            return false;
        };

        if (isDirichlet()) {
            b[j] = g(nodes[j].coordinates);
        }
    }
    return b;
}

std::vector<double> buildF(const std::vector<Node>& nodes, const std::vector<Element>& elements, const std::string& interpolationMethod,
                           const std::vector<int>& dirichletBoundaryNodes) {
    std::vector<double> f(nodes.size());
    std::vector<double> fe(3);

    for (const auto& element : elements) {
        // VTK_TRIANGLE
        if (element.type == 5) {
            const Node& node0 = element.nodes[0];
            const double x0 = node0.coordinates[0];
            const double y0 = node0.coordinates[1];
            const Node& node1 = element.nodes[1];
            const double x1 = node1.coordinates[0];
            const double y1 = node1.coordinates[1];
            const Node& node2 = element.nodes[2];
            const double x2 = node2.coordinates[0];
            const double y2 = node2.coordinates[1];

            // LINEAR INTERPOLATION
            if (interpolationMethod == "linear") {
                double D = x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1);
                double S = std::abs(D) / 2;
                std::cout << "S: " << S << std::endl;

                fe[0] = F * S / 3;
                fe[1] = F * S / 3;
                fe[2] = F * S / 3;

                f[node0.id] += fe[0];
                f[node1.id] += fe[1];
                f[node2.id] += fe[2];

                fs::path outputDir = OUTPUT_DIR;
                std::string feFile = (outputDir / "ISO/ISO_fe.txt").string();
                outputVector(feFile, fe);
            }
        }
    }

    for (int j = 0; j < f.size(); ++j) {
        auto isDirichlet = [j, &dirichletBoundaryNodes]() {
            for (const int& dirichletNode : dirichletBoundaryNodes) {
                if (j == dirichletNode) {
                    return true;
                }
            }
            return false;
        };

        if (isDirichlet()) {
            f[j] = 0;
        }
    }
    return f;
}