#include <filesystem>
#include <iostream>

#include "../model/Domain2.h"
#include "../src/outputMatrix.h"

namespace fs = std::filesystem;

// Simple functions for e1, e2, g, q
double e1Func(std::vector<double> x) { return x[0]; }
double e2Func(std::vector<double> x) { return x[1]; }
double gFunc(std::vector<double> x) { return 5; }
double qFunc(std::vector<double> x) { return 5; }

int main() {
    std::function<double(std::vector<double>)> e1 = e1Func;
    std::function<double(std::vector<double>)> e2 = e2Func;
    std::function<double(std::vector<double>)> g = gFunc;
    std::function<double(std::vector<double>)> q = qFunc;

    Domain2d domain(e1, e2);

    std::array<std::string, 4> boundaryConditions = {"neumann", "neumann", "dirichlet", "neumann"};

    domain.setBoundary(0.0, 1.0, 0.0, 1.0, g, q, boundaryConditions);
    domain.setNodes(30, 30);  // 3x3 grid
    domain.setElements(30, 30, "linear");

    fs::path outputDir = OUTPUT_DIR;

    std::string statusFilePath = (outputDir / "domain_status.txt").string();
    std::string vtkFilePath = (outputDir / "domain.vtk").string();

    std::cout << "Attempting to write status file to: " << statusFilePath << std::endl;
    std::cout << "Attempting to write vtk file to: " << vtkFilePath << std::endl;

    domain.printStatus(statusFilePath);
    domain.outputVTK(vtkFilePath);

    std::string uFilePath = (outputDir / "u.txt").string();

    outputVector(uFilePath, domain.u);

    std::cout << "Domain setup completed and files generated." << std::endl;
    return 0;
}
