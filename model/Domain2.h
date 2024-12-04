#ifndef DOMAIN_H
#define DOMAIN_H

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "Element.h"
#include "Node.h"

struct Domain2d {
    // 2 basis parameter functions. Cartesian coordinates should be transformed into 2d parameter space.
    const std::function<double(std::vector<double>)> e1;
    const std::function<double(std::vector<double>)> e2;

    // Boundary Condition function for Dirichlet BC (x1,x2)
    std::function<double(std::vector<double>)> g;
    // Boundary Condition function for Neumann BC (x1,x2)
    std::function<double(std::vector<double>)> q;

    // Boundary values in (e1,e2)
    double e1Min;
    double e1Max;
    double e2Min;
    double e2Max;

    // boundary assignment. supports "dirichlet" and "neumann". to be updated...
    std::array<std::string, 4> boundaryConditions;  // {bottom, right, top, left} in (e1,e2) space.

    // array of nodes the domain holds which is described in Cartesian coordinates
    std::vector<Node> nodes;

    // solution.
    std::vector<double> u;

    // these store node ids on the boundary
    std::vector<int> dirichletBoundaryNodes;
    std::vector<std::vector<int>> neumannBoundaryNodes;

    // array of elements the domain holds
    std::vector<Element> elements;

    std::string interpolationMethod;

    // constructor
    Domain2d(const std::function<double(std::vector<double>)>& e1, const std::function<double(std::vector<double>)>& e2) : e1(e1), e2(e2) {}

    // putting boundary.
    void setBoundary(double e1Min, double e1Max, double e2Min, double e2Max, const std::function<double(std::vector<double>)>& g,
                     const std::function<double(std::vector<double>)>& q, std::array<std::string, 4> boundaryConditions);
    // setting nodes. distributes equally for now (it might be updated soon).
    void setNodes(int N, int M);
    // setting elements. should be updated with Delaunay triangulation.
    void setElements(int N, int M, std::string interpolationMethod);
    // legacy vtk
    void outputVTK(std::string fileName) const;
    void printStatus(std::string fileName) const;

    // subroutines
   private:
    Node findNode(int nodeId) const;
};

#endif  // DOMAIN_H
