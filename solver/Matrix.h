#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

static const int MAX_ITERATIONS = 10000;
static const double EPSILON = 1e-7;

typedef struct Matrix {
   public:
    Matrix(int rows, int cols);
    Matrix(const std::vector<std::vector<double>>& elements);

    int getRows() const;
    int getCols() const;

    Matrix transpose() const;
    Matrix multiply(const Matrix& other) const;
    std::vector<double> solve(const std::vector<double>& b) const;
    std::vector<double> CG(const std::vector<double>& b) const;

    std::vector<double>& operator[](int index);
    const std::vector<double>& operator[](int index) const;

    void print(const std::string& fileName) const;

    Matrix removeDirichletRowsAndCols(const std::vector<int>& dirichletBoundaryNodes) const;

   private:
    int rows;
    int cols;
    std::vector<std::vector<double>> elements;
} Matrix;

#endif  // MATRIX_H