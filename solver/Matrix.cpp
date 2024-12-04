#include "Matrix.h"

#include <algorithm>
#include <iomanip>
#include <stdexcept>

// Constructor to initialize matrix with given dimensions and zero values
Matrix::Matrix(int rows, int cols) : rows(rows), cols(cols), elements(rows, std::vector<double>(cols, 0.0)) {}

// Constructor to initialize matrix with given elements
Matrix::Matrix(const std::vector<std::vector<double>>& elements) : rows(elements.size()), cols(elements[0].size()), elements(elements) {}

// Get the number of rows
int Matrix::getRows() const { return rows; }

// Get the number of columns
int Matrix::getCols() const { return cols; }

// Transpose the matrix
Matrix Matrix::transpose() const {
    Matrix result(cols, rows);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[j][i] = elements[i][j];
        }
    }
    return result;
}

// Matrix multiplication
Matrix Matrix::multiply(const Matrix& other) const {
    if (cols != other.getRows()) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
    }
    Matrix result(rows, other.getCols());
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.getCols(); ++j) {
            for (int k = 0; k < cols; ++k) {
                result[i][j] += elements[i][k] * other[k][j];
            }
        }
    }
    return result;
}

// Solving Ax = b using Gaussian elimination
std::vector<double> Matrix::solve(const std::vector<double>& b) const {
    if (rows != cols || rows != b.size()) {
        throw std::invalid_argument("Matrix dimensions do not match for solving.");
    }

    // Create augmented matrix
    std::vector<std::vector<double>> augmented(rows, std::vector<double>(cols + 1));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            augmented[i][j] = elements[i][j];
        }
        augmented[i][cols] = b[i];
    }

    // Perform Gaussian elimination
    for (int i = 0; i < rows; ++i) {
        // Search for maximum in this column
        double maxEl = std::abs(augmented[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < rows; ++k) {
            if (std::abs(augmented[k][i]) > maxEl) {
                maxEl = std::abs(augmented[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row
        for (int k = i; k < cols + 1; ++k) {
            std::swap(augmented[maxRow][k], augmented[i][k]);
        }

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < rows; ++k) {
            double c = -augmented[k][i] / augmented[i][i];
            for (int j = i; j < cols + 1; ++j) {
                if (i == j) {
                    augmented[k][j] = 0;
                } else {
                    augmented[k][j] += c * augmented[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix
    std::vector<double> x(rows);
    for (int i = rows - 1; i >= 0; --i) {
        x[i] = augmented[i][cols] / augmented[i][i];
        for (int k = i - 1; k >= 0; --k) {
            augmented[k][cols] -= augmented[k][i] * x[i];
        }
    }
    return x;
}

// Solve Ax = b using the Conjugate Gradient method
std::vector<double> Matrix::CG(const std::vector<double>& b) const {
    if (rows != cols || rows != b.size()) {
        throw std::invalid_argument("Matrix dimensions do not match for solving.");
    }

    std::vector<double> x(rows, 0.0);  // Initial guess: x = 0
    std::vector<double> r = b;         // Initial residual: r = b - A*x (x=0 -> r=b)
    std::vector<double> p = r;         // Initial direction: p = r
    std::vector<double> Ap(rows);

    double rsold = 0.0;
    for (int i = 0; i < rows; ++i) {
        rsold += r[i] * r[i];
    }

    for (int k = 0; k < MAX_ITERATIONS; ++k) {
        // Ap = A * p
        for (int i = 0; i < rows; ++i) {
            Ap[i] = 0.0;
            for (int j = 0; j < cols; ++j) {
                Ap[i] += elements[i][j] * p[j];
            }
        }

        double alpha = 0.0;
        double pdotAp = 0.0;
        for (int i = 0; i < rows; ++i) {
            pdotAp += p[i] * Ap[i];
        }
        alpha = rsold / pdotAp;

        for (int i = 0; i < rows; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rsnew = 0.0;
        for (int i = 0; i < rows; ++i) {
            rsnew += r[i] * r[i];
        }

        std::printf("CG iteration %d: Residual = %e\n", k, std::sqrt(rsnew));

        if (std::sqrt(rsnew) < EPSILON) {
            break;
        }

        for (int i = 0; i < rows; ++i) {
            p[i] = r[i] + (rsnew / rsold) * p[i];
        }

        rsold = rsnew;
    }

    return x;
}

// Access element at index
std::vector<double>& Matrix::operator[](int index) { return elements[index]; }

const std::vector<double>& Matrix::operator[](int index) const { return elements[index]; }

// Print the matrix
void Matrix::print(const std::string& fileName) const {
    std::ofstream outFile(fileName);
    if (!outFile) {
        throw std::runtime_error("Cannot open file for writing: " + fileName);
    }
    for (const auto& row : elements) {
        for (const auto& el : row) {
            outFile << std::setw(10) << el << " ";
        }
        outFile << std::endl;
    }
}

Matrix Matrix::removeDirichletRowsAndCols(const std::vector<int>& dirichletBoundaryNodes) const {
    int newSize = rows - dirichletBoundaryNodes.size();
    Matrix result(newSize, newSize);
    int newRow = 0;

    for (int i = 0; i < rows; ++i) {
        if (std::find(dirichletBoundaryNodes.begin(), dirichletBoundaryNodes.end(), i) != dirichletBoundaryNodes.end()) {
            continue;
        }
        int newCol = 0;
        for (int j = 0; j < cols; ++j) {
            if (std::find(dirichletBoundaryNodes.begin(), dirichletBoundaryNodes.end(), j) != dirichletBoundaryNodes.end()) {
                continue;
            }
            result[newRow][newCol] = elements[i][j];
            ++newCol;
        }
        ++newRow;
    }

    return result;
}
