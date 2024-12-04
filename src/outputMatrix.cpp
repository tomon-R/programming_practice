// outputMatrix.cpp

#include "outputMatrix.h"

#include <fstream>
#include <iomanip>
#include <iostream>

template <typename T>
void outputMatrix(const std::string& fileName, const std::vector<std::vector<T>>& matrix, int precision) {
    std::ofstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing: " << fileName << std::endl;
        return;
    }

    file << std::setprecision(precision) << std::fixed;

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            file << std::setw(precision + 3) << matrix[i][j];
            if (j < matrix[i].size() - 1) {
                file << ", ";
            }
        }
        file << std::endl;
    }
    file << std::endl;

    file.close();
}

// Overloaded function to output a Matrix object to a file
void outputMatrix(const std::string& fileName, const Matrix& matrix, int precision) {
    std::ofstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing: " << fileName << std::endl;
        return;
    }

    file << std::setprecision(precision) << std::fixed;

    for (int i = 0; i < matrix.getRows(); ++i) {
        for (int j = 0; j < matrix.getCols(); ++j) {
            file << std::setw(precision + 3) << matrix[i][j];
            if (j < matrix.getCols() - 1) {
                file << ", ";
            }
        }
        file << std::endl;
    }

    file.close();
}

template <typename T>
void outputVector(const std::string& fileName, const std::vector<T>& vector, int precision) {
    std::ofstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing: " << fileName << std::endl;
        return;
    }

    file << std::setprecision(precision) << std::fixed;

    if (!vector.empty()) {
        for (size_t i = 0; i < vector.size() - 1; ++i) {
            file << std::setw(precision + 3) << vector[i] << ", " << std::endl;
        }
        file << vector.back();
    }
    file << std::endl;

    file.close();
}

// Explicit instantiation for types that will be used
template void outputMatrix(const std::string&, const std::vector<std::vector<int>>&, int);
template void outputMatrix(const std::string&, const std::vector<std::vector<double>>&, int);
template void outputVector(const std::string&, const std::vector<int>&, int);
template void outputVector(const std::string&, const std::vector<double>&, int);
