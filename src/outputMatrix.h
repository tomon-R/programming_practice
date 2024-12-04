// outputMatrix.h

#ifndef OUTPUT_MATRIX_H
#define OUTPUT_MATRIX_H

#include <string>
#include <vector>

#include "../solver/Matrix.h"

// Function to output a matrix to a file
template <typename T>
void outputMatrix(const std::string& fileName, const std::vector<std::vector<T>>& matrix, int precision = 4);

// Overloaded function to output a Matrix object to a file
void outputMatrix(const std::string& fileName, const Matrix& matrix, int precision = 4);

// Function to output a vector to a file
template <typename T>
void outputVector(const std::string& fileName, const std::vector<T>& vector, int precision = 4);

#endif  // OUTPUT_MATRIX_H
