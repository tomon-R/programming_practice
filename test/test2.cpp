#include <filesystem>
#include <iostream>
#include <vector>

#include "../src/outputMatrix.h"

namespace fs = std::filesystem;

int main() {
    // Example data: a matrix and a vector
    std::vector<std::vector<double>> matrixData = {{1.23456789, 2.3456789, 3.456789}, {4.56789, 5.6789, 6.789}, {7.89, 8.9, 9.055555}};
    Matrix matrix(matrixData);

    std::vector<double> vector = {1.23456789, 2.3456789, 3.456789};

    // Use the output directory defined in CMakeLists.txt
    fs::path outputDir = OUTPUT_DIR;

    std::string matrixFile = (outputDir / "matrix_output.txt").string();
    std::string vectorFile = (outputDir / "vector_output.txt").string();

    // Test outputMatrix and outputVector functions
    outputMatrix(matrixFile, matrix);  // Uses default precision of 4
    outputVector(vectorFile, vector);  // Uses default precision of 4

    std::cout << "Data written to files: " << matrixFile << " and " << vectorFile << std::endl;

    return 0;
}
