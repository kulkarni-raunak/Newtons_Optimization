#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <stdexcept>
#include "helper_function.h"

std::vector<float> multiple_mat_vec(std::vector<std::vector<float>> A, std::vector<float> B)
{
    std::vector<float> C(B.size());
    for (int i = 0; i < A.size(); i++){
        float C_temp = 0.0; 
        for (int j = 0; j < B.size(); j++){
            C_temp = C_temp + (A[i][j] * B[j]); 
        }
        C.push_back(C_temp);
    }
    return C;
}

std::vector<float> negate_vector(std::vector<float> A)
{
    std::vector<float> C;
    for (int i = 0; i < A.size(); i++){
            C.at(i) = -A.at(i); 
    }
    return C;
}

float getDeterminant(const std::vector<std::vector<float>> vect) {
    if(vect.size() != vect[0].size()) {
        throw std::runtime_error("Matrix is not quadratic");
    } 
    int dimension = vect.size();

    if(dimension == 0) {
        return 1;
    }

    if(dimension == 1) {
        return vect[0][0];
    }

    //Formula for 2x2-matrix
    if(dimension == 2) {
        return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
    }

    float result = 0;
    int sign = 1;
    for(int i = 0; i < dimension; i++) {

        //Submatrix
        std::vector<std::vector<float>> subVect(dimension - 1, std::vector<float> (dimension - 1));
        for(int m = 1; m < dimension; m++) {
            int z = 0;
            for(int n = 0; n < dimension; n++) {
                if(n != i) {
                    subVect[m-1][z] = vect[m][n];
                    z++;
                }
            }
        }

        //recursive call
        result = result + sign * vect[0][i] * getDeterminant(subVect);
        sign = -sign;
    }

    return result;
}

std::vector<std::vector<float>> getTranspose(const std::vector<std::vector<float>> matrix1) {

    //Transpose-matrix: height = width(matrix), width = height(matrix)
    std::vector<std::vector<float>> solution(matrix1[0].size(), std::vector<float> (matrix1.size()));

    //Filling solution-matrix
    for(size_t i = 0; i < matrix1.size(); i++) {
        for(size_t j = 0; j < matrix1[0].size(); j++) {
            solution[j][i] = matrix1[i][j];
        }
    }
    return solution;
}

std::vector<std::vector<float>> getCofactor(const std::vector<std::vector<float>> vect) {
    if(vect.size() != vect[0].size()) {
        throw std::runtime_error("Matrix is not quadratic");
    } 

    std::vector<std::vector<float>> solution(vect.size(), std::vector<float> (vect.size()));
    std::vector<std::vector<float>> subVect(vect.size() - 1, std::vector<float> (vect.size() - 1));

    for(std::size_t i = 0; i < vect.size(); i++) {
        for(std::size_t j = 0; j < vect[0].size(); j++) {

            int p = 0;
            for(size_t x = 0; x < vect.size(); x++) {
                if(x == i) {
                    continue;
                }
                int q = 0;

                for(size_t y = 0; y < vect.size(); y++) {
                    if(y == j) {
                        continue;
                    }

                    subVect[p][q] = vect[x][y];
                    q++;
                }
                p++;
            }
            solution[i][j] = pow(-1, i + j) * getDeterminant(subVect);
        }
    }
    return solution;
}

std::vector<std::vector<float>> getInverse(const std::vector<std::vector<float>> vect) {
    if(getDeterminant(vect) == 0) {
        throw std::runtime_error("Determinant is 0");
    } 
    float d = 1.0/getDeterminant(vect);
    std::vector<std::vector<float>> solution(vect.size(), std::vector<float> (vect.size()));

    for(size_t i = 0; i < vect.size(); i++) {
        for(size_t j = 0; j < vect.size(); j++) {
            solution[i][j] = vect[i][j] * d; 
        }
    }

    return getTranspose(getCofactor(solution));
}

void printMatrix(const std::vector<std::vector<float>> vect) {
    for(std::size_t i = 0; i < vect.size(); i++) {
        for(std::size_t j = 0; j < vect[0].size(); j++) {
            std::cout << std::setw(8) << vect[i][j] << " ";
        }
        std::cout << "\n";
    }
}

int main() {
    std::vector<std::vector<float>> matrix(3, std::vector<float> (3));
    matrix = {
        {1,2,3},
        {4,5,6},
        {7,8,8}
    };

    printMatrix(getInverse(matrix));
    return 0;
}