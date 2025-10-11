#include <vector>
#include <cmath>

double norm1(const std::vector<std::vector<double>>& A) {
    double maximum = 0;
    size_t n = A.size();
    size_t m = A[0].size();
    
    for (size_t line = 0; line < n; line++) {
        double sum = 0;
        for (size_t j = 0; j < m; j++) {
            sum += std::abs(A[line][j]);
        }
        if (sum > maximum) {
            maximum = sum;
        }
    }
    return maximum;
}