#include <vector>
#include <cmath>
//функция нахождения нормы матрицы
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
// Функция нахождения вектора F в методе Ньютона
std::vector<double> find_F(const std::vector<double>& w0) {
    double F1 = w0[0] * w0[0] + w0[1] * w0[1] - 1.0;
    double F2 = w0[1] - std::tan(w0[0]);
    return {F1, F2};
}

// Функция нахождения обратной матрицы Якоби
std::vector<std::vector<double>> inverse_Jacobi_matrix(const std::vector<double>& w0) {
    double cos_w0 = std::cos(w0[0]);
    double det_J = 2.0 * (w0[0] + w0[1] / cos_w0);
    
    std::vector<std::vector<double>> inverse_J(2, std::vector<double>(2));
    
    // Вычисление элекментов обратной матрицы
    inverse_J[0][0] = (1.0 / det_J) * 1.0;
    inverse_J[0][1] = (1.0 / det_J) * (-2.0 * w0[1]);
    inverse_J[1][0] = (1.0 / det_J) * (1.0 / (cos_w0 * cos_w0));
    inverse_J[1][1] = (1.0 / det_J) * (2.0 * w0[0]);
    
    return inverse_J;
}