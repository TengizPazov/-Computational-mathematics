#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// Функция нахождения нормы вектора
double vector_norm(const std::vector<double>& v) {
    double sum = 0.0;
    for (double val : v) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

// Функция умножения матрицы на вектор
std::vector<double> matrix_vector_mult(const std::vector<std::vector<double>>& A, 
                                      const std::vector<double>& v) {
    size_t n = A.size();
    std::vector<double> result(n, 0.0);
    
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < v.size(); j++) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}

// Функция вычитания векторов
std::vector<double> vector_subtract(const std::vector<double>& a, 
                                   const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
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
    double cos_sq = cos_w0 * cos_w0;
    double det_J = 2.0 * w0[0] + 2.0 * w0[1] / cos_sq;
    
    std::vector<std::vector<double>> inverse_J(2, std::vector<double>(2));
    
    // Вычисление элементов обратной матрицы
    inverse_J[0][0] = (1.0 / det_J) * 1.0;
    inverse_J[0][1] = (1.0 / det_J) * (-2.0 * w0[1]);
    inverse_J[1][0] = (1.0 / det_J) * (-1.0 / cos_sq);
    inverse_J[1][1] = (1.0 / det_J) * (2.0 * w0[0]);
    
    return inverse_J;
}

int main() {
    // Начальное приближение
    std::vector<double> w0 = {-0.7, -0.7};
    
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Начальное приближение: (" << w0[0] << ", " << w0[1] << ")" << std::endl;
    
    // Первая итерация
    std::vector<double> F_now = find_F(w0);
    std::vector<std::vector<double>> inv_J = inverse_Jacobi_matrix(w0);
    std::vector<double> delta = matrix_vector_mult(inv_J, F_now);
    std::vector<double> w = vector_subtract(w0, delta);
    
    int iteration = 1;
    
    // Итерационный процесс
    while (vector_norm(F_now) >= 1e-6) {
        w0 = w;
        F_now = find_F(w0);
        inv_J = inverse_Jacobi_matrix(w0);
        delta = matrix_vector_mult(inv_J, F_now);
        w = vector_subtract(w0, delta);
        iteration++;
    }   
    return 0;
}