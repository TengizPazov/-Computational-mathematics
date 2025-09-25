#include <vector>
#include <iostream>

std::vector<double> dividedDiff(const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
    
    // Создаем матрицу n x n
    std::vector<std::vector<double>> coef(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        coef[i][0] = y[i];
    }
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            coef[i][j] = (coef[i+1][j-1] - coef[i][j-1]) / (x[i+j] - x[i]);
        }
    }
    std::vector<double> result;
    for (int j = 0; j < n; j++) {
        result.push_back(coef[0][j]);
    }
    return result;
}