#include <vector>
#include <iostream>
#include <cmath>

//Функция для вычисления разделенных разностей
std::vector<double> dividedDiff(const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
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

//Функция для вычисления значения полинома Ньютона
double newtonPoly(const std::vector<double>& coef, const std::vector<double>& x_data, double x) {
    int n = x_data.size() - 1;
    double p = coef[n];
    
    for (int k = 1; k <= n; k++) {
        p = coef[n-k] + (x - x_data[n-k]) * p;
    }
    
    return p;
}

//Алгоритм Томаса
std::vector<double> thomas_algorithm(const std::vector<double>& h, const std::vector<double>& u) {
    int N = h.size();
    int n_unknown = N - 1;
    std::vector<double> xi(n_unknown, 0.0);
    std::vector<double> eta(n_unknown, 0.0);
    std::vector<double> c(N + 1, 0.0);
    
    if (n_unknown >= 1) {
        double beta1 = 2.0;
        double gamma1 = h[1] / (h[0] + h[1]);
        
        xi[0] = -gamma1 / beta1;
        eta[0] = (6.0 * u[0]) / beta1;
    }
    
    for (int i = 1; i < n_unknown - 1; ++i) {
        double alpha = h[i] / (h[i] + h[i + 1]);
        double beta = 2.0;
        double gamma = h[i + 1] / (h[i] + h[i + 1]);
        double denominator = beta + alpha * xi[i - 1];
        
        xi[i] = -gamma / denominator;
        eta[i] = (6.0 * u[i] - alpha * eta[i - 1]) / denominator;
    }

    if (n_unknown >= 2) {
        int i = n_unknown - 1;
        double alpha = h[N - 1] / (h[N - 2] + h[N - 1]);
        double beta = 2.0;
        double denominator = beta + alpha * xi[i - 1];
        eta[i] = (6.0 * u[i] - alpha * eta[i - 1]) / denominator;
    }
    
    if (n_unknown >= 1) {
        c[N - 1] = eta[n_unknown - 1];
        for (int i = n_unknown - 2; i >= 0; --i) {
            c[i + 1] = xi[i] * c[i + 2] + eta[i];
        }
    }
    
    
    c[0] = 0.0;
    c[N] = 0.0;
    
    return c;
}
struct SplineCoefficients {
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> d;
    std::vector<double> h;
};

//Функция для вычисления всех коэффициентов кубического сплайна
SplineCoefficients cubic_spline_coefficients(const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
    std::vector<double> h(n - 1);
    for (int i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }
    
    std::vector<double> u(n - 2);
    for (int i = 0; i < n - 2; i++) {
        double h_prev = x[i + 1] - x[i];
        double h_curr = x[i + 2] - x[i + 1];
        u[i] = ((y[i + 2] - y[i + 1]) / h_curr - (y[i + 1] - y[i]) / h_prev) / (h_prev + h_curr);
    }
    
    std::vector<double> c_vec = thomas_algorithm(h, u);
    
    std::vector<double> a(n - 1);
    std::vector<double> b(n - 1);
    std::vector<double> d(n - 1);
    std::vector<double> c_result(n - 1);
    
    for (int i = 0; i < n - 1; i++) {
        a[i] = y[i];
        d[i] = (c_vec[i + 1] - c_vec[i]) / (3 * h[i]);
        b[i] = (y[i + 1] - y[i]) / h[i] - (h[i] / 3) * (2 * c_vec[i] + c_vec[i + 1]);
        c_result[i] = c_vec[i];
    }
    
    return {a, b, c_result, d, h};
}

//Функция для вычисления значения сплайна в точке
double evaluate_spline(const std::vector<double>& x, const std::vector<double>& y, 
                      const std::vector<double>& a, const std::vector<double>& b, 
                      const std::vector<double>& c, const std::vector<double>& d, 
                      const std::vector<double>& h, double x_eval) {
    int n = x.size();
    
    for (int i = 0; i < n - 1; i++) {
        if (x[i] <= x_eval && x_eval <= x[i + 1]) {
            double dx = x_eval - x[i];
            return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        }
    }
    
    if (x_eval < x[0]) {
        double dx = x_eval - x[0];
        return a[0] + b[0] * dx + c[0] * dx * dx + d[0] * dx * dx * dx;
    } else {
        int i = n - 2;
        double dx = x_eval - x[i];
        return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
    }
}

int main() {
    // Исходные данные
    std::vector<double> x = {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000};
    std::vector<double> y = {92228496, 106021537, 123202624, 132164569, 151325798, 
                             179323175, 203211926, 226545805, 248709873, 281421906};
    
    //Вычисление коэффициентов разделенных разностей
    std::vector<double> divided_diff_coef = dividedDiff(x, y);
    
    std::cout << "Коэффициенты разделенных разностей (полином Ньютона):" << std::endl;
    for (int i = 0; i < divided_diff_coef.size(); i++) {
        std::cout << "a[" << i << "] = " << divided_diff_coef[i] << std::endl;
    }
    std::cout << std::endl;
    
    //Вычисление коэффициентов кубического сплайна
    SplineCoefficients spline_coef = cubic_spline_coefficients(x, y);
    
    std::cout << "Коэффициенты кубического сплайна:" << std::endl;
    for (int i = 0; i < spline_coef.a.size(); i++) {
        std::cout << "Отрезок [" << x[i] << ", " << x[i+1] << "]:" << std::endl;
        std::cout << "  a[" << i << "] = " << spline_coef.a[i] << std::endl;
        std::cout << "  b[" << i << "] = " << spline_coef.b[i] << std::endl;
        std::cout << "  c[" << i << "] = " << spline_coef.c[i] << std::endl;
        std::cout << "  d[" << i << "] = " << spline_coef.d[i] << std::endl;
    }
    std::cout << std::endl;
    double test_x = 2010;
    double real_2010 = 308745538;
    
    double newton_value = newtonPoly(divided_diff_coef, x, test_x);
    
    double spline_value = evaluate_spline(x, y, spline_coef.a, spline_coef.b, 
                                         spline_coef.c, spline_coef.d, spline_coef.h, test_x);
    
    std::cout << "Экстраполяция на 2010 год:" << std::endl;
    std::cout << "Полином Ньютона: " << newton_value << std::endl;
    std::cout << "Кубический сплайн: " << spline_value << std::endl;
    std::cout << "Реальное значение: " << real_2010 << std::endl;
    std::cout << std::endl;
    
    std::cout << "Ошибки экстраполяции:" << std::endl;
    std::cout << "Полином Ньютона: " << std::abs(newton_value - real_2010) << std::endl;
    std::cout << "Кубический сплайн: " << std::abs(spline_value - real_2010) << std::endl;
    
    return 0;
}